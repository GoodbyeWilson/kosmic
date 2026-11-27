# Vectorised per-study DESeq2 fit for the case-control permutation (ADR-019).
#
# 'fit_vec' follows PyDESeq2 0.5.4's 'DeseqDataSet.deseq2()' and the Wald test
# of 'DeseqStats' step by step, as array code over genes: size factors, rough
# and moment dispersions, the fitted means, gene-wise Cox-Reid dispersions, the
# parametric trend, the prior variance, MAP dispersions, the log fold changes
# by IRLS, Cook's distances with outlier replacement and refit, the Wald test
# and the Cook's p-value filter. It exists because the exact case-control
# permutation (ADR-018) re-fits every study on every shuffle and PyDESeq2 does
# that gene by gene, about 300 s per permutation for five studies; this takes
# about 30 s and returns the same fold changes, standard errors and p-values
# (tests/test_vec_deseq2.py compares them with 'fit_deseq2_wald').
#
# Where it differs from PyDESeq2 is in the one-dimensional and IRLS solvers:
#   * dispersions and, for a gene whose IRLS diverges, fold changes are found
#     with scipy's own L-BFGS-B, as PyDESeq2 does, but run on all such genes
#     in lockstep with the objective evaluated for all of them at once
#     ('vec_lbfgsb.py'). A gene on which L-BFGS-B does not converge gets
#     PyDESeq2's own grid search.
#   * the IRLS iterations for all genes are solved together.
# It does not compute shrunken fold changes or the independent filter, which
# change neither the unshrunk estimates nor the p-values the pooling reads.
#
# 'fit_study' takes the same inputs as 'kosmic.de.de_analysis.fit_deseq2_wald'
# and applies the same donor selection, expression filter and covariate rule,
# so a study fitted here is the study's own DE on the labels given.

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln, polygamma
from scipy.stats import f as f_dist
from scipy.stats import norm, trim_mean

from kosmic.meta_analysis.vec_lbfgsb import lbfgsb_batch, lbfgsb_lockstep

LN2 = np.log(2.0)
MIN_MU = 0.5
MIN_DISP = 1e-8
BETA_TOL = 1e-8
MIN_REPLICATES = 7


# Likelihood pieces

def nb_nll(y, mu, alpha):
    """Negative NB log-likelihood per gene, as 'pydeseq2.utils.nb_nll' (alpha per gene)."""
    a_inv = 1.0 / alpha
    logbinom = gammaln(y + a_inv) - gammaln(y + 1.0) - gammaln(a_inv)
    return (a_inv * np.log(alpha) - logbinom
            + (y + a_inv) * np.log(mu + a_inv) - y * np.log(mu)).sum(axis=0)


def _xtwx(X, W):
    return np.einsum('ni,ng,nj->gij', X, W, X)


def fit_alpha(y, X, mu, alpha_init, min_disp, max_disp, prior_var=None, cr_reg=True):
    """Cox-Reid (and, with 'prior_var', prior-regularised) NB dispersion for every gene.

    This is 'pydeseq2.utils.fit_alpha_mle' for all genes: scipy's L-BFGS-B on
    log(alpha), started at log(alpha_init) with bounds [log(min_disp), log(max_disp)],
    run in lockstep over genes. The loss and gradient are evaluated for all genes at
    once with PyDESeq2's own formulas in PyDESeq2's operation order (genes along rows,
    donors along the contiguous axis), which reproduces its values bit for bit. That
    matters because for a gene started at min_disp the gradient is a difference of
    terms of order 1e8, and whether it is under scipy's tolerance depends on rounding.
    A gene on which L-BFGS-B does not succeed gets PyDESeq2's grid search.

    With 'prior_var' the prior is centred on log(alpha_init), as PyDESeq2 centres it
    on the fitted trend value.
    """
    n = y.shape[0]
    lo, hi = np.log(min_disp), np.log(max_disp)
    Y = np.ascontiguousarray(y.T)
    MU = np.ascontiguousarray(mu.T)
    XT = X.T[None]

    def evaluate(idx, x):
        la = x[:, 0]
        alpha = np.exp(la)
        a_inv = 1.0 / alpha
        Yk, Mk = Y[idx], MU[idx]
        A, AL = a_inv[:, None], alpha[:, None]
        logbinom = gammaln(Yk + A) - gammaln(Yk + 1) - gammaln(A)
        nll = (n * a_inv) * np.log(alpha) + (
            -logbinom + (Yk + A) * np.log(A + Mk) - Yk * np.log(Mk)).sum(axis=1)
        ll = (a_inv ** 2) * (polygamma(0, A) - polygamma(0, Yk + A) + np.log(1 + Mk * AL)
                             + (Yk - Mk) / (Mk + A)).sum(axis=1)
        reg = 0.0
        reg_grad = 0.0
        if cr_reg:
            W = Mk / (1 + Mk * AL)
            M = np.matmul(XT * W[:, None, :], X)
            dM = np.matmul(XT * (-(W ** 2))[:, None, :], X)
            reg = reg + 0.5 * np.linalg.slogdet(M)[1]
            reg_grad = reg_grad + 0.5 * (np.linalg.inv(M) * dM).sum(axis=(1, 2)) * alpha
        if prior_var is not None:
            centre = np.log(alpha_init[idx])
            reg = reg + (la - centre) ** 2 / (2 * prior_var)
            reg_grad = reg_grad + (la - centre) / prior_var
        return nll + reg, (alpha * (-ll) + reg_grad)[:, None]

    xs, ok = lbfgsb_lockstep(evaluate, np.log(alpha_init)[:, None], lo, hi)
    out = xs[:, 0].copy()
    if (~ok).any():
        from pydeseq2.utils import grid_fit_alpha
        for k in np.flatnonzero(~ok):
            out[k] = grid_fit_alpha(y[:, k], X, mu[:, k], alpha_init[k], min_disp, max_disp,
                                    prior_var, cr_reg, prior_var is not None)
    return np.exp(out)


# IRLS

def irls(y, s, X, disp, min_mu=MIN_MU, beta_tol=BETA_TOL, maxiter=250):
    """NB GLM by IRLS for every gene ('pydeseq2.utils.irls_solver', vectorised).

    Returns the coefficients (G, p), the unthresholded fitted means (n, G), the hat
    matrix diagonal (n, G) and which genes diverged and fell back to the bounded solver.
    """
    n, n_genes = y.shape
    p = X.shape[1]
    ridge = 1e-6 * np.eye(p)[None]
    full_rank = np.linalg.matrix_rank(X) == p
    if full_rank:
        beta = np.linalg.lstsq(X, np.log(y / s[:, None] + 0.1), rcond=None)[0].T
    else:
        beta = np.zeros((n_genes, p))
        beta[:, 0] = np.log(y / s[:, None]).mean(axis=0)
    mu = np.maximum(s[:, None] * np.exp(X @ beta.T), min_mu)
    beta0 = beta.copy()
    dev = np.full(n_genes, 1000.0)
    active = np.ones(n_genes, dtype=bool)
    diverged = np.zeros(n_genes, dtype=bool)
    for _ in range(maxiter):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        ya, ma, da = y[:, idx], mu[:, idx], disp[idx]
        W = ma / (1.0 + ma * da[None, :])
        z = np.log(ma / s[:, None]) + (ya - ma) / ma
        H = _xtwx(X, W) + ridge
        rhs = np.einsum('ni,ng->gi', X, W * z)
        bh = np.linalg.solve(H, rhs[..., None])[..., 0]
        bad = (np.abs(bh) > 30).any(axis=1)
        if bad.any():
            bh[bad] = np.clip(bh[bad], -30, 30)
            diverged[idx[bad]] = True
        beta[idx] = bh
        mu_new = np.maximum(s[:, None] * np.exp(X @ bh.T), min_mu)
        mu[:, idx] = mu_new
        dev_new = -2.0 * nb_nll(ya, mu_new, da)
        ratio = np.abs(dev_new - dev[idx]) / (np.abs(dev_new) + 0.1)
        dev[idx] = dev_new
        active[idx[(ratio <= beta_tol) | bad]] = False
    # Genes whose IRLS diverged or did not converge: PyDESeq2's own solver (restart from
    # the initial coefficients, bounded L-BFGS-B, grid search for a two-column design).
    diverged |= active
    div = np.flatnonzero(diverged)
    if div.size:
        b, ok = lbfgsb_batch(y[:, div], s, X, disp[div], beta0[div], min_mu=min_mu)
        if p <= 2:                              # PyDESeq2's last resort, two-column design
            from pydeseq2.utils import grid_fit_beta
            for k in np.flatnonzero(~ok):
                b[k] = grid_fit_beta(y[:, div[k]], s, X, disp[div[k]])
        beta[div] = b
    W = np.maximum(s[:, None] * np.exp(X @ beta.T), min_mu)      # thresholded, for H
    Wt = W / (1.0 + W * disp[None, :])
    M = _xtwx(X, Wt) + ridge
    Minv = np.linalg.inv(M)
    hat = np.einsum('ni,gij,nj->ng', X, Minv, X) * Wt            # diag of W^1/2 X (X'WX)^-1 X' W^1/2
    mu_out = s[:, None] * np.exp(X @ beta.T)                      # unthresholded, as in the R code
    return beta, mu_out, hat, diverged


# Dispersion trend

def _trend_glm(cov, tgt):
    X1 = np.column_stack([np.ones(len(cov)), cov])

    def loss(c):
        m = X1 @ c
        return np.nanmean(tgt / m + np.log(m), axis=0)

    def grad(c):
        m = X1 @ c
        return -np.nanmean(((tgt / m - 1)[:, None] * X1) / m[:, None], axis=0)
    res = minimize(loss, x0=np.array([1.0, 1.0]), jac=grad, method='L-BFGS-B',
                   bounds=[(1e-12, np.inf)])
    return res.x, X1 @ res.x, res.success


def fit_trend(genewise, normed_means):
    """Parametric a0 + a1/mu trend with PyDESeq2's iteration; the mean trend on failure."""
    ok = np.isfinite(1.0 / normed_means)
    cov, tgt = (1.0 / normed_means)[ok], genewise[ok]
    old, coeffs = np.array([0.1, 0.1]), np.array([1.0, 1.0])
    while (coeffs > 1e-10).all() and (np.log(np.abs(coeffs / old)) ** 2).sum() >= 1e-6:
        old = coeffs
        coeffs, pred, converged = _trend_glm(cov, tgt)
        if not converged or (coeffs <= 1e-10).any():
            m = trim_mean(genewise[genewise > 10 * MIN_DISP], proportiontocut=0.001)
            return np.full_like(genewise, m), ('mean', m)
        ratio = tgt / pred
        keep = ~((ratio < 1e-4) | (ratio >= 15))
        cov, tgt = cov[keep], tgt[keep]
    return coeffs[0] + coeffs[1] / normed_means, ('parametric', coeffs)


def _trend_at(kind, normed_means):
    if kind[0] == 'mean':
        return np.full_like(normed_means, kind[1])
    return kind[1][0] + kind[1][1] / normed_means


def _mad(x):
    c = np.median(x)
    return np.median(np.abs(x - c)) / norm.ppf(0.75)


# Cook's distance helpers

def _group_ids(X):
    _, inv, counts = np.unique(X, axis=0, return_inverse=True, return_counts=True)
    return inv.ravel(), counts


def _n_or_more(X, k):
    inv, counts = _group_ids(X)
    return counts[inv] >= k


def _trimmed_mean(x, trim):
    s = np.sort(x, axis=0)
    n = x.shape[0]
    nt = int(np.floor(n * trim))
    return s[nt:n - nt].mean(axis=0)


def _robust_mom_disp(normed, X):
    three = _n_or_more(X, 3)
    inv, _counts = _group_ids(X)
    if three.any():
        cells = inv[three]
        sub = normed[three]
        ns = {c: int((cells == c).sum()) for c in np.unique(cells)}
        trimratio = (1 / 3, 1 / 4, 1 / 8)

        def tf(x):
            return 2 if x >= 23.5 else 1 if x >= 3.5 else 0
        sq = np.zeros_like(sub)
        for c in np.unique(cells):
            sel = cells == c
            cm = _trimmed_mean(sub[sel], trimratio[tf(ns[c])])
            sq[sel] = sub[sel] - cm[None, :]
        sq **= 2
        est = []
        for c in np.unique(cells):
            sel = cells == c
            scale = [2.04, 1.86, 1.51][tf(ns[c])]
            est.append(scale * _trimmed_mean(sq[sel], trimratio[tf(ns[c])]))
        v = np.max(est, axis=0)
    else:
        rm = _trimmed_mean(normed, 0.125)
        v = 1.51 * _trimmed_mean((normed - rm) ** 2, 0.125)
    m = normed.mean(axis=0)
    return np.maximum((v - m) / m ** 2, 0.04)


# The fit

def _dispersions(y, s, normed, X):
    """Moment dispersions, the fitted means and the gene-wise Cox-Reid MLE."""
    n = y.shape[0]
    p = X.shape[1]
    max_disp = max(10.0, n)
    beta_ls = np.linalg.lstsq(X, normed, rcond=None)[0]
    y_hat = np.maximum(X @ beta_ls, 1.0)
    rde = np.maximum((((normed - y_hat) ** 2 - y_hat) / ((n - p) * y_hat ** 2)).sum(axis=0), 0.0)
    s_mean_inv = (1.0 / s).mean()
    mu_ = normed.mean(axis=0)
    sig = normed.var(axis=0, ddof=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        mde = np.nan_to_num((sig - s_mean_inv * mu_) / mu_ ** 2)
    mom = np.clip(np.minimum(rde, mde), MIN_DISP, max_disp)
    if _group_ids(X)[1].size == p:
        mu_hat = np.maximum(s[:, None] * (X @ beta_ls), MIN_MU)
    else:
        _, mu_hat, _, _ = irls(y, s, X, mom)
    genewise = np.clip(fit_alpha(y, X, mu_hat, mom, MIN_DISP, max_disp), MIN_DISP, max_disp)
    return {'mom': mom, 'mu_hat': mu_hat, 'max_disp': max_disp, 'genewise': genewise}


def design_matrix(roles, covariates):
    """The design PyDESeq2 builds for ``~ covariates + condition``.

    Parameters
    ----------
    roles : array-like of str
        'disease' / 'control' per donor; the last column is the disease indicator.
    covariates : DataFrame or None
        The usable categorical covariates, one row per donor. Each is coded as
        indicators for every level but the first in sorted order, as PyDESeq2 does.
    """
    n = len(roles)
    cols = [np.ones(n)]
    if covariates is not None:
        for c in covariates.columns:
            v = covariates[c].astype(str).to_numpy()
            for lv in sorted(set(v))[1:]:
                cols.append((v == lv).astype(float))
    cols.append((np.asarray(roles) == 'disease').astype(float))
    return np.column_stack(cols)


def fit_vec(counts, roles, covariates=None, cooks_filter=True):
    """The per-study fit for one labelling, for all genes at once.

    Parameters
    ----------
    counts : ndarray (n, G)
        Integer counts of the donors in the contrast, after the gene filter, with no
        all-zero gene.
    roles : array-like of str
        'disease' / 'control' per donor.
    covariates : DataFrame, optional
        Usable categorical covariates, one row per donor.
    cooks_filter : bool
        Give a Cook's outlier gene no p-value (NaN), as PyDESeq2's ``cooks_filter``.

    Returns
    -------
    (ndarray, ndarray, ndarray)
        log2 fold change (disease over control), its standard error and the Wald
        p-value for each gene.

    Raises
    ------
    NotImplementedError
        When every gene has a zero count: PyDESeq2 then switches to iterative size
        factors, which this does not implement.
    """
    from pydeseq2.preprocessing import deseq2_norm_fit, deseq2_norm_transform
    n, n_genes = counts.shape
    y = counts.astype(float)
    X = design_matrix(roles, covariates)
    p = X.shape[1]

    if (y == 0).any(axis=0).all():
        raise NotImplementedError('every gene has a zero: PyDESeq2 switches to iterative '
                                  'size factors')
    logmeans, filtered = deseq2_norm_fit(y)
    normed, sf = deseq2_norm_transform(y, logmeans, filtered)
    s = np.asarray(sf, dtype=float)
    normed = np.asarray(normed, dtype=float)
    normed_means = normed.mean(axis=0)

    d = _dispersions(y, s, normed, X)
    genewise, mu_hat, max_disp = d['genewise'], d['mu_hat'], d['max_disp']
    fitted, kind = fit_trend(genewise, normed_means)
    res = np.log(genewise) - np.log(fitted)
    above = genewise >= 100 * MIN_DISP
    sq_logres = _mad(res[above]) ** 2
    prior_var = max(sq_logres - polygamma(1, (n - p) / 2.0), 0.25)
    map_disp = np.clip(fit_alpha(y, X, mu_hat, fitted, MIN_DISP, max_disp, prior_var=prior_var),
                       MIN_DISP, max_disp)
    outlier_genes = np.log(genewise) > np.log(fitted) + 2 * np.sqrt(sq_logres)
    disp = np.where(outlier_genes, genewise, map_disp)

    beta, mu_lfc, hat, _div = irls(y, s, X, disp)

    # Cook's distances, outlier replacement and refit
    disp_rob = _robust_mom_disp(normed, X)
    V = mu_lfc ** 2 * disp_rob[None, :] + mu_lfc
    cooks = (y - mu_lfc) ** 2 / V / p
    cooks = cooks * hat / (1.0 - hat) ** 2
    cutoff = f_dist.ppf(0.99, p, n - p)
    replaceable = _n_or_more(X, MIN_REPLICATES)
    refitted = np.zeros(n_genes, dtype=bool)
    new_zero = np.zeros(n_genes, dtype=bool)
    if replaceable.any():
        idx = cooks > cutoff
        replaced_genes = idx.any(axis=0)
        if replaced_genes.any():
            rg = np.flatnonzero(replaced_genes)
            sub_y = y[:, rg].copy()
            tm = _trimmed_mean(sub_y / s[:, None], 0.2)
            repl = (tm[None, :] * s[:, None]).astype(int).astype(float)
            m = replaceable[:, None] & idx[:, rg]
            sub_y[m] = repl[m]
            new_zero_sub = (sub_y == 0).all(axis=0)
            new_zero[rg[new_zero_sub]] = True
            keep_sub = ~new_zero_sub
            refitted[rg[keep_sub]] = True
            if keep_sub.any():
                ry = sub_y[:, keep_sub]
                rg2 = rg[keep_sub]
                rnormed = ry / s[:, None]
                dd = _dispersions(ry, s, rnormed, X)
                r_gw = dd['genewise']
                r_means = rnormed.mean(axis=0)
                r_fit = _trend_at(kind, r_means)
                r_map = np.clip(fit_alpha(ry, X, dd['mu_hat'], r_fit, MIN_DISP, max_disp,
                                          prior_var=prior_var), MIN_DISP, max_disp)
                r_out = np.log(r_gw) > np.log(r_fit) + 2 * np.sqrt(sq_logres)
                r_disp = np.where(r_out, r_gw, r_map)
                r_beta, _, _, _ = irls(ry, s, X, r_disp)
                beta[rg2] = r_beta
                disp[rg2] = r_disp

    # Wald test
    mu_w = s[:, None] * np.exp(X @ beta.T)
    W = mu_w / (1.0 + mu_w * disp[None, :])
    M = _xtwx(X, W)
    Hinv = np.linalg.inv(M + 1e-6 * np.eye(p)[None])
    Hc = Hinv[:, :, -1]
    se = np.sqrt(np.einsum('gi,gij,gj->g', Hc, M, Hc))
    lfc = np.where(new_zero, 0.0, beta[:, -1])
    stat = lfc / se
    pval = 2.0 * norm.sf(np.abs(stat))
    se = np.where(new_zero, 0.0, se)
    pval = np.where(new_zero, 1.0, pval)

    # Cook's p-value filter
    use_for_max = _n_or_more(X, 3)
    cks = cooks
    if refitted.any():
        cks = cooks.copy()
        cks[np.ix_(replaceable, refitted)] = 0.0
    out_mask = (cks[use_for_max, :] > cutoff).any(axis=0)
    if out_mask.any():
        pos = cooks[:, out_mask].argmax(axis=0)
        cnt = (y[:, out_mask] > y[:, out_mask][pos, np.arange(len(pos))]).sum(axis=0)
        out_mask[out_mask] = cnt < 3
    if cooks_filter:
        pval = np.where(out_mask, np.nan, pval)
    return lfc / LN2, se / LN2, pval


def fit_study(expression, sample_df, gene_names, covariates=None,
              filter_min_count=None, filter_min_samples=None, study_col=None):
    """One study's DE on the labels in 'sample_df', as 'fit_deseq2_wald' would make it.

    Takes the inputs of ``kosmic.de.de_analysis.fit_deseq2_wald`` and applies the same
    donor selection ('exclude' dropped), rounding, expression filter, zero-gene drop and
    usable-covariate rule, then fits with 'fit_vec'.

    Parameters
    ----------
    expression : ndarray (n_samples, n_genes)
        Pseudobulk sum counts.
    sample_df : DataFrame
        Must hold 'role' and the columns named in 'covariates'.
    gene_names : list of str
    covariates : list of str, optional

    Returns
    -------
    (ndarray, ndarray, ndarray) or None
        log2 fold change, standard error and p-value over 'gene_names' (NaN for a gene
        not fitted); None when fewer than two donors, or no gene, remain.
    """
    from kosmic import DE_FILTER_MIN_COUNT, DE_FILTER_MIN_SAMPLES
    from kosmic.de.de_analysis import (
        _study_values, filter_by_expression, usable_covariates)

    if filter_min_count is None:
        filter_min_count = DE_FILTER_MIN_COUNT
    if filter_min_samples is None:
        filter_min_samples = DE_FILTER_MIN_SAMPLES
    keep = sample_df['role'].isin(['disease', 'control']).values
    sample_df = sample_df[keep].reset_index(drop=True)
    counts = np.asarray(expression)[keep].round().astype(int)
    if len(sample_df) < 2:
        return None
    roles = sample_df['role'].astype(str).values
    keep_genes = filter_by_expression(
        counts, roles, min_count=filter_min_count, min_samples=filter_min_samples,
        studies=_study_values(sample_df, study_col))
    gene_idx = np.flatnonzero(keep_genes)
    if gene_idx.size == 0:
        return None
    counts = counts[:, gene_idx]
    nonzero = counts.sum(axis=0) > 0
    gene_idx, counts = gene_idx[nonzero], counts[:, nonzero]
    if gene_idx.size == 0:
        return None
    used = usable_covariates(sample_df, covariates, roles)
    lfc, se, pval = fit_vec(counts, roles, sample_df[used] if used else None)
    n_genes = len(gene_names)
    out = tuple(np.full(n_genes, np.nan) for _ in range(3))
    out[0][gene_idx], out[1][gene_idx], out[2][gene_idx] = lfc, se, pval
    return out
