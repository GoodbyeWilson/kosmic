# CC-permutation accelerator for two-group NB GLM DE.
#
# This is *not* a from-scratch DESeq2 reimplementation and *not* a drop-in
# replacement for PyDESeq2. It is a two-tier hybrid built solely to make
# the consensus CC permutation in 'kosmic/meta_analysis/cc_permutation.py'
# tractable:
#
# * 'precompute_deseq2' runs the full PyDESeq2 pipeline once per
#   dataset (Cox-Reid-adjusted dispersion MLE, parametric trend with
#   mean fallback, EB shrinkage, Cook's filtering). The expensive
#   per-gene work happens here, exactly once.
# * 'fast_deseq2_perm' is the inner loop. It reuses the cached
#   dispersions and runs IRLS + Wald entirely as broadcast numpy across
#   the whole gene matrix -- no Python per-gene loop. ~3-5 ms per call.
#
# The DE tab ('kosmic/gui/de_analysis/') does *not* use this file; it calls
# plain PyDESeq2 via 'kosmic/de/de_analysis._run_deseq2'. The file's only
# production caller is 'cc_permutation.py', hence its location here.
#
# See 'dev/simulations/validate_deseq2_fast.py' for the head-to-head
# validation against full PyDESeq2.

import numpy as np
import pandas as pd
from scipy.stats import norm

_IRLS_MAX_ITER = 50
_IRLS_TOL = 1e-8


# Size Factors (median-of-ratios)

def _compute_size_factors(counts: np.ndarray) -> np.ndarray:
    """Median-of-ratios size factor estimation (DESeq2 method).

    Parameters
    ----------
    counts : ndarray (n_samples, n_genes)
        Raw integer counts.

    Returns
    -------
    size_factors : ndarray (n_samples,)
    """
    with np.errstate(divide='ignore'):
        log_counts = np.log(counts.astype(np.float64))
    log_counts[~np.isfinite(log_counts)] = np.nan

    log_geo_means = np.nanmean(log_counts, axis=0)

    finite_mask = np.isfinite(log_geo_means)
    if finite_mask.sum() == 0:
        totals = counts.sum(axis=1).astype(np.float64)
        return totals / np.exp(np.mean(np.log(totals)))

    log_ratios = log_counts[:, finite_mask] - log_geo_means[finite_mask]
    size_factors = np.exp(np.nanmedian(log_ratios, axis=1))
    size_factors = np.maximum(size_factors, 1e-8)
    return size_factors



# LFC via IRLS (vectorised NB GLM for 2-group design)

def _fit_lfc_irls(counts, size_factors, dispersions, condition_mask,
                  max_iter=_IRLS_MAX_ITER, tol=_IRLS_TOL):
    """Fit log-fold changes using IRLS for two-group NB GLM.

    The design matrix is [1, x] where x=0 for control, x=1 for disease.
    The 2x2 weighted normal equations are solved analytically.

    Parameters
    ----------
    counts : ndarray (n_samples, n_genes)
    size_factors : ndarray (n_samples,)
    dispersions : ndarray (n_genes,)
    condition_mask : ndarray of bool (n_samples,)
        True for disease samples, False for control.

    Returns
    -------
    intercept : ndarray (n_genes,)
    lfc : ndarray (n_genes,)
        Log fold change (natural log).
    mu : ndarray (n_samples, n_genes)
    converged : ndarray of bool (n_genes,)
    """
    n_samples, n_genes = counts.shape
    y = counts.astype(np.float64)
    s = size_factors[:, np.newaxis]
    alpha = dispersions[np.newaxis, :]

    ctrl = ~condition_mask
    dis = condition_mask

    # Initialize: log(group mean / size_factor)
    with np.errstate(divide='ignore', invalid='ignore'):
        ctrl_mean = (y[ctrl] / s[ctrl]).mean(axis=0)
        dis_mean = (y[dis] / s[dis]).mean(axis=0)
    ctrl_mean = np.maximum(ctrl_mean, 1e-8)
    dis_mean = np.maximum(dis_mean, 1e-8)

    beta0 = np.log(ctrl_mean)
    beta1 = np.log(dis_mean) - np.log(ctrl_mean)

    converged_flag = np.zeros(n_genes, dtype=bool)

    for iteration in range(max_iter):
        # mu = s * exp(beta0 + beta1 * x)
        eta = np.zeros_like(y)
        eta[ctrl] = beta0[np.newaxis, :]
        eta[dis] = (beta0 + beta1)[np.newaxis, :]
        mu = s * np.exp(eta)
        mu = np.clip(mu, 1e-10, 1e15)

        # NB weights: w = mu / (1 + alpha * mu)
        w = mu / (1.0 + alpha * mu)

        # Working residuals
        working_resid = (y - mu) / mu

        # Components of X'WX
        w_ctrl = w[ctrl]
        w_dis = w[dis]
        S00 = w_ctrl.sum(axis=0) + w_dis.sum(axis=0)
        S01 = w_dis.sum(axis=0)
        S11 = w_dis.sum(axis=0)

        # X'W(z) where z = eta + (y-mu)/mu
        wz_ctrl = (w_ctrl * (eta[ctrl] + working_resid[ctrl])).sum(axis=0)
        wz_dis = (w_dis * (eta[dis] + working_resid[dis])).sum(axis=0)
        rhs0 = wz_ctrl + wz_dis
        rhs1 = wz_dis

        # Solve 2x2
        det = S00 * S11 - S01 * S01
        det = np.where(np.abs(det) < 1e-20, 1e-20, det)

        beta0_new = (S11 * rhs0 - S01 * rhs1) / det
        beta1_new = (S00 * rhs1 - S01 * rhs0) / det

        delta = np.maximum(np.abs(beta0_new - beta0), np.abs(beta1_new - beta1))
        newly_converged = delta < tol
        converged_flag = converged_flag | newly_converged

        beta0 = beta0_new
        beta1 = beta1_new

        if converged_flag.all():
            break

    # Final mu
    eta = np.zeros_like(y)
    eta[ctrl] = beta0[np.newaxis, :]
    eta[dis] = (beta0 + beta1)[np.newaxis, :]
    mu = s * np.exp(eta)
    mu = np.clip(mu, 1e-10, 1e15)

    return beta0, beta1, mu, converged_flag



# Wald Test

def _wald_test(counts, mu, dispersions, condition_mask, lfc):
    """Wald test for the LFC coefficient.

    Parameters
    ----------
    counts : ndarray (n_samples, n_genes)
    mu : ndarray (n_samples, n_genes)
    dispersions : ndarray (n_genes,)
    condition_mask : ndarray of bool (n_samples,)
    lfc : ndarray (n_genes,)

    Returns
    -------
    se : ndarray (n_genes,)
    stat : ndarray (n_genes,)
    pvals : ndarray (n_genes,)
    """
    alpha = dispersions[np.newaxis, :]
    w = mu / (1.0 + alpha * mu)

    S00 = w.sum(axis=0)
    S01 = w[condition_mask].sum(axis=0)
    S11 = S01  # x^2 = x for binary

    det = S00 * S11 - S01**2
    det = np.where(np.abs(det) < 1e-20, 1e-20, det)

    # (X'WX)^{-1}[1,1] = S00 / det
    var_lfc = np.maximum(S00 / det, 1e-20)
    se = np.sqrt(var_lfc)

    with np.errstate(divide='ignore', invalid='ignore'):
        stat = np.where(se > 1e-15, lfc / se, 0.0)

    pvals = 2.0 * norm.sf(np.abs(stat))
    pvals = np.clip(pvals, 1e-300, 1.0)

    return se, stat, pvals



# Pre-compute (PyDESeq2 one-off) + Fast Permutation (vectorised inner loop)

def precompute_deseq2(expression, is_disease_full, two_group_full):
    """Pre-compute size factors and dispersions for CC permutation.

    Uses PyDESeq2's full DESeq2 pipeline for the one-time per-dataset
    fit (MoM init, Cox-Reid-adjusted MLE via L-BFGS-B, parametric-trend
    with mean fallback, EB shrinkage, outlier replacement, Cook's
    filtering). The result is cached and handed to 'fast_deseq2_perm',
    which runs the vectorised NB GLM IRLS + Wald per permutation
    (~3-5 ms per call).

    Parameters
    ----------
    expression : ndarray (n_all_samples, n_genes)
        Sum counts (pseudobulk).
    is_disease_full : ndarray of bool, shape (n_all_samples,)
        True iff the sample's role is 'disease' on the original
        labels. Pre-resolved by the caller.
    two_group_full : ndarray of bool, shape (n_all_samples,)
        True iff the sample is in disease ∪ control on the original
        labels. Excluded samples are False.

    Returns
    -------
    dict or None
        Returns 'None' if data is insufficient. Otherwise a dict with::

            counts          : ndarray (n_all_samples, n_kept_genes)
                              float64 counts on the post-filter gene axis,
                              over ALL samples (perms reshuffle labels
                              across the full sample set, not just the
                              two-group subset).
            size_factors    : ndarray (n_all_samples,)
                              median-of-ratios on the post-filter genes.
            dispersions     : ndarray (n_kept_genes,)
                              PyDESeq2 MAP dispersions.
            gene_mask       : ndarray of bool (n_genes,)
                              Which genes survived the pre-fit filter
                              AND PyDESeq2's all-zero filter.
            n_genes_total   : int
                              Original gene-axis length, for output array
                              sizing in 'fast_deseq2_perm'.
            sample_mask     : ndarray of bool (n_all_samples,)
                              Which samples are in {disease, control};
                              other groups are excluded.
    """
    from pydeseq2.dds import DeseqDataSet

    n_genes = expression.shape[1]

    if two_group_full.sum() < 2:
        return None

    counts_two_group = expression[two_group_full].astype(np.int64)
    is_disease = is_disease_full[two_group_full]

    if is_disease.sum() < 1 or (~is_disease).sum() < 1:
        return None

    # Gene filter (same criterion PyDESeq2 applies internally + our
    # min_expressing_samples rule).
    gene_mask = (((counts_two_group > 0).sum(axis=0) >= 2)
                 & (counts_two_group.sum(axis=0) > 0))
    if gene_mask.sum() == 0:
        return None

    counts_f = counts_two_group[:, gene_mask]

    # Delegate dispersion fitting to PyDESeq2.  Synthesise generic
    # 'disease' / 'control' labels from the resolved role mask --
    # PyDESeq2 still needs string labels for its design matrix; the
    # role model just guarantees we never expose lab-specific labels.
    synth_conds = np.where(is_disease, 'disease', 'control')
    sample_ids = [f'S{i}' for i in range(len(synth_conds))]
    gene_ids = [f'G{i}' for i in range(counts_f.shape[1])]
    counts_df = pd.DataFrame(counts_f, index=sample_ids, columns=gene_ids)
    metadata_df = pd.DataFrame(
        {'condition': pd.Categorical(
            synth_conds, categories=['control', 'disease'])},
        index=sample_ids)

    dds = DeseqDataSet(
        counts=counts_df, metadata=metadata_df,
        design='~condition', quiet=True, refit_cooks=True)
    dds.deseq2()

    # PyDESeq2 drops all-zero genes internally; map non_zero back to
    # the pre-filter gene indices so the final gene_mask refers to the
    # original gene axis of 'expression'.
    non_zero = dds.var['non_zero'].values
    py_dispersions = np.asarray(dds.var['dispersions'])[non_zero]
    kept_gene_indices = np.where(gene_mask)[0][non_zero]
    final_gene_mask = np.zeros(n_genes, dtype=bool)
    final_gene_mask[kept_gene_indices] = True

    # Size factors and counts over ALL samples (permutations may reshuffle
    # labels across all samples, not just the two-group subset).
    full_counts_f = expression[:, final_gene_mask].astype(np.float64)
    full_size_factors = _compute_size_factors(full_counts_f)

    return {
        'counts': full_counts_f,
        'size_factors': full_size_factors,
        'dispersions': py_dispersions,
        'gene_mask': final_gene_mask,
        'n_genes_total': n_genes,
        'sample_mask': two_group_full,  # which samples are disease/control
    }


def fast_deseq2_perm(precomp, is_disease, two_group_mask):
    """Run NB GLM + Wald test using pre-computed dispersions.

    Inner loop for CC permutation -- re-fits the GLM coefficients and
    runs the Wald test using the dispersions cached by
    'precompute_deseq2'. Genes that PyDESeq2 would handle via
    independent filtering or Cook's outlier replacement (zero counts in
    one group, IRLS non-convergence) are returned as NaN instead of
    nonsense β; the consensus CC permutation loop drops NaN entries via
    its 'np.isfinite' filter, mirroring PyDESeq2's behaviour of
    excluding such genes from the FDR pool.

    Parameters
    ----------
    precomp : dict from 'precompute_deseq2'
    is_disease : ndarray of bool, shape (n_all_samples,)
        True iff this sample's role is 'disease' in this
        permutation. Pre-resolved by the caller; engines never see
        literal disease/control labels.
    two_group_mask : ndarray of bool, shape (n_all_samples,)
        True iff this sample is in disease ∪ control in this
        permutation. 'False' for excluded samples.

    Returns
    -------
    logfcs, pvals, ses : ndarray (n_genes_total,)
        NaN for genes the engine cannot fit reliably.
    """
    counts_f = precomp['counts']       # (all_samples, kept_genes)
    size_factors = precomp['size_factors']
    dispersions = precomp['dispersions']
    gene_mask = precomp['gene_mask']
    n_genes = precomp['n_genes_total']

    nan_out = (np.full(n_genes, np.nan),) * 3

    if two_group_mask.sum() < 2:
        return nan_out

    counts_sub = counts_f[two_group_mask]
    sf_sub = size_factors[two_group_mask]
    is_disease = is_disease[two_group_mask]

    if is_disease.sum() < 1 or (~is_disease).sum() < 1:
        return nan_out

    # Pre-filter ill-conditioned genes (β → ±∞ in NB GLM): all-zero
    # counts in either group. PyDESeq2 handles these via independent
    # filtering; our fast path marks them missing.
    ctrl_zero = counts_sub[~is_disease].sum(axis=0) == 0
    dis_zero = counts_sub[is_disease].sum(axis=0) == 0
    bad = ctrl_zero | dis_zero

    # 10 IRLS iterations: brings unconverged fraction on well-behaved
    # genes from ~70% (3 iter) to <1% on simulated and real pseudobulk.
    _, lfc_nat, mu_final, converged = _fit_lfc_irls(
        counts_sub, sf_sub, dispersions, is_disease, max_iter=10)

    se_nat, _, pvals_f = _wald_test(
        counts_sub, mu_final, dispersions, is_disease, lfc_nat)

    bad = bad | (~converged) | ~np.isfinite(lfc_nat) | ~np.isfinite(pvals_f)
    lfc_nat = np.where(bad, np.nan, lfc_nat)
    se_nat = np.where(bad, np.nan, se_nat)
    pvals_f = np.where(bad, np.nan, pvals_f)

    logfcs = np.full(n_genes, np.nan)
    pvals = np.full(n_genes, np.nan)
    ses = np.full(n_genes, np.nan)

    idx = np.where(gene_mask)[0]
    logfcs[idx] = lfc_nat / np.log(2)
    pvals[idx] = pvals_f
    ses[idx] = se_nat / np.log(2)

    return logfcs, pvals, ses
