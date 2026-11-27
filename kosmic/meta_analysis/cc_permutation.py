# Case-control permutation calibration for meta-analysis p-values.
#
# Public entry points:
#     'consensus_cc_permutation' -- gene-level, one observed DE per method
#     'pathway_cc_permutation'   -- pathway-level equivalent
#     'load_pseudobulk'          -- read a pseudobulk CSV written by the DE tab
#     'load_pseudobulk_covariates' -- the DE run's covariates from that CSV
#     'fast_ttest_de'            -- vectorised Welch engine (used internally
#                                     + by the benchmark layer)
#
# Pipeline (gene-level): permute disease/control labels per dataset, re-run
# DE for the permuted labels, pool across studies with the chosen pooling
# method, and build an empirical null of '-log10(p)' to calibrate the
# observed analytical p-values with the permutation p-value
# '(1 + n_extreme) / (1 + n_null)' ('permutation_p').
#
# DE engine per permutation is selectable via 'de_method' (ADR-018):
#     'deseq2_exact'     -- the study's own DE ('fit_deseq2_wald') re-run on
#                           every permutation, with its covariates, on roles
#                           shuffled within the covariate groups. Exact, and
#                           slow: about 30-40 s per study per core. Each
#                           permutation is saved so a stopped run resumes.
#     'deseq2_vec'       -- the same model and the same shuffles as the exact
#                           engine, fitted by array code over genes
#                           ('kosmic/meta_analysis/vec_deseq2.fit_study', ADR-019):
#                           the same fold changes, standard errors and p-values
#                           (tests/test_vec_deseq2.py) five to ten times
#                           faster. Needs scipy's private L-BFGS-B interface, which
#                           'vec_lbfgsb.available()' checks.
#     'deseq2_fast' ('deseq2') -- an APPROXIMATION of it: a two-group NB GLM
#                           Wald test at PyDESeq2 dispersions fitted once on
#                           the observed labels (precomputed per dataset; see
#                           'kosmic/meta_analysis/deseq2_fast.precompute_deseq2').
#                           Its null is heavier-tailed than the exact one, so
#                           it calls fewer genes; see
#                           'dev/validation/cc_exact_vs_fast.py'.
#     'welch_cpm' / 'welch_cpm_eb' / 'welch_raw' -- Welch t-test variants
#                           using 'fast_ttest_de' here
#
# Pooling is dispatched via '_MATRIX_POOLERS'; every key corresponds to a
# method exposed by the app UI (see ``kosmic/gui/meta_analysis/dialogs/
# ma_settings_dialog.py::METHOD_FAMILIES`` + the Top50/HKSJ modifiers on
# the consensus page).

from __future__ import annotations
from kosmic import MIN_STUDIES, CC_N_PERMS
from kosmic.paths import PSEUDOBULK_COVARIATE_PREFIX

import numpy as np
import pandas as pd


def load_pseudobulk(filepath: str):
    """Load a pseudobulk CSV saved by MetaExportWorker.

    Returns
    -------
    expression : ndarray (n_samples, n_genes)
        Sum counts per sample.
    conditions : ndarray of str (n_samples,)
    gene_names : list of str
    sample_names : list of str
    n_cells : ndarray of int (n_samples,) or None
        Number of cells per sample, if the column is present.
    roles : ndarray of str (n_samples,) or None
        Per-sample role assignment ('disease' / 'control' / 'exclude').
        'None' if the pseudobulk CSV has no 'role' column, in which case
        downstream pooling raises a clear error prompting the user to set
        roles on the Inspect tab and re-run DE.
    """
    df = pd.read_csv(str(filepath), index_col=0)

    conditions = df['condition'].values
    n_cells = None
    drop_cols = ['condition']
    if 'n_cells' in df.columns:
        n_cells = df['n_cells'].values.astype(int)
        drop_cols.append('n_cells')
    if 'total_counts' in df.columns:
        # Per-donor depth written by 'pseudobulk_frame'; metadata, not a gene.
        drop_cols.append('total_counts')
    roles = None
    if 'role' in df.columns:
        roles = df['role'].astype(str).values
        drop_cols.append('role')
    # Covariates of the DE run (ADR-018) sit between 'role' and the genes;
    # they are read by 'load_pseudobulk_covariates', never as genes.
    drop_cols.extend(c for c in df.columns
                     if str(c).startswith(PSEUDOBULK_COVARIATE_PREFIX))
    expression = df.drop(columns=drop_cols).values.astype(float)
    gene_names = [c for c in df.columns if c not in drop_cols]
    sample_names = list(df.index)
    return expression, conditions, gene_names, sample_names, n_cells, roles


def load_pseudobulk_covariates(filepath: str):
    """Read the sample-level covariates a pseudobulk CSV carries.

    Returns
    -------
    pandas.DataFrame or None
        One row per donor (index: the donor id, in the file's order), one
        column per covariate, named without the 'cov:' prefix. 'None' for a
        file with no covariate columns -- one written before ADR-018, or
        for a run whose design had none.
    """
    header = pd.read_csv(str(filepath), nrows=0)
    cov_cols = [c for c in header.columns
                if str(c).startswith(PSEUDOBULK_COVARIATE_PREFIX)]
    if not cov_cols:
        return None
    df = pd.read_csv(str(filepath), usecols=[header.columns[0], *cov_cols],
                     index_col=0, dtype=str)
    df.columns = [c[len(PSEUDOBULK_COVARIATE_PREFIX):] for c in df.columns]
    return df

def fast_ttest_de(expression, is_disease, two_group_mask,
                  eb_moderate=False, normalize='cpm', size_factors=None):
    """Run Welch's t-test per gene on pseudobulk expression.

    Applies normalization and log2 transformation per sample before
    running the t-test.  LogFC and SE are computed directly on the
    log2 scale (no delta method).

    Parameters
    ----------
    expression : ndarray (n_samples, n_genes)
        Sum counts per sample (pseudobulk), or mean counts if
        normalize='log2'.
    is_disease : ndarray of bool, shape (n_samples,)
        True iff this sample's role is 'disease'. Pre-resolved by
        the caller; engines never see literal labels.
    two_group_mask : ndarray of bool, shape (n_samples,)
        True iff this sample is in disease ∪ control. 'False' for
        excluded samples.
    eb_moderate : bool
        If True, apply empirical Bayes variance moderation (limma-style)
        to shrink per-gene variances toward a common prior.
    normalize : str
        'cpm' (default): CPM normalize then log2(CPM+1).
        'log2': skip CPM, just log2(x+1) -- for mean-count input.
        'deseq2': use pre-computed size factors (must pass size_factors).
    size_factors : ndarray (n_samples,), optional
        Pre-computed DESeq2-style size factors.  Required when
        normalize='deseq2'.

    Returns
    -------
    logfcs : ndarray (n_genes,)
    pvals : ndarray (n_genes,)
    ses : ndarray (n_genes,)
    """
    from scipy.stats import t as t_dist

    if normalize == 'log2':
        log2_cpm = np.log2(expression + 1)
    elif normalize == 'deseq2' and size_factors is not None:
        normalized = expression / size_factors[:, None]
        log2_cpm = np.log2(normalized + 1)
    else:
        lib_sizes = expression.sum(axis=1, keepdims=True)
        lib_sizes = np.clip(lib_sizes, 1, None)
        cpm = expression / lib_sizes * 1e6
        log2_cpm = np.log2(cpm + 1)

    disease_mask = is_disease & two_group_mask
    control_mask = (~is_disease) & two_group_mask

    disease_expr = log2_cpm[disease_mask]
    control_expr = log2_cpm[control_mask]

    n_d = int(disease_mask.sum())
    n_c = int(control_mask.sum())

    if n_d < 2 or n_c < 2:
        n_genes = expression.shape[1]
        return np.zeros(n_genes), np.ones(n_genes), np.full(n_genes, np.nan)

    mean_d = disease_expr.mean(axis=0)
    mean_c = control_expr.mean(axis=0)
    var_d = disease_expr.var(axis=0, ddof=1)
    var_c = control_expr.var(axis=0, ddof=1)

    logfcs = mean_d - mean_c  # log2-CPM difference

    if eb_moderate:
        # Limma-style squeezeVar moderation.
        df_gene = n_d + n_c - 2
        s2 = ((n_d - 1) * var_d + (n_c - 1) * var_c) / max(df_gene, 1)

        s2_pos = s2[s2 > 0]
        if len(s2_pos) >= 3:
            from scipy.special import polygamma
            log_s2 = np.log(s2_pos)
            s0_sq = np.exp(np.median(log_s2))
            var_log_s2 = np.var(log_s2, ddof=1)
            mean_trigamma = polygamma(1, df_gene / 2)
            if np.isscalar(mean_trigamma):
                mean_trigamma = float(mean_trigamma)
            else:
                mean_trigamma = float(np.mean(mean_trigamma))
            trigamma_d0_half = max(var_log_s2 - mean_trigamma, 0.01)
            d0 = min(max(2.0 / trigamma_d0_half, 1.0), 50.0)

            s2_mod = (df_gene * s2 + d0 * s0_sq) / (df_gene + d0)
            df_mod = df_gene + d0

            se_mod = np.sqrt(np.clip(s2_mod * (1.0 / n_d + 1.0 / n_c), 1e-20, None))
            t_stat = logfcs / se_mod
            pvals = 2 * t_dist.sf(np.abs(t_stat), df_mod)
            pvals = np.clip(pvals, 1e-300, 1.0)
            return logfcs, pvals, se_mod

    se_sq = var_d / n_d + var_c / n_c
    ses = np.sqrt(np.clip(se_sq, 1e-20, None))

    t_stat = logfcs / ses

    # Welch-Satterthwaite df.
    num = se_sq ** 2
    denom = (var_d / n_d) ** 2 / (n_d - 1) + (var_c / n_c) ** 2 / (n_c - 1)
    denom = np.clip(denom, 1e-20, None)
    df = num / denom

    pvals = 2 * t_dist.sf(np.abs(t_stat), df)
    pvals = np.clip(pvals, 1e-300, 1.0)

    return logfcs, pvals, ses


def fast_pathway_de(expression, is_disease, two_group_mask,
                    gene_names, pathway_gene_sets):
    """Run a Welch t-test per pathway on pseudobulk expression.

    For each pathway, computes the mean log2-CPM of its member genes per
    sample, then runs a t-test of disease vs control on those per-sample
    pathway scores.

    Parameters
    ----------
    expression : ndarray (n_samples, n_genes)
        Sum counts per sample (pseudobulk).
    is_disease : ndarray of bool, shape (n_samples,)
        True iff this sample's role is 'disease'. Pre-resolved by
        the caller.
    two_group_mask : ndarray of bool, shape (n_samples,)
        True iff this sample is in disease ∪ control.
    gene_names : list of str
        Column names corresponding to expression columns.
    pathway_gene_sets : dict[str, list[str]]
        {pathway_name: [gene_symbols]}.  Genes not present in
        gene_names are dropped.

    Returns
    -------
    pathway_names : list of str
        Pathway names that had >= 1 gene in the data.
    logfcs : ndarray (n_pathways,)
    pvals : ndarray (n_pathways,)
    ses : ndarray (n_pathways,)
    """
    from scipy.stats import t as t_dist

    lib_sizes = expression.sum(axis=1, keepdims=True)
    lib_sizes = np.clip(lib_sizes, 1, None)
    cpm = expression / lib_sizes * 1e6
    log2_cpm = np.log2(cpm + 1)

    gene_idx = {g.upper(): i for i, g in enumerate(gene_names)}

    pathway_names = []
    pathway_scores = []  # list of (n_samples,) arrays
    for pw_name, genes in pathway_gene_sets.items():
        cols = [gene_idx[g.upper()] for g in genes if g.upper() in gene_idx]
        if not cols:
            continue
        pw_score = log2_cpm[:, cols].mean(axis=1)
        pathway_names.append(pw_name)
        pathway_scores.append(pw_score)

    if not pathway_names:
        return [], np.array([]), np.array([]), np.array([])

    score_matrix = np.column_stack(pathway_scores)

    disease_mask = is_disease & two_group_mask
    control_mask = (~is_disease) & two_group_mask

    n_d = int(disease_mask.sum())
    n_c = int(control_mask.sum())

    if n_d < 2 or n_c < 2:
        n_pw = len(pathway_names)
        return (pathway_names, np.zeros(n_pw),
                np.ones(n_pw), np.full(n_pw, np.nan))

    disease_scores = score_matrix[disease_mask]
    control_scores = score_matrix[control_mask]

    mean_d = disease_scores.mean(axis=0)
    mean_c = control_scores.mean(axis=0)
    var_d = disease_scores.var(axis=0, ddof=1)
    var_c = control_scores.var(axis=0, ddof=1)

    logfcs = mean_d - mean_c
    se_sq = var_d / n_d + var_c / n_c
    ses = np.sqrt(np.clip(se_sq, 1e-20, None))
    t_stat = logfcs / ses

    num = se_sq ** 2
    denom = (var_d / n_d) ** 2 / max(n_d - 1, 1) + \
            (var_c / n_c) ** 2 / max(n_c - 1, 1)
    denom = np.clip(denom, 1e-20, None)
    df = num / denom

    pvals = 2 * t_dist.sf(np.abs(t_stat), df)
    pvals = np.clip(pvals, 1e-300, 1.0)

    return pathway_names, logfcs, pvals, ses


# Fast matrix-level pooling for CC permutation

# Pathway-level calibration. Not wired into the GUI yet (the pathway
# meta-analysis page reports analytical p-values); kept and tested so
# hypothesis mode can get the same calibration the gene page has.
def pathway_cc_permutation(pb_data, observed_dfs, pathway_gene_sets,
                           pooling_method_fns,
                           n_perms=CC_N_PERMS, seed=42, progress_cb=None):
    """Run case/control permutation calibration at the pathway level.

    For each permutation, shuffles condition labels in each study,
    re-computes per-pathway scores (mean log2-CPM of pathway genes),
    runs a per-study t-test, and pools across studies with each
    requested method.  After all permutations, calibrates the observed
    pooled p-values against the empirical null per pathway per method.

    Parameters
    ----------
    pb_data : list[dict]
        One entry per study with keys: 'expression' (n_samples, n_genes),
        'conditions' (n_samples,), 'gene_names' (list of str),
        'role' (ndarray of str: 'disease' / 'control' / 'exclude' per
        sample), 'dataset_name' (str).
    observed_dfs : dict[str, pd.Series]
        {method_key: per-pathway observed p-value Series indexed by
        pathway name}.  Used to know which pathways need calibration.
    pathway_gene_sets : dict[str, list[str]]
        {pathway_name: [gene_symbols]}.
    pooling_method_fns : dict[str, callable]
        {method_key: pool_fn}.  Each pool_fn takes a list of per-study
        result DataFrames (with names, logfoldchanges, se, pvals)
        and returns a per-pathway p-value Series indexed by pathway name.
    n_perms : int
    seed : int
    progress_cb : callable, optional
        Called with (current_perm, total_perms).

    Returns
    -------
    calibrated : dict[str, pd.Series]
        {method_key: calibrated p-value Series indexed by pathway name}.
    """
    rng = np.random.default_rng(seed)

    # Stable pathways = intersection of all per-method observed Series.
    all_pathway_names = None
    for series in observed_dfs.values():
        names = set(series.index)
        if all_pathway_names is None:
            all_pathway_names = names
        else:
            all_pathway_names &= names
    if not all_pathway_names:
        return {k: pd.Series(dtype=float) for k in observed_dfs}
    all_pathway_names = sorted(all_pathway_names)

    # -log10(p) is the calibration test statistic.
    obs_stats = {}
    for k, series in observed_dfs.items():
        s = series.reindex(all_pathway_names).astype(float)
        obs_stats[k] = -np.log10(np.clip(s.values, 1e-300, 1.0))

    null_stats = {k: np.empty((n_perms, len(all_pathway_names)))
                  for k in observed_dfs}

    # Pre-resolve role masks; per-perm code index-permutes them.
    is_disease_full_per_study = []
    two_group_full_per_study = []
    for study in pb_data:
        role = study.get('role')
        if role is None:
            raise ValueError(
                f"Pathway CC perm: study '{study.get('dataset_name', '?')}' "
                "has no 'role' column. Re-confirm the disease/control "
                "assignment in the Inspect tab and regenerate the pseudobulk."
            )
        role = np.asarray(role)
        is_dis = (role == 'disease')
        is_ctl = (role == 'control')
        is_disease_full_per_study.append(is_dis)
        two_group_full_per_study.append(is_dis | is_ctl)

    for perm_i in range(n_perms):
        # Permute role masks (label-shuffle is equivalent to mask-shuffle).
        per_study_dfs = []
        for study_i, study in enumerate(pb_data):
            n_samp = len(study['conditions'])
            perm_idx = rng.permutation(n_samp)
            is_disease_perm = is_disease_full_per_study[study_i][perm_idx]
            two_group_perm = two_group_full_per_study[study_i][perm_idx]
            pw_names, pw_logfcs, pw_pvals, pw_ses = fast_pathway_de(
                study['expression'], is_disease_perm, two_group_perm,
                study['gene_names'], pathway_gene_sets,
            )
            if not pw_names:
                continue
            per_study_dfs.append(pd.DataFrame({
                'names': pw_names,
                'logfoldchanges': pw_logfcs,
                'se': pw_ses,
                'pvals': pw_pvals,
                'pvals_adj': pw_pvals,
                'dataset': study['dataset_name'],
            }))

        if len(per_study_dfs) < 2:
            for k in observed_dfs:
                null_stats[k][perm_i, :] = np.nan
            continue

        for k, pool_fn in pooling_method_fns.items():
            try:
                null_pvals = pool_fn(per_study_dfs)
                null_pvals = null_pvals.reindex(all_pathway_names).astype(float)
                null_stats[k][perm_i, :] = -np.log10(
                    np.clip(null_pvals.values, 1e-300, 1.0))
            except (ValueError, np.linalg.LinAlgError, KeyError, ZeroDivisionError):
                null_stats[k][perm_i, :] = np.nan

        if progress_cb is not None and (perm_i + 1) % max(n_perms // 50, 1) == 0:
            progress_cb(perm_i + 1, n_perms)

    if progress_cb is not None:
        progress_cb(n_perms, n_perms)

    # p_calibrated = (1 + sum(null >= obs)) / (1 + n_valid_perms)
    calibrated = {}
    for k in observed_dfs:
        nulls = null_stats[k]
        obs = obs_stats[k]
        cal = np.ones(len(all_pathway_names))
        for i in range(len(all_pathway_names)):
            col = nulls[:, i]
            valid = np.isfinite(col)
            if valid.sum() == 0 or not np.isfinite(obs[i]):
                cal[i] = 1.0
                continue
            n_ge = (col[valid] >= obs[i]).sum()
            cal[i] = (1 + n_ge) / (1 + valid.sum())
        calibrated[k] = pd.Series(
            np.clip(cal, 1e-300, 1.0),
            index=all_pathway_names)

    return calibrated


def _apply_top_proportion_mask(var_matrix, valid_mask, proportion_top):
    """Reduce valid_mask to the top-(proportion_top * k) most-precise studies
    per gene (smallest SE). Mirrors the selection used by 'dl_fast' /
    'reml_fast' so CC-perm nulls match the observed analytical pool.
    """
    if proportion_top >= 1.0:
        return valid_mask
    se_matrix = np.where(valid_mask, np.sqrt(var_matrix), np.inf)
    rank = se_matrix.argsort(axis=1).argsort(axis=1)
    k_arr = valid_mask.sum(axis=1).astype(int)
    k_top = np.maximum(np.ceil(k_arr * proportion_top).astype(int), 2)
    k_top = np.minimum(k_top, k_arr)
    return valid_mask & (rank < k_top[:, None])


def _pool_dl_pvals(effects, var_matrix, valid_mask,
                   proportion_top=1.0, hksj=False):
    """DerSimonian-Laird pooling, per-gene 2-sided p-values.
    Supports 'proportion_top < 1.0' (Top-p% most precise studies) and
    'hksj=True' (Hartung-Knapp-Sidik-Jonkman t-based variance)."""
    from scipy.stats import norm, t as t_dist
    from kosmic.meta_analysis.pooling.pooling_core import _vectorized_dl, _hksj_se
    mask = _apply_top_proportion_mask(var_matrix, valid_mask, proportion_top)
    pooled, se_pooled, _, _, weights_re, _ = _vectorized_dl(
        effects, var_matrix, mask)
    if hksj:
        eff = np.where(mask, effects, 0.0)
        se = _hksj_se(pooled, eff, weights_re, mask)
        k = mask.sum(axis=1).astype(float)
        df = np.maximum(k - 1, 1)
        t_stat = np.where(se > 0, pooled / se, 0.0)
        pvals = 2 * t_dist.sf(np.abs(t_stat), df)
    else:
        z = np.where(se_pooled > 0, pooled / se_pooled, 0.0)
        pvals = 2 * norm.sf(np.abs(z))
    return np.clip(pvals, 1e-300, 1.0)


def _pool_reml_pvals(effects, var_matrix, valid_mask,
                     proportion_top=1.0, hksj=False):
    """REML pooling, per-gene 2-sided p-values.
    Supports 'proportion_top < 1.0' and 'hksj=True' (matches DL's API)."""
    from scipy.stats import norm, t as t_dist
    from kosmic.meta_analysis.pooling.reml import _vectorized_reml
    from kosmic.meta_analysis.pooling.pooling_core import _hksj_se
    mask = _apply_top_proportion_mask(var_matrix, valid_mask, proportion_top)
    pooled, se_pooled, _, _, weights_re = _vectorized_reml(
        effects, var_matrix, mask)
    if hksj:
        eff = np.where(mask, effects, 0.0)
        se = _hksj_se(pooled, eff, weights_re, mask)
        k = mask.sum(axis=1).astype(float)
        df = np.maximum(k - 1, 1)
        t_stat = np.where(se > 0, pooled / se, 0.0)
        pvals = 2 * t_dist.sf(np.abs(t_stat), df)
    else:
        z = np.where(se_pooled > 0, pooled / se_pooled, 0.0)
        pvals = 2 * norm.sf(np.abs(z))
    return np.clip(pvals, 1e-300, 1.0)


# --- Rank-based matrix poolers ---

def _irwin_hall_cdf_vec(x_arr, n):
    """Vectorized Irwin-Hall CDF for an array of x values at fixed n.

    Parameters
    ----------
    x_arr : ndarray (m,)
    n : int

    Returns
    -------
    cdf : ndarray (m,)
    """
    from math import comb, factorial
    x = np.asarray(x_arr, dtype=np.float64)
    result = np.zeros_like(x)
    result[x >= n] = 1.0

    mask = (x > 0) & (x < n)
    if not mask.any():
        return result

    xm = x[mask]
    total = np.zeros_like(xm)
    fn = factorial(n)
    for k_val in range(n + 1):
        diff = xm - k_val
        contrib = np.where(diff > 0, diff, 0.0) ** n
        sign = (-1) ** k_val * comb(n, k_val)
        total += sign * contrib
    result[mask] = total / fn
    return np.clip(result, 0.0, 1.0)


def _pool_sumrank_pvals(pval_matrix, effects_matrix, valid_mask, proportion_top=1.0):
    """SumRank pooling at the matrix level (no DataFrames).

    First ranks genes within each study (fractional rank 0-1),
    then sums ranks across studies and applies Irwin-Hall CDF.

    Parameters
    ----------
    pval_matrix : ndarray (n_genes, n_datasets)
    effects_matrix : ndarray (n_genes, n_datasets)  -- logFCs for direction
    valid_mask : ndarray of bool (n_genes, n_datasets)
    proportion_top : float

    Returns
    -------
    pvals : ndarray (n_genes,)
    """
    from scipy.stats import rankdata

    n_genes, n_ds = pval_matrix.shape
    pvals_out = np.ones(n_genes)

    rank_up = np.full((n_genes, n_ds), np.nan)
    rank_down = np.full((n_genes, n_ds), np.nan)

    for ds_i in range(n_ds):
        vmask = valid_mask[:, ds_i]
        n_valid = vmask.sum()
        if n_valid < 2:
            continue

        p_col = pval_matrix[vmask, ds_i]
        lfc_col = effects_matrix[vmask, ds_i]

        # One-sided directional p: small p + positive lfc -> small up-rank.
        p_up = np.where(lfc_col > 0, p_col / 2, 1.0 - p_col / 2)
        p_down = np.where(lfc_col <= 0, p_col / 2, 1.0 - p_col / 2)

        r_up = rankdata(p_up) / (n_valid + 1)
        r_down = rankdata(p_down) / (n_valid + 1)

        rank_up[vmask, ds_i] = r_up
        rank_down[vmask, ds_i] = r_down

    k_per_gene = valid_mask.sum(axis=1)

    # Group by k so the Irwin-Hall null can be vectorised per group.
    for k in range(2, n_ds + 1):
        gmask = k_per_gene == k
        if not gmask.any():
            continue

        k_eff = k
        if proportion_top < 1.0:
            k_eff = max(int(np.ceil(k * proportion_top)), 1)

        ru_group = np.where(np.isnan(rank_up[gmask]), 1.0, rank_up[gmask])
        rd_group = np.where(np.isnan(rank_down[gmask]), 1.0, rank_down[gmask])

        if k_eff < k:
            ru_sorted = np.sort(ru_group, axis=1)[:, :k_eff]
            rd_sorted = np.sort(rd_group, axis=1)[:, :k_eff]
            su = ru_sorted.sum(axis=1)
            sd = rd_sorted.sum(axis=1)
        else:
            su = ru_group.sum(axis=1)
            sd = rd_group.sum(axis=1)

        su = np.minimum(su, k_eff / 2.0)
        sd = np.minimum(sd, k_eff / 2.0)

        if proportion_top < 1.0:
            from kosmic.meta_analysis.pooling.sumrank import _top_r_sum_cdf
            p_up_arr = np.array([_top_r_sum_cdf(s, k, k_eff) for s in su])
            p_down_arr = np.array([_top_r_sum_cdf(s, k, k_eff) for s in sd])
            pvals_out[gmask] = np.minimum(2.0 * np.minimum(p_up_arr, p_down_arr), 1.0)
        else:
            p_up = _irwin_hall_cdf_vec(su, k_eff)
            p_down = _irwin_hall_cdf_vec(sd, k_eff)
            pvals_out[gmask] = np.minimum(2.0 * np.minimum(p_up, p_down), 1.0)

    return np.clip(pvals_out, 1e-300, 1.0)


def _pool_sumrank_full(pval_matrix, effects_matrix, valid_mask):
    return _pool_sumrank_pvals(pval_matrix, effects_matrix, valid_mask, 1.0)


def _pool_sumrank_top50(pval_matrix, effects_matrix, valid_mask):
    return _pool_sumrank_pvals(pval_matrix, effects_matrix, valid_mask, 0.5)




def _pool_gwop_pvals(pval_matrix, effects_matrix, var_matrix, valid_mask,
                     strict=False):
    """gwOP pooling at the matrix level -- per-gene 1/SE precision weights.

    Weights are binned to 1dp for PIT-lookup efficiency.
    """
    from scipy import stats as scipy_stats
    from scipy.stats import beta as beta_dist

    n_genes, n_ds = pval_matrix.shape
    pvals_out = np.ones(n_genes)

    # Per-gene 1/SE weights, normalised to sum to k per gene.
    se_matrix = np.sqrt(np.where(valid_mask & (var_matrix > 0), var_matrix, np.nan))
    weight_matrix = np.where(np.isfinite(se_matrix) & (se_matrix > 0),
                             1.0 / se_matrix, 1.0)
    k_per_gene = valid_mask.sum(axis=1)
    w_sums = np.where(valid_mask, weight_matrix, 0.0).sum(axis=1)
    for g in range(n_genes):
        if w_sums[g] > 0 and k_per_gene[g] > 0:
            weight_matrix[g] = weight_matrix[g] * k_per_gene[g] / w_sums[g]

    w_binned = np.round(weight_matrix, 1)

    _null_grid = np.linspace(0, 1, 5001)
    weight_cdfs = {}
    unique_w = set(w_binned[valid_mask].ravel())
    for w in unique_w:
        w_key = round(float(w), 1)
        if w_key not in weight_cdfs:
            z_orig = scipy_stats.norm.ppf(1 - _null_grid / 2)
            z_w = z_orig * np.sqrt(max(w_key, 0.1))
            p_w = np.clip(2 * scipy_stats.norm.sf(np.abs(z_w)), 1e-300, 1.0)
            sorted_idx = np.argsort(p_w)
            weight_cdfs[w_key] = (p_w[sorted_idx],
                                   np.linspace(0, 1, len(p_w)))

    # Directional gwOP: split p_up/p_down by sign of lfc and pool each
    # side separately; final p = 2 * min(p_up, p_down). Discordant genes
    # cannot achieve significance.
    p_safe = np.where(valid_mask & (pval_matrix > 0) & np.isfinite(pval_matrix),
                      pval_matrix, np.nan)
    lfc_pos = effects_matrix > 0
    p_up_arr = np.where(lfc_pos, p_safe / 2.0, 1.0 - p_safe / 2.0)
    p_down_arr = np.where(~lfc_pos, p_safe / 2.0, 1.0 - p_safe / 2.0)

    w_eff = np.where(valid_mask, weight_matrix, 1.0)

    def _pool_dir(p_dir_arr):
        z_m = scipy_stats.norm.ppf(1 - p_dir_arr / 2)
        z_w = z_m * np.sqrt(w_eff)
        p_w = np.clip(2 * scipy_stats.norm.sf(np.abs(z_w)), 1e-300, 1.0)
        u_mat = np.full((n_genes, n_ds), np.nan)
        for w_key, (p_sorted, cdf_vals) in weight_cdfs.items():
            mask = valid_mask & (np.abs(w_binned - w_key) < 0.05)
            if mask.any():
                u_mat[mask] = np.clip(
                    np.interp(p_w[mask], p_sorted, cdf_vals),
                    1e-15, 1 - 1e-15)
        u_sorted = np.sort(u_mat, axis=1)
        p_dir_out = np.ones(n_genes)
        for k in range(2, n_ds + 1):
            gmask = k_per_gene == k
            if not gmask.any():
                continue
            if strict:
                # Unanimous direction: r = k (largest-order uniform must be small).
                r = k
            else:
                r = int(np.ceil(k / 2))
            u_r = u_sorted[gmask, r - 1]
            p_dir_out[gmask] = beta_dist.cdf(u_r, r, k - r + 1)
        return p_dir_out

    p_up_out = _pool_dir(p_up_arr)
    p_down_out = _pool_dir(p_down_arr)
    pvals_out = np.clip(2.0 * np.minimum(p_up_out, p_down_out), 1e-16, 1.0)

    return pvals_out


def _pool_fisher_pvals(effects, var_matrix, valid_mask, pval_matrix=None, **kw):
    """Fisher's combined probability test (CC-permutation form).

    'pval_matrix' is supplied by the dispatcher under the 'pvals' calling
    convention. If None, p is derived from 'effects / sqrt(var_matrix)'.
    """
    from scipy.stats import chi2
    if pval_matrix is None:
        with np.errstate(divide='ignore', invalid='ignore'):
            z = np.where(var_matrix > 0,
                         effects / np.sqrt(var_matrix), 0.0)
        from scipy.stats import norm
        pval_matrix = 2 * norm.sf(np.abs(z))
    p_valid = valid_mask & np.isfinite(pval_matrix) & (pval_matrix > 0)
    k_arr = p_valid.sum(axis=1).astype(int)
    p_safe = np.where(p_valid, np.clip(pval_matrix, 1e-300, 1.0), 1.0)
    log_p = np.where(p_valid, np.log(p_safe), 0.0)
    fisher_stat = -2.0 * log_p.sum(axis=1)
    df = 2 * k_arr
    pvals = np.where(k_arr > 0, chi2.sf(fisher_stat, df), 1.0)
    return np.clip(pvals, 1e-300, 1.0)


def _pool_stouffer_pvals(effects, var_matrix, valid_mask, pval_matrix=None,
                         weighted=False, **kw):
    """Stouffer's Z combined for CC permutation."""
    from scipy.stats import norm
    if pval_matrix is None:
        with np.errstate(divide='ignore', invalid='ignore'):
            z = np.where(var_matrix > 0,
                         effects / np.sqrt(var_matrix), 0.0)
        pval_matrix = 2 * norm.sf(np.abs(z))
    p_valid = (valid_mask & np.isfinite(pval_matrix)
               & (pval_matrix > 0) & (pval_matrix < 1))
    k_arr = p_valid.sum(axis=1).astype(int)
    p_safe = np.where(p_valid,
                      np.clip(pval_matrix, 1e-300, 1 - 1e-12), 0.5)
    z_mag = norm.isf(p_safe / 2.0)
    signs = np.where(effects > 0, 1.0,
                     np.where(effects < 0, -1.0, 0.0))
    z_signed = np.where(p_valid, z_mag * signs, 0.0)
    if weighted:
        se_safe = np.where(valid_mask & (var_matrix > 0),
                           np.sqrt(var_matrix), 1.0)
        w = np.where(p_valid, 1.0 / se_safe, 0.0)
        num = (w * z_signed).sum(axis=1)
        denom = np.sqrt((w ** 2).sum(axis=1))
        pooled_z = np.where(denom > 0, num / denom, 0.0)
    else:
        pooled_z = np.where(
            k_arr > 0,
            z_signed.sum(axis=1)
            / np.sqrt(np.maximum(k_arr, 1).astype(float)),
            0.0)
    pvals = np.where(k_arr > 0, 2.0 * norm.sf(np.abs(pooled_z)), 1.0)
    return np.clip(pvals, 1e-300, 1.0)


# Dispatch: pooling_key -> (pool_fn, calling_convention).
# Rule methods (pooling_dispatch.RULE_METHODS, vote counting) have no
# p-value to calibrate and are not listed.
# Keys mirror the method options in ma_settings_dialog.METHOD_FAMILIES
# combined with the top50 / hksj modifier checkboxes.
_MATRIX_POOLERS = {
    # Effect-size families (DL/REML), 4 variants each: base / top50 / hksj / top50_hksj.
    'dl':              (lambda e, v, m: _pool_dl_pvals(e, v, m), True),
    'dl_top50':        (lambda e, v, m: _pool_dl_pvals(e, v, m, proportion_top=0.5), True),
    'dl_hksj':         (lambda e, v, m: _pool_dl_pvals(e, v, m, hksj=True), True),
    'dl_top50_hksj':   (lambda e, v, m: _pool_dl_pvals(e, v, m, proportion_top=0.5, hksj=True), True),
    'reml':            (lambda e, v, m: _pool_reml_pvals(e, v, m), True),
    'reml_top50':      (lambda e, v, m: _pool_reml_pvals(e, v, m, proportion_top=0.5), True),
    'reml_hksj':       (lambda e, v, m: _pool_reml_pvals(e, v, m, hksj=True), True),
    'reml_top50_hksj': (lambda e, v, m: _pool_reml_pvals(e, v, m, proportion_top=0.5, hksj=True), True),

    # Rank-based
    'sumrank':       (_pool_sumrank_full, 'rank'),
    'sumrank_top50': (_pool_sumrank_top50, 'rank'),
    'gwop':          (_pool_gwop_pvals, 'weighted_rank'),

    # P-value combining
    'fisher':            (_pool_fisher_pvals, 'pvals'),
    'stouffer':          (_pool_stouffer_pvals, 'pvals'),
    'stouffer_weighted': (lambda e, v, m, pval_matrix=None, **kw:
                          _pool_stouffer_pvals(e, v, m, pval_matrix=pval_matrix,
                                               weighted=True),
                          'pvals'),
}



# Worker-process state, populated once by '_worker_init'. Reading from
# globals avoids re-pickling pseudobulk / precomputed / perm-indices on
# every per-perm submit.
_W_STATE = None


# ---------------------------------------------------------------------------
# Exact engine (ADR-018)
#
# 'deseq2_exact' re-runs the study's own DE on every permutation: the same
# function the DE workspace calls ('fit_deseq2_wald'), with the design the
# run used, on labels shuffled within the covariate groups. It is the
# reference the fast engine is an approximation of; it costs about 30-40 s
# per study per permutation on one core, against well under a second.

CC_ENGINE_FAST = 'deseq2_fast'
CC_ENGINE_EXACT = 'deseq2_exact'
CC_ENGINE_VEC = 'deseq2_vec'

#: The engines that re-fit each study's own model on every permutation, with its
#: covariates and the covariate-stratified shuffle. They share the per-study state,
#: the saved permutations and the refusals of the exact engine.
_REFIT_ENGINES = (CC_ENGINE_EXACT, CC_ENGINE_VEC)


def normalise_de_method(de_method):
    """The engine a 'de_method' selector names.

    'deseq2' and 'deseq2_fast' both name the fast engine. Before ADR-018
    only 'deseq2' did, and the function's own default, 'deseq2_fast', fell
    through to the Welch t-test. 'deseq2_exact', 'deseq2_vec' and the Welch
    selectors are returned as they are.
    """
    return CC_ENGINE_FAST if de_method in ('deseq2', CC_ENGINE_FAST) else de_method


class CCPermutationCancelled(RuntimeError):
    """The caller's 'cancelled' check turned true during the run. Every
    permutation finished before that is kept on disk and is reused when the
    run is started again with the same settings."""


def shuffle_roles(role, strata, rng):
    """Permute the disease/control roles of one study.

    Only donors with a role of 'disease' or 'control' take part: an excluded
    donor is not in the analysis, so it is neither moved into the contrast
    nor given a role. With 'strata' (one integer per donor) the roles are
    permuted within each stratum, so the covariate-by-condition table of the
    study is the same on every permutation; without it, among all donors in
    the contrast.
    """
    role = np.asarray(role)
    out = role.copy()
    in_model = np.isin(role, ('disease', 'control'))
    if strata is None:
        groups = [np.flatnonzero(in_model)]
    else:
        strata = np.asarray(strata)
        groups = [np.flatnonzero(in_model & (strata == s))
                  for s in np.unique(strata[in_model])]
    for g in groups:
        if len(g) > 1:
            out[g] = role[g][rng.permutation(len(g))]
    return out


def exact_study_state(pb, seed, study_index, stratify=True):
    """Everything the exact engine needs about one study, built once.

    Raises ValueError, naming the study, when the permutation cannot
    reproduce the study's DE: its run adjusted for covariates its pseudobulk
    file does not carry, or no donor can change role.
    """
    import math
    from kosmic.de.de_analysis import usable_covariates

    name = pb.get('dataset_name') or f'study {study_index + 1}'
    role = pb.get('role')
    if role is None:
        raise ValueError(f"Exact CC perm: study '{name}' has no 'role' column.")
    role = np.asarray(role).astype(str)
    missing = list(pb.get('missing_covariates') or ())
    if missing:
        raise ValueError(
            f"Exact CC perm: the DE run of study '{name}' adjusted for "
            f"{', '.join(missing)}, but its pseudobulk file does not carry "
            "them. Repeat that study's DE run so the file is written with "
            "its covariates (ADR-018).")
    samples = pb.get('samples')
    samples = ([str(s) for s in samples] if samples is not None
               else [f'{name}_{i}' for i in range(len(role))])

    sample_df = pd.DataFrame({'sample': samples, 'role': role})
    covariates = []
    cov = pb.get('covariates')
    if cov is not None:
        cov = cov.reindex(samples)
        for c in cov.columns:
            sample_df[c] = cov[c].astype(str).to_numpy()
            covariates.append(str(c))

    in_model = np.isin(role, ('disease', 'control'))
    n_dis = int((role == 'disease').sum())
    n_ctl = int((role == 'control').sum())
    if n_dis < 1 or n_ctl < 1:
        raise ValueError(
            f"Exact CC perm: study '{name}' needs at least one disease and "
            "one control donor.")
    usable = (usable_covariates(sample_df[in_model].reset_index(drop=True),
                                covariates, role[in_model])
              if covariates else [])

    strata, n_labelings = None, math.comb(n_dis + n_ctl, n_dis)
    if stratify and usable:
        keys = sample_df[usable].astype(str).agg('|'.join, axis=1).to_numpy()
        strata = pd.factorize(keys)[0]
        n_labelings = 1
        movable = False
        for s in np.unique(strata[in_model]):
            sel = in_model & (strata == s)
            k = int((role[sel] == 'disease').sum())
            n_labelings *= math.comb(int(sel.sum()), k)
            movable = movable or (0 < k < int(sel.sum()))
        if not movable:
            raise ValueError(
                f"Exact CC perm: in study '{name}' no donor can change role "
                f"within the {', '.join(usable)} groups, so the permutation "
                "has nothing to shuffle.")

    return {
        'name': name, 'study_index': study_index, 'seed': seed,
        'expr': np.asarray(pb['expression'], dtype=float),
        'genes': list(pb['gene_names']),
        'gene_pos': {g: i for i, g in enumerate(pb['gene_names'])},
        'sample_df': sample_df, 'covariates': covariates, 'usable': usable,
        'role': role, 'strata': strata, 'n_labelings': int(n_labelings),
    }


def _exact_de(state, perm_i):
    """The study's DE on the labels of permutation 'perm_i'.

    The shuffle is seeded by (run seed, permutation, study), so a run that
    is stopped and started again, or extended, uses the same shuffles.
    Returns (log2 fold change, p-value, standard error) over the study's
    genes -- the unshrunk estimates the pooling reads -- and whether the fit
    ran. A fit that fails leaves the study missing from that permutation.
    """
    from kosmic.de.de_analysis import fit_deseq2_wald

    n = len(state['genes'])
    lfc, pvals, ses = (np.full(n, np.nan) for _ in range(3))
    rng = np.random.default_rng(
        [state['seed'], perm_i, state['study_index']])
    sample_df = state['sample_df'].copy()
    sample_df['role'] = shuffle_roles(state['role'], state['strata'], rng)
    # None = the study's usual expression filter. 0 = its CPM rule off: the
    # pseudobulk came from a DE run whose genes already passed it, and with
    # only those genes the library sizes are smaller, so the rule would drop
    # more. The summed-count rule is unchanged by a shuffle.
    fmc = state.get('filter_min_count')
    fmc_kw = {} if fmc is None else {'filter_min_count': fmc}
    if state.get('engine') == CC_ENGINE_VEC:
        # The same model and shuffle as below, fitted by array code (ADR-019).
        from kosmic.meta_analysis.vec_deseq2 import fit_study
        try:
            fit = fit_study(state['expr'], sample_df, state['genes'],
                            covariates=state['covariates'], **fmc_kw)
        except Exception:  # noqa: BLE001 -- one failed fit must not end a run
            return lfc, pvals, ses, False
        if fit is None:
            return lfc, pvals, ses, False
        lfc, ses, pvals = fit
        return lfc, pvals, ses, True
    try:
        fit = fit_deseq2_wald(state['expr'], sample_df, state['genes'],
                              covariates=state['covariates'],
                              n_cpus=1, quiet=True, **fmc_kw)
    except Exception:  # noqa: BLE001 -- one failed fit must not end a run
        return lfc, pvals, ses, False
    if fit is None:
        return lfc, pvals, ses, False
    res = fit.ds.results_df
    pos = np.array([state['gene_pos'].get(g, -1) for g in res.index])
    ok = pos >= 0
    lfc[pos[ok]] = res['log2FoldChange'].to_numpy(float)[ok]
    ses[pos[ok]] = res['lfcSE'].to_numpy(float)[ok]
    pvals[pos[ok]] = res['pvalue'].to_numpy(float)[ok]
    return lfc, pvals, ses, True


def _worker_init(pseudobulk_data, precomputed,
                 is_disease_full_per_ds, two_group_full_per_ds,
                 ds_gene_indices, perm_indices,
                 n_genes, n_ds, min_studies,
                 de_method, eb_moderate, method_keys,
                 universe_mask=None):
    """Push static-per-run data into the worker process's globals.

    Used as 'ProcessPoolExecutor(initializer=...)' (also called inline
    in the single-worker path so '_consensus_cc_perm_one' can reuse it).

    'is_disease_full_per_ds' / 'two_group_full_per_ds' are
    pre-resolved per-dataset boolean masks (full-sample indexed); per-perm
    code index-permutes them, engines never see literal labels.
    """
    global _W_STATE
    method_specs = []
    for key in method_keys:
        if key not in _MATRIX_POOLERS:
            continue
        pool_fn, kind = _MATRIX_POOLERS[key]
        if kind == 'rank':
            mode = 'rank'
        elif kind == 'weighted_rank':
            mode = 'weighted_rank'
        elif kind == 'pvals':
            mode = 'pvals'
        else:
            mode = 'effect'
        method_specs.append((key, pool_fn, mode))

    _W_STATE = {
        'pseudobulk_data': pseudobulk_data,
        'precomputed': precomputed,
        'is_disease_full_per_ds': is_disease_full_per_ds,
        'two_group_full_per_ds': two_group_full_per_ds,
        'ds_gene_indices': ds_gene_indices,
        'perm_indices': perm_indices,
        'n_genes': n_genes,
        'n_ds': n_ds,
        'min_studies': min_studies,
        'eb_moderate': eb_moderate,
        'method_specs': method_specs,
        # Boolean over the unified gene index, or None: when set, the null
        # keeps only these genes (the observed table's gene universe).
        'universe_mask': universe_mask,
    }


def _cc_perm_core(perm_i):
    """Run a single consensus CC permutation. Reads config from '_W_STATE'.

    Returns
    -------
    (dict, int)
        The dict maps method_key -> 1-D ndarray of -log10(p) for genes
        passing the min_studies filter on this permutation (empty if none
        pass). The int is the number of study fits that failed (exact
        engine only; 0 otherwise).
    """
    s = _W_STATE
    precomputed = s['precomputed']
    is_disease_full_per_ds = s['is_disease_full_per_ds']
    two_group_full_per_ds = s['two_group_full_per_ds']
    ds_gene_indices = s['ds_gene_indices']
    perm_indices = s['perm_indices']
    n_genes = s['n_genes']
    n_ds = s['n_ds']
    min_studies = s['min_studies']
    eb_moderate = s['eb_moderate']
    method_specs = s['method_specs']

    effects_matrix = np.full((n_genes, n_ds), np.nan)
    var_matrix = np.full((n_genes, n_ds), np.nan)
    pval_matrix = np.full((n_genes, n_ds), np.nan)

    failed_fits = 0
    for ds_i in range(n_ds):
        pc = precomputed[ds_i]
        if 'exact' in pc:
            logfcs, pvals, ses, fit_ok = _exact_de(pc['exact'], perm_i)
            failed_fits += int(not fit_ok)
        else:
            perm_idx = perm_indices[ds_i][perm_i]
            is_disease_perm = is_disease_full_per_ds[ds_i][perm_idx]
            two_group_perm = two_group_full_per_ds[ds_i][perm_idx]

            if 'deseq2_precomp' in pc:
                from kosmic.meta_analysis.deseq2_fast import fast_deseq2_perm
                logfcs, pvals, ses = fast_deseq2_perm(
                    pc['deseq2_precomp'], is_disease_perm, two_group_perm)
            else:
                logfcs, pvals, ses = fast_ttest_de(
                    pc['expr'], is_disease_perm, two_group_perm,
                    eb_moderate=eb_moderate,
                    normalize=pc['norm'],
                    size_factors=pc['size_factors'])

        gi = ds_gene_indices[ds_i]
        # A gene with no p-value (a Cook's outlier) is not a valid estimate
        # for any method, so the null and the observed table agree on it.
        valid = (np.isfinite(logfcs) & np.isfinite(ses) & (ses > 0)
                 & np.isfinite(pvals))
        effects_matrix[gi[valid], ds_i] = logfcs[valid]
        var_matrix[gi[valid], ds_i] = ses[valid] ** 2
        pval_matrix[gi[valid], ds_i] = pvals[valid]

    universe_mask = s.get('universe_mask')
    if universe_mask is not None:
        # Studies are still fitted on every gene (size factors, dispersion
        # trend); only the genes the observed table tests are pooled, so
        # rank-based methods rank among the same genes as the observed run.
        outside = ~universe_mask
        effects_matrix[outside] = np.nan
        var_matrix[outside] = np.nan
        pval_matrix[outside] = np.nan
    valid_mask = (np.isfinite(effects_matrix) &
                  np.isfinite(var_matrix) & (var_matrix > 0))
    gene_mask = valid_mask.sum(axis=1) >= min_studies

    one_result = {}
    if not gene_mask.any():
        return one_result, failed_fits

    for key, pool_fn, mode in method_specs:
        if mode == 'weighted_rank':
            pp = pool_fn(pval_matrix, effects_matrix, var_matrix,
                         valid_mask)
        elif mode == 'rank':
            pp = pool_fn(pval_matrix, effects_matrix, valid_mask)
        elif mode == 'pvals':
            pp = pool_fn(effects_matrix, var_matrix, valid_mask,
                         pval_matrix=pval_matrix)
        else:
            pp = pool_fn(effects_matrix, var_matrix, valid_mask)

        p = pp[gene_mask]
        vp = np.isfinite(p) & (p > 0)
        if vp.any():
            one_result[key] = -np.log10(p[vp])

    return one_result, failed_fits


def _consensus_cc_perm_one(perm_i):
    """One permutation's null values, method_key -> 1-D ndarray of
    -log10(p) (see '_cc_perm_core')."""
    return _cc_perm_core(perm_i)[0]


def _cc_perm_task(perm_i):
    """Executor entry point: the permutation index with its result."""
    one, failed_fits = _cc_perm_core(perm_i)
    return perm_i, one, failed_fits


def permutation_p(null_values, observed):
    """Permutation p-value of each observed statistic against a null sample.

    '(1 + number of null values >= observed) / (1 + number of null values)':
    the observed statistic counts as one draw of its own null, so the
    p-value is never zero, and a tie with a null value counts against the
    observed one.

    Parameters
    ----------
    null_values : array-like of float
        The null statistics (here -log10 of the pooled p-values), all genes
        and permutations together.
    observed : array-like of float
    """
    null_sorted = np.sort(np.asarray(null_values, dtype=float))
    observed = np.asarray(observed, dtype=float)
    n_ge = null_sorted.size - np.searchsorted(null_sorted, observed, side='left')
    return (1.0 + n_ge) / (null_sorted.size + 1.0)


def _fingerprint(state):
    """A short hash of one study's donors, roles, genes and counts, so saved
    permutations are reused only for the data they were made from."""
    import hashlib
    h = hashlib.sha1()
    h.update('|'.join(state['sample_df']['sample']).encode())
    h.update('|'.join(state['role']).encode())
    h.update('|'.join(state['genes']).encode())
    h.update(repr(float(state['expr'].sum())).encode())
    return h.hexdigest()[:12]


def _run_exact_permutations(init_args, n_perms, method_keys, n_workers,
                            resume_dir, signature, cancelled,
                            progress_cb, status_cb):
    """Run (or resume) the exact permutations; return the null chunks.

    Each finished permutation is written to 'resume_dir/perm_<n>.npz', so a
    stopped run keeps its work. A saved permutation is reused when the
    manifest matches this run's settings and the file holds every requested
    method; otherwise the folder is cleared and the run starts again.

    Returns
    -------
    (dict, dict)
        method_key -> list of 1-D arrays of -log10(p), one per permutation;
        and a record of the run (permutations done, resumed, failed study
        fits, seconds per permutation).
    """
    import json
    import time
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from pathlib import Path

    say = status_cb or (lambda msg: None)
    label = 'Exact' if signature['engine'] == CC_ENGINE_EXACT else 'Vectorised'
    null_chunks = {k: [] for k in method_keys}
    done = set()
    failed_total = 0
    folder = Path(resume_dir) if resume_dir else None

    def take(one):
        for k in method_keys:
            if k in one:
                null_chunks[k].append(one[k])

    if folder is not None:
        folder.mkdir(parents=True, exist_ok=True)
        manifest = folder / 'manifest.json'
        stale = False
        if manifest.exists():
            try:
                stale = json.loads(manifest.read_text()) != signature
            except (OSError, ValueError):
                stale = True
        if stale:
            for f in folder.glob('perm_*.npz'):
                f.unlink()
            say("Saved permutations were made with other settings or data; "
                "starting the permutations again.")
        manifest.write_text(json.dumps(signature, indent=1, sort_keys=True))
        for i in range(n_perms):
            f = folder / f'perm_{i:05d}.npz'
            if not f.exists():
                continue
            try:
                with np.load(f) as z:
                    if not set(method_keys) <= set(z['_methods'].tolist()):
                        continue
                    take({k: z[k] for k in method_keys if k in z.files})
                    failed_total += int(z['_failed'])
            except (OSError, ValueError, KeyError):
                continue          # unreadable (e.g. cut off mid-write): redo it
            done.add(i)
    n_resumed = len(done)
    todo = [i for i in range(n_perms) if i not in done]
    if n_resumed:
        say(f"{label} CC permutation: {n_resumed} of {n_perms} permutations "
            "already saved; continuing.")

    started = time.time()
    n_now = 0

    def finished(i, one, failed):
        nonlocal failed_total, n_now
        failed_total += failed
        n_now += 1
        done.add(i)
        take(one)
        if folder is not None:
            np.savez_compressed(
                folder / f'perm_{i:05d}.npz',
                _methods=np.array(sorted(method_keys)), _failed=failed, **one)
        if failed == signature['n_studies'] and n_now == 1:
            raise RuntimeError(
                f"{label} CC permutation: every study fit failed on the first "
                "permutation; the design cannot be fitted on shuffled labels.")
        per = (time.time() - started) / n_now
        left = per * (len(todo) - n_now) / max(1, min(n_workers, len(todo)))
        if progress_cb:
            progress_cb(len(done), n_perms)
        say(f"{label} CC permutation {len(done)}/{n_perms}; about "
            f"{left / 60:.0f} min left")

    if todo:
        workers = max(1, min(n_workers, len(todo)))
        if workers <= 1:
            _worker_init(*init_args)
            for i in todo:
                if cancelled and cancelled():
                    raise CCPermutationCancelled(
                        f"stopped after {len(done)} of {n_perms} permutations")
                finished(*_cc_perm_task(i))
        else:
            cost = ("about 30-40 s per study per permutation on one core"
                    if signature['engine'] == CC_ENGINE_EXACT
                    else "about 6-14 s per study per permutation on one core")
            say(f"{label} CC permutation: {len(todo)} permutations on "
                f"{workers} processes ({cost}).")
            executor = ProcessPoolExecutor(
                max_workers=workers, initializer=_worker_init,
                initargs=init_args)
            stopped = False
            try:
                futures = [executor.submit(_cc_perm_task, i) for i in todo]
                for fut in as_completed(futures):
                    finished(*fut.result())
                    if cancelled and cancelled():
                        stopped = True
                        break
            finally:
                executor.shutdown(wait=not stopped, cancel_futures=True)
                if stopped:
                    for proc in list(getattr(executor, '_processes', {}).values()):
                        try:
                            proc.terminate()
                        except Exception:  # noqa: BLE001 -- best effort
                            pass
            if stopped:
                raise CCPermutationCancelled(
                    f"stopped after {len(done)} of {n_perms} permutations")

    return null_chunks, {
        'n_done': len(done), 'n_resumed': n_resumed,
        'failed_fits': failed_total,
        'seconds_per_permutation': ((time.time() - started) / n_now
                                    if n_now else None),
        'n_workers': n_workers,
    }


# Consensus CC permutation -- shared DE, three poolers

def consensus_cc_permutation(pseudobulk_data, method_meta_dfs,
                             n_perms=CC_N_PERMS, seed=42,
                             de_method='deseq2_fast',
                             min_studies=MIN_STUDIES,
                             progress_cb=None,
                             status_cb=None,
                             stratify=True,
                             resume_dir=None,
                             cancelled=None,
                             n_workers=None,
                             gene_universe=None,
                             refilter=True):
    """CC permutation calibration for consensus: shares DE across all poolers.

    Parameters
    ----------
    pseudobulk_data : list of dict
        From 'load_pseudobulk_set'. Each dict must include a per-sample
        'role' ndarray ('disease' / 'control' / 'exclude'). The exact
        engine also reads 'samples' and 'covariates' (ADR-018).
    method_meta_dfs : dict[str, DataFrame]
        Observed analytical results per method key (e.g. 'dl', 'sumrank_top50', 'wop').
    n_perms, seed, min_studies : standard.
    de_method : str
        The permutation DE engine. 'deseq2_fast' (also spelled 'deseq2')
        is the fast approximation: dispersions fitted once on the observed
        labels, condition-only model. 'deseq2_exact' re-runs the study's own
        DE on every permutation. 'deseq2_vec' fits the same model on the same
        shuffles by array code, five to ten times faster (ADR-019); it raises
        RuntimeError when scipy's private L-BFGS-B interface is not usable.
        'welch_*' are the Welch variants.
    progress_cb : callable(current, total), optional
    status_cb : callable(message), optional
        Setup-phase status (DE precompute, gene-index union, pool startup) --
        the first perm doesn't complete for ~30-60s.
    stratify : bool
        Exact engine: permute within the covariate groups of the design.
    resume_dir : path, optional
        Exact engine: folder where each finished permutation is saved, so
        a stopped run continues instead of starting again.
    cancelled : callable() -> bool, optional
        Exact engine: polled after each permutation; when true the run
        stops with CCPermutationCancelled and keeps what it saved.
    n_workers : int, optional
        Exact engine: processes to use (default: cores - 2).
    gene_universe : iterable of str, optional
        The genes the observed results test. When given, each permutation
        still fits every study on all its genes (size factors and the
        dispersion trend as in the real fit), but only these genes enter
        the null distribution, and only these genes are calibrated (rows
        outside it are dropped from the returned frames). Default: every
        gene pooled from at least *min_studies* studies.
    refilter : bool
        Exact and vectorised engines: True (default) re-applies the DE
        expression filter to the pseudobulk on every permutation. False
        turns off its CPM rule (min_count 0), for pseudobulk saved by a DE
        run whose genes already passed the filter: with only those genes
        saved, library sizes are smaller than in the real fit and the CPM
        rule drops some genes the real fit tested. The summed-count rule
        stays; a shuffle cannot change it. Other engines raise ValueError
        when False.

    Returns
    -------
    dict[str, DataFrame]
        Calibrated results per method key. Each frame's
        'attrs["cc_permutation"]' records how the calibration was made.
    """
    rng = np.random.default_rng(seed)

    # Role masks come from the pseudobulk's 'role' column (populated upstream
    # by load_pseudobulk_set; this function does not read settings.json).
    is_disease_full_per_ds = []
    two_group_full_per_ds = []
    for pb in pseudobulk_data:
        role = pb.get('role')
        if role is None:
            raise ValueError(
                f"Consensus CC perm: study "
                f"'{pb.get('dataset_name', '?')}' has no 'role' column. "
                "Re-confirm the disease/control assignment in the "
                "Inspect tab and regenerate the pseudobulk."
            )
        role = np.asarray(role)
        is_dis = (role == 'disease')
        is_ctl = (role == 'control')
        is_disease_full_per_ds.append(is_dis)
        two_group_full_per_ds.append(is_dis | is_ctl)

    obs_neglogp = {}
    for key, mdf in method_meta_dfs.items():
        pval_col = 'pvals_pooled' if 'pvals_pooled' in mdf.columns else 'pval'
        p_obs = mdf[pval_col].to_numpy(float)
        usable = np.isfinite(p_obs) & (p_obs > 0)
        obs_neglogp[key] = dict(zip(mdf['names'].to_numpy()[usable],
                                    -np.log10(p_obs[usable])))

    # 'deseq2_fast' -> PyDESeq2 dispersions once per dataset +
    # fast_deseq2_perm NB GLM Wald per perm. 'deseq2_exact' -> the study's
    # own DE on every permutation. 'welch_*' -> Welch t-test variants.
    engine = normalise_de_method(de_method)
    exact = engine in _REFIT_ENGINES
    if not refilter and not exact:
        raise ValueError(
            "refilter=False needs the exact or vectorised engine; the "
            f"'{de_method}' engine has no expression filter to skip.")
    if engine == CC_ENGINE_VEC:
        from kosmic.meta_analysis.vec_lbfgsb import available
        if not available():
            raise RuntimeError(
                "The vectorised CC permutation engine needs scipy's private "
                "L-BFGS-B interface, which this scipy version does not provide "
                "in the form it expects. Use the exact engine instead.")
    n_ds = len(pseudobulk_data)
    precomputed = []
    exact_states = []
    eb_moderate = engine in ('welch_ttest_eb', 'welch_cpm_eb')
    if status_cb:
        if engine == CC_ENGINE_EXACT:
            import pydeseq2
            status_cb(
                f"CC perm DE engine: exact -- each study's own DE (pydeseq2 "
                f"v{pydeseq2.__version__}) re-run on every permutation")
        elif engine == CC_ENGINE_VEC:
            import pydeseq2
            status_cb(
                f"CC perm DE engine: vectorised -- each study's own DE model "
                f"(pydeseq2 v{pydeseq2.__version__}) fitted by array code on "
                "every permutation")
        elif engine == CC_ENGINE_FAST:
            import pydeseq2
            status_cb(
                f"CC perm DE engine: fast, an approximation (pydeseq2 "
                f"v{pydeseq2.__version__} dispersions cached, NB GLM + Wald "
                "per permutation)")
        else:
            status_cb(
                f"CC perm DE engine: Welch t-test ({de_method}, "
                f"eb_moderate={eb_moderate})")
    for pb_i, pb in enumerate(pseudobulk_data):
        expr = pb['expression']
        if exact:
            state = exact_study_state(pb, seed, pb_i, stratify=stratify)
            state['engine'] = engine
            if not refilter:
                state['filter_min_count'] = 0
            exact_states.append(state)
            precomputed.append({'exact': state, 'expr': None})
        elif engine == CC_ENGINE_FAST:
            ds_name = pb.get('dataset_name') or pb.get('study_id') or f'dataset {pb_i + 1}'
            n_samp = expr.shape[0] if hasattr(expr, 'shape') else len(pb['conditions'])
            if status_cb:
                status_cb(
                    f"Pre-computing DE: {pb_i + 1}/{n_ds} "
                    f"({ds_name}, {n_samp} samples)…")
            from kosmic.meta_analysis.deseq2_fast import precompute_deseq2
            pc = precompute_deseq2(
                expr,
                is_disease_full_per_ds[pb_i],
                two_group_full_per_ds[pb_i])
            if pc is not None:
                precomputed.append({
                    'norm': 'deseq2', 'size_factors': None,
                    'expr': expr, 'deseq2_precomp': pc,
                })
            else:
                # PyDESeq2 fit failed -- fall back to Welch on CPM.
                precomputed.append({'norm': 'cpm', 'size_factors': None, 'expr': expr})
        else:
            precomputed.append({'norm': 'cpm', 'size_factors': None, 'expr': expr})

    if status_cb:
        status_cb("Unifying gene index across studies…")
    gene_to_idx = {}
    for pb in pseudobulk_data:
        for g in pb['gene_names']:
            if g not in gene_to_idx:
                gene_to_idx[g] = len(gene_to_idx)
    n_genes = len(gene_to_idx)

    ds_gene_indices = [
        np.array([gene_to_idx[g] for g in pb['gene_names']])
        for pb in pseudobulk_data
    ]

    universe_mask = None
    universe = None
    if gene_universe is not None:
        universe = {str(g) for g in gene_universe}
        universe_mask = np.fromiter((g in universe for g in gene_to_idx),
                                    dtype=bool, count=n_genes)
        if not universe_mask.any():
            raise ValueError("gene_universe shares no gene with the pseudobulk.")
        obs_neglogp = {k: {g: v for g, v in d.items() if g in universe}
                       for k, d in obs_neglogp.items()}

    # The exact engine draws its own shuffles, seeded per permutation and
    # study, so a resumed or extended run reproduces them.
    if exact:
        perm_indices = None
    else:
        if status_cb:
            status_cb(f"Generating {n_perms} permutation index arrays…")
        perm_indices = [
            np.array([rng.permutation(len(pb['conditions']))
                      for _ in range(n_perms)])
            for pb in pseudobulk_data
        ]

    # ProcessPoolExecutor initializer ships the static state (pseudobulk,
    # precomputed dispersions, perm indices) once per worker. Each submit
    # carries only a perm index, so the parent gets a progress callback
    # after every individual permutation.
    import os
    from concurrent.futures import ProcessPoolExecutor, as_completed

    method_keys = list(method_meta_dfs.keys())
    null_chunks = {k: [] for k in method_keys}
    run_info = {}

    default_workers = max(1, (os.cpu_count() or 1) - 2)
    if n_workers is None:
        n_workers = default_workers
        if n_perms < 20 and not exact:
            n_workers = 1

    init_args = (None if exact else pseudobulk_data, precomputed,
                 is_disease_full_per_ds, two_group_full_per_ds,
                 ds_gene_indices, perm_indices,
                 n_genes, n_ds, min_studies,
                 de_method, eb_moderate, method_keys, universe_mask)

    if exact:
        import pydeseq2
        signature = {
            'format': 1, 'engine': engine, 'seed': int(seed),
            'stratify': bool(stratify), 'min_studies': int(min_studies),
            'pydeseq2': pydeseq2.__version__, 'n_studies': n_ds,
            'studies': [{'name': s['name'], 'donors': len(s['role']),
                         'genes': len(s['genes']),
                         'covariates': s['covariates'], 'usable': s['usable'],
                         'fingerprint': _fingerprint(s)}
                        for s in exact_states],
        }
        if engine == CC_ENGINE_VEC:
            import scipy
            signature['scipy'] = scipy.__version__
        # Recorded only when set, so runs with the defaults keep the
        # signature their saved permutations were written under.
        if universe is not None:
            import hashlib
            signature['gene_universe'] = {
                'genes': int(universe_mask.sum()),
                'sha1': hashlib.sha1('\n'.join(sorted(universe)).encode()).hexdigest()}
        if not refilter:
            signature['refilter'] = False
        null_chunks, run_info = _run_exact_permutations(
            init_args, n_perms, method_keys, n_workers, resume_dir,
            signature, cancelled, progress_cb, status_cb)
    elif n_workers <= 1:
        # Single-process: populate the worker globals locally and loop.
        _worker_init(*init_args)
        progress_step = max(1, n_perms // 50)
        for perm_i in range(n_perms):
            one = _consensus_cc_perm_one(perm_i)
            for k, vec in one.items():
                null_chunks[k].append(vec)
            if progress_cb and ((perm_i + 1) % progress_step == 0
                                or perm_i + 1 == n_perms):
                progress_cb(perm_i + 1, n_perms)
    else:
        if status_cb:
            status_cb(
                f"Spawning {n_workers} worker processes "
                f"(shipping cached DE state)…")
        with ProcessPoolExecutor(
                max_workers=n_workers,
                initializer=_worker_init,
                initargs=init_args) as executor:
            futures = {executor.submit(_consensus_cc_perm_one, i): i
                       for i in range(n_perms)}
            completed = 0
            failed = 0
            for future in as_completed(futures):
                try:
                    one = future.result()
                    for k, vec in one.items():
                        null_chunks[k].append(vec)
                except Exception:
                    failed += 1
                completed += 1
                if progress_cb:
                    progress_cb(completed, n_perms)

        if failed and status_cb:
            status_cb(
                f"CC perm: {failed}/{n_perms} permutations failed "
                "(worker process death likely; null distribution "
                "built from surviving permutations)")

    calibrated_dfs = {}
    for key, mdf in method_meta_dfs.items():
        chunks = null_chunks.get(key, [])
        if not chunks:
            calibrated_dfs[key] = mdf
            continue

        null_arr = np.sort(np.concatenate(chunks))
        n_null = null_arr.size
        obs = obs_neglogp[key]

        pval_col = 'pvals_pooled' if 'pvals_pooled' in mdf.columns else 'pval'
        result = mdf.copy()
        n_outside = 0
        if universe is not None:
            inside = result['names'].astype(str).isin(universe).to_numpy()
            n_outside = int((~inside).sum())
            result = result[inside]
        result['pvals_analytical'] = result[pval_col].copy()

        # The permutation p-value: (1 + null values at least as extreme as
        # the observed one) / (1 + number of null values). The observed
        # statistic counts as one draw of its own null, so a calibrated
        # p-value is never zero and never below 1 / (n + 1).
        obs_vals = np.array([obs.get(g, np.nan) for g in result['names']],
                            dtype=float)
        have = np.isfinite(obs_vals)
        calibrated = np.ones(len(result))
        calibrated[have] = permutation_p(null_arr, obs_vals[have])

        result['pvals_pooled'] = calibrated
        result = result.sort_values('pvals_pooled').reset_index(drop=True)
        result.attrs['cc_permutation'] = {
            'engine': engine, 'n_permutations_requested': int(n_perms),
            'n_permutations': int(run_info.get('n_done', len(chunks))),
            'seed': int(seed), 'stratified': bool(exact and stratify),
            'null_values': int(n_null),
            'smallest_calibrated_p': 1.0 / (n_null + 1.0),
            'approximation': engine == CC_ENGINE_FAST,
            'covariates': ({s['name']: s['usable'] for s in exact_states}
                           if exact else {}),
            'fewest_distinct_labelings': (
                min(s['n_labelings'] for s in exact_states)
                if exact else None),
            'failed_study_fits': int(run_info.get('failed_fits', 0)),
            'seconds_per_permutation': run_info.get('seconds_per_permutation'),
            'resumed_permutations': int(run_info.get('n_resumed', 0)),
            'gene_universe': (None if universe is None
                              else int(universe_mask.sum())),
            'rows_outside_universe_dropped': n_outside,
            'refilter': bool(refilter),
        }
        calibrated_dfs[key] = result

    return calibrated_dfs
