# Case-control permutation calibration for meta-analysis p-values.
#
# Public entry points:
#     'consensus_cc_permutation' -- gene-level, one observed DE per method
#     'pathway_cc_permutation'   -- pathway-level equivalent
#     'load_pseudobulk'          -- read a pseudobulk CSV written by the DE tab
#     'fast_ttest_de'            -- vectorised Welch engine (used internally
#                                     + by the benchmark layer)
#
# Pipeline (gene-level): permute disease/control labels per dataset, re-run
# DE for the permuted labels, pool across studies with the chosen pooling
# method, and build an empirical null of '-log10(p)' to calibrate the
# observed analytical p-values.
#
# DE engine per permutation is selectable via 'de_method':
#     'deseq2' (default) -- NB GLM Wald test at PyDESeq2 dispersions
#                           (precomputed once per dataset; see
#                           'kosmic/meta_analysis/deseq2_fast.precompute_deseq2')
#     'welch_cpm' / 'welch_cpm_eb' / 'welch_raw' -- Welch t-test variants
#                           using 'fast_ttest_de' here
#
# Pooling is dispatched via '_MATRIX_POOLERS'; every key corresponds to a
# method exposed by the app UI (see ``kosmic/gui/meta_analysis/dialogs/
# ma_settings_dialog.py::METHOD_FAMILIES`` + the Top50/HKSJ modifiers on
# the consensus page).

from __future__ import annotations
from kosmic import MIN_STUDIES, CC_N_PERMS

import numpy as np
import pandas as pd
from scipy.stats import ecdf as scipy_ecdf


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
    roles = None
    if 'role' in df.columns:
        roles = df['role'].astype(str).values
        drop_cols.append('role')
    expression = df.drop(columns=drop_cols).values.astype(float)
    gene_names = [c for c in df.columns if c not in drop_cols]
    sample_names = list(df.index)
    return expression, conditions, gene_names, sample_names, n_cells, roles

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


def _worker_init(pseudobulk_data, precomputed,
                 is_disease_full_per_ds, two_group_full_per_ds,
                 ds_gene_indices, perm_indices,
                 n_genes, n_ds, min_studies,
                 de_method, eb_moderate, method_keys):
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
    }


def _consensus_cc_perm_one(perm_i):
    """Run a single consensus CC permutation. Reads config from '_W_STATE'.

    Returns
    -------
    dict
        method_key -> 1-D ndarray of -log10(p) for genes passing the
        min_studies filter on this permutation. Empty if none pass.
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

    for ds_i in range(n_ds):
        pc = precomputed[ds_i]
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
        valid = np.isfinite(logfcs) & np.isfinite(ses) & (ses > 0)
        effects_matrix[gi[valid], ds_i] = logfcs[valid]
        var_matrix[gi[valid], ds_i] = ses[valid] ** 2
        pval_matrix[gi[valid], ds_i] = pvals[valid]

    valid_mask = (np.isfinite(effects_matrix) &
                  np.isfinite(var_matrix) & (var_matrix > 0))
    gene_mask = valid_mask.sum(axis=1) >= min_studies

    one_result = {}
    if not gene_mask.any():
        return one_result

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

    return one_result


# Consensus CC permutation -- shared DE, three poolers

def consensus_cc_permutation(pseudobulk_data, method_meta_dfs,
                             n_perms=CC_N_PERMS, seed=42,
                             de_method='deseq2_fast',
                             min_studies=MIN_STUDIES,
                             progress_cb=None,
                             status_cb=None):
    """CC permutation calibration for consensus: shares DE across all poolers.

    Parameters
    ----------
    pseudobulk_data : list of dict
        From 'load_pseudobulk_set'. Each dict must include a per-sample
        'role' ndarray ('disease' / 'control' / 'exclude').
    method_meta_dfs : dict[str, DataFrame]
        Observed analytical results per method key (e.g. 'dl', 'sumrank_top50', 'wop').
    n_perms, seed, de_method, min_studies : standard.
    progress_cb : callable(current, total), optional
    status_cb : callable(message), optional
        Setup-phase status (DE precompute, gene-index union, pool startup) --
        the first perm doesn't complete for ~30-60s.

    Returns
    -------
    dict[str, DataFrame]
        Calibrated results per method key.
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
        obs = {}
        for _, row in mdf.iterrows():
            p = row[pval_col]
            if np.isfinite(p) and p > 0:
                obs[row['names']] = -np.log10(p)
        obs_neglogp[key] = obs

    # 'deseq2' -> PyDESeq2 dispersions once per dataset + fast_deseq2_perm
    # NB GLM Wald per perm. 'welch_*' -> Welch t-test variants.
    n_ds = len(pseudobulk_data)
    precomputed = []
    eb_moderate = de_method in ('welch_ttest_eb', 'welch_cpm_eb')
    if status_cb:
        if de_method == 'deseq2':
            import pydeseq2
            status_cb(
                f"CC perm DE engine: pydeseq2 v{pydeseq2.__version__} "
                "(DESeq2 dispersions cached, NB GLM + Wald per permutation)")
        else:
            status_cb(
                f"CC perm DE engine: Welch t-test ({de_method}, "
                f"eb_moderate={eb_moderate})")
    for pb_i, pb in enumerate(pseudobulk_data):
        expr = pb['expression']
        if de_method == 'deseq2':
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

    if status_cb:
        status_cb(f"Generating {n_perms} permutation index arrays…")
    perm_indices = [
        np.array([rng.permutation(len(pb['conditions'])) for _ in range(n_perms)])
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

    n_workers = max(1, (os.cpu_count() or 1) - 2)
    if n_perms < 20:
        n_workers = 1

    init_args = (pseudobulk_data, precomputed,
                 is_disease_full_per_ds, two_group_full_per_ds,
                 ds_gene_indices, perm_indices,
                 n_genes, n_ds, min_studies,
                 de_method, eb_moderate, method_keys)

    if n_workers <= 1:
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

        null_arr = np.concatenate(chunks)
        null_ecdf = scipy_ecdf(null_arr).cdf
        obs = obs_neglogp[key]

        pval_col = 'pvals_pooled' if 'pvals_pooled' in mdf.columns else 'pval'
        result = mdf.copy()
        result['pvals_analytical'] = result[pval_col].copy()

        calibrated = np.ones(len(result))
        for i, row in result.iterrows():
            gene = row['names']
            if gene in obs:
                cal_p = 1.0 - null_ecdf.evaluate(obs[gene])
                calibrated[i] = max(float(cal_p), 1.0 / (len(null_arr) + 1))

        result['pvals_pooled'] = calibrated
        calibrated_dfs[key] = result.sort_values('pvals_pooled').reset_index(drop=True)

    return calibrated_dfs
