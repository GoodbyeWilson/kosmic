# Shared vectorised helpers for gene-level meta-analysis pooling.
#
# Per-method entry points live in their own files ('dl.py', 'reml.py',
# 'sumrank.py', 'rop.py', 'fisher.py', 'stouffer.py') and import
# the four shared cores defined here:
#
# - '_build_effect_matrices': list of DataFrames → gene × study matrices
# - '_vectorized_dl': DerSimonian-Laird τ² + pooled estimate
# - '_hksj_se': Hartung-Knapp-Sidik-Jonkman SE correction
# - '_build_result_df': assemble the standard output DataFrame
#
# '_vectorized_reml' lives in 'reml.py' (only REML uses it; it calls
# '_vectorized_dl' from here as a warm-start).

import numpy as np
import pandas as pd


# Matrix construction helpers

def _build_effect_matrices(datasets):
    """Build (n_genes, n_datasets) matrices of effects, variances, and metadata.

    Parameters
    ----------
    datasets : list of DataFrames
        Each with columns: names, logfoldchanges, se, dataset (optional),
        pathways (optional), pvals/pvals_adj (optional).

    Returns
    -------
    gene_names : list of str
        Sorted unique gene names.
    effects_matrix : ndarray (n_genes, n_datasets)
        Log fold changes; NaN where gene absent or invalid.
    var_matrix : ndarray (n_genes, n_datasets)
        Variances (SE^2); NaN where gene absent or invalid.
    valid_mask : ndarray bool (n_genes, n_datasets)
        True where both effect and variance are finite and SE > 0.
    pathway_map : dict
        gene name -> pathway string (first non-null found).
    study_info : dict
        gene name -> list of per-study dicts for study_effects column.
    dataset_names : list of str
        Dataset name per column.
    """
    # Collect all gene names and build index
    gene_to_idx = {}
    for df in datasets:
        for g in df['names'].unique():
            if g not in gene_to_idx:
                gene_to_idx[g] = len(gene_to_idx)

    n_genes = len(gene_to_idx)
    n_ds = len(datasets)

    effects_matrix = np.full((n_genes, n_ds), np.nan)
    var_matrix = np.full((n_genes, n_ds), np.nan)
    pval_matrix = np.full((n_genes, n_ds), np.nan)
    pathway_map = {}
    study_info = {}  # gene -> list of per-study dicts
    dataset_names = []

    for ds_i, df in enumerate(datasets):
        # Determine dataset name
        if 'dataset' in df.columns and len(df) > 0:
            ds_name = df['dataset'].iloc[0]
        else:
            ds_name = f'Study_{ds_i + 1}'
        dataset_names.append(ds_name)

        pval_col = None
        for col in ('pvals_adj', 'pvals'):
            if col in df.columns:
                pval_col = col
                break

        # Vectorized extraction per dataset
        names = df['names'].values
        logfcs = df['logfoldchanges'].values.astype(float)
        ses = df['se'].values.astype(float) if 'se' in df.columns else np.full(len(df), np.nan)
        pvals = df[pval_col].values.astype(float) if pval_col else np.ones(len(df))

        # Per-row dataset names
        if 'dataset' in df.columns:
            row_ds_names = df['dataset'].values
        else:
            row_ds_names = np.full(len(df), ds_name, dtype=object)

        has_pathways = 'pathways' in df.columns
        if has_pathways:
            pathways_col = df['pathways'].values

        # Vectorized: map gene names to matrix indices
        gene_indices = np.array([gene_to_idx[g] for g in names])

        # Valid mask for this dataset
        valid = np.isfinite(logfcs) & np.isfinite(ses) & (ses > 0)
        valid_idx = np.where(valid)[0]

        # Bulk assign into matrices
        gi_valid = gene_indices[valid_idx]
        effects_matrix[gi_valid, ds_i] = logfcs[valid_idx]
        var_matrix[gi_valid, ds_i] = ses[valid_idx] ** 2
        pval_matrix[gi_valid, ds_i] = pvals[valid_idx]

        # Build study_info for forest plots (needs per-entry dicts)
        for j in valid_idx:
            gene = names[j]
            if gene not in study_info:
                study_info[gene] = []
            study_info[gene].append({
                'logfc': float(logfcs[j]),
                'se': float(ses[j]),
                'dataset': str(row_ds_names[j]),
                'pval': float(pvals[j]),
            })

        # Pathway map (first non-null per gene)
        if has_pathways:
            for j in range(len(df)):
                gene = names[j]
                if gene not in pathway_map and pd.notna(pathways_col[j]):
                    pathway_map[gene] = pathways_col[j]

    # Build ordered gene_names list
    gene_names = [''] * n_genes
    for g, idx in gene_to_idx.items():
        gene_names[idx] = g

    valid_mask = np.isfinite(effects_matrix) & np.isfinite(var_matrix) & (var_matrix > 0)

    return (gene_names, effects_matrix, var_matrix, valid_mask,
            pathway_map, study_info, dataset_names, pval_matrix)


# Vectorized DerSimonian-Laird core

def _vectorized_dl(effects_matrix, var_matrix, valid_mask):
    """Vectorized DerSimonian-Laird for all genes at once.

    Parameters
    ----------
    effects_matrix : ndarray (n_genes, n_datasets), NaN where absent
    var_matrix : ndarray (n_genes, n_datasets), NaN where absent
    valid_mask : ndarray bool (n_genes, n_datasets)

    Returns
    -------
    pooled : ndarray (n_genes,)
    se_pooled : ndarray (n_genes,)
    tau2 : ndarray (n_genes,)
    I2 : ndarray (n_genes,)  -- fractional [0, 1], NOT percentage
    weights_re : ndarray (n_genes, n_datasets)
    Q : ndarray (n_genes,)  -- heterogeneity statistic
    """
    # Replace invalid entries: 0 effect, inf variance -> zero weight
    eff = np.where(valid_mask, effects_matrix, 0.0)
    var = np.where(valid_mask, var_matrix, np.inf)

    weights_fe = np.where(valid_mask, 1.0 / var, 0.0)
    sum_w = weights_fe.sum(axis=1)

    # Fixed-effects pooled estimate
    pooled_fe = np.where(sum_w > 0,
                         (weights_fe * eff).sum(axis=1) / sum_w,
                         0.0)

    # Q statistic
    resid = np.where(valid_mask, eff - pooled_fe[:, None], 0.0)
    Q = (weights_fe * resid ** 2).sum(axis=1)

    k = valid_mask.sum(axis=1).astype(float)
    df = k - 1

    # C for tau^2 estimation
    sum_w2 = (weights_fe ** 2).sum(axis=1)
    C = sum_w - np.where(sum_w > 0, sum_w2 / sum_w, 0.0)

    # tau^2 -- gate the divide at the operation, not at the surrounding
    # np.where, so we don't trip a RuntimeWarning on the C == 0 rows.
    tau2 = np.zeros_like(Q)
    np.divide(Q - df, C, out=tau2, where=(C > 0) & (df > 0))
    np.maximum(0, tau2, out=tau2)

    # Random-effects weights and pooled estimate
    weights_re = np.where(valid_mask,
                          1.0 / (var + tau2[:, None]),
                          0.0)
    sum_w_re = weights_re.sum(axis=1)

    pooled = np.where(sum_w_re > 0,
                      (weights_re * eff).sum(axis=1) / sum_w_re,
                      pooled_fe)
    se_pooled = np.where(sum_w_re > 0,
                         np.sqrt(1.0 / sum_w_re),
                         0.0)

    # I^2 (fractional, 0 to 1) -- gate the divide as above.
    I2 = np.zeros_like(Q)
    np.divide(Q - df, Q, out=I2, where=Q > 0)
    np.maximum(0, I2, out=I2)

    return pooled, se_pooled, tau2, I2, weights_re, Q


def _hksj_se(pooled, effects_matrix, weights_re, valid_mask):
    """Hartung-Knapp-Sidik-Jonkman variance correction.

    Replaces the standard DL/REML SE with a corrected version that
    accounts for uncertainty in the between-study variance estimate.
    Should be used with a t-distribution (k-1 df) instead of normal.

    Parameters
    ----------
    pooled : ndarray (n_genes,)
        Pooled effect estimates.
    effects_matrix : ndarray (n_genes, n_datasets)
        Per-study effects (0 where absent).
    weights_re : ndarray (n_genes, n_datasets)
        Random-effects weights (0 where absent).
    valid_mask : ndarray bool (n_genes, n_datasets)

    Returns
    -------
    se_hksj : ndarray (n_genes,)
        HKSJ-corrected standard errors.
    """
    eff = np.where(valid_mask, effects_matrix, 0.0)
    k = valid_mask.sum(axis=1).astype(float)
    df = k - 1

    # Weighted residual sum of squares
    resid = np.where(valid_mask, eff - pooled[:, None], 0.0)
    q_hk = np.where(df > 0,
                     (weights_re * resid ** 2).sum(axis=1) / df,
                     0.0)

    # Truncation: q_hk must be >= 1 to prevent the HKSJ SE from being
    # smaller than the standard DL/REML SE. Without this, genes with
    # very consistent effects across studies get artificially small SEs,
    # leading to anti-conservative p-values. This is the standard
    # correction used by R's metafor (test="knha") and Stata.
    q_hk = np.maximum(q_hk, 1.0)

    # HKSJ SE: sqrt(q_hk / sum(w_i))
    sum_w = weights_re.sum(axis=1)
    se_hksj = np.where(sum_w > 0,
                       np.sqrt(q_hk / sum_w),
                       0.0)

    return se_hksj


# Helper: build output DataFrame from vectorized results

def _build_result_df(gene_names, pooled, se_pooled, I2, k_arr,
                     valid_mask, pathway_map, study_info,
                     min_studies, hksj=False, pooling_method=None,
                     k_filter=None):
    """Build the standard output DataFrame from vectorized arrays.

    Parameters
    ----------
    gene_names : list of str
    pooled, se_pooled, I2 : ndarray (n_genes,)
    k_arr : ndarray (n_genes,) int
        Number of studies actually pooled for this gene. Drives the
        'n_studies' output column and the HKSJ t-test degrees of freedom.
        Under Top 50% this is the post-selection count.
    valid_mask : ndarray bool (n_genes, n_datasets)
    pathway_map : dict
    study_info : dict gene -> list of dicts
    min_studies : int
    hksj : bool
        If True, use t-distribution with k-1 df instead of normal.
        se_pooled should already be HKSJ-corrected.
    k_filter : ndarray (n_genes,) int, optional
        If provided, used for the 'k >= min_studies' admission filter
        instead of 'k_arr'. Top 50% paths pass the raw per-gene count
        here so admission is decided on "gene appeared in ≥ N studies"
        rather than on the post-selection pool size.

    Returns
    -------
    DataFrame with standard columns.
    """
    from scipy import stats as scipy_stats

    k_for_filter = k_arr if k_filter is None else k_filter
    mask = k_for_filter >= min_studies
    indices = np.where(mask)[0]

    _STANDARD_COLS = [
        'names', 'logfoldchanges', 'se', 'z_stat', 'ci_lower',
        'ci_upper', 'pvals_pooled', 'n_studies', 'pathways',
        'study_effects', 'heterogeneity_i2', 'pooling_method']
    if len(indices) == 0:
        return pd.DataFrame(columns=_STANDARD_COLS)

    se_idx = se_pooled[indices]
    pooled_idx = pooled[indices]

    if hksj:
        # Hartung-Knapp: t-distribution with k-1 degrees of freedom
        df_t = k_arr[indices].astype(float) - 1
        z_stat = np.where(se_idx > 0, pooled_idx / se_idx, 0.0)
        # Two-sided p-value from t-distribution
        pvals = np.where(df_t > 0,
                         2 * scipy_stats.t.sf(np.abs(z_stat), df_t),
                         1.0)
        # CI using t critical values per gene (varies with k)
        t_crit = np.where(df_t > 0,
                          scipy_stats.t.ppf(0.975, df_t),
                          1.96)
        ci_lower = pooled_idx - t_crit * se_idx
        ci_upper = pooled_idx + t_crit * se_idx
    else:
        # Standard: normal distribution
        z_stat = np.where(se_idx > 0, pooled_idx / se_idx, 0.0)
        pvals = np.clip(2 * scipy_stats.norm.sf(np.abs(z_stat)),
                        1e-300, 1.0)
        ci_lower = pooled_idx - 1.96 * se_idx
        ci_upper = pooled_idx + 1.96 * se_idx

    results = []
    for j, i in enumerate(indices):
        gene = gene_names[i]
        results.append({
            'names': gene,
            'logfoldchanges': float(pooled[i]),
            'se': float(se_pooled[i]),
            'z_stat': float(z_stat[j]),
            'ci_lower': float(ci_lower[j]),
            'ci_upper': float(ci_upper[j]),
            'pvals_pooled': float(pvals[j]),
            'n_studies': int(k_arr[i]),
            'pathways': pathway_map.get(gene),
            'study_effects': study_info.get(gene, []),
            'heterogeneity_i2': float(I2[i]),
            'pooling_method': pooling_method or 'dl',
        })

    return pd.DataFrame(results)

