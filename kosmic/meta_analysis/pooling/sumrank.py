# SumRank meta-analysis — rank-sum method with Irwin-Hall null.
#
# Ranks each per-study p-value (signed by logFC direction) within its
# dataset, sums the per-gene ranks across datasets, and compares the
# observed sum against the analytical Irwin-Hall distribution. Supports
# a Top-r proportion variant (sum only the r smallest ranks per gene).
#
# For proper permutation calibration (case/control shuffling), the
# SumRank matrix pooler is registered in
# 'kosmic.meta_analysis.cc_permutation' and invoked via the GUI's
# "CC permutation" calibration option.
from __future__ import annotations
from kosmic import MIN_STUDIES

import numpy as np
import pandas as pd


# Rank matrices (gene × study)

def _build_rank_matrices(datasets):
    """Build (n_genes, n_datasets) rank matrices using vectorized pandas."""
    gene_to_idx = {}
    for df in datasets:
        for g in df['names'].unique():
            if g not in gene_to_idx:
                gene_to_idx[g] = len(gene_to_idx)

    n_genes = len(gene_to_idx)
    n_ds = len(datasets)
    rank_up = np.full((n_genes, n_ds), np.nan)
    rank_down = np.full((n_genes, n_ds), np.nan)
    logfc_matrix = np.full((n_genes, n_ds), np.nan)
    pathway_map = {}
    study_info = {}
    dataset_names = []

    for ds_i, df in enumerate(datasets):
        df = df.copy()
        pval_col = 'pvals' if 'pvals' in df.columns else 'pvals_adj'
        if pval_col not in df.columns:
            dataset_names.append(f'Study_{ds_i + 1}')
            continue

        ds_name = (df['dataset'].iloc[0] if 'dataset' in df.columns and len(df) > 0
                   else f'Study_{ds_i + 1}')
        dataset_names.append(ds_name)

        pvals = df[pval_col].clip(1e-300, 1.0)
        signs = np.sign(df['logfoldchanges'].fillna(0))
        signs[signs == 0] = 1
        signed_score = signs * (-np.log10(pvals))

        n = len(df)
        r_up = (signed_score.rank(ascending=False).values - 1) / max(n - 1, 1)
        r_down = (signed_score.rank(ascending=True).values - 1) / max(n - 1, 1)

        names_arr = df['names'].values
        lfc_arr = df['logfoldchanges'].values
        se_arr = df['se'].values if 'se' in df.columns else np.full(n, np.nan)

        gene_indices = np.array([gene_to_idx[g] for g in names_arr])
        rank_up[gene_indices, ds_i] = r_up
        rank_down[gene_indices, ds_i] = r_down
        logfc_matrix[gene_indices, ds_i] = lfc_arr

        for j in range(n):
            gene = names_arr[j]
            if gene not in study_info:
                study_info[gene] = []
            study_info[gene].append({
                'logfc': float(lfc_arr[j]),
                'se': float(se_arr[j]) if np.isfinite(se_arr[j]) else 0.0,
                'dataset': ds_name,
                'signed_rank': float((1.0 - 2.0 * r_up[j]) if r_up[j] < r_down[j]
                                     else -(1.0 - 2.0 * r_down[j])),
            })

        if 'pathways' in df.columns:
            pw_col = df['pathways'].values
            for j in range(n):
                gene = names_arr[j]
                if gene not in pathway_map and pd.notna(pw_col[j]):
                    pathway_map[gene] = pw_col[j]

    gene_names = [''] * n_genes
    for g, idx in gene_to_idx.items():
        gene_names[idx] = g

    n_studies = np.sum(~np.isnan(rank_up), axis=1).astype(int)
    return (gene_names, rank_up, rank_down, logfc_matrix,
            n_studies, pathway_map, study_info, dataset_names)


# Null distributions (Irwin-Hall + top-r Monte Carlo)

# Cache for top-r order statistic null lookup tables.
# Key: (k, r) -> (sorted_sums, None) for fast np.searchsorted.
_top_r_lookup_cache = {}


def _build_top_r_lookup(k, r, n_sim=5_000_000):
    """Build an empirical CDF lookup for sum of r smallest from k Uniforms.

    Uses Monte Carlo simulation cached for reuse. 5M samples gives
    minimum p-value resolution of 2e-7 (sufficient for BH FDR).
    Verified against the exact integral formula to <0.001 accuracy.
    """
    rng = np.random.default_rng(42 + k * 100 + r)
    draws = rng.random((n_sim, k))
    draws.sort(axis=1)
    sums = draws[:, :r].sum(axis=1)
    sums.sort()
    return sums, None  # cdf computed on the fly via searchsorted


def _top_r_sum_cdf(x, k, r):
    """CDF of the sum of the r smallest order statistics from k Uniform(0,1).

    When r == k this is the Irwin-Hall distribution. When r < k, uses
    a cached empirical CDF (5M Monte Carlo samples, <0.001 accuracy).
    """
    if r >= k:
        return _irwin_hall_cdf_fast(x, k)
    if x <= 0:
        return 0.0
    if x >= r:
        return 1.0
    if r == 1:
        return 1.0 - (1.0 - min(x, 1.0)) ** k

    key = (k, r)
    if key not in _top_r_lookup_cache:
        _top_r_lookup_cache[key] = _build_top_r_lookup(k, r)

    sums, cdf = _top_r_lookup_cache[key]
    idx = np.searchsorted(sums, x, side='right')
    return float(idx / len(sums))


def _irwin_hall_cdf_fast(x, n):
    """CDF of Irwin-Hall distribution."""
    from math import factorial, comb
    if x <= 0:
        return 0.0
    if x >= n:
        return 1.0
    total = 0.0
    for k in range(int(np.floor(x)) + 1):
        total += (-1) ** k * comb(n, k) * (x - k) ** n
    return total / factorial(n)


# Public gene-level entry point (analytical Irwin-Hall)

def sumrank_fast(datasets, min_studies=MIN_STUDIES, proportion_top=1.0,
                 progress_callback=None):
    """Vectorized SumRank with Irwin-Hall analytical p-value.

    Parameters
    ----------
    datasets : list of DataFrames
    min_studies : int
        Admission criterion on raw per-gene study count.
    proportion_top : float
        Fraction of studies to use per gene (1.0 = all, 0.5 = best half).
        When < 1.0, for each gene sorts ranks and takes the k_top smallest,
        then applies the top-r null with k_top.
    progress_callback : callable, optional
    """
    if progress_callback:
        progress_callback(0, 1)

    (gene_names, rank_up, rank_down, logfc_matrix,
     n_studies, pathway_map, study_info, _) = _build_rank_matrices(datasets)

    gene_mask = n_studies >= min_studies
    ru = np.where(np.isnan(rank_up), 0.0, rank_up)
    rd = np.where(np.isnan(rank_down), 0.0, rank_down)

    lfc_filled = np.where(np.isnan(logfc_matrix), 0.0, logfc_matrix)
    lfc_count = np.sum(~np.isnan(logfc_matrix), axis=1).clip(1)
    mean_logfc = lfc_filled.sum(axis=1) / lfc_count

    # Precompute cached lookups for top-r null (triggers Monte Carlo once per (k, r))
    if proportion_top < 1.0:
        unique_k = set(n_studies[gene_mask])
        for k_val in unique_k:
            k_top = max(int(np.ceil(k_val * proportion_top)), 1)
            if k_top < k_val:
                _top_r_sum_cdf(0.5, int(k_val), k_top)  # triggers cache build

    results = []
    n_to_process = gene_mask.sum()
    processed = 0

    for i in np.where(gene_mask)[0]:
        k = n_studies[i]

        if proportion_top < 1.0:
            k_top = max(int(np.ceil(k * proportion_top)), 1)
            ru_nz = np.sort(ru[i][ru[i] > 0])
            rd_nz = np.sort(rd[i][rd[i] > 0])
            s_up = ru_nz[:k_top].sum() if len(ru_nz) >= k_top else ru[i].sum()
            s_down = rd_nz[:k_top].sum() if len(rd_nz) >= k_top else rd[i].sum()
            k_eff = k_top
        else:
            s_up = ru[i].sum()
            s_down = rd[i].sum()
            k_eff = k

        su_cap = min(s_up, k_eff / 2.0)
        sd_cap = min(s_down, k_eff / 2.0)
        if proportion_top < 1.0:
            p_up = _top_r_sum_cdf(su_cap, k, k_eff)
            p_down = _top_r_sum_cdf(sd_cap, k, k_eff)
        else:
            p_up = _irwin_hall_cdf_fast(su_cap, k_eff)
            p_down = _irwin_hall_cdf_fast(sd_cap, k_eff)

        if p_up <= p_down:
            pval = min(2 * p_up, 1.0)
            direction = 'up'
            rank_sum = s_up
            mean_signed_rank = 1.0 - (2.0 * s_up / k_eff)
        else:
            pval = min(2 * p_down, 1.0)
            direction = 'down'
            rank_sum = s_down
            mean_signed_rank = -(1.0 - (2.0 * s_down / k_eff))

        gene = gene_names[i]
        results.append({
            'names': gene, 'rank_sum': rank_sum,
            'pvals_pooled': pval, 'logfoldchanges': mean_logfc[i],
            'direction': direction,
            'mean_signed_rank': mean_signed_rank, 'n_studies': int(k),
            'pathways': pathway_map.get(gene, ''),
            'study_effects': study_info.get(gene, []),
            'pooling_method': 'sumrank',
        })
        processed += 1
        if progress_callback and (processed % 500 == 0 or processed == n_to_process):
            progress_callback(processed, n_to_process)

    result_df = pd.DataFrame(results)
    if not result_df.empty:
        result_df = result_df.sort_values('pvals_pooled').reset_index(drop=True)
        # Pseudo-Z from p-value so the volcano Z-axis works
        from scipy.stats import norm as _norm
        p_arr = result_df['pvals_pooled'].values.astype(float)
        lfc_arr = result_df['logfoldchanges'].values.astype(float)
        raw_z = _norm.isf(np.clip(p_arr, 1e-300, 1.0) / 2.0)
        result_df['z_stat'] = np.sign(lfc_arr) * np.clip(raw_z, 0, 38)
    return result_df

