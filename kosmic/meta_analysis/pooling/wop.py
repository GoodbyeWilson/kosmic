# Weighted Ordered P-value (wOP) meta-analysis.
#
# Li & Tseng (2011). For each study i, assigns a precision-based weight wᵢ
# (median 1/SE across genes in that study), normalised so Σwᵢ = k.
# Transforms each per-study p-value via the weighted z-inflation
# 'z_w = Φ⁻¹(1 − p/2)·√w', then through a probability integral transform
# to Uniform(0,1) using the weight-specific null CDF. The r-th order
# statistic of the transformed values follows Beta(r, k−r+1) under H0.
#
# The analytical null is exact up to the PIT lookup resolution
# (10k-point grid per weight). Single public entry point 'wop_fast'.
from __future__ import annotations
from kosmic import MIN_STUDIES

from typing import Dict

import numpy as np
import pandas as pd


def wop_fast(datasets, r=None, min_studies=MIN_STUDIES, proportion_top=1.0,
             progress_callback=None):
    """Weighted Ordered P-value with analytical Beta null.

    For each study weight w, the null CDF of the weighted p-value is
    known (monotone transform of Uniform(0,1)). We transform each
    weighted p-value to Uniform via its weight-specific null CDF, then
    the r-th order statistic follows Beta(r, k−r+1).

    Parameters
    ----------
    datasets : list of DataFrames
        Each must have 'names' and 'pvals' (or 'pvals_adj') columns.
    r : int or None
        Which order statistic. 'None' → ⌈k/2⌉ (majority rule).
    min_studies : int
        Minimum datasets a gene must appear in.
    proportion_top : float
        Fraction of studies to keep per gene (1.0 = all, 0.5 = best half).
        'k_top = max(⌈k·proportion_top⌉, 1)'; r clamped to '[1, k_top]'.
    progress_callback : callable, optional

    Returns
    -------
    DataFrame with columns: names, pvals_pooled, logfoldchanges, direction,
        n_studies, r_used, study_effects, pathways, pooling_method.
    """
    from scipy import stats as scipy_stats
    from scipy.stats import beta as beta_dist

    # Per-study weights: median precision (1/SE) across genes in the study.
    dataset_weights = []
    for df in datasets:
        if 'se' in df.columns:
            valid_se = df['se'][(df['se'] > 0) & df['se'].notna()]
            if len(valid_se) > 0:
                dataset_weights.append(1.0 / valid_se.median())
            else:
                dataset_weights.append(1.0)
        else:
            dataset_weights.append(1.0)

    # Normalise so the weights sum to k.
    w_arr = np.array(dataset_weights)
    w_arr = w_arr * len(w_arr) / w_arr.sum()

    # Precompute null CDF per unique weight. Under H₀ each p ~ U(0,1),
    # so the CDF of p_w = 2·Φ(|Φ⁻¹(1−p/2)|·√w) is a monotone transform
    # of Uniform and can be sampled on a 10k-point grid once per weight.
    # Grid endpoints avoided (norm.ppf(1) = inf).
    _null_grid = np.linspace(1e-6, 1 - 1e-6, 10001)

    def _null_cdf_for_weight(w):
        """Return (p_w_sorted, cdf_vals) arrays such that
        'np.interp(observed_p_w, p_w_sorted, cdf_vals)' ≈ F_w(observed).
        """
        # Clamp w so sqrt(w) is well-defined.
        w_safe = max(float(w), 1e-3)
        z_orig = scipy_stats.norm.ppf(1 - _null_grid / 2)
        z_w = z_orig * np.sqrt(w_safe)
        p_w = np.clip(2 * scipy_stats.norm.sf(np.abs(z_w)), 1e-300, 1.0)
        sorted_idx = np.argsort(p_w)
        p_w_sorted = p_w[sorted_idx]
        cdf_vals = np.linspace(0, 1, len(p_w_sorted))
        return p_w_sorted, cdf_vals

    # Cache per unique study weight.
    weight_cdfs = {}
    for w in w_arr:
        w_key = round(float(w), 8)
        if w_key not in weight_cdfs:
            weight_cdfs[w_key] = _null_cdf_for_weight(w)

    # Collect per-gene data.
    gene_data: Dict[str, list] = {}
    pathway_map: Dict[str, str] = {}

    for ds_i, df in enumerate(datasets):
        pval_col = 'pvals_adj' if 'pvals_adj' in df.columns else 'pvals'
        if pval_col not in df.columns:
            continue
        for _, row in df.iterrows():
            gene = row['names']
            pv = row[pval_col]
            if not np.isfinite(pv) or pv <= 0:
                continue
            if gene not in gene_data:
                gene_data[gene] = []
            gene_data[gene].append({
                'pval': pv,
                'logfc': row.get('logfoldchanges', 0.0),
                'se': row.get('se', np.nan),
                'dataset': row.get('dataset', f'Study {ds_i + 1}'),
                'weight': w_arr[ds_i],
            })
            if gene not in pathway_map and 'pathways' in row.index:
                pw = row.get('pathways')
                if pd.notna(pw) and pw:
                    pathway_map[gene] = pw

    gene_items = list(gene_data.items())
    n_total = len(gene_items)

    if progress_callback is not None:
        progress_callback(0, n_total)

    # Group genes by k so the Beta CDF can be vectorised within each group.
    genes_by_k: Dict[int, list] = {}
    for gene, entries in gene_items:
        k = len(entries)
        if k < min_studies:
            continue
        genes_by_k.setdefault(k, []).append((gene, entries))

    results = []
    processed = 0

    for k, gene_group in genes_by_k.items():
        n_genes_k = len(gene_group)

        pval_arr = np.full((n_genes_k, k), np.nan)
        weight_arr = np.full((n_genes_k, k), 1.0)
        logfc_arr = np.full((n_genes_k, k), np.nan)

        for g_i, (_gene, entries) in enumerate(gene_group):
            for e_i, e in enumerate(entries):
                pval_arr[g_i, e_i] = e['pval']
                weight_arr[g_i, e_i] = e['weight']
                logfc_arr[g_i, e_i] = e['logfc']

        # Vectorised z-transform + weight inflation.
        z_matrix = scipy_stats.norm.ppf(1 - pval_arr / 2)
        z_w_matrix = z_matrix * np.sqrt(weight_arr)
        p_w_matrix = np.clip(
            2 * scipy_stats.norm.sf(np.abs(z_w_matrix)), 1e-300, 1.0)

        # PIT transform per unique weight (weights differ per gene-study
        # pair because different genes may be present in different studies).
        uniform_matrix = np.full_like(p_w_matrix, np.nan)
        unique_weights = set(round(float(w), 8)
                             for w in weight_arr.ravel()
                             if np.isfinite(w))
        for w_key in unique_weights:
            if w_key not in weight_cdfs:
                continue
            p_sorted, cdf_vals = weight_cdfs[w_key]
            mask = np.abs(weight_arr - w_key) < 1e-6
            if mask.any():
                uniform_matrix[mask] = np.clip(
                    np.interp(p_w_matrix[mask], p_sorted, cdf_vals),
                    1e-15, 1 - 1e-15)

        uniform_sorted = np.sort(uniform_matrix, axis=1)

        if proportion_top < 1.0:
            k_top = max(int(np.ceil(k * proportion_top)), 1)
        else:
            k_top = k

        r_used = r if r is not None else int(np.ceil(k_top / 2))
        r_used = min(r_used, k_top)

        # r-th order statistic across all genes in the group.
        u_r = uniform_sorted[:, r_used - 1]
        meta_pvals = beta_dist.cdf(u_r, r_used, k_top - r_used + 1)

        # Mean logFC + direction for display.
        lfc_valid = np.where(np.isfinite(logfc_arr), logfc_arr, 0.0)
        lfc_count = np.sum(np.isfinite(logfc_arr), axis=1).clip(1)
        mean_lfcs = lfc_valid.sum(axis=1) / lfc_count

        for g_i, (gene, entries) in enumerate(gene_group):
            study_effects = [
                {'dataset': e['dataset'], 'logfc': e['logfc'],
                 'se': e['se'], 'pval': e['pval'],
                 'weight': float(w_arr[e_i])}
                for e_i, e in enumerate(entries)
            ]
            results.append({
                'names': gene,
                'logfoldchanges': float(mean_lfcs[g_i]),
                'pvals_pooled': float(meta_pvals[g_i]),
                'pooling_method': 'wop',
                'direction': 'up' if mean_lfcs[g_i] > 0 else 'down',
                'n_studies': k,
                'r_used': r_used,
                'pathways': pathway_map.get(gene, ''),
                'study_effects': study_effects,
            })

        processed += n_genes_k
        if progress_callback is not None:
            progress_callback(processed, n_total)

    result_df = pd.DataFrame(results)
    if not result_df.empty:
        result_df = result_df.sort_values('pvals_pooled').reset_index(drop=True)
    return result_df
