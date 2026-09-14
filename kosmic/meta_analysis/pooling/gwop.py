# Gene-weighted Ordered P-value (gwOP) meta-analysis with directional pooling.
#
# Variant of wOP where the precision weight is **per gene-study pair** (1/SE
# of that specific gene's estimate in that specific study) rather than a
# single study-level median. Combined with directional pooling: each study's
# two-sided p is split into one-sided 'p_up' and 'p_down' by logFC sign,
# each direction is pooled via the wOP analytical null (Beta(r, k−r+1) on the
# r-th order statistic of PIT-transformed weighted p-values), and the final
# pooled p is '2·min(p_up, p_down)' clipped to '[0, 1]'. A gene that
# flips direction across studies cannot achieve a small combined p in either
# direction, which halves the wrong-direction rate in cross-study replication
# benchmarks relative to undirected wOP/rOP.
#
# Single public entry point 'gwop_fast'.
from __future__ import annotations
from kosmic import MIN_STUDIES

from typing import Dict

import numpy as np
import pandas as pd


def gwop_fast(datasets, r=None, min_studies=MIN_STUDIES, proportion_top=1.0,
              progress_callback=None, strict=False):
    """Gene-weighted Ordered P-value with directional pooling.

    Per-gene weighting (1/SE of each gene's estimate in each study)
    rather than wOP's study-median weight. Each study's two-sided p-value
    is split into one-sided p_up / p_down by logFC sign; each direction
    is pooled via the wOP analytical null; final pooled p = 2·min(p_up,
    p_down). Flips across studies cannot produce a small combined p.

    Parameters
    ----------
    datasets : list of DataFrames
        Each must have 'names', 'pvals' (or 'pvals_adj'), and 'se'.
    r : int or None
        Order statistic index. 'None' → ⌈k_top/2⌉ (majority rule).
        Ignored if 'strict=True'.
    min_studies : int
        Minimum datasets a gene must appear in.
    proportion_top : float
        Fraction of studies to keep per gene (1.0 = all, 0.5 = best half).
    strict : bool, default False
        If True, require unanimous directional agreement (r = k_top). Most
        conservative rule; the largest order statistic being small requires
        every study to show evidence in the same direction.
    progress_callback : callable, optional

    Returns
    -------
    DataFrame with columns: names, pvals_pooled, logfoldchanges, direction,
        n_studies, r_used, study_effects, pathways, pooling_method,
        mean_signed_rank, z_stat.
    """
    from scipy import stats as scipy_stats
    from scipy.stats import beta as beta_dist
    from scipy.stats import norm as _norm

    # Null-CDF cache shared across all genes. Each gene's weights may
    # differ (per-gene 1/SE), so a grid binned at 1 decimal place keeps
    # the cache small while remaining numerically accurate (wOP PIT is
    # smooth in w).
    _null_grid = np.linspace(1e-6, 1 - 1e-6, 10001)
    _weight_cdf_cache = {}

    def _null_cdf_for_weight(w):
        w_key = round(float(w), 1)
        if w_key not in _weight_cdf_cache:
            w_safe = max(float(w), 1e-3)
            z_orig = scipy_stats.norm.ppf(1 - _null_grid / 2)
            z_w = z_orig * np.sqrt(w_safe)
            p_w = np.clip(2 * scipy_stats.norm.sf(np.abs(z_w)), 1e-300, 1.0)
            sorted_idx = np.argsort(p_w)
            _weight_cdf_cache[w_key] = (
                p_w[sorted_idx], np.linspace(0, 1, len(p_w)))
        return _weight_cdf_cache[w_key]

    # Collect per-gene data with per-gene SE as weight.
    gene_data: Dict[str, list] = {}
    pathway_map: Dict[str, str] = {}

    for ds_i, df in enumerate(datasets):
        pval_col = 'pvals_adj' if 'pvals_adj' in df.columns else 'pvals'
        if pval_col not in df.columns:
            continue
        for _, row in df.iterrows():
            gene = row['names']
            pv = row[pval_col]
            se = row.get('se', np.nan)
            if not np.isfinite(pv) or pv <= 0:
                continue
            # Floor per-study p at 1/n_genes. No single study has more
            # precision than its own multiple-testing burden; without this,
            # combined-p methods produce absurdly small pooled values (1e-30)
            # from multiplying tail artifacts.
            pv = max(pv, 1.0 / max(len(df), 1))
            # Per-gene weight: 1/SE (or 1.0 if SE missing).
            w = 1.0 / se if np.isfinite(se) and se > 0 else 1.0
            if gene not in gene_data:
                gene_data[gene] = []
            gene_data[gene].append({
                'pval': pv,
                'logfc': row.get('logfoldchanges', 0.0),
                'se': se,
                'dataset': row.get('dataset', f'Study {ds_i + 1}'),
                'weight': w,
            })
            if gene not in pathway_map and 'pathways' in row.index:
                pw = row.get('pathways')
                if pd.notna(pw) and pw:
                    pathway_map[gene] = pw

    gene_items = list(gene_data.items())
    n_total = len(gene_items)

    if progress_callback is not None:
        progress_callback(0, n_total)

    # Group genes by k for vectorised Beta CDF.
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
            raw_weights = np.array([e['weight'] for e in entries])
            w_sum = raw_weights.sum()
            if w_sum > 0:
                raw_weights = raw_weights * k / w_sum
            for e_i, e in enumerate(entries):
                pval_arr[g_i, e_i] = e['pval']
                weight_arr[g_i, e_i] = raw_weights[e_i]
                logfc_arr[g_i, e_i] = e['logfc']

        # 'weight_arr' and 'k' are bound via default args so each
        # iteration's closure owns its own values (silences B023; the
        # closure is only called within this iteration anyway).
        def _pool_directional(p_directional_arr,
                              _w=weight_arr, _k=k):
            z_m = scipy_stats.norm.ppf(1 - p_directional_arr / 2)
            z_w = z_m * np.sqrt(_w)
            p_w = np.clip(
                2 * scipy_stats.norm.sf(np.abs(z_w)), 1e-300, 1.0)
            u_mat = np.full_like(p_w, np.nan)
            w_b = np.round(_w, 1)
            uniq = set(w_b[np.isfinite(w_b)].ravel())
            for wk in uniq:
                p_sorted, cdf_vals = _null_cdf_for_weight(wk)
                m = np.abs(w_b - wk) < 0.05
                if m.any():
                    u_mat[m] = np.clip(
                        np.interp(p_w[m], p_sorted, cdf_vals),
                        1e-15, 1 - 1e-15)
            u_sorted = np.sort(u_mat, axis=1)
            if proportion_top < 1.0:
                kt = max(int(np.ceil(_k * proportion_top)), 1)
            else:
                kt = _k
            u_top = u_sorted[:, :kt]
            if strict:
                # Unanimous directional agreement: r = k_top. Largest
                # order statistic must be small, requiring every study
                # to show evidence in the same direction.
                ru = kt
            else:
                ru = r if r is not None else int(np.ceil(kt / 2))
            ru = min(ru, kt)
            return beta_dist.cdf(u_top[:, ru - 1], ru, kt - ru + 1), kt, ru

        # Split two-sided p into one-sided p_up / p_down by logFC sign.
        pos = logfc_arr > 0
        p_up_arr = np.where(pos, pval_arr / 2, 1.0 - pval_arr / 2)
        p_down_arr = np.where(~pos, pval_arr / 2, 1.0 - pval_arr / 2)

        p_up, _k_top, r_used = _pool_directional(p_up_arr)
        p_down, _, _ = _pool_directional(p_down_arr)
        meta_pvals = np.clip(2.0 * np.minimum(p_up, p_down), 0.0, 1.0)

        # IV-weighted pooled logFC for display (fixed-effects-style;
        # good enough for display since DL and REML agree on the point
        # estimate). gwOP pools p-values, not effects.
        se_arr = np.array([[e['se'] for e in gene_group[g][1]]
                           for g in range(n_genes_k)], dtype=float)
        se_valid = np.isfinite(se_arr) & (se_arr > 0)
        var_arr = np.where(se_valid, se_arr ** 2, np.inf)
        lfc_mask = np.isfinite(logfc_arr) & se_valid
        iv_w = np.where(lfc_mask, 1.0 / var_arr, 0.0)
        iv_sum = iv_w.sum(axis=1)
        mean_lfcs = np.where(
            iv_sum > 0,
            (iv_w * np.where(lfc_mask, logfc_arr, 0.0)).sum(axis=1) / iv_sum,
            np.where(np.isfinite(logfc_arr), logfc_arr, 0.0).sum(axis=1)
            / np.sum(np.isfinite(logfc_arr), axis=1).clip(1)
        )

        for g_i, (gene, entries) in enumerate(gene_group):
            raw_w = np.array([e['weight'] for e in entries])
            w_sum = raw_w.sum()
            if w_sum > 0:
                norm_w = raw_w * k / w_sum
            else:
                norm_w = np.ones(k)
            study_effects = [
                {'dataset': e['dataset'], 'logfc': e['logfc'],
                 'se': e['se'], 'pval': e['pval'],
                 'weight': float(norm_w[e_i])}
                for e_i, e in enumerate(entries)
            ]
            # Signed consistency score: sign(logFC) · (1 − p) per study,
            # averaged. Range [−1, +1]. ±1 = all studies strongly agree;
            # 0 = no evidence or conflicting direction.
            signed_evidence = [
                np.sign(e['logfc']) * (1.0 - e['pval'])
                for e in entries
                if np.isfinite(e['logfc']) and np.isfinite(e['pval'])
            ]
            msr = (float(np.mean(signed_evidence))
                   if signed_evidence else 0.0)

            results.append({
                'names': gene,
                'logfoldchanges': float(mean_lfcs[g_i]),
                'mean_signed_rank': msr,
                'pvals_pooled': float(meta_pvals[g_i]),
                'pooling_method': 'gwop',
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
        # Pseudo-Z from p-value so the volcano Z-axis works.
        p_arr = result_df['pvals_pooled'].values.astype(float)
        lfc_arr = result_df['logfoldchanges'].values.astype(float)
        raw_z = _norm.isf(np.clip(p_arr, 1e-300, 1.0) / 2.0)
        result_df['z_stat'] = np.sign(lfc_arr) * np.clip(raw_z, 0, 38)
    return result_df
