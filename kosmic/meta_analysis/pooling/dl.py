# DerSimonian-Laird random-effects meta-analysis — unified entry point.
#
# Single 'dl_fast' covering plain DL, DL on top-50% most-precise studies
# per gene, DL with HKSJ correction, and any combination. Vectorised over
# genes; per-row callers (single gene or pathway, K studies) build a 1-row
# DataFrame per study and call 'dl_fast' directly — see
# 'kosmic.visualisation.forest_plot' and
# 'gui.meta_analysis.pages.pathway_ma_page' for examples.
# Cross-checked against statsmodels in 'tests/meta_analysis/test_dl_correctness.py'.
from __future__ import annotations
from kosmic import MIN_STUDIES

import numpy as np

from .pooling_core import (
    _build_effect_matrices,
    _vectorized_dl,
    _hksj_se,
    _build_result_df,
)


def dl_fast(datasets, min_studies=MIN_STUDIES, proportion_top=1.0, hksj=False,
            method='random', progress_callback=None):
    """DerSimonian-Laird random-effects meta-analysis.

    Parameters
    ----------
    datasets : list of DataFrames
        Each with columns: names, logfoldchanges, se, pathways (optional).
    min_studies : int
        Minimum datasets a gene must appear in post-filter.
    proportion_top : float, default 1.0
        Fraction of studies to keep per gene, selected by lowest SE
        (most precise first). 1.0 = keep all; 0.5 = top-50%. Clamped
        to a floor of 2 studies per gene.
    hksj : bool, default False
        If True, apply Hartung-Knapp-Sidik-Jonkman correction: replace
        the pooled SE with the HKSJ-corrected form, use t-distribution
        with k-1 df for p-values. More conservative under few studies.
    method : {'random', 'fixed'}, default 'random'
        'random' for DL random-effects (the usual case); 'fixed' for
        inverse-variance-weighted fixed-effects (τ²=0 by construction).
    progress_callback : callable, optional
        Called with (current, total).

    Returns
    -------
    DataFrame
        Core columns: names, logfoldchanges, se, z_stat, ci_lower, ci_upper,
        pvals_pooled, n_studies, pathways, study_effects, heterogeneity_i2,
        pooling_method.
        Random-effects reporting columns: tau_squared, Q, p_heterogeneity,
        pi_lower, pi_upper (95% prediction interval). Fixed-effects fixes
        tau_squared to 0 and reports the fixed-effects heterogeneity Q.
    """
    from scipy.stats import chi2

    if progress_callback:
        progress_callback(0, 1)

    (gene_names, effects_matrix, var_matrix, valid_mask,
     pathway_map, study_info, _dataset_names, _) = _build_effect_matrices(datasets)

    k_arr = valid_mask.sum(axis=1).astype(int)

    # Optional: reduce to the top-k most-precise studies per gene
    if proportion_top < 1.0:
        se_matrix = np.where(valid_mask, np.sqrt(var_matrix), np.inf)
        # rank[i, j] = 0 means study j has the smallest SE for gene i
        rank = se_matrix.argsort(axis=1).argsort(axis=1)
        k_top = np.maximum(np.ceil(k_arr * proportion_top).astype(int), 2)
        k_top = np.minimum(k_top, k_arr)  # never exceed available k
        analysis_mask = valid_mask & (rank < k_top[:, None])
        k_analysis = analysis_mask.sum(axis=1).astype(int)
    else:
        analysis_mask = valid_mask
        k_analysis = k_arr

    # Pool using DL (random) or plain inverse-variance (fixed)
    if method == 'random':
        pooled, se_pooled, tau2, I2, weights_re, Q = _vectorized_dl(
            effects_matrix, var_matrix, analysis_mask)
    elif method == 'fixed':
        eff = np.where(analysis_mask, effects_matrix, 0.0)
        var = np.where(analysis_mask, var_matrix, np.inf)
        weights_fe = np.where(analysis_mask, 1.0 / var, 0.0)
        sum_w = weights_fe.sum(axis=1)
        pooled = np.where(sum_w > 0,
                          (weights_fe * eff).sum(axis=1) / sum_w, 0.0)
        se_pooled = np.where(sum_w > 0, np.sqrt(1.0 / sum_w), 0.0)
        resid = np.where(analysis_mask, eff - pooled[:, None], 0.0)
        Q = (weights_fe * resid ** 2).sum(axis=1)
        df_arr = k_analysis.astype(float) - 1
        with np.errstate(invalid='ignore', divide='ignore'):
            I2 = np.where(Q > 0, np.maximum(0, (Q - df_arr) / Q), 0.0)
        tau2 = np.zeros_like(pooled)  # fixed-effects: no between-study variance
        weights_re = weights_fe  # identical under fixed-effects for HKSJ path
    else:
        raise ValueError(f"method must be 'random' or 'fixed', got {method!r}")

    # Optional HKSJ correction on the pooled SE + t-distribution p-value
    if hksj:
        eff_for_hksj = np.where(analysis_mask, effects_matrix, 0.0)
        se_pooled = _hksj_se(pooled, eff_for_hksj, weights_re, analysis_mask)

    # Label for the pooling_method column (informational; not parsed)
    parts = ['dl' if method == 'random' else 'dl_fixed']
    if proportion_top < 1.0:
        parts.append(f'top{int(round(proportion_top * 100))}')
    if hksj:
        parts.append('hksj')
    pooling_method = '_'.join(parts)

    result = _build_result_df(
        gene_names, pooled, se_pooled, I2, k_analysis,
        analysis_mask, pathway_map, study_info, min_studies,
        hksj=hksj, pooling_method=pooling_method,
        k_filter=k_arr)

    # Append random-effects reporting columns.
    # Order of operations: the admission filter in _build_result_df drops
    # sub-min_studies genes, so we must align Q/tau2 to the surviving rows.
    if not result.empty:
        gene_to_idx = {g: i for i, g in enumerate(gene_names)}
        row_idx = np.array([gene_to_idx[g] for g in result['names']])
        tau2_row = tau2[row_idx]
        Q_row = Q[row_idx]
        k_row = k_analysis[row_idx].astype(float)
        df_row = np.maximum(k_row - 1, 0)
        se_row = result['se'].values.astype(float)

        with np.errstate(invalid='ignore', divide='ignore'):
            p_het = np.where(df_row > 0,
                             1.0 - chi2.cdf(Q_row, df_row),
                             1.0)
        pi_se = np.sqrt(se_row ** 2 + tau2_row)
        # HKSJ uses t-distribution for the diamond CI; the prediction
        # interval still uses the Z multiplier by convention.
        pi_lower = result['logfoldchanges'].values - 1.96 * pi_se
        pi_upper = result['logfoldchanges'].values + 1.96 * pi_se

        result['tau_squared'] = tau2_row
        result['Q'] = Q_row
        result['p_heterogeneity'] = p_het
        result['pi_lower'] = pi_lower
        result['pi_upper'] = pi_upper

    if progress_callback:
        progress_callback(1, 1)

    return result
