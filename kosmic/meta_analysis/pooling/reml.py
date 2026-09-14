# REML random-effects meta-analysis — unified entry point.
#
# Single 'reml_fast' covering plain REML and REML with HKSJ correction,
# vectorised over genes. Per-row callers (single pathway / gene, K studies)
# build a 1-row DataFrame per study and call 'reml_fast' directly — see
# 'gui.meta_analysis.pages.pathway_ma_page._pool_advanced'.
#
# Matches R's 'metafor::rma(method="REML")' for the τ² estimate to ~1e-7
# (Fisher-scoring iteration reaches the same maximum of the REML log-likelihood
# that metafor does; cross-checked against 'scipy.optimize.minimize_scalar'
# on the REML log-likelihood and against the published 'dat.bcg' values in
# 'tests/meta_analysis/test_reml_correctness.py'). With 'hksj=True' matches
# 'metafor::rma(method="REML", test="knha")'.
from __future__ import annotations
from kosmic import MIN_STUDIES

import numpy as np

from .pooling_core import (
    _build_effect_matrices,
    _vectorized_dl,
    _hksj_se,
    _build_result_df,
)


# Vectorised REML τ² (gene × study matrices)

def _vectorized_reml(effects_matrix, var_matrix, valid_mask,
                     max_iter=100, tol=1e-6):
    """Fisher-scoring REML τ² across all genes at once.

    Starts from the DL estimate, iterates until convergence. Returns
    '(pooled, se_pooled, tau2, I2, weights_re)'.
    """
    _, _, tau2, _, _, _ = _vectorized_dl(effects_matrix, var_matrix, valid_mask)

    eff = np.where(valid_mask, effects_matrix, 0.0)
    var = np.where(valid_mask, var_matrix, np.inf)

    for _ in range(max_iter):
        total_var = np.where(valid_mask, var + tau2[:, None], np.inf)
        w = np.where(valid_mask, 1.0 / total_var, 0.0)
        W = w.sum(axis=1)

        pooled_iter = np.where(W > 0, (w * eff).sum(axis=1) / W, 0.0)
        resid = np.where(valid_mask, eff - pooled_iter[:, None], 0.0)

        w2 = w ** 2
        w3 = w ** 3
        sum_w2 = w2.sum(axis=1)
        sum_w3 = w3.sum(axis=1)

        num = (w2 * resid ** 2).sum(axis=1) - (W - np.where(W > 0, sum_w2 / W, 0.0))
        denom = sum_w2 - np.where(W > 0, sum_w3 / W, 0.0)

        safe_denom = np.where(np.abs(denom) > 1e-15, denom, 1.0)
        tau2_new = np.maximum(0, tau2 + num / safe_denom)
        tau2_new = np.where(np.abs(denom) > 1e-15, tau2_new, tau2)

        if np.all(np.abs(tau2_new - tau2) < tol):
            tau2 = tau2_new
            break
        tau2 = tau2_new

    total_var_final = np.where(valid_mask, var + tau2[:, None], np.inf)
    weights_re = np.where(valid_mask, 1.0 / total_var_final, 0.0)
    sum_w_re = weights_re.sum(axis=1)

    pooled = np.where(sum_w_re > 0,
                      (weights_re * eff).sum(axis=1) / sum_w_re, 0.0)
    se_pooled = np.where(sum_w_re > 0,
                         np.sqrt(1.0 / sum_w_re), 0.0)

    weights_fe = np.where(valid_mask, 1.0 / var, 0.0)
    sum_w_fe = weights_fe.sum(axis=1)
    pooled_fe = np.where(sum_w_fe > 0,
                         (weights_fe * eff).sum(axis=1) / sum_w_fe, 0.0)
    resid_fe = np.where(valid_mask, eff - pooled_fe[:, None], 0.0)
    Q = (weights_fe * resid_fe ** 2).sum(axis=1)
    k = valid_mask.sum(axis=1).astype(float)
    df = k - 1
    I2 = np.where(Q > 0, np.maximum(0, (Q - df) / Q), 0.0)

    return pooled, se_pooled, tau2, I2, weights_re


# Public gene-level entry point

def reml_fast(datasets, min_studies=MIN_STUDIES, proportion_top=1.0, hksj=False,
              progress_callback=None):
    """REML random-effects meta-analysis.

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
        the pooled SE with the HKSJ-corrected form and use t-distribution
        with k−1 df for p-values.
    progress_callback : callable, optional
        Called with (current, total).

    Returns
    -------
    DataFrame
        Core columns: names, logfoldchanges, se, z_stat, ci_lower, ci_upper,
        pvals_pooled, n_studies, pathways, study_effects, heterogeneity_i2,
        pooling_method.
        Random-effects reporting columns: tau_squared, Q, p_heterogeneity,
        pi_lower, pi_upper (95% prediction interval).
    """
    from scipy.stats import chi2 as _chi2

    if progress_callback:
        progress_callback(0, 1)

    (gene_names, effects_matrix, var_matrix, valid_mask,
     pathway_map, study_info, _dataset_names, _) = _build_effect_matrices(datasets)

    k_arr = valid_mask.sum(axis=1).astype(int)

    if proportion_top < 1.0:
        se_matrix = np.where(valid_mask, np.sqrt(var_matrix), np.inf)
        rank = se_matrix.argsort(axis=1).argsort(axis=1)
        k_top = np.maximum(np.ceil(k_arr * proportion_top).astype(int), 2)
        k_top = np.minimum(k_top, k_arr)
        analysis_mask = valid_mask & (rank < k_top[:, None])
        k_analysis = analysis_mask.sum(axis=1).astype(int)
    else:
        analysis_mask = valid_mask
        k_analysis = k_arr

    pooled, se_pooled, tau2, I2, weights_re = _vectorized_reml(
        effects_matrix, var_matrix, analysis_mask)

    # Q statistic (fixed-effects Cochran's Q on the analysis subset)
    eff = np.where(analysis_mask, effects_matrix, 0.0)
    var = np.where(analysis_mask, var_matrix, np.inf)
    w_fe = np.where(analysis_mask, 1.0 / var, 0.0)
    sum_w_fe = w_fe.sum(axis=1)
    pooled_fe = np.where(sum_w_fe > 0,
                         (w_fe * eff).sum(axis=1) / sum_w_fe, 0.0)
    resid_fe = np.where(analysis_mask, eff - pooled_fe[:, None], 0.0)
    Q_arr = (w_fe * resid_fe ** 2).sum(axis=1)

    if hksj:
        se_pooled = _hksj_se(pooled, eff, weights_re, analysis_mask)

    parts = ['reml']
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

    # Append RE reporting columns (same shape as dl_fast).
    if not result.empty:
        gene_to_idx = {g: i for i, g in enumerate(gene_names)}
        row_idx = np.array([gene_to_idx[g] for g in result['names']])
        tau2_row = tau2[row_idx]
        Q_row = Q_arr[row_idx]
        k_row = k_analysis[row_idx].astype(float)
        df_row = np.maximum(k_row - 1, 0)
        se_row = result['se'].values.astype(float)

        with np.errstate(invalid='ignore', divide='ignore'):
            p_het = np.where(df_row > 0,
                             1.0 - _chi2.cdf(Q_row, df_row),
                             1.0)
        pi_se = np.sqrt(se_row ** 2 + tau2_row)
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


