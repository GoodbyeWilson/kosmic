# Stouffer's Z-score combined-probability method (1949).
#
# Per-study signed Z::
#
#     zᵢ = Φ⁻¹(1 − pᵢ/2) · sign(logFCᵢ)
#
# Pooled (unweighted)::
#
#     Z = Σᵢ zᵢ / √k    ∼ N(0, 1) under H0
#
# Pooled (precision-weighted, wᵢ = 1/SEᵢ)::
#
#     Z = Σᵢ wᵢ·zᵢ / √(Σᵢ wᵢ²)
#
# Direction-aware (uses the sign of logFC). The √k boost gives it access
# to subthreshold-consistent signal — a property of any √k-boosted
# direction-aware combined-Z method rather than something specific to
# DerSimonian-Laird's variance-weighted pooling.
from __future__ import annotations
from kosmic import MIN_STUDIES

import numpy as np

from .pooling_core import (
    _build_effect_matrices,
    _vectorized_dl,
    _build_result_df,
)


def stouffer_fast(datasets, min_studies=MIN_STUDIES, weighted=False,
                  progress_callback=None):
    """Stouffer's Z meta-analysis (vectorised).

    Parameters
    ----------
    datasets : list of DataFrames
        Each with columns: names, logfoldchanges, se, pvals (or pvals_adj).
    min_studies : int
        Admission criterion on the raw per-gene study count.
    weighted : bool, default False
        If True, weight per-study z by 1/SE (precision weights).
    progress_callback : callable, optional

    Returns
    -------
    DataFrame
        Standard meta-analysis columns. Pooled logFC / SE / I² come from
        DL inverse-variance weighting for the volcano and CI display.
        'z_stat' is the actual Stouffer pooled Z (the real test
        statistic); 'pvals_pooled' is the two-sided normal tail of
        '|Z|'.
    """
    from scipy.stats import norm

    if progress_callback:
        progress_callback(0, 1)

    (gene_names, effects_matrix, var_matrix, valid_mask,
     pathway_map, study_info, _dataset_names, pval_matrix) = \
        _build_effect_matrices(datasets)

    p_valid = (valid_mask & np.isfinite(pval_matrix)
               & (pval_matrix > 0) & (pval_matrix < 1))
    k_arr = p_valid.sum(axis=1).astype(int)

    # Floor per-study p at 1/n_genes (same rationale as Fisher).
    n_genes = len(gene_names)
    p_floor = 1.0 / max(n_genes, 1)
    p_safe = np.where(p_valid,
                      np.clip(pval_matrix, p_floor, 1 - 1e-12), 0.5)
    z_mag = norm.isf(p_safe / 2.0)
    signs = np.where(effects_matrix > 0, 1.0,
                     np.where(effects_matrix < 0, -1.0, 0.0))
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
            0.0,
        )

    pvals = np.where(k_arr > 0,
                     2.0 * norm.sf(np.abs(pooled_z)), 1.0)
    pvals = np.clip(pvals, 1e-300, 1.0)

    # IV-weighted pooled logFC + SE + I² for display.
    mean_lfc, se_pooled, _tau2, I2, _wre, _Q = _vectorized_dl(
        effects_matrix, var_matrix, valid_mask)

    pooling_method = 'stouffer_weighted' if weighted else 'stouffer'
    result = _build_result_df(
        gene_names, mean_lfc, se_pooled, I2, k_arr,
        valid_mask, pathway_map, study_info, min_studies,
        pooling_method=pooling_method)
    if not result.empty:
        pval_by_gene = dict(zip(gene_names, pvals))
        result['pvals_pooled'] = result['names'].map(pval_by_gene).astype(float)
        # Override z_stat with the actual Stouffer pooled Z (the real
        # test statistic). _build_result_df computed effect/se = 0/0 = 0
        # because Stouffer doesn't have a per-gene SE.
        z_by_gene = dict(zip(gene_names, pooled_z))
        result['z_stat'] = result['names'].map(z_by_gene).astype(float)

    if progress_callback:
        progress_callback(1, 1)

    return result
