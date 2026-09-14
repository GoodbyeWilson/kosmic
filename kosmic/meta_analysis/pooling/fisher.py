# Fisher's combined probability method (1925).
#
# Combines per-study p-values via
#
#     X = −2 · Σᵢ log(pᵢ)  ~ χ²(2·k) under H0
#
# Direction-blind: a gene up in 4 studies and down in 4 can still produce
# a small pooled p if per-study p-values are small. Kept as the oldest
# and most cited combined-p procedure; included in the consensus dispatch
# alongside the direction-aware alternatives.
from __future__ import annotations
from kosmic import MIN_STUDIES

import numpy as np

from .pooling_core import (
    _build_effect_matrices,
    _vectorized_dl,
    _build_result_df,
)


def fisher_fast(datasets, min_studies=MIN_STUDIES, progress_callback=None):
    """Fisher's combined-probability meta-analysis (vectorised).

    Parameters
    ----------
    datasets : list of DataFrames
        Each with columns: names, logfoldchanges, se, pvals (or pvals_adj).
    min_studies : int
        Admission criterion on the raw per-gene study count.
    progress_callback : callable, optional

    Returns
    -------
    DataFrame
        Standard meta-analysis columns. Pooled logFC / SE / I² come from
        DL inverse-variance weighting for the volcano and CI display
        (Fisher has no native effect estimate). 'pvals_pooled' is the
        χ² combined-p; 'z_stat' is a signed pseudo-Z from the p-value
        for the volcano axis.
    """
    from scipy.stats import chi2, norm as _norm

    if progress_callback:
        progress_callback(0, 1)

    (gene_names, effects_matrix, var_matrix, valid_mask,
     pathway_map, study_info, _dataset_names, pval_matrix) = \
        _build_effect_matrices(datasets)

    p_valid = valid_mask & np.isfinite(pval_matrix) & (pval_matrix > 0)
    k_arr = p_valid.sum(axis=1).astype(int)

    # Floor per-study p at 1/n_genes: no single study has more precision
    # than its own multiple-testing burden. Without this, Fisher's
    # combined-p produces absurdly small pooled values.
    n_genes = len(gene_names)
    p_floor = 1.0 / max(n_genes, 1)
    p_safe = np.where(p_valid, np.clip(pval_matrix, p_floor, 1.0), 1.0)
    log_p = np.where(p_valid, np.log(p_safe), 0.0)
    fisher_stat = -2.0 * log_p.sum(axis=1)
    df = 2 * k_arr

    pvals = np.where(k_arr > 0, chi2.sf(fisher_stat, df), 1.0)
    pvals = np.clip(pvals, 1e-300, 1.0)

    # IV-weighted pooled logFC + SE + I² for display.
    mean_lfc, se_pooled, _tau2, I2, _wre, _Q = _vectorized_dl(
        effects_matrix, var_matrix, valid_mask)

    result = _build_result_df(
        gene_names, mean_lfc, se_pooled, I2, k_arr,
        valid_mask, pathway_map, study_info, min_studies,
        pooling_method='fisher')
    if not result.empty:
        pval_by_gene = dict(zip(gene_names, pvals))
        result['pvals_pooled'] = result['names'].map(pval_by_gene).astype(float)
        # Signed pseudo-Z from p for the volcano Z-axis.
        p_arr = result['pvals_pooled'].values.astype(float)
        lfc_arr = result['logfoldchanges'].values.astype(float)
        raw_z = _norm.isf(np.clip(p_arr, 1e-300, 1.0) / 2.0)
        result['z_stat'] = np.sign(lfc_arr) * np.clip(raw_z, 0, 38)

    if progress_callback:
        progress_callback(1, 1)

    return result
