# rOP (rth Ordered P-value) meta-analysis, Song & Tseng 2014.
#
# For each gene with k valid per-study p-values:
#
#     * sort the p-values ascending,
#     * take the r-th smallest (default r = ⌈k/2⌉),
#     * pooled p-value = Beta-CDF(p_r, r, k − r + 1) under the uniform null.
#
# The majority-rule default (r = ⌈k/2⌉) means a gene must be significant
# in at least half the studies for a small pooled p. Direction awareness
# comes from the signed mean logFC on the output, not from the statistic
# itself (rOP operates on two-sided per-study p-values).
from __future__ import annotations
from kosmic import MIN_STUDIES

import numpy as np

from .pooling_core import (
    _build_effect_matrices,
    _vectorized_dl,
    _build_result_df,
)


def rop_fast(datasets, r=None, min_studies=MIN_STUDIES, progress_callback=None):
    """Vectorised rOP meta-analysis.

    Parameters
    ----------
    datasets : list of DataFrames
        Each with columns: names, logfoldchanges, se, pvals (or pvals_adj).
    r : int or None, default None
        Order statistic index (1-based). 'None' picks the majority rule
        'r = ⌈k/2⌉' per gene. Clamped to '[1, k]'.
    min_studies : int
        Admission criterion on the raw per-gene study count.
    progress_callback : callable, optional

    Returns
    -------
    DataFrame
        Standard meta-analysis columns plus an 'r_used' column.
        Pooled logFC / SE / I² come from DL inverse-variance weighting
        for the volcano and CI display (rOP has no native effect
        estimate; matches Fisher/Stouffer's convention here). 'z_stat'
        is a signed pseudo-Z from the pooled p for the volcano axis.
    """
    from scipy.stats import beta as beta_dist, norm as _norm

    if progress_callback:
        progress_callback(0, 1)

    (gene_names, effects_matrix, var_matrix, valid_mask,
     pathway_map, study_info, _dataset_names, pval_matrix) = \
        _build_effect_matrices(datasets)

    # Valid p-value mask: finite p > 0 from a valid study
    p_valid = valid_mask & np.isfinite(pval_matrix) & (pval_matrix > 0)
    k_arr = p_valid.sum(axis=1).astype(int)

    # Sort p-values per gene ascending; invalid entries pushed to the end
    p_for_sort = np.where(p_valid, pval_matrix, np.inf)
    p_sorted = np.sort(p_for_sort, axis=1)

    # Per-gene r: default = ⌈k/2⌉; clamp to [1, k]
    if r is None:
        r_per_gene = np.ceil(k_arr / 2).astype(int)
    else:
        r_per_gene = np.minimum(np.full_like(k_arr, r), k_arr)
    r_per_gene = np.maximum(r_per_gene, 1)

    # r-th order statistic per gene
    n_ds = p_sorted.shape[1] if p_sorted.size else 1
    r_idx = np.clip(r_per_gene - 1, 0, max(n_ds - 1, 0))
    p_r = np.take_along_axis(p_sorted, r_idx[:, None], axis=1).ravel()

    # Meta p-value under Beta(r, k-r+1) null
    a = r_per_gene.astype(float)
    b = np.maximum(k_arr - r_per_gene + 1, 1).astype(float)
    with np.errstate(invalid='ignore'):
        meta_p = beta_dist.cdf(p_r, a, b)
    meta_p = np.where(k_arr > 0, meta_p, 1.0)
    meta_p = np.clip(meta_p, 1e-300, 1.0)

    # IV-weighted pooled logFC/SE/I² (display only)
    mean_lfc, se_pooled, _tau2, I2, _wre, _Q = _vectorized_dl(
        effects_matrix, var_matrix, valid_mask)

    result = _build_result_df(
        gene_names, mean_lfc, se_pooled, I2, k_arr,
        valid_mask, pathway_map, study_info, min_studies,
        pooling_method='rop')

    if not result.empty:
        pval_by_gene = dict(zip(gene_names, meta_p))
        r_by_gene = dict(zip(gene_names, r_per_gene))
        result['pvals_pooled'] = result['names'].map(pval_by_gene).astype(float)
        result['r_used'] = result['names'].map(r_by_gene).astype(int)
        # Signed pseudo-Z from p for the volcano axis.
        p_arr = result['pvals_pooled'].values.astype(float)
        lfc_arr = result['logfoldchanges'].values.astype(float)
        raw_z = _norm.isf(np.clip(p_arr, 1e-300, 1.0) / 2.0)
        result['z_stat'] = np.sign(lfc_arr) * np.clip(raw_z, 0, 38)
        result = result.sort_values('pvals_pooled').reset_index(drop=True)

    if progress_callback:
        progress_callback(1, 1)

    return result
