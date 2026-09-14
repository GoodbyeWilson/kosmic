# Pooling-method dispatch + consensus merge helper.
#
# Two functions:
#   'get_analytical_pool_fn(key, min_studies)' returns a callable that
#   runs one pooling method. Used by both single-method and consensus
#   gene-level MA ('gene_ma_page'), pathway MA ('pathway_ma_page'),
#   methods comparison ('methods_comparison_page'), and LOO validation.
#
#   'merge_consensus_methods(method_keys, method_dfs)' merges N pooled
#   DataFrames into a consensus frame (per-method FDRs + max-FDR consensus).
#   Only used in consensus mode.

from __future__ import annotations

from typing import Callable, Dict, List

import numpy as np
import pandas as pd


def get_analytical_pool_fn(pooling_key: str,
                           min_studies: int) -> Callable:
    """Return 'fn(datasets, min_studies=...)' for the given pooling key.

    Keys must match the UI's method list (see METHOD_FAMILIES) plus the
    Top50 / HKSJ modifier combinations exposed by the consensus sidebar.
    """
    from kosmic.meta_analysis.pooling.dl import dl_fast
    from kosmic.meta_analysis.pooling.reml import reml_fast
    from kosmic.meta_analysis.pooling.sumrank import sumrank_fast
    from kosmic.meta_analysis.pooling.fisher import fisher_fast
    from kosmic.meta_analysis.pooling.stouffer import stouffer_fast
    from kosmic.meta_analysis.pooling.gwop import gwop_fast

    _ms = min_studies
    # Each lambda accepts an optional min_studies override; defaults to
    # the value captured when the dispatcher was constructed.
    dispatch = {
        'dl':               lambda ds, min_studies=_ms: dl_fast(ds, min_studies=min_studies),
        'dl_top50':         lambda ds, min_studies=_ms: dl_fast(ds, min_studies=min_studies, proportion_top=0.5),
        'dl_hksj':          lambda ds, min_studies=_ms: dl_fast(ds, min_studies=min_studies, hksj=True),
        'dl_top50_hksj':    lambda ds, min_studies=_ms: dl_fast(ds, min_studies=min_studies, proportion_top=0.5, hksj=True),
        'reml':             lambda ds, min_studies=_ms: reml_fast(ds, min_studies=min_studies),
        'reml_top50':       lambda ds, min_studies=_ms: reml_fast(ds, min_studies=min_studies, proportion_top=0.5),
        'reml_hksj':        lambda ds, min_studies=_ms: reml_fast(ds, min_studies=min_studies, hksj=True),
        'reml_top50_hksj':  lambda ds, min_studies=_ms: reml_fast(ds, min_studies=min_studies, proportion_top=0.5, hksj=True),
        'sumrank':          lambda ds, min_studies=_ms: sumrank_fast(ds, min_studies=min_studies),
        'sumrank_top50':    lambda ds, min_studies=_ms: sumrank_fast(ds, min_studies=min_studies, proportion_top=0.5),
        'gwop':             lambda ds, min_studies=_ms: gwop_fast(ds, min_studies=min_studies),
        'fisher':           lambda ds, min_studies=_ms: fisher_fast(ds, min_studies=min_studies),
        'stouffer':         lambda ds, min_studies=_ms: stouffer_fast(ds, min_studies=min_studies),
        'stouffer_weighted': lambda ds, min_studies=_ms: stouffer_fast(ds, min_studies=min_studies, weighted=True),
    }
    if pooling_key not in dispatch:
        raise KeyError(
            f"Unknown pooling_key {pooling_key!r}. "
            f"Known keys: {sorted(dispatch.keys())}")
    return dispatch[pooling_key]


def merge_consensus_methods(method_keys: List[str],
                            method_dfs: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Merge N pooling-method DataFrames into a consensus DataFrame.

    Per-method raw p-values are stored as 'pvals_<key>' and per-method
    BH-FDRs as 'fdr_<key>'. The consensus FDR is the max of the
    per-method FDRs, so a gene passes consensus iff every individual
    method passes its FDR threshold.

    The base method's columns (logfoldchanges, ci_lower, ci_upper,
    heterogeneity_i2, n_studies, etc.) are preserved on the merged
    output. The base method is the first 'dl'-prefixed key if any
    (DerSimonian-Laird and its HKSJ variants populate the standard
    effect-size columns); otherwise the first key in 'method_keys'.

    Parameters
    ----------
    method_keys : list of str
        Method keys in the order they were pooled.
    method_dfs : dict
        Mapping '{key: pooled_DataFrame}'. Each pooled DataFrame must
        have a 'names' column and a 'pvals_pooled' column.

    Returns
    -------
    DataFrame
        With columns: names, logfoldchanges, ..., pvals_<key> (per
        method), fdr_<key> (per method), pvals_pooled (worst case),
        fdr (consensus = max per-method FDR), pooling_method='consensus'.
    """
    from kosmic.numerical import bh_fdr

    if not method_keys:
        raise ValueError("method_keys must not be empty")
    if not method_dfs:
        raise ValueError("method_dfs must not be empty")

    # Per-method input validation -- a None result or a missing
    # 'names' column is a dispatch bug upstream (e.g. a pooling
    # function that returned 'pd.DataFrame()' instead of a
    # schema-preserving empty frame). Surface it with a descriptive
    # error instead of letting it blow up later on a KeyError.
    for key in method_keys:
        mdf = method_dfs.get(key)
        if mdf is None:
            raise ValueError(
                f"Pooling method '{key}' returned None. "
                f"This is a dispatch bug in get_analytical_pool_fn.")
        if 'names' not in mdf.columns:
            raise ValueError(
                f"Pooling method '{key}' returned a DataFrame without "
                f"a 'names' column (cols: {list(mdf.columns)}, "
                f"len={len(mdf)}). Check the pooling function for "
                f"this method.")

    base_key = method_keys[0]
    for k in method_keys:
        if k.startswith('dl'):
            base_key = k
            break

    base_df = method_dfs[base_key].set_index('names')
    consensus = base_df.copy()

    for key in method_keys:
        mdf = method_dfs[key].set_index('names')
        pvals = mdf['pvals_pooled'].reindex(
            consensus.index, fill_value=1.0).values.astype(float)
        # Clamp exact zeros to 1e-300 so FDR ordering survives underflow.
        pvals = np.where(pvals == 0.0, 1e-300, pvals)
        # Coerce infs/NaNs to 1.0 so they do not corrupt the BH step.
        pvals = np.where(np.isfinite(pvals) & (pvals > 0), pvals, 1.0)
        consensus[f'pvals_{key}'] = pvals

        consensus[f'fdr_{key}'] = bh_fdr(pvals)

    pval_cols = [f'pvals_{k}' for k in method_keys]
    consensus['pvals_pooled'] = consensus[pval_cols].max(axis=1)

    fdr_cols = [f'fdr_{k}' for k in method_keys]
    consensus['fdr'] = consensus[fdr_cols].max(axis=1)

    consensus['pooling_method'] = 'consensus'
    return consensus.reset_index()
