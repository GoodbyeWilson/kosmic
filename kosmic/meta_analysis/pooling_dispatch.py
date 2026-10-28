# Pooling-method dispatch + consensus merge helper.
#
# Two functions:
#   'get_analytical_pool_fn(key, min_studies)' returns a callable that
#   runs one pooling method. Used by both single-method and consensus
#   gene-level MA ('gene_ma_page'), pathway MA ('pathway_ma_page'),
#   methods comparison ('methods_comparison_page'), and LOO validation.
#   Every result carries the per-study direction columns added by
#   'kosmic.meta_analysis.direction.annotate_direction' (n_up, n_down,
#   direction_conflict), whichever method produced it.
#
#   'merge_consensus_methods(method_keys, method_dfs)' merges N pooled
#   DataFrames into a consensus frame (per-method FDRs + max-FDR consensus).
#   Only used in consensus mode.
#
#   'describe_pooling_method(key, n_studies)' states in words what a key
#   runs, for the methods text; 'hksj_smallest_p(n_studies)' is the
#   smallest p-value the HKSJ t-test can give with that many studies.

from __future__ import annotations

import math
from typing import Callable, Dict, List, Optional

import numpy as np
import pandas as pd


#: Methods that apply a rule rather than a test: their table carries a
#: boolean '<rule column>' instead of a p-value to adjust, and CC
#: permutation does not calibrate them.
RULE_METHODS = {'vote_count': 'vote_pass'}

_METHOD_NAMES = {
    'dl': 'DerSimonian-Laird random effects',
    'reml': 'REML random effects',
    'sumrank': 'SumRank',
    'gwop': 'gwOP',
    'fisher': "Fisher's method",
    'stouffer': "Stouffer's Z",
    'stouffer_weighted': "weighted Stouffer's Z",
    'vote_count': 'vote counting',
}


def describe_pooling_method(key: str, n_studies: Optional[int] = None) -> str:
    """The pooling a key runs, in words, for the methods text.

    The Top 50% modifier does different things by method: SumRank sums
    each gene's best ceil(k/2) study ranks and tests that sum against the
    top-r null; DL and REML pool each gene's ceil(k/2) most precise
    studies (smallest standard error, at least 2).
    """
    hksj = key.endswith('_hksj')
    base = key.removesuffix('_hksj')
    top50 = base.endswith('_top50')
    base = base.removesuffix('_top50')
    text = _METHOD_NAMES.get(base, base)
    if top50 and base == 'sumrank':
        r = (f", r = {math.ceil(n_studies / 2)} of k = {n_studies}"
             if n_studies else "")
        text = f"SumRank over the best ceil(k/2) study ranks per gene (top-r sum{r})"
    elif top50:
        text += (" on the ceil(k/2) most precise studies per gene "
                 "(smallest standard error, at least 2)")
    if hksj:
        text += (" with the Hartung-Knapp-Sidik-Jonkman adjustment "
                 "(t-test on k - 1 degrees of freedom)")
    return text


def hksj_smallest_p(n_studies: int, t_max: float = 20.0) -> float:
    """Two-sided p-value of the HKSJ t-test at 't_max' with k - 1 df.

    With few studies even the strongest gene cannot reach the p-value a
    genome-wide FDR needs: at k = 5 a t of 20 gives p of about 4e-5,
    while Benjamini-Hochberg over 11,000 genes needs about 4e-6.
    """
    from scipy.stats import t
    if n_studies < 2:
        return 1.0
    return float(2 * t.sf(t_max, df=n_studies - 1))


def get_analytical_pool_fn(pooling_key: str,
                           min_studies: int,
                           min_votes: Optional[int] = None) -> Callable:
    """Return 'fn(datasets, min_studies=...)' for the given pooling key.

    Keys must match the UI's method list (see METHOD_FAMILIES) plus the
    Top50 / HKSJ modifier combinations exposed by the consensus sidebar.
    'min_votes' is the vote-counting threshold (default: every study).
    The returned table has the columns added by 'annotate_direction'.
    """
    from kosmic.meta_analysis.direction import annotate_direction
    from kosmic.meta_analysis.pooling.dl import dl_fast
    from kosmic.meta_analysis.pooling.reml import reml_fast
    from kosmic.meta_analysis.pooling.sumrank import sumrank_fast
    from kosmic.meta_analysis.pooling.fisher import fisher_fast
    from kosmic.meta_analysis.pooling.stouffer import stouffer_fast
    from kosmic.meta_analysis.pooling.gwop import gwop_fast
    from kosmic.meta_analysis.pooling.vote_count import vote_count_fast

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
        'vote_count':       lambda ds, min_studies=_ms: vote_count_fast(ds, min_studies=min_studies, min_votes=min_votes),
    }
    if pooling_key not in dispatch:
        raise KeyError(
            f"Unknown pooling_key {pooling_key!r}. "
            f"Known keys: {sorted(dispatch.keys())}")
    pool = dispatch[pooling_key]

    def pool_with_direction(ds, min_studies=_ms):
        return annotate_direction(pool(ds, min_studies=min_studies), ds)

    return pool_with_direction


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
        rule_col = RULE_METHODS.get(key)
        if rule_col and rule_col in mdf.columns:
            # A rule has no p-value to adjust: passing is FDR 0, failing 1.
            passes = mdf[rule_col].reindex(
                consensus.index, fill_value=False).astype(bool).values
            consensus[f'pvals_{key}'] = np.where(passes, 1e-300, 1.0)
            consensus[f'fdr_{key}'] = np.where(passes, 0.0, 1.0)
            for col in ('votes', 'vote_fraction', 'vote_conflict', 'vote_pass',
                        'min_votes', 'n_tested'):
                if col in mdf.columns and col not in consensus.columns:
                    consensus[col] = mdf[col].reindex(consensus.index).values
            continue
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
