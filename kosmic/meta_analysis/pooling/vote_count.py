# Vote counting: how many studies call a gene significant in one direction.
#
# For each gene, a study votes 'up' when its adjusted p-value is below the
# study-level FDR threshold and its log2 fold change is positive, 'down'
# when the p-value is below the threshold and the fold change negative.
# The gene's votes are the larger of the two counts, out of the studies in
# which it could be tested (present, with a finite effect and standard
# error, and an adjusted p-value). Studies without a 'pvals_adj' column
# are left out of the votes and of the denominator, as in
# 'kosmic.meta_analysis.direction'.
#
# Vote counting is a rule, not a test: a gene passes when it has at least
# 'min_votes' votes and no study calls it significant in the opposite
# direction. There is no p-value; the pooled table carries 'vote_pass',
# and 'merge_consensus_methods' turns that into an FDR of 0 (passes) or 1
# (does not), so the consensus, the significant-gene counts and the Venn
# diagram work as for the other methods. CC permutation calibrates the
# other methods only.
#
# 'vote_levels' groups genes by the number of studies agreeing (0/k ...
# k/k, and 'conflict'), for the Results filter and the Study Agreement
# figure.
#
# Vote counting ignores effect sizes and standard errors and gives no
# pooled estimate or confidence interval; it is a robustness check beside
# DL / REML, not a replacement. The pooled log2 fold change column holds
# the inverse-variance mean, as for Fisher's method, for the plots.
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from kosmic import DEFAULT_FDR, MIN_STUDIES

from .pooling_core import _build_effect_matrices, _build_result_df, _vectorized_dl

CONFLICT = 'conflict'


def count_votes(effects: np.ndarray, padj: np.ndarray, testable: np.ndarray,
                fdr: float = DEFAULT_FDR) -> dict:
    """Per-gene votes from (n_genes, n_studies) matrices.

    Parameters
    ----------
    effects : ndarray
        Per-study log2 fold changes.
    padj : ndarray
        Per-study adjusted p-values.
    testable : ndarray of bool
        True where the gene could be tested in the study.
    fdr : float
        Study-level FDR threshold for a vote.

    Returns
    -------
    dict of ndarray
        'n_tested', 'n_up', 'n_down', 'votes' (the larger count) and
        'conflict' (significant in both directions).
    """
    sig = testable & np.isfinite(padj) & (padj < fdr)
    n_up = (sig & (effects > 0)).sum(axis=1)
    n_down = (sig & (effects < 0)).sum(axis=1)
    return {'n_tested': testable.sum(axis=1), 'n_up': n_up, 'n_down': n_down,
            'votes': np.maximum(n_up, n_down),
            'conflict': (n_up > 0) & (n_down > 0)}


def vote_count_fast(datasets, min_studies=MIN_STUDIES, min_votes: Optional[int] = None,
                    fdr=DEFAULT_FDR, progress_callback=None):
    """Vote-counting meta-analysis over per-study DE tables.

    Parameters
    ----------
    datasets : list of DataFrames
        Each with columns: names, logfoldchanges, se, pvals_adj.
    min_studies : int
        Admission criterion on the number of studies a gene was tested in.
    min_votes : int, optional
        Votes a gene needs to pass; default every study ('len(datasets)').
    fdr : float
        Study-level FDR threshold for a vote.
    progress_callback : callable, optional

    Returns
    -------
    DataFrame
        Standard meta-analysis columns plus 'n_tested', 'votes',
        'vote_level' ('v/k' or 'conflict', as in 'vote_levels'),
        'vote_fraction', 'vote_conflict', 'vote_pass' and 'min_votes'.
        'pvals_pooled' is 1e-300 where the gene passes and 1 where it does
        not;
        'z_stat' is the votes, signed by the direction with more votes.
    """
    if progress_callback:
        progress_callback(0, 1)
    if min_votes is None:
        min_votes = len(datasets)

    (gene_names, effects_matrix, var_matrix, valid_mask,
     pathway_map, study_info, _dataset_names, _) = _build_effect_matrices(datasets)

    # _build_effect_matrices falls back to raw p-values for a study without
    # 'pvals_adj'; such a study cannot vote, so it is not testable here.
    has_adj = np.array(['pvals_adj' in df.columns for df in datasets], dtype=bool)
    padj = np.full(effects_matrix.shape, np.nan)
    gene_idx = {g: i for i, g in enumerate(gene_names)}
    for j, df in enumerate(datasets):
        if not has_adj[j]:
            continue
        rows = df.drop_duplicates('names', keep='last')
        idx = rows['names'].map(gene_idx).to_numpy()
        padj[idx, j] = rows['pvals_adj'].to_numpy(float)
    testable = valid_mask & has_adj[None, :] & np.isfinite(padj)

    res = count_votes(effects_matrix, padj, testable, fdr)
    passes = (res['votes'] >= min_votes) & ~res['conflict']

    mean_lfc, se_pooled, _tau2, I2, _w, _Q = _vectorized_dl(
        effects_matrix, var_matrix, testable)
    result = _build_result_df(
        gene_names, mean_lfc, se_pooled, I2, res['n_tested'],
        testable, pathway_map, study_info, min_studies,
        pooling_method='vote_count')
    if not result.empty:
        rows = result['names'].map(gene_idx).to_numpy()
        sign = np.where(res['n_up'][rows] >= res['n_down'][rows], 1.0, -1.0)
        result['n_tested'] = res['n_tested'][rows]
        result['votes'] = res['votes'][rows]
        result['vote_level'] = vote_level(
            res['n_up'][rows], res['n_down'][rows], len(datasets)).to_numpy()
        result['vote_fraction'] = res['votes'][rows] / np.maximum(res['n_tested'][rows], 1)
        result['vote_conflict'] = res['conflict'][rows]
        result['vote_pass'] = passes[rows]
        result['min_votes'] = int(min_votes)
        # Not a p-value: 1e-300 (not 0, which BH treats as invalid) marks
        # a pass, so a caller that adjusts it gets FDR ~0 / 1.
        result['pvals_pooled'] = np.where(passes[rows], 1e-300, 1.0)
        result['z_stat'] = sign * res['votes'][rows]

    if progress_callback:
        progress_callback(1, 1)
    return result


def describe_selection(method_keys, meta_df: pd.DataFrame,
                       n_studies: int, fdr: float = DEFAULT_FDR) -> tuple:
    """How a run picks its genes, in words: (noun, criterion).

    A run of one method calls genes significant; several methods give
    consensus genes; vote counting alone is a rule that genes pass, with
    no FDR. Used for the counts on the Gene Meta-Analysis page.

    Returns
    -------
    tuple of (str, str)
        For example ('consensus genes', 'FDR < 0.05') or
        ('genes pass', 'at least 5 of 5 studies agreeing').
    """
    keys = list(method_keys or [])
    others = [k for k in keys if k != 'vote_count']
    rule = None
    if 'vote_count' in keys and 'min_votes' in meta_df.columns and len(meta_df):
        n = int(meta_df['min_votes'].iloc[0])
        rule = f"at least {n} of {n_studies} studies agreeing"
    if rule and not others:
        return 'genes pass', rule
    criterion = f"FDR < {fdr:g}" + (f" and {rule}" if rule else "")
    return ('consensus genes' if len(keys) > 1 else 'significant genes'), criterion


def vote_level(n_up, n_down, n_studies: int) -> pd.Series:
    """'v/k' per gene (v = votes in the agreeing direction), or 'conflict'."""
    up = np.asarray(n_up, dtype=int)
    down = np.asarray(n_down, dtype=int)
    levels = [f'{v}/{n_studies}' for v in np.maximum(up, down)]
    return pd.Series(np.where((up > 0) & (down > 0), CONFLICT, levels))


def vote_levels(meta_df: pd.DataFrame, n_studies: int,
                significant=None) -> pd.DataFrame:
    """Genes per agreement level, split by direction.

    Parameters
    ----------
    meta_df : DataFrame
        Pooled table with the 'n_up' / 'n_down' columns that
        'annotate_direction' adds.
    n_studies : int
        Studies in the run (k).
    significant : array of bool, optional
        Genes a pooling method calls significant (aligned with
        'meta_df'), counted per level in a 'significant' column.

    Returns
    -------
    DataFrame
        One row per level, '0/k' ... 'k/k' then 'conflict', with columns
        'up', 'down', 'total' and, with 'significant', 'significant'.
    """
    up = meta_df['n_up'].to_numpy(int)
    down = meta_df['n_down'].to_numpy(int)
    level = vote_level(up, down, n_studies).to_numpy()
    order = [f'{v}/{n_studies}' for v in range(n_studies + 1)] + [CONFLICT]
    rows = []
    for lv in order:
        at = level == lv
        rows.append({
            'level': lv,
            'up': int((at & (up > down)).sum()),
            'down': int((at & (down > up)).sum()),
            'total': int(at.sum()),
            **({'significant': int((at & np.asarray(significant, bool)).sum())}
               if significant is not None else {}),
        })
    return pd.DataFrame(rows).set_index('level')
