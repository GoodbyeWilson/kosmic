"""One long table for every method's result (ADR-017 §5).

Each adapter takes one method's own output and returns rows with the
columns of :data:`COLUMNS`:

- ``gene``, ``method``;
- ``effect``: the log2 fold change (merged DESeq2's, the pooled effect;
  SumRank's mean of the studies' log2 fold changes); NA where undefined;
- ``direction``: -1, 0 or 1;
- ``stat``: the method's statistic (DESeq2 Wald statistic, pooled Z,
  SumRank rank sum, vote counting's signed number of votes);
- ``pvalue``, ``padj``: raw and BH-adjusted p-values (NA for vote counting);
- ``tested``: whether the method produced a result for the gene (an
  adjusted p-value, or for vote counting a vote result);
- ``called``: the method's own call under the default rule
  (``padj < DEFAULT_FDR``, or ``vote_pass``), which
  :mod:`kosmic.comparison.agreement` re-applies for the rule chosen;
- ``rank``: 1 for the strongest tested gene; by p-value, ties by larger
  |effect| (vote counting: by more votes, then larger |effect|). Untested
  genes have no rank.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

COLUMNS = ['gene', 'method', 'effect', 'direction', 'stat', 'pvalue', 'padj', 'tested',
           'called', 'rank']

#: Method id -> display name.
METHOD_LABELS = {'pb_merged': 'Merged pseudobulk', 'ma_reml': 'REML',
                 'ma_sumrank': 'SumRank', 'ma_stouffer': 'Stouffer',
                 'ma_vote': 'Vote counting'}


def _num(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.nan, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors='coerce')


def _finish(out: pd.DataFrame, rank_by: list, ascending: list) -> pd.DataFrame:
    """Drop duplicate genes, add ``rank`` over tested genes, order columns."""
    out = out.dropna(subset=['gene'])
    out['gene'] = out['gene'].astype(str)
    out = out.sort_values(['tested', *rank_by], ascending=[False, *ascending],
                          na_position='last').drop_duplicates('gene')
    rank = pd.Series(pd.NA, index=out.index, dtype='Int64')
    tested = out['tested'].to_numpy(bool)
    rank[tested] = np.arange(1, int(tested.sum()) + 1)
    out['rank'] = rank
    out['direction'] = out['direction'].fillna(0).astype(int)
    out['called'] = out['called'].fillna(False).astype(bool) & out['tested']
    return out[COLUMNS].reset_index(drop=True)


def _p_based(df, method, effect, stat, pvalue, padj, alpha, direction=None):
    out = pd.DataFrame({'gene': df['names'], 'method': method, 'effect': effect,
                        'stat': stat, 'pvalue': pvalue, 'padj': padj})
    out['direction'] = np.sign(out['effect']) if direction is None else direction
    out['tested'] = np.isfinite(out['padj']) & np.isfinite(out['pvalue'])
    out['called'] = out['padj'] < alpha
    out['abs_effect'] = -out['effect'].abs()
    return _finish(out, ['pvalue', 'abs_effect'], [True, True])


def merged_to_standard(df: pd.DataFrame, alpha: float = None) -> pd.DataFrame:
    """Merged pseudobulk: a DESeq2 table (``names, logfoldchanges, pvals,
    pvals_adj``, ``statistic``). Genes the DESeq2 filters removed have no
    adjusted p-value and are not tested."""
    from kosmic import DEFAULT_FDR
    return _p_based(df, 'pb_merged', _num(df, 'logfoldchanges'), _num(df, 'statistic'),
                    _num(df, 'pvals'), _num(df, 'pvals_adj'),
                    DEFAULT_FDR if alpha is None else alpha)


def _pooled_p(df, method, alpha):
    from kosmic import DEFAULT_FDR
    from kosmic.numerical import bh_fdr
    pvalue = _num(df, 'pvals_pooled')
    padj = _num(df, 'fdr') if 'fdr' in df.columns else pd.Series(
        bh_fdr(pvalue.values), index=df.index)
    padj = padj.where(np.isfinite(pvalue))
    return pvalue, padj, DEFAULT_FDR if alpha is None else alpha


def reml_to_standard(df: pd.DataFrame, alpha: float = None) -> pd.DataFrame:
    """REML: the pooled effect and its Z (``z_stat``)."""
    pvalue, padj, a = _pooled_p(df, 'ma_reml', alpha)
    return _p_based(df, 'ma_reml', _num(df, 'logfoldchanges'), _num(df, 'z_stat'),
                    pvalue, padj, a)


def stouffer_to_standard(df: pd.DataFrame, alpha: float = None) -> pd.DataFrame:
    """Stouffer: the inverse-variance pooled effect and the Stouffer Z."""
    pvalue, padj, a = _pooled_p(df, 'ma_stouffer', alpha)
    return _p_based(df, 'ma_stouffer', _num(df, 'logfoldchanges'), _num(df, 'z_stat'),
                    pvalue, padj, a)


def sumrank_to_standard(df: pd.DataFrame, alpha: float = None) -> pd.DataFrame:
    """SumRank: the rank sum as statistic, the mean log2 fold change as effect,
    and its own ``direction`` where it gives one."""
    pvalue, padj, a = _pooled_p(df, 'ma_sumrank', alpha)
    effect = _num(df, 'logfoldchanges')
    direction = np.sign(effect)
    if 'direction' in df.columns:
        named = df['direction'].astype(str).str.lower().map({'up': 1.0, 'down': -1.0})
        direction = named.where(named.notna(), direction)
    return _p_based(df, 'ma_sumrank', effect, _num(df, 'rank_sum'), pvalue, padj, a,
                    direction=direction)


def vote_to_standard(df: pd.DataFrame, alpha: float = None) -> pd.DataFrame:
    """Vote counting: called where ``vote_pass``; the statistic is the signed
    number of votes (``z_stat``); no p-value."""
    stat = _num(df, 'z_stat')
    effect = _num(df, 'logfoldchanges')
    out = pd.DataFrame({'gene': df['names'], 'method': 'ma_vote', 'effect': effect,
                        'stat': stat, 'pvalue': np.nan, 'padj': np.nan})
    out['direction'] = np.where(stat.fillna(0) != 0, np.sign(stat.fillna(0)), np.sign(effect))
    tested = _num(df, 'n_tested') if 'n_tested' in df.columns else _num(df, 'n_studies')
    out['tested'] = tested.notna() & (tested > 0)
    passed = df['vote_pass'] if 'vote_pass' in df.columns else pd.Series(False, index=df.index)
    out['called'] = passed.astype(str).str.lower().isin(['true', '1'])
    out['votes'] = -stat.abs()
    out['abs_effect'] = -effect.abs()
    return _finish(out, ['votes', 'abs_effect'], [True, True])


ADAPTERS = {'pb_merged': merged_to_standard, 'ma_reml': reml_to_standard,
            'ma_sumrank': sumrank_to_standard, 'ma_stouffer': stouffer_to_standard,
            'ma_vote': vote_to_standard}


def to_standard(method: str, df: pd.DataFrame, alpha: float = None) -> pd.DataFrame:
    """The standard rows for one method's output."""
    return ADAPTERS[method](df, alpha)
