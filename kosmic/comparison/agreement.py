"""Which genes each method calls, and how the methods agree (ADR-017 §5).

The comparison is made on the genes tested by every method compared
(:func:`shared_universe`). A rule decides which of them each method calls
(:class:`Rule`):

- ``significance`` (default): an adjusted p-value below ``alpha``; for
  vote counting, its own ``vote_pass``;
- ``top_n``: the ``n`` strongest genes by the method's rank, re-ranked
  within the shared genes, so every method calls the same number;
- ``significance_and_direction``: the significance rule, keeping only
  genes whose direction is the same in every method compared.

:func:`pairwise` reports, for each pair of methods, the genes both call,
the Jaccard index, the fraction of the shared calls with the same
direction, and the Spearman correlation of the two methods' ranks over
the shared genes. The last two are reported under every rule.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import combinations

import numpy as np
import pandas as pd

RULES = ('significance', 'top_n', 'significance_and_direction')

RULE_LABELS = {'significance': 'Adjusted p-value below the threshold (vote counting: vote pass)',
               'top_n': 'The top N genes by rank',
               'significance_and_direction': 'Significant, with the same direction in every '
                                             'method'}


@dataclass(frozen=True)
class Rule:
    """A call rule and its thresholds."""

    rule: str = 'significance'
    alpha: float = 0.05
    top_n: int = 100

    def __post_init__(self):
        if self.rule not in RULES:
            raise ValueError(f"Unknown rule '{self.rule}'; expected one of {RULES}.")

    def as_dict(self) -> dict:
        return asdict(self)

    @property
    def text(self) -> str:
        if self.rule == 'top_n':
            return f"top {self.top_n} genes by rank"
        text = f"adjusted p < {self.alpha:g} (vote counting: vote pass)"
        return text + (", same direction in every method"
                       if self.rule == 'significance_and_direction' else '')


def methods_in(results: pd.DataFrame) -> list:
    """The methods in a standard table, in their first-seen order."""
    return list(dict.fromkeys(results['method']))


def shared_universe(results: pd.DataFrame, methods=None) -> set:
    """The genes every one of ``methods`` tested."""
    methods = list(methods or methods_in(results))
    tested = results[results['tested'] & results['method'].isin(methods)]
    counts = tested.groupby('gene')['method'].nunique()
    return set(counts.index[counts == len(methods)])


def apply_rule(results: pd.DataFrame, rule: Rule, methods=None) -> pd.DataFrame:
    """``results`` for ``methods`` with ``in_universe`` and ``significant``.

    Every gene of the input is kept, so the table can be re-cut under
    another rule; ``significant`` is False outside the shared genes.
    """
    methods = list(methods or methods_in(results))
    out = results[results['method'].isin(methods)].copy()
    universe = shared_universe(out, methods)
    out['in_universe'] = out['gene'].isin(universe)
    if rule.rule == 'top_n':
        inside = out[out['in_universe']]
        local = inside.groupby('method')['rank'].rank(method='first')
        sig = pd.Series(False, index=out.index)
        sig[local.index] = local <= rule.top_n
    else:
        vote = out['method'] == 'ma_vote'
        sig = (out['padj'] < rule.alpha) | (vote & out['called'])
        sig = sig & out['in_universe']
        if rule.rule == 'significance_and_direction':
            dirs = out[out['in_universe']].groupby('gene')['direction'].agg(
                lambda d: d.nunique() == 1 and d.iloc[0] != 0)
            agree = set(dirs.index[dirs.astype(bool)])
            sig = sig & out['gene'].isin(agree)
    out['significant'] = sig.fillna(False).astype(bool)
    return out


def call_sets(table: pd.DataFrame) -> dict:
    """``{method: set of genes}`` called in a table from :func:`apply_rule`."""
    return {m: set(table.loc[(table['method'] == m) & table['significant'], 'gene'])
            for m in methods_in(table)}


def pairwise(table: pd.DataFrame) -> pd.DataFrame:
    """One row per pair of methods in a table from :func:`apply_rule`.

    Columns: ``method_a, method_b, n_a, n_b, n_both, jaccard,
    n_same_direction, direction_concordance, spearman``. Concordance is
    over the genes both call (NaN when there are none); the Spearman
    correlation is of the two methods' ranks over the shared genes.
    """
    from scipy.stats import spearmanr

    inside = table[table['in_universe']]
    sets = call_sets(table)
    rank = inside.pivot_table(index='gene', columns='method', values='rank', aggfunc='first')
    direction = inside.pivot_table(index='gene', columns='method', values='direction',
                                   aggfunc='first')
    rows = []
    for a, b in combinations(methods_in(table), 2):
        both = sets[a] & sets[b]
        union = sets[a] | sets[b]
        same = int((direction.loc[list(both), a] == direction.loc[list(both), b]).sum()) \
            if both else 0
        r = rank[[a, b]].dropna().astype(float) if {a, b} <= set(rank.columns) else None
        rho = float(spearmanr(r[a], r[b]).statistic) if r is not None and len(r) > 2 \
            else np.nan
        rows.append({'method_a': a, 'method_b': b, 'n_a': len(sets[a]), 'n_b': len(sets[b]),
                     'n_both': len(both),
                     'jaccard': len(both) / len(union) if union else np.nan,
                     'n_same_direction': same,
                     'direction_concordance': same / len(both) if both else np.nan,
                     'spearman': rho})
    return pd.DataFrame(rows, columns=['method_a', 'method_b', 'n_a', 'n_b', 'n_both',
                                       'jaccard', 'n_same_direction',
                                       'direction_concordance', 'spearman'])


def summary(table: pd.DataFrame) -> pd.DataFrame:
    """Per method: genes in the shared universe, genes called, and how many
    of them go up and down."""
    inside = table[table['in_universe']]
    rows = []
    for m in methods_in(table):
        t = inside[inside['method'] == m]
        called = t[t['significant']]
        rows.append({'method': m, 'n_tested': len(t), 'n_called': len(called),
                     'n_up': int((called['direction'] > 0).sum()),
                     'n_down': int((called['direction'] < 0).sum())})
    return pd.DataFrame(rows, columns=['method', 'n_tested', 'n_called', 'n_up', 'n_down'])


def venn_regions(sets: dict, universe=None) -> pd.DataFrame:
    """Every called gene with the methods that call it, for a Venn diagram.

    ``sets`` maps each method to its called genes (two or three methods, in
    drawing order). Returns one row per gene called by at least one of them:
    ``gene``, ``region`` (the methods that call it, joined with ' + ') and a
    boolean column per method. With ``universe``, a gene outside it is left
    out.
    """
    methods = list(sets)
    genes = sorted(set().union(*sets.values()))
    if universe is not None:
        genes = [g for g in genes if g in universe]
    rows = []
    for g in genes:
        inside = [m for m in methods if g in sets[m]]
        rows.append({'gene': g, 'region': ' + '.join(inside),
                     **{m: g in sets[m] for m in methods}})
    return pd.DataFrame(rows, columns=['gene', 'region', *methods])


def region_counts(sets: dict) -> dict:
    """``{tuple of membership flags: number of genes}`` for every region of
    the methods in ``sets``, (1, 0), (0, 1), (1, 1) and so on, the empty
    region left out."""
    from itertools import product
    methods = list(sets)
    every = set().union(*sets.values())
    out = {}
    for flags in product((0, 1), repeat=len(methods)):
        if not any(flags):
            continue
        out[flags] = sum(1 for g in every
                         if all((g in sets[m]) == bool(f) for m, f in zip(methods, flags)))
    return out
