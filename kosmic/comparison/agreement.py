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

Where each method's FDR is corrected is part of the rule (``fdr_over``):

- ``compared`` (default): every p-value method is corrected again, with
  Benjamini-Hochberg, over the genes compared, from its raw p-values, so
  every method's FDR has the same denominator. A p-value method has tested
  a gene when it gave a raw p-value; DESeq2's independent filter withholds
  only the adjusted p-value, and the correction is made again here, as
  the DE workspace's committed-gene FDR does. Merged pseudobulk fitted
  genome-wide is then corrected over the same genes as the meta-analysis
  without refitting: the fit, and so the raw p-values, do not depend on it.
- ``own``: each method's adjusted p-value as the method gave it, and the
  genes compared are those with one in every method.

:func:`gene_characteristics` and :func:`method_groups` describe the genes
each group of methods calls, from the studies' own results and the REML
fit: whether the studies agree in direction, how heterogeneous they are,
whether the 95% prediction interval excludes 0, and whether any single
study calls the gene.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import combinations

import numpy as np
import pandas as pd

RULES = ('significance', 'top_n', 'significance_and_direction')
FDR_SCOPES = ('compared', 'own')

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
    fdr_over: str = 'compared'

    def __post_init__(self):
        if self.rule not in RULES:
            raise ValueError(f"Unknown rule '{self.rule}'; expected one of {RULES}.")
        if self.fdr_over not in FDR_SCOPES:
            raise ValueError(f"Unknown FDR scope '{self.fdr_over}'; expected one of "
                             f"{FDR_SCOPES}.")

    def as_dict(self) -> dict:
        return asdict(self)

    @property
    def text(self) -> str:
        if self.rule == 'top_n':
            return f"top {self.top_n} genes by rank"
        text = f"adjusted p < {self.alpha:g} (vote counting: vote pass)"
        text += (", same direction in every method"
                 if self.rule == 'significance_and_direction' else '')
        return text + ("; each method's FDR corrected over the genes compared"
                       if self.fdr_over == 'compared' else
                       "; each method's FDR as the method gave it")


def methods_in(results: pd.DataFrame) -> list:
    """The methods in a standard table, in their first-seen order."""
    return list(dict.fromkeys(results['method']))


def _testable(results: pd.DataFrame, fdr_over: str = 'own') -> pd.Series:
    """Whether each row's method tested its gene: under ``own``, the
    standard table's ``tested``; under ``compared``, a raw p-value for a
    p-value method (vote counting keeps ``tested``)."""
    tested = results['tested'].astype(bool)
    if fdr_over != 'compared':
        return tested
    p = pd.to_numeric(results['pvalue'], errors='coerce')
    vote = results['method'] == 'ma_vote'
    return (vote & tested) | (~vote & np.isfinite(p))


def shared_universe(results: pd.DataFrame, methods=None, fdr_over: str = 'own') -> set:
    """The genes every one of ``methods`` tested."""
    methods = list(methods or methods_in(results))
    tested = results[_testable(results, fdr_over) & results['method'].isin(methods)]
    counts = tested.groupby('gene')['method'].nunique()
    return set(counts.index[counts == len(methods)])


def apply_rule(results: pd.DataFrame, rule: Rule, methods=None) -> pd.DataFrame:
    """``results`` for ``methods`` with ``in_universe``, ``padj_used`` and
    ``significant``.

    Every gene of the input is kept, so the table can be re-cut under
    another rule; ``significant`` is False outside the shared genes.
    ``padj_used`` is the adjusted p-value the call uses: under
    ``fdr_over='compared'`` the method's raw p-values corrected over the
    shared genes (NA outside them), else the method's own.
    """
    from kosmic.numerical import bh_fdr

    methods = list(methods or methods_in(results))
    out = results[results['method'].isin(methods)].copy()
    universe = shared_universe(out, methods, rule.fdr_over)
    out['in_universe'] = out['gene'].isin(universe)
    vote = out['method'] == 'ma_vote'
    if rule.fdr_over == 'compared':
        used = pd.Series(np.nan, index=out.index, dtype=float)
        for m in methods:
            rows = out.index[(out['method'] == m) & out['in_universe'] & ~vote]
            if len(rows):
                used[rows] = bh_fdr(pd.to_numeric(out.loc[rows, 'pvalue']).to_numpy(float))
        out['padj_used'] = used
    else:
        out['padj_used'] = pd.to_numeric(out['padj'], errors='coerce')
    if rule.rule == 'top_n':
        # Within the shared genes, in the method's own order: p-value, then
        # larger |effect| (vote counting: its rank, by votes).
        inside = out[out['in_universe']].copy()
        inside['_k1'] = np.where(inside['method'] == 'ma_vote',
                                 pd.to_numeric(inside['rank'], errors='coerce'),
                                 pd.to_numeric(inside['pvalue'], errors='coerce'))
        inside['_k2'] = -pd.to_numeric(inside['effect'], errors='coerce').abs()
        inside = inside.sort_values(['_k1', '_k2'], na_position='last', kind='mergesort')
        local = inside.groupby('method').cumcount() + 1
        sig = pd.Series(False, index=out.index)
        sig[local.index] = local <= rule.top_n
    else:
        sig = (out['padj_used'] < rule.alpha) | (vote & out['called'])
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

    ``sets`` maps each method to its called genes (two to four methods, in
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


def gene_characteristics(study_tables, reml: pd.DataFrame = None,
                         alpha: float = None) -> pd.DataFrame:
    """Per gene, what the studies and the REML fit say about it.

    ``study_tables`` is ``[(name, DataFrame)]``, each study's own DE table
    (``names``, ``logfoldchanges``, ``pvals_adj``); ``reml`` the REML pooled
    table, or None. Columns:

    - ``n_studies``: studies reporting a log2 fold change;
    - ``same_direction``: reported by every study, all with the same sign;
    - ``n_studies_significant``: studies whose own FDR is below ``alpha``;
    - ``no_single_study``: none of them;
    - ``reml_effect``, ``i2``, ``pi_lower``, ``pi_upper`` from REML, and
      ``pi_excludes_0``: the 95% prediction interval (pooled +/- t(k-2)
      sqrt(SE^2 + tau^2)) lies wholly above or below 0. NA without REML.
    """
    from kosmic import DEFAULT_FDR
    alpha = DEFAULT_FDR if alpha is None else alpha
    lfc, sig = {}, {}
    for name, df in study_tables:
        d = df.drop_duplicates('names').set_index('names')
        lfc[name] = pd.to_numeric(d['logfoldchanges'], errors='coerce')
        padj = (pd.to_numeric(d['pvals_adj'], errors='coerce') if 'pvals_adj' in d
                else pd.Series(np.nan, index=d.index))
        sig[name] = padj < alpha
    L = pd.DataFrame(lfc)
    S = pd.DataFrame(sig).reindex(L.index).fillna(False).astype(bool)
    out = pd.DataFrame(index=L.index)
    out.index.name = 'gene'
    out['n_studies'] = L.notna().sum(axis=1)
    signs = np.sign(L)
    out['same_direction'] = (L.notna().all(axis=1) & (signs.nunique(axis=1) == 1)
                             & (signs.abs().min(axis=1) > 0))
    out['n_studies_significant'] = S.sum(axis=1)
    out['no_single_study'] = out['n_studies_significant'] == 0
    cols = {'reml_effect': 'logfoldchanges', 'i2': 'heterogeneity_i2',
            'pi_lower': 'pi_lower', 'pi_upper': 'pi_upper'}
    if reml is not None and len(reml):
        r = reml.drop_duplicates('names').set_index('names')
        for new, old in cols.items():
            out[new] = (pd.to_numeric(r[old], errors='coerce').reindex(out.index)
                        if old in r else np.nan)
        known = out['pi_lower'].notna() & out['pi_upper'].notna()
        excl = (out['pi_lower'] > 0) | (out['pi_upper'] < 0)
        out['pi_excludes_0'] = excl.astype(object).where(known, np.nan)
    else:
        for new in cols:
            out[new] = np.nan
        out['pi_excludes_0'] = np.nan
    return out


GROUP_COLUMNS = ['methods', 'n_methods', 'genes', 'same_direction',
                       'median_abs_reml_effect', 'median_i2', 'pi_excludes_0',
                       'no_single_study']


def method_groups(table: pd.DataFrame, genes: pd.DataFrame = None) -> pd.DataFrame:
    """One row per group of methods that call the same genes.

    For every shared gene called by at least one method in ``table`` (from
    :func:`apply_rule`), the group is exactly the methods calling it
    (``methods``, joined with ' + ' in the table's method order). Per group:
    the number of genes and, from ``genes`` (:func:`gene_characteristics`),
    the fraction with the same direction in every study, the median |REML
    log2 fold change|, the median I^2, the fraction whose prediction
    interval excludes 0 and the fraction no single study calls.
    """
    order = methods_in(table)
    called = table[table['in_universe'] & table['significant']]
    groups = called.groupby('gene')['method'].agg(
        lambda ms: ' + '.join(m for m in order if m in set(ms)))
    g = genes.reindex(groups.index) if genes is not None else None

    def frac(col, idx):
        if g is None or col not in g:
            return np.nan
        v = g.loc[idx, col].dropna().astype(float)
        return float(v.mean()) if len(v) else np.nan

    def med(col, idx, absolute=False):
        if g is None or col not in g:
            return np.nan
        v = pd.to_numeric(g.loc[idx, col], errors='coerce').dropna()
        v = v.abs() if absolute else v
        return float(v.median()) if len(v) else np.nan

    rows = []
    for label, idx in groups.groupby(groups).groups.items():
        rows.append({'methods': label, 'n_methods': label.count(' + ') + 1,
                     'genes': len(idx), 'same_direction': frac('same_direction', idx),
                     'median_abs_reml_effect': med('reml_effect', idx, absolute=True),
                     'median_i2': med('i2', idx),
                     'pi_excludes_0': frac('pi_excludes_0', idx),
                     'no_single_study': frac('no_single_study', idx)})
    out = pd.DataFrame(rows, columns=GROUP_COLUMNS)
    return out.sort_values(['n_methods', 'genes'], ascending=[False, False]
                           ).reset_index(drop=True)
