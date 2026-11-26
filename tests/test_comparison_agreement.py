"""Agreement rules, pairwise measures and the saved comparison (ADR-017 §5, §6)."""
import json

import numpy as np
import pandas as pd
import pytest

from kosmic.comparison import bundle as B
from kosmic.comparison.agreement import (
    Rule, apply_rule, call_sets, pairwise, shared_universe, summary,
)
from kosmic.comparison.standard import COLUMNS


def _rows(method, genes, padj, direction, tested=None, called=None):
    n = len(genes)
    df = pd.DataFrame({'gene': genes, 'method': method, 'effect': np.array(direction) * 1.0,
                       'direction': direction, 'stat': 0.0,
                       'pvalue': padj, 'padj': padj,
                       'tested': [True] * n if tested is None else tested,
                       'called': [False] * n if called is None else called})
    df['rank'] = pd.Series(pd.Series(padj).rank(method='first'), dtype='Int64')
    return df[COLUMNS]


@pytest.fixture
def table():
    """Merged calls A B C D, REML calls A B E; vote counting passes A and E.

    Gene F was not tested by REML, so it is outside the shared genes.
    """
    g = ['A', 'B', 'C', 'D', 'E', 'F']
    merged = _rows('pb_merged', g, [0.001, 0.01, 0.02, 0.03, 0.5, 0.001], [1, 1, -1, 1, 1, 1])
    reml = _rows('ma_reml', g, [0.002, 0.02, 0.4, 0.6, 0.03, np.nan], [1, -1, -1, 1, 1, 0],
                 tested=[True, True, True, True, True, False])
    vote = _rows('ma_vote', g, [np.nan] * 6, [1, 1, -1, 1, 1, 1],
                 called=[True, False, False, False, True, False])
    vote['rank'] = pd.Series([1, 3, 4, 5, 2, 6], dtype='Int64')
    return pd.concat([merged, reml, vote], ignore_index=True)


def test_the_shared_genes_are_those_every_method_tested(table):
    assert shared_universe(table) == {'A', 'B', 'C', 'D', 'E'}
    assert shared_universe(table, ['pb_merged', 'ma_vote']) == set('ABCDEF')


def test_the_significance_rule_uses_padj_and_vote_pass(table):
    sets = call_sets(apply_rule(table, Rule()))
    # F is called by merged but outside the shared genes.
    assert sets == {'pb_merged': {'A', 'B', 'C', 'D'}, 'ma_reml': {'A', 'B', 'E'},
                    'ma_vote': {'A', 'E'}}


def test_the_top_n_rule_calls_n_genes_per_method_within_the_shared_genes(table):
    sets = call_sets(apply_rule(table, Rule('top_n', top_n=2)))
    assert sets == {'pb_merged': {'A', 'B'}, 'ma_reml': {'A', 'B'}, 'ma_vote': {'A', 'E'}}


def test_the_direction_rule_drops_genes_whose_methods_disagree(table):
    sets = call_sets(apply_rule(table, Rule('significance_and_direction')))
    # B goes up in merged and down in REML.
    assert sets['pb_merged'] == {'A', 'C', 'D'} and sets['ma_reml'] == {'A', 'E'}


def test_pairwise_measures_on_known_sets(table):
    pw = pairwise(apply_rule(table, Rule())).set_index(['method_a', 'method_b'])
    row = pw.loc[('pb_merged', 'ma_reml')]
    assert (row['n_a'], row['n_b'], row['n_both']) == (4, 3, 2)
    assert row['jaccard'] == pytest.approx(2 / 5)
    assert row['n_same_direction'] == 1 and row['direction_concordance'] == pytest.approx(0.5)
    # Ranks over A-E: merged 1 2 3 4 5 by padj, REML 1 2 4 5 3.
    assert row['spearman'] == pytest.approx(0.7)
    assert pw.loc[('ma_reml', 'ma_vote'), 'n_both'] == 2


def test_no_shared_calls_gives_no_concordance():
    t = pd.concat([_rows('pb_merged', ['A', 'B'], [0.01, 0.9], [1, 1]),
                   _rows('ma_reml', ['A', 'B'], [0.9, 0.01], [1, 1])], ignore_index=True)
    row = pairwise(apply_rule(t, Rule())).iloc[0]
    assert row['n_both'] == 0 and np.isnan(row['direction_concordance'])


def test_the_summary_counts_calls_by_direction(table):
    s = summary(apply_rule(table, Rule())).set_index('method')
    assert s.loc['pb_merged'].tolist() == [5, 4, 3, 1]


def test_an_unknown_rule_is_refused():
    with pytest.raises(ValueError):
        Rule('majority')


def test_a_saved_comparison_reloads_and_can_be_recut(table, tmp_path):
    b = B.ComparisonBundle.build(table, Rule(), {'design_formula': '~sex + study + condition'})
    folder = b.save(tmp_path, 'decontX_counts/Ventricular_Cardiomyocyte')
    assert folder == tmp_path / 'methods_comparison/decontX_counts/Ventricular_Cardiomyocyte'
    params = json.loads((folder / 'params.json').read_text())
    assert params['design_formula'] == '~sex + study + condition'
    assert params['n_genes_compared'] == 5 and params['n_called']['pb_merged'] == 4
    assert params['rule']['rule'] == 'significance'
    record = json.loads((folder / 'provenance.json').read_text())
    assert [s['stage'] for s in record['stages']] == ['methods_comparison']

    back = B.load(folder)
    assert back.sets == b.sets and back.rule == b.rule
    pd.testing.assert_frame_equal(back.pairwise, b.pairwise)
    assert back.recut(Rule('top_n', top_n=1)).sets['pb_merged'] == {'A'}
    assert B.list_saved(tmp_path) == [('decontX_counts/Ventricular_Cardiomyocyte', folder)]

    # A second save replaces the stage rather than adding one.
    b.save(tmp_path, 'decontX_counts/Ventricular_Cardiomyocyte')
    record = json.loads((folder / 'provenance.json').read_text())
    assert len(record['stages']) == 1


def test_an_unreadable_folder_is_not_a_comparison(tmp_path):
    assert B.load(tmp_path) is None and B.list_saved(tmp_path) == []
