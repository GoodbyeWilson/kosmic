"""Vote counting: a gene passes when at least N studies call it significant
in the same direction and none in the opposite one (#30)."""
from __future__ import annotations

import pandas as pd
import pytest

from kosmic.meta_analysis.pooling.vote_count import vote_count_fast, vote_levels
from kosmic.meta_analysis.pooling_dispatch import (
    get_analytical_pool_fn, merge_consensus_methods,
)


def _study(name, rows):
    """rows: {gene: (log2fc, padj)}."""
    return pd.DataFrame({
        'names': list(rows),
        'logfoldchanges': [v[0] for v in rows.values()],
        'se': 0.3,
        'pvals': [v[1] / 10 for v in rows.values()],
        'pvals_adj': [v[1] for v in rows.values()],
        'dataset': name,
    })


@pytest.fixture
def studies():
    return [
        _study('A', {'UP': (1.0, 0.001), 'MIX': (1.0, 0.001), 'NONE': (0.1, 0.6)}),
        _study('B', {'UP': (1.2, 0.002), 'MIX': (-1.0, 0.001), 'NONE': (-0.1, 0.7)}),
        _study('C', {'UP': (0.9, 0.010), 'MIX': (1.1, 0.002), 'NONE': (0.2, 0.5)}),
        _study('D', {'MIX': (0.5, 0.3), 'NONE': (0.0, 0.9)}),   # UP not measured
    ]


def _row(df, gene):
    return df.set_index('names').loc[gene]


def test_votes_count_significant_studies_in_one_direction(studies):
    res = vote_count_fast(studies, min_studies=2, min_votes=3)
    up, mix = _row(res, 'UP'), _row(res, 'MIX')
    assert (up['votes'], up['n_tested']) == (3, 3)        # absent from D: not in the denominator
    assert (mix['votes'], mix['n_tested']) == (2, 4)      # 2 up, 1 down, 1 not significant
    assert _row(res, 'NONE')['votes'] == 0


def test_a_gene_passes_with_n_votes_and_no_conflict(studies):
    res = vote_count_fast(studies, min_studies=2, min_votes=2)
    assert bool(_row(res, 'UP')['vote_pass'])
    assert not _row(res, 'MIX')['vote_pass']              # 2 up but 1 down: a conflict
    assert bool(_row(res, 'MIX')['vote_conflict'])
    assert not _row(res, 'NONE')['vote_pass']
    # Every gene keeps its level, as the Study Agreement figure counts it.
    assert [_row(res, g)['vote_level'] for g in ('UP', 'MIX', 'NONE')] == [
        '3/4', 'conflict', '0/4']
    assert not _row(vote_count_fast(studies, min_studies=2), 'UP')['vote_pass']  # default: all 4


def test_the_rule_is_an_fdr_of_0_or_1_in_the_consensus(studies):
    vote = get_analytical_pool_fn('vote_count', 2, min_votes=3)(studies)
    reml = get_analytical_pool_fn('reml', 2)(studies)
    merged = merge_consensus_methods(['reml', 'vote_count'],
                                     {'reml': reml, 'vote_count': vote}).set_index('names')
    assert merged.loc['UP', 'fdr_vote_count'] == 0.0
    assert merged.loc['MIX', 'fdr_vote_count'] == 1.0
    assert merged.loc['MIX', 'fdr'] == 1.0                # fails the rule, fails consensus
    assert merged.loc['UP', 'votes'] == 3                 # carried for the Results table
    assert merged.loc['UP', 'vote_level'] == '3/4'


def test_a_study_without_adjusted_p_values_cannot_vote(studies):
    studies[0] = studies[0].drop(columns='pvals_adj')
    up = _row(vote_count_fast(studies, min_studies=2), 'UP')
    assert (up['votes'], up['n_tested']) == (2, 2)


def test_min_studies_admits_on_the_tested_count(studies):
    res = vote_count_fast(studies, min_studies=4)
    assert 'UP' not in set(res['names']) and 'MIX' in set(res['names'])


def test_levels_count_genes_by_studies_agreeing(studies):
    pooled = get_analytical_pool_fn('reml', 2)(studies)
    levels = vote_levels(pooled, 4, significant=pooled['names'].eq('UP'))
    assert list(levels.index) == ['0/4', '1/4', '2/4', '3/4', '4/4', 'conflict']
    assert levels.loc['3/4', ['up', 'down', 'total', 'significant']].tolist() == [1, 0, 1, 1]
    assert levels.loc['conflict', 'total'] == 1 and levels.loc['0/4', 'total'] == 1


def test_cc_permutation_leaves_the_rule_alone():
    from kosmic.meta_analysis.cc_permutation import _MATRIX_POOLERS
    from kosmic.meta_analysis.pooling_dispatch import RULE_METHODS
    assert not set(RULE_METHODS) & set(_MATRIX_POOLERS)


def test_a_saved_method_selection_does_not_hide_vote_counting(tmp_path):
    from PyQt6.QtCore import QSettings

    from kosmic.gui.meta_analysis.dialogs.ma_settings_dialog import (
        load_enabled_methods, save_enabled_methods)
    settings = QSettings(str(tmp_path / 's.ini'), QSettings.Format.IniFormat)
    settings.setValue('meta_analysis/enabled_methods', ['dl', 'fisher'])   # saved before
    assert load_enabled_methods(settings) == {'dl', 'fisher', 'vote_count'}
    save_enabled_methods(settings, {'dl'})            # now unticked on purpose
    assert load_enabled_methods(settings) == {'dl'}


def test_the_selection_is_worded_by_the_methods_run(studies):
    from kosmic.meta_analysis.pooling.vote_count import describe_selection
    vote = get_analytical_pool_fn('vote_count', 2, min_votes=3)(studies)
    assert describe_selection(['vote_count'], vote, 4) == (
        'genes pass', 'at least 3 of 4 studies agreeing')
    assert describe_selection(['reml'], vote, 4) == ('significant genes', 'FDR < 0.05')
    assert describe_selection(['reml', 'vote_count'], vote, 4) == (
        'consensus genes', 'FDR < 0.05 and at least 3 of 4 studies agreeing')
