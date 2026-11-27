"""The Methods Comparison's FDR over the genes compared, the meta-analysis
gene-set settings it pools with, and the genes each group of methods calls
(kosmic.comparison.agreement, .ma, .bundle)."""
import numpy as np
import pandas as pd
import pytest

from kosmic.comparison import bundle as B
from kosmic.comparison.agreement import (
    Rule, apply_rule, call_sets, gene_characteristics, method_groups, shared_universe,
)
from kosmic.comparison.ma import MASettings, _settings_from, pool
from kosmic.comparison.standard import COLUMNS
from kosmic.numerical import bh_fdr


def _rows(method, genes, p, direction, padj=None, tested=None, called=None):
    n = len(genes)
    df = pd.DataFrame({'gene': genes, 'method': method,
                       'effect': np.array(direction, dtype=float), 'direction': direction,
                       'stat': 0.0, 'pvalue': p, 'padj': p if padj is None else padj,
                       'tested': [True] * n if tested is None else tested,
                       'called': [False] * n if called is None else called})
    df['rank'] = pd.Series(pd.Series(p).rank(method='first'), dtype='Int64')
    return df[COLUMNS]


@pytest.fixture
def table():
    """Merged tested A-F (F's adjusted p withheld by DESeq2's independent
    filter, raw p kept) and G; REML tested A-F, not G."""
    g = list('ABCDEFG')
    merged = _rows('pb_merged', g, [.001, .004, .01, .2, .5, .002, .001],
                   [1, 1, -1, 1, 1, 1, 1],
                   padj=[.004, .01, .02, .3, .6, np.nan, .004],
                   tested=[True, True, True, True, True, False, True])
    reml = _rows('ma_reml', g, [.002, .01, .3, .4, .02, .003, np.nan],
                 [1, -1, -1, 1, 1, 1, 0],
                 tested=[True] * 6 + [False])
    return pd.concat([merged, reml], ignore_index=True)


def test_compared_genes_are_those_every_method_gave_a_raw_p_value(table):
    assert shared_universe(table) == set('ABCDE')                      # own: padj needed
    assert shared_universe(table, fdr_over='compared') == set('ABCDEF')


def test_each_method_is_corrected_again_over_the_genes_compared(table):
    out = apply_rule(table, Rule())
    assert Rule().fdr_over == 'compared'
    m = out[(out.method == 'pb_merged') & out.in_universe].set_index('gene')
    expected = bh_fdr(np.array([.001, .004, .01, .2, .5, .002]))
    assert m.loc[list('ABCDEF'), 'padj_used'].to_numpy() == pytest.approx(expected)
    assert out.loc[out.gene == 'G', 'padj_used'].isna().all()          # outside
    assert call_sets(out)['pb_merged'] == set('ABCF')


def test_own_keeps_each_methods_adjusted_p_values(table):
    out = apply_rule(table, Rule(fdr_over='own'))
    m = out[(out.method == 'pb_merged') & out.in_universe].set_index('gene')
    assert list(m['padj_used']) == [.004, .01, .02, .3, .6]
    assert call_sets(out)['pb_merged'] == set('ABC')


def test_an_unknown_fdr_scope_is_refused():
    with pytest.raises(ValueError):
        Rule(fdr_over='genome')


def test_the_rule_text_says_where_fdr_is_corrected():
    assert 'over the genes compared' in Rule().text
    assert 'as the method gave it' in Rule(fdr_over='own').text


# --- per gene, and by the methods calling it --------------------------------

def _study(lfc, padj):
    return pd.DataFrame({'names': list('ABC'), 'logfoldchanges': lfc, 'pvals_adj': padj})


STUDIES = [('S1', _study([1.0, -0.5, 0.2], [0.01, 0.5, 0.9])),
           ('S2', _study([0.8, 0.4, 0.1], [0.20, 0.6, 0.9])),
           ('S3', _study([1.2, -0.1, np.nan], [0.03, 0.7, np.nan]))]
REML = pd.DataFrame({'names': list('ABC'), 'logfoldchanges': [1.0, -0.1, 0.15],
                     'heterogeneity_i2': [0.0, 0.6, 0.1],
                     'pi_lower': [0.2, -1.5, -0.4], 'pi_upper': [1.8, 1.3, 0.7]})


def test_gene_characteristics_from_the_studies_and_reml():
    g = gene_characteristics(STUDIES, REML)
    assert list(g['same_direction']) == [True, False, False]   # C missing in S3
    assert list(g['n_studies_significant']) == [2, 0, 0]
    assert list(g['no_single_study']) == [False, True, True]
    assert list(g['pi_excludes_0']) == [True, False, False]
    assert g.at['B', 'i2'] == 0.6


def test_without_reml_the_reml_columns_are_empty():
    g = gene_characteristics(STUDIES, None)
    assert g['pi_excludes_0'].isna().all() and g['i2'].isna().all()


def test_method_groups_summarise_the_genes_each_group_calls():
    t = pd.concat([_rows('pb_merged', list('ABC'), [.001, .001, .9], [1, -1, 1]),
                   _rows('ma_reml', list('ABC'), [.001, .9, .9], [1, -1, 1])],
                  ignore_index=True)
    groups = method_groups(apply_rule(t, Rule(fdr_over='own')),
                           gene_characteristics(STUDIES, REML)).set_index('methods')
    assert list(groups.index) == ['pb_merged + ma_reml', 'pb_merged']
    assert groups.at['pb_merged + ma_reml', 'genes'] == 1
    assert groups.at['pb_merged + ma_reml', 'pi_excludes_0'] == 1.0
    assert groups.at['pb_merged', 'no_single_study'] == 1.0
    assert groups.at['pb_merged', 'median_abs_reml_effect'] == pytest.approx(0.1)


def test_a_saved_comparison_keeps_its_genes_and_groups(tmp_path):
    t = pd.concat([_rows('pb_merged', list('ABC'), [.001, .001, .9], [1, -1, 1]),
                   _rows('ma_reml', list('ABC'), [.001, .9, .9], [1, -1, 1])],
                  ignore_index=True)
    bundle = B.ComparisonBundle.build(t, Rule(), {}, genes=gene_characteristics(STUDIES, REML))
    folder = bundle.save(tmp_path, 'decontX_counts/Endothelial_Cell')
    assert (folder / B.GROUPS_FILE).is_file() and (folder / B.GENES_FILE).is_file()
    again = B.load(folder)
    assert again.rule.fdr_over == 'compared'
    pd.testing.assert_frame_equal(again.method_groups, bundle.method_groups,
                                  check_dtype=False)
    assert again.recut(Rule(fdr_over='own')).genes is not None


def test_a_comparison_saved_before_the_fdr_scope_reads_as_own(tmp_path):
    t = _rows('pb_merged', list('AB'), [.01, .9], [1, 1])
    folder = B.ComparisonBundle.build(t, Rule(fdr_over='own'), {}).save(tmp_path, 'c/x')
    import json
    params = json.loads((folder / B.PARAMS_FILE).read_text())
    params['rule'].pop('fdr_over')
    (folder / B.PARAMS_FILE).write_text(json.dumps(params))
    assert B.load(folder).rule.fdr_over == 'own'


# --- the meta-analysis gene set ----------------------------------------------

def _de(names, seed):
    rng = np.random.default_rng(seed)
    n = len(names)
    return pd.DataFrame({'names': names, 'logfoldchanges': rng.normal(0, 1, n),
                         'se': rng.uniform(.1, .3, n), 'pvals': rng.uniform(.001, .9, n),
                         'pct_disease': .5, 'pct_control': .5})


def test_the_saved_runs_gene_set_settings_are_read():
    s = _settings_from({'min_studies': 3, 'current_symbols_only': True,
                        'drop_inconsistent_merges': True}, 'run.csv')
    assert s.current_symbols_only and s.drop_inconsistent_merges
    old = _settings_from({'min_studies': 3}, 'old.csv')
    assert not old.current_symbols_only and not old.drop_inconsistent_merges


def test_pooling_uses_the_gene_set_settings():
    names = ['TTN', 'NPPA', 'HNRNPU', 'ENSG00000000001']
    tables = [(f'S{i}', _de(names, i)) for i in range(3)]
    plain = MASettings(min_studies=3, filter_active=False, filter_threshold=None,
                       filter_stat=None)
    assert set(pool(tables, 'reml', plain)['names']) == set(names)
    rules = MASettings(min_studies=3, filter_active=False, filter_threshold=None,
                       filter_stat=None, current_symbols_only=True,
                       drop_inconsistent_merges=True)
    assert set(pool(tables, 'reml', rules, exclude={'HNRNPU'})['names']) == {'TTN', 'NPPA'}
