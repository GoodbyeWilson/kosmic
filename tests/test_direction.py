"""Per-study direction counts on pooled tables (issue #62).

'annotate_direction' counts, per gene, the studies with a significant
positive and a significant negative effect, and flags genes where both
occur. These tests build studies with known calls and check the counts,
the handling of a study without adjusted p-values, and that the columns
survive the pooling dispatcher and the consensus merge.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from kosmic.meta_analysis.direction import (
    DIRECTION_COLUMNS,
    annotate_direction,
    studies_without_adjusted_p,
)
from kosmic.meta_analysis.pooling_dispatch import (
    get_analytical_pool_fn,
    merge_consensus_methods,
)


def _study(name, rows, adjusted=True):
    """rows: list of (gene, log2fc, p). 'p' goes to pvals_adj or pvals."""
    df = pd.DataFrame(rows, columns=['names', 'logfoldchanges', 'p'])
    df['se'] = 0.2
    df['dataset'] = name
    return df.rename(columns={'p': 'pvals_adj' if adjusted else 'pvals'})


def _three_studies():
    # CONFLICT: up in A, down in B.       UP: up in A and B.
    # WEAK: opposite signs, none significant.  ABSENT: only in C.
    a = _study('A', [('CONFLICT', 1.0, 0.001), ('UP', 0.8, 0.01),
                     ('WEAK', 0.5, 0.40)])
    b = _study('B', [('CONFLICT', -1.2, 0.002), ('UP', 0.9, 0.03),
                     ('WEAK', -0.5, 0.30)])
    c = _study('C', [('CONFLICT', 0.1, 0.90), ('UP', 0.7, 0.20),
                     ('ABSENT', -2.0, 0.001)])
    return [a, b, c]


def _meta(names):
    return pd.DataFrame({'names': names, 'pvals_pooled': 0.01})


def test_counts_and_conflict_flag():
    out = annotate_direction(
        _meta(['CONFLICT', 'UP', 'WEAK', 'ABSENT', 'MISSING']),
        _three_studies(), fdr=0.05).set_index('names')

    assert out.loc['CONFLICT', ['n_up', 'n_down']].tolist() == [1, 1]
    assert bool(out.loc['CONFLICT', 'direction_conflict'])

    assert out.loc['UP', ['n_up', 'n_down']].tolist() == [2, 0]
    assert not out.loc['UP', 'direction_conflict']

    # Opposite signs without significance is not a conflict.
    assert out.loc['WEAK', ['n_up', 'n_down']].tolist() == [0, 0]
    assert not out.loc['WEAK', 'direction_conflict']

    assert out.loc['ABSENT', ['n_up', 'n_down']].tolist() == [0, 1]
    assert out.loc['MISSING', ['n_up', 'n_down']].tolist() == [0, 0]


def test_threshold_is_strict_and_configurable():
    studies = [_study('A', [('G', 1.0, 0.05)]),
               _study('B', [('G', -1.0, 0.01)])]
    at_005 = annotate_direction(_meta(['G']), studies, fdr=0.05)
    assert at_005[['n_up', 'n_down']].iloc[0].tolist() == [0, 1]
    at_010 = annotate_direction(_meta(['G']), studies, fdr=0.10)
    assert bool(at_010['direction_conflict'].iloc[0])


def test_study_without_adjusted_p_is_not_counted():
    studies = [_study('A', [('G', 1.0, 0.001)]),
               _study('RAW', [('G', -1.0, 0.001)], adjusted=False)]
    out = annotate_direction(_meta(['G']), studies)
    assert out[['n_up', 'n_down']].iloc[0].tolist() == [1, 0]
    assert not out['direction_conflict'].iloc[0]
    assert studies_without_adjusted_p(studies) == ['RAW']


def test_unnamed_study_is_reported_by_position():
    df = _study('X', [('G', 1.0, 0.1)], adjusted=False).drop(columns='dataset')
    assert studies_without_adjusted_p([_three_studies()[0], df]) == ['Study_2']


def test_input_not_modified_and_empty_table():
    meta = _meta(['CONFLICT'])
    annotate_direction(meta, _three_studies())
    assert 'n_up' not in meta.columns

    empty = annotate_direction(_meta([]), _three_studies())
    assert list(empty.columns[-3:]) == ['n_up', 'n_down', 'direction_conflict']
    assert len(empty) == 0


def _pooling_fixture():
    """30 genes x 3 studies; G0 is strongly up in S1 and down in S2."""
    rng = np.random.default_rng(0)
    genes = [f'G{i}' for i in range(30)]
    studies = []
    for s in range(3):
        lfc = rng.normal(0, 0.3, 30)
        padj = rng.uniform(0.2, 1.0, 30)
        if s == 0:
            lfc[0], padj[0] = 2.0, 1e-6
        elif s == 1:
            lfc[0], padj[0] = -2.0, 1e-6
        studies.append(pd.DataFrame({
            'names': genes, 'logfoldchanges': lfc, 'se': 0.2,
            'pvals': padj / 10, 'pvals_adj': padj, 'dataset': f'S{s + 1}',
        }))
    return studies


def test_dispatcher_output_carries_direction_columns():
    studies = _pooling_fixture()
    for key in ('dl', 'fisher', 'sumrank', 'gwop', 'stouffer'):
        df = get_analytical_pool_fn(key, 2)(studies).set_index('names')
        assert {'n_up', 'n_down', 'direction_conflict'} <= set(df.columns), key
        assert bool(df.loc['G0', 'direction_conflict']), key
        assert int(df['direction_conflict'].sum()) == 1, key


def test_consensus_merge_keeps_direction_columns():
    studies = _pooling_fixture()
    keys = ['dl', 'fisher']
    dfs = {k: get_analytical_pool_fn(k, 2)(studies) for k in keys}
    merged = merge_consensus_methods(keys, dfs).set_index('names')
    assert bool(merged.loc['G0', 'direction_conflict'])
    assert merged.loc['G0', ['n_up', 'n_down']].tolist() == [1, 1]


def test_count_significant_conflicts():
    from kosmic.meta_analysis.direction import count_significant_conflicts
    df = pd.DataFrame({
        'fdr': [0.01, 0.02, 0.20, 0.01],
        'direction_conflict': [True, False, True, False],
    })
    assert count_significant_conflicts(df, fdr=0.05) == (3, 1)
    assert count_significant_conflicts(df.drop(columns='direction_conflict')) == (3, None)
    assert count_significant_conflicts(df.drop(columns='fdr')) == (0, None)


def test_no_countable_study_adds_no_columns():
    """Absent counts must not read as 'no conflict' (DESeq2+VIF pathway route)."""
    studies = [_study('A', [('G', 1.0, 0.001)], adjusted=False),
               _study('B', [('G', -1.0, 0.001)], adjusted=False)]
    out = annotate_direction(_meta(['G']), studies)
    assert not set(DIRECTION_COLUMNS) & set(out.columns)


def test_unadjusted_warning():
    from kosmic.meta_analysis.direction import unadjusted_warning
    adj = _study('A', [('G', 1.0, 0.001)])
    raw = _study('RAW', [('G', 1.0, 0.001)], adjusted=False)
    assert unadjusted_warning([adj]) is None
    assert 'RAW' in unadjusted_warning([adj, raw])
    assert 'not computed' in unadjusted_warning([raw])


def test_add_wald_pvalues_matches_hand_calculation():
    from scipy.stats import norm

    from kosmic.meta_analysis.direction import add_wald_pvalues
    from kosmic.numerical import bh_fdr

    df = pd.DataFrame({'names': ['P1', 'P2', 'P3', 'P4'],
                       'logfoldchanges': [1.0, -0.6, 0.1, 0.5],
                       'se': [0.25, 0.2, 0.3, np.nan]})
    out = add_wald_pvalues(df)
    expected = 2 * norm.sf(np.abs([4.0, 3.0, 0.1 / 0.3]))
    np.testing.assert_allclose(out['pvals'].iloc[:3], expected)
    # BH over the three valid rows only; the SE-less row gets 1.
    np.testing.assert_allclose(out['pvals_adj'].iloc[:3], bh_fdr(expected))
    assert out['pvals'].iloc[3] == 1.0 and out['pvals_adj'].iloc[3] == 1.0
    assert 'pvals' not in df.columns


def test_wald_pvalues_do_not_change_effect_size_pooling():
    """DL and REML pool the effects; the added p-values must not alter them."""
    from kosmic.meta_analysis.direction import add_wald_pvalues

    bare = [df.drop(columns=['pvals', 'pvals_adj']) for df in _pooling_fixture()]
    with_p = [add_wald_pvalues(df) for df in bare]
    cols = ['names', 'logfoldchanges', 'se', 'pvals_pooled']
    for key in ('dl', 'reml'):
        a = get_analytical_pool_fn(key, 2)(bare)[cols]
        b = get_analytical_pool_fn(key, 2)(with_p)[cols]
        pd.testing.assert_frame_equal(a, b)
        assert 'direction_conflict' not in get_analytical_pool_fn(key, 2)(bare)
