"""Methods that combine study-level p-values use the raw p-values (ADR-018).

Adjusted p-values are not uniform under the null, so Fisher, Stouffer, gwOP
and WOP combine 'pvals'. 'pvals_adj' stands in only for a table that has no
raw column. Each test gives a study table both columns, with adjusted values
that differ plainly from the raw ones, and checks which one was used.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from kosmic.meta_analysis.pooling.fisher import fisher_fast
from kosmic.meta_analysis.pooling.gwop import gwop_fast
from kosmic.meta_analysis.pooling.stouffer import stouffer_fast
from kosmic.meta_analysis.pooling.wop import wop_fast


def _studies(with_raw=True, with_adj=True, k=3, n_genes=40, seed=3):
    rng = np.random.default_rng(seed)
    studies = []
    for s in range(k):
        p = rng.uniform(0.001, 0.9, n_genes)
        df = pd.DataFrame({
            'names': [f'G{i:02d}' for i in range(n_genes)],
            'logfoldchanges': rng.normal(0, 0.5, n_genes),
            'se': rng.uniform(0.15, 0.35, n_genes),
            'dataset': f'S{s + 1}',
        })
        if with_raw:
            df['pvals'] = p
        if with_adj:
            df['pvals_adj'] = np.minimum(p * 20, 1.0)   # plainly not the raw p
        studies.append(df)
    return studies


METHODS = [
    pytest.param(fisher_fast, id='fisher'),
    pytest.param(stouffer_fast, id='stouffer'),
    pytest.param(gwop_fast, id='gwop'),
    pytest.param(wop_fast, id='wop'),
]


def _pooled(fn, studies):
    return fn(studies, min_studies=2).set_index('names')['pvals_pooled']


@pytest.mark.parametrize('fn', METHODS)
def test_both_columns_give_the_raw_result(fn):
    both = _pooled(fn, _studies(True, True))
    raw_only = _pooled(fn, _studies(True, False))
    pd.testing.assert_series_equal(both, raw_only)


@pytest.mark.parametrize('fn', METHODS)
def test_the_adjusted_column_is_not_what_was_pooled(fn):
    both = _pooled(fn, _studies(True, True))
    adj_only = _pooled(fn, _studies(False, True))
    common = both.index.intersection(adj_only.index)
    assert not np.allclose(both[common], adj_only[common])


@pytest.mark.parametrize('fn', METHODS)
def test_a_table_without_raw_p_values_falls_back_to_adjusted(fn):
    out = _pooled(fn, _studies(False, True))
    assert len(out) > 0 and np.isfinite(out).all()
