"""Tests for kosmic.scrna.inspect.gene_group: gene-group expression summary.

Sentinel: known per-gene means / percent-expressing and the per-group
dot-plot matrices, plus missing-gene handling.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy.sparse import csr_matrix

from kosmic.scrna.inspect.gene_group import summarise_gene_group


def _toy_adata(sparse=False):
    # 6 cells x 3 genes; first 3 cells group 'x', last 3 group 'y'.
    X = np.array([
        [1, 0, 1],
        [2, 0, 0],
        [3, 0, 1],
        [0, 4, 0],
        [0, 4, 0],
        [0, 4, 2],
    ], dtype=np.float32)
    obs = pd.DataFrame({'grp': ['x', 'x', 'x', 'y', 'y', 'y']})
    var = pd.DataFrame(index=['A', 'B', 'C'])
    return AnnData(X=csr_matrix(X) if sparse else X, obs=obs, var=var)


@pytest.mark.parametrize('sparse', [False, True])
def test_per_gene_stats_and_missing(sparse):
    a = _toy_adata(sparse)
    s = summarise_gene_group(a, ['A', 'B', 'C', 'ZZZ', 'A'], group_col='grp')

    # Duplicates collapsed, absent gene routed to 'missing'.
    assert s.present == ['A', 'B', 'C']
    assert s.missing == ['ZZZ']

    # Global per-gene means: A=(1+2+3)/6=1.0, B=12/6=2.0, C=4/6≈0.667
    assert s.per_gene.loc['A', 'mean_expr'] == pytest.approx(1.0)
    assert s.per_gene.loc['B', 'mean_expr'] == pytest.approx(2.0)
    # Percent expressing (value > 0): A 3/6, B 3/6, C 3/6
    assert s.per_gene.loc['A', 'pct_expressing'] == pytest.approx(0.5)
    assert s.per_gene.loc['C', 'pct_expressing'] == pytest.approx(0.5)


def test_per_group_matrices_and_peak():
    a = _toy_adata()
    s = summarise_gene_group(a, ['A', 'B', 'C'], group_col='grp')

    assert s.groups == ['x', 'y']
    assert s.mean_matrix.shape == (3, 2)
    # A is only in group x, B only in group y.
    assert s.mean_matrix.loc['A', 'x'] == pytest.approx(2.0)
    assert s.mean_matrix.loc['A', 'y'] == pytest.approx(0.0)
    assert s.mean_matrix.loc['B', 'y'] == pytest.approx(4.0)
    assert s.pct_matrix.loc['B', 'y'] == pytest.approx(1.0)
    # Peak group = highest mean expression.
    assert s.per_gene.loc['A', 'peak_group'] == 'x'
    assert s.per_gene.loc['B', 'peak_group'] == 'y'


def test_no_grouping_leaves_matrices_empty():
    a = _toy_adata()
    s = summarise_gene_group(a, ['A', 'B'])
    assert s.groups == []
    assert s.mean_matrix.empty
    assert 'peak_group' not in s.per_gene.columns


def test_all_genes_missing_returns_empty():
    a = _toy_adata()
    s = summarise_gene_group(a, ['NOPE1', 'NOPE2'], group_col='grp')
    assert s.present == []
    assert s.missing == ['NOPE1', 'NOPE2']
    assert s.per_gene.empty
