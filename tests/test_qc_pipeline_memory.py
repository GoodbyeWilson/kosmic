"""run_qc_pipeline filters in place, once, and makes the counts layer last.

The result must be what the earlier view-and-copy version produced:
the same cells kept for every threshold, and the counts layer equal to
the surviving cells' original counts.
"""
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from kosmic.scrna.qc.filter import cell_filter_mask, run_qc_pipeline


def _study(n=300, g=60, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.poisson(3, size=(n, g)).astype(np.float32)
    X[:20, :] = 0                                   # empty-ish cells
    X[:20, :3] = 1                                  # 3 genes each
    names = [f'MT-G{j}' if j < 5 else f'G{j}' for j in range(g)]
    X[280:, :5] *= 40                               # high-MT cells at the end
    a = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(index=[f'c{i}' for i in range(n)]),
                   var=pd.DataFrame(index=names))
    return a


def _reference(a, p):
    """The previous implementation, step for step."""
    import scanpy as sc
    from kosmic.scrna.qc.filter import detect_mitochondrial_genes
    a = a.copy()
    detect_mitochondrial_genes(a, 'human')
    sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    a.layers['counts'] = a.X.copy()
    if p['min_genes'] > 0:
        sc.pp.filter_cells(a, min_genes=p['min_genes'])
    if p['max_genes'] > 0:
        sc.pp.filter_cells(a, max_genes=p['max_genes'])
    if p['min_counts'] > 0:
        sc.pp.filter_cells(a, min_counts=p['min_counts'])
    if p['max_counts'] > 0:
        a = a[a.obs['total_counts'] <= p['max_counts'], :].copy()
    sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    if p['max_mt'] > 0:
        a = a[a.obs['pct_counts_mt'] < p['max_mt'], :].copy()
    return a


def test_same_cells_and_counts_as_the_view_and_copy_version():
    p = {'min_genes': 10, 'max_genes': 0, 'min_counts': 50, 'max_counts': 400, 'max_mt': 30, 'min_cells': 0}
    a = _study()
    ref = _reference(a, p)
    out, stats = run_qc_pipeline(a, p)
    assert out is a                                   # in place, same object
    assert list(out.obs_names) == list(ref.obs_names)
    assert (out.layers['counts'] != ref.layers['counts']).nnz == 0
    assert stats['n_cells_after'] == ref.n_obs and stats['cells_removed'] == 300 - ref.n_obs
    assert 20 <= stats['cells_removed'] < 300


def test_nothing_removed_means_no_subset_and_a_counts_layer():
    a = _study()
    X_id = id(a.X)
    out, stats = run_qc_pipeline(a, {'min_genes': 1, 'max_genes': 0, 'min_counts': 0, 'max_counts': 0, 'max_mt': 0})
    assert out is a and id(a.X) == X_id and stats['cells_removed'] == 0
    assert 'counts' in a.layers and (a.layers['counts'] != a.X).nnz == 0


def test_cell_filter_mask_thresholds():
    obs = pd.DataFrame({'n_genes_by_counts': [5, 50, 500], 'total_counts': [10, 100, 1000],
                        'pct_counts_mt': [1.0, 20.0, 5.0]})
    assert cell_filter_mask(obs).tolist() == [True, True, True]
    assert cell_filter_mask(obs, min_genes=10).tolist() == [False, True, True]
    assert cell_filter_mask(obs, max_genes=100).tolist() == [True, True, False]
    assert cell_filter_mask(obs, min_counts=50, max_counts=500).tolist() == [False, True, False]
    assert cell_filter_mask(obs, max_mt=10).tolist() == [True, False, True]
    assert cell_filter_mask(obs.drop(columns=['pct_counts_mt']), max_mt=10).tolist() == [True, True, True]
