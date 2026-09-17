"""The embedding and clustering steps work on one matrix for the included
cells (kosmic.scrna.cluster.workset), never on a copy of the whole study.
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from kosmic.scrna.cluster.workset import apply_results, working_subset

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


def _study(n=60, g=40, excluded=None, normalised=True, seed=0):
    rng = np.random.default_rng(seed)
    counts = sp.csr_matrix(rng.poisson(2, size=(n, g)).astype(np.float32))
    a = ad.AnnData(X=counts.copy(),
                   obs=pd.DataFrame(index=[f'c{i}' for i in range(n)]),
                   var=pd.DataFrame(index=[f'G{j}' for j in range(g)]))
    a.obs['sample'] = pd.Categorical(['s1', 's2', 's3'] * (n // 3))
    a.obs['condition'] = pd.Categorical(['dcm', 'nf'] * (n // 2))
    if excluded is not None:
        role = np.array(['disease'] * n, dtype=object)
        role[excluded] = 'exclude'
        a.obs['_role'] = pd.Categorical(role)
    if normalised:
        a.layers['counts'] = counts
        import scanpy as sc
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    return a


def test_shared_matrix_when_nothing_excluded_and_only_read():
    a = _study()
    w = working_subset(a)
    assert w.X is a.X                       # no copy
    assert not w.layers and w.raw is None
    assert w.obs is not a.obs and w.var is not a.var


def test_private_matrix_on_request():
    a = _study()
    w = working_subset(a, copy_matrix=True)
    assert w.X is not a.X and (w.X != a.X).nnz == 0


def test_masked_subset_carries_rows_and_slots():
    a = _study(excluded=np.arange(0, 60, 4))
    mask = (a.obs['_role'] != 'exclude').to_numpy()
    a.obsm['X_pca'] = np.arange(60 * 3, dtype=np.float32).reshape(60, 3)
    a.uns['pca'] = {'k': 1}
    graph = sp.random(60, 60, density=0.1, format='csr', random_state=0)
    a.obsp['connectivities'] = graph
    w = working_subset(a, mask, obsm_keys=('X_pca', 'missing'),
                       uns_keys=('pca',), obsp_keys=('connectivities',))
    assert w.n_obs == 45 and list(w.obs_names) == list(a.obs_names[mask])
    assert list(w.obs['_role'].cat.categories) == ['disease']   # unused category gone
    assert (w.X != a.X[mask]).nnz == 0
    assert np.array_equal(w.obsm['X_pca'], a.obsm['X_pca'][mask])
    assert 'missing' not in w.obsm and w.uns['pca'] is a.uns['pca']
    idx = np.flatnonzero(mask)
    assert (w.obsp['connectivities'] != graph[idx][:, idx]).nnz == 0


def test_apply_results_scatters_with_mask_and_assigns_without():
    a = _study(excluded=[0, 1])
    mask = (a.obs['_role'] != 'exclude').to_numpy()
    w = working_subset(a, mask)
    w.obsm['X_umap'] = np.ones((w.n_obs, 2), dtype=np.float32)
    w.obs['leiden'] = pd.Categorical(['0'] * w.n_obs)
    w.var['highly_variable'] = [True] * 10 + [False] * 30
    w.uns['leiden'] = {'params': 1}
    w.obsp['connectivities'] = sp.eye(w.n_obs, format='csr')
    apply_results(a, w, mask, obsm_keys=('X_umap',), obs_cols=('leiden',),
                  var_cols=('highly_variable',), uns_keys=('leiden',),
                  obsp_keys=('connectivities',), drop_uns=('neighbors',))
    assert np.isnan(a.obsm['X_umap'][0]).all() and (a.obsm['X_umap'][2] == 1).all()
    assert pd.isna(a.obs['leiden'].iloc[0]) and a.obs['leiden'].iloc[2] == '0'
    assert a.var['highly_variable'].sum() == 10
    assert 'connectivities' not in a.obsp          # meaningless off the subset
    assert a.uns['leiden'] == {'params': 1}

    b = _study()
    w = working_subset(b)
    w.obsm['X_umap'] = np.ones((60, 2), dtype=np.float32)
    w.obs['leiden'] = pd.Categorical(['1'] * 60)
    w.obsp['connectivities'] = sp.eye(60, format='csr')
    apply_results(b, w, None, obsm_keys=('X_umap',), obs_cols=('leiden',),
                  obsp_keys=('connectivities',))
    assert (b.obsm['X_umap'] == 1).all() and (b.obs['leiden'] == '1').all()
    assert 'connectivities' in b.obsp


@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.mark.parametrize('excluded', [None, list(range(0, 60, 3))])
def test_pca_and_cluster_workers_leave_the_study_matrices_alone(app, tmp_path, excluded):
    """End to end through both workers: results land on the full study,
    X and the counts layer are untouched, no second copy is made."""
    from kosmic.gui.scrna.tabs.cluster_tab import ClusteringWorker, PCAHarmonyWorker
    a = _study(excluded=excluded)
    X_before = a.X.copy()
    counts_before = a.layers['counts'].copy()
    out = tmp_path / 'study.h5ad'

    w = PCAHarmonyWorker(a, n_hvg=20, n_pcs=5, output_path=str(out))
    full, msg = w._run()
    assert full is a
    assert (a.X != X_before).nnz == 0 and (a.layers['counts'] != counts_before).nnz == 0
    assert 'X_pca' in a.obsm and a.obsm['X_pca'].shape == (60, 5)
    assert a.var['highly_variable'].sum() == 20
    if excluded:
        assert np.isnan(a.obsm['X_pca'][0]).all() and not np.isnan(a.obsm['X_pca'][1]).any()

    w = ClusteringWorker(a, {'n_pcs': 5, 'n_neighbors': 5, 'resolution': 0.5}, str(out))
    full, msg = w._run()
    assert full is a
    assert (a.X != X_before).nnz == 0 and (a.layers['counts'] != counts_before).nnz == 0
    assert 'leiden' in a.obs and 'X_umap' in a.obsm and 'rank_genes_groups' in a.uns
    if excluded:
        assert pd.isna(a.obs['leiden'].iloc[0]) and not pd.isna(a.obs['leiden'].iloc[1])
        assert 'connectivities' not in a.obsp
    else:
        assert a.obs['leiden'].notna().all() and 'connectivities' in a.obsp
    saved = ad.read_h5ad(out)
    assert saved.n_obs == 60 and 'leiden' in saved.obs and list(saved.layers) == ['counts']


def test_pca_worker_on_counts_normalises_only_its_private_copy(app):
    from kosmic.gui.scrna.tabs.cluster_tab import PCAHarmonyWorker
    a = _study(normalised=False)           # fresh import: counts in X, no layer
    X_before = a.X.copy()
    full, _ = PCAHarmonyWorker(a, n_hvg=20, n_pcs=5)._run()
    assert full is a and (a.X != X_before).nnz == 0 and not a.layers
    assert 'X_pca' in a.obsm


def test_marker_genes_rank_on_a_per_cluster_subsample():
    """Above the cap the Wilcoxon ranking runs on a random subset per
    cluster; the study itself is untouched and the result says so."""
    from kosmic.scrna.annotate.marker_genes import compute_cluster_marker_genes
    a = _study(n=600, g=40)
    a.obs['leiden'] = pd.Categorical(['0'] * 400 + ['1'] * 150 + ['2'] * 50)
    X = a.X.toarray()
    X[:400, :5] += 3.0                                     # cluster 0 markers G0..G4
    a.X = sp.csr_matrix(X)
    X_before = a.X.copy()
    compute_cluster_marker_genes(a, max_cells_per_cluster=100)
    r = a.uns['rank_genes_groups']
    assert r['params']['max_cells_per_cluster'] == 100 and r['params']['n_cells_used'] == 250
    assert set(pd.DataFrame(r['names'])['0'][:5]) == {'G0', 'G1', 'G2', 'G3', 'G4'}
    assert (a.X != X_before).nnz == 0 and a.n_obs == 600
    compute_cluster_marker_genes(a, max_cells_per_cluster=0)   # every cell
    assert 'n_cells_used' not in a.uns['rank_genes_groups']['params']
