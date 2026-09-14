"""Neighbour-graph reuse.

Building the kNN graph is the most expensive step in a clustering run
(~35s per 40k cells, minutes on a combined master) and it does not depend
on the Leiden resolution. Rebuilding it to try a different resolution was
pure waste, so these pin when it may be skipped -- and, more importantly,
when it must not be.
"""
import numpy as np
import pandas as pd
import anndata as ad

from kosmic.scrna.cluster.neighbors import (
    compute_neighbors,
    neighbors_are_current,
)


def _adata(n=200, n_pcs=10, seed=0):
    rng = np.random.default_rng(seed)
    a = ad.AnnData(
        X=np.zeros((n, 2), dtype=np.float32),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]))
    a.obsm["X_pca"] = rng.normal(0, 1, (n, n_pcs)).astype(np.float32)
    return a


def test_a_fresh_object_has_no_graph():
    assert not neighbors_are_current(_adata(), 15, 10, "X_pca")


def test_matching_parameters_count_as_current():
    a = _adata()
    compute_neighbors(a, n_neighbors=15, n_pcs=10)
    assert neighbors_are_current(a, 15, 10, "X_pca")


def test_a_different_neighbour_count_is_not_current():
    a = _adata()
    compute_neighbors(a, n_neighbors=15, n_pcs=10)
    assert not neighbors_are_current(a, 30, 10, "X_pca")


def test_a_different_pc_count_is_not_current():
    a = _adata()
    compute_neighbors(a, n_neighbors=15, n_pcs=10)
    assert not neighbors_are_current(a, 15, 5, "X_pca")


def test_switching_representation_is_not_current():
    """A graph built before Harmony ran must not be reused after."""
    a = _adata()
    compute_neighbors(a, n_neighbors=15, n_pcs=10)
    assert not neighbors_are_current(a, 15, 10, "X_pca_harmony")


def test_reuse_leaves_the_graph_untouched():
    a = _adata()
    compute_neighbors(a, n_neighbors=15, n_pcs=10)
    before = a.obsp["connectivities"].copy()
    compute_neighbors(a, n_neighbors=15, n_pcs=10)
    assert (a.obsp["connectivities"] != before).nnz == 0


def test_force_rebuilds_even_when_current():
    a = _adata()
    compute_neighbors(a, n_neighbors=15, n_pcs=10)
    a.obsp["connectivities"] *= 0          # sabotage it
    compute_neighbors(a, n_neighbors=15, n_pcs=10, force=True)
    assert a.obsp["connectivities"].nnz > 0


def test_a_graph_recorded_without_obsp_is_not_trusted():
    """uns can survive an operation that dropped the matrices."""
    a = _adata()
    compute_neighbors(a, n_neighbors=15, n_pcs=10)
    del a.obsp["connectivities"]
    assert not neighbors_are_current(a, 15, 10, "X_pca")


def test_harmony_representation_is_preferred_when_present():
    a = _adata()
    a.obsm["X_pca_harmony"] = a.obsm["X_pca"].copy()
    compute_neighbors(a, n_neighbors=15, n_pcs=10)
    assert a.uns["neighbors"]["params"]["use_rep"] == "X_pca_harmony"
