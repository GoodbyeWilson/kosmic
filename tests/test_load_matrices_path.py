"""File-level candidate discovery and targeted reads.

The in-memory functions are fine for a loaded object, but the import path has
to choose a matrix *before* loading one: a CellxGene h5ad can hold the same
1.35 billion non-zeros as both a normalised X and raw counts, and reading both
to pick one runs out of memory. These sentinels cover the file-level route.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from kosmic.scrna.load.matrices import (
    LOG_NORMALISED,
    RAW_COUNTS,
    describe_count_matrices,
    describe_count_matrices_path,
    read_counts,
)

N_CELLS, N_GENES = 30, 10


def _counts(seed=0):
    rng = np.random.default_rng(seed)
    return csr_matrix(rng.integers(0, 40, size=(N_CELLS, N_GENES)).astype(np.float32))


def _lognorm(counts):
    dense = counts.toarray()
    dense = dense / np.maximum(dense.sum(axis=1, keepdims=True), 1) * 1e4
    return csr_matrix(np.log1p(dense).astype(np.float32))


def _frame(n, prefix):
    return pd.DataFrame(index=[f"{prefix}{i}" for i in range(n)])


@pytest.fixture
def cellxgene_like(tmp_path):
    """Normalised X, raw counts in .raw -- the layout that broke the import."""
    counts = _counts()
    a = ad.AnnData(X=_lognorm(counts), obs=_frame(N_CELLS, "cell"),
                   var=_frame(N_GENES, "gene"))
    a.obs["donor"] = ["d%d" % (i % 3) for i in range(N_CELLS)]
    a.raw = ad.AnnData(X=counts, obs=_frame(N_CELLS, "cell"),
                       var=_frame(N_GENES, "gene"))
    path = tmp_path / "cellxgene.h5ad"
    a.write_h5ad(path)
    return path, counts


def test_path_description_matches_in_memory(cellxgene_like):
    path, _ = cellxgene_like
    from_file = {d["slot"]: d["verdict"] for d in describe_count_matrices_path(path)}
    in_memory = {d["slot"]: d["verdict"] for d in describe_count_matrices(ad.read_h5ad(path))}
    assert from_file == in_memory
    assert from_file["X"] == LOG_NORMALISED
    assert from_file["raw"] == RAW_COUNTS


def test_path_description_sees_layers(tmp_path):
    a = ad.AnnData(X=_counts(0), obs=_frame(N_CELLS, "cell"),
                   var=_frame(N_GENES, "gene"))
    a.layers["cellranger_raw"] = _counts(1)
    path = tmp_path / "layered.h5ad"
    a.write_h5ad(path)
    slots = {d["slot"] for d in describe_count_matrices_path(path)}
    assert slots == {"X", "layers:cellranger_raw"}


def test_read_counts_returns_the_raw_matrix(cellxgene_like):
    path, counts = cellxgene_like
    got = read_counts(path, "raw")
    assert np.allclose(got.X.toarray(), counts.toarray())
    assert "donor" in got.obs.columns          # obs comes from the top level
    assert describe_count_matrices(got)[0]["verdict"] == RAW_COUNTS


def test_read_counts_x_and_layer(tmp_path):
    a = ad.AnnData(X=_counts(0), obs=_frame(N_CELLS, "cell"),
                   var=_frame(N_GENES, "gene"))
    a.layers["cellranger_raw"] = _counts(1)
    path = tmp_path / "layered.h5ad"
    a.write_h5ad(path)

    assert np.allclose(read_counts(path, "X").X.toarray(), a.X.toarray())
    assert np.allclose(read_counts(path, "layers:cellranger_raw").X.toarray(),
                       a.layers["cellranger_raw"].toarray())


@pytest.mark.parametrize("slot", ["raw", "layers:missing", "nonsense"])
def test_read_counts_rejects_absent_slot(tmp_path, slot):
    a = ad.AnnData(X=_counts(), obs=_frame(N_CELLS, "cell"),
                   var=_frame(N_GENES, "gene"))
    path = tmp_path / "plain.h5ad"
    a.write_h5ad(path)
    with pytest.raises(ValueError):
        read_counts(path, slot)
