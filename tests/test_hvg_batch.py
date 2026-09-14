"""Batch-aware highly-variable-gene selection.

Pooling every cell before ranking variance lets genes that differ between
studies score as "variable", so the gene set carries the batch effect into
PCA -- before Harmony gets a chance to remove it. Selecting within batch
and ranking by how many batches call a gene variable is the usual guard.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from kosmic.scrna.cluster.hvg import find_hvg


def _adata(n_cells=300, n_genes=60, n_batches=3, seed=0):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(4, (n_cells, n_genes)).astype(np.float32)
    batch = np.repeat([f"B{i}" for i in range(n_batches)],
                      n_cells // n_batches)
    # A few genes are wildly different in one batch and flat elsewhere:
    # pooled they look highly variable, per-batch only one batch calls them.
    counts[batch == "B0", :3] *= 40
    # find_hvg documents log-transformed input; flavor='seurat' bins on
    # mean expression and cannot bin raw counts.
    X = np.log1p(counts / counts.sum(1, keepdims=True) * 1e4).astype(np.float32)
    return ad.AnnData(
        X=X,
        var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]),
        obs=pd.DataFrame({"batch": batch, "one_level": "x"},
                         index=[f"c{i}" for i in range(n_cells)]))


def test_batch_key_changes_the_selection():
    pooled = _adata()
    per_batch = _adata()
    find_hvg(pooled, n_hvg=10)
    find_hvg(per_batch, n_hvg=10, batch_key="batch")
    assert (set(pooled.var_names[pooled.var["highly_variable"]])
            != set(per_batch.var_names[per_batch.var["highly_variable"]]))


def test_batch_key_reports_how_many_batches_call_each_gene():
    a = _adata()
    find_hvg(a, n_hvg=10, batch_key="batch")
    assert "highly_variable_nbatches" in a.var.columns


def test_a_single_level_key_is_ignored_not_an_error():
    """scanpy raises on a one-level batch key; a study with one sample
    must still get HVGs."""
    a = _adata()
    find_hvg(a, n_hvg=10, batch_key="one_level")
    assert int(a.var["highly_variable"].sum()) > 0


def test_a_missing_column_is_ignored():
    a = _adata()
    find_hvg(a, n_hvg=10, batch_key="not_a_column")
    assert int(a.var["highly_variable"].sum()) > 0


def test_a_per_cell_unique_key_is_ignored():
    """A key with one level per cell would make every gene its own batch."""
    a = _adata()
    a.obs["cell_id"] = list(a.obs_names)
    find_hvg(a, n_hvg=10, batch_key="cell_id")
    assert int(a.var["highly_variable"].sum()) > 0


def test_no_batch_key_is_unchanged():
    a, b = _adata(), _adata()
    find_hvg(a, n_hvg=10)
    find_hvg(b, n_hvg=10, batch_key=None)
    assert list(a.var["highly_variable"]) == list(b.var["highly_variable"])


@pytest.mark.parametrize("n_hvg", [5, 20])
def test_raw_survives_the_call(n_hvg):
    a = _adata()
    a.raw = a
    find_hvg(a, n_hvg=n_hvg, batch_key="batch")
    assert a.raw is not None


# --- flavour selection ------------------------------------------------------

def _counts(n_cells=300, n_genes=60, seed=1):
    rng = np.random.default_rng(seed)
    X = rng.poisson(4, (n_cells, n_genes)).astype(np.float32)
    return ad.AnnData(
        X=X,
        var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]),
        obs=pd.DataFrame({"batch": np.repeat(["B0", "B1", "B2"], n_cells // 3)},
                         index=[f"c{i}" for i in range(n_cells)]))


def test_auto_picks_seurat_v3_on_counts():
    """seurat_v3 models the mean-variance curve of counts; it is the right
    selector when X has not been normalised yet."""
    a = _counts()
    find_hvg(a, n_hvg=10)
    assert a.uns["hvg"]["flavor"] == "seurat_v3"


def test_auto_falls_back_to_seurat_on_log_data():
    """seurat_v3 gives nonsense on log values, so the flavour has to
    follow the matrix rather than be fixed."""
    a = _adata()          # already log-normalised
    find_hvg(a, n_hvg=10)
    assert a.uns["hvg"]["flavor"] == "seurat"


def test_an_explicit_flavour_is_respected():
    a = _counts()
    find_hvg(a, n_hvg=10, flavor="seurat")
    assert a.uns["hvg"]["flavor"] == "seurat"


def test_seurat_v3_combines_with_batch_key():
    a = _counts()
    find_hvg(a, n_hvg=10, batch_key="batch")
    assert a.uns["hvg"]["flavor"] == "seurat_v3"
    assert "highly_variable_nbatches" in a.var.columns
