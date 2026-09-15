"""The import's second test for raw counts: row sums against the per-cell
totals the depositor recorded.

Integer values alone do not prove a matrix is raw counts -- SCTransform's
corrected counts are integers with every cell rescaled to a common depth,
and GSE292067 was imported that way. The three signatures pinned here
are the ones measured on the DCM project's files: raw counts equal the
recorded totals (Koenig, Youness), a gene-filtered or ambient-corrected
copy sits at or below them for every cell (Reichart via CELLxGENE,
Chaffin's CellBender X), and SCT-corrected counts scatter above and
below (Guo).
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from kosmic.scrna.load.matrices import (
    DEPTH_BELOW, DEPTH_MATCHES, DEPTH_RESCALED, RAW_COUNTS, UNCLEAR,
    best_counts_slot, compare_with_recorded_totals, describe_count_matrices,
    describe_depth,
)

N_CELLS, N_GENES = 200, 50


def _counts(seed=0):
    rng = np.random.default_rng(seed)
    # Wide depth range so a rescaling is visible.
    depth = rng.integers(500, 8000, size=N_CELLS)
    probs = rng.dirichlet(np.ones(N_GENES) * 0.3)
    return np.vstack([rng.multinomial(d, probs) for d in depth]).astype(np.float32)


def _adata(X, total=None, column="nCount_RNA"):
    obs = pd.DataFrame(index=[f"c{i}" for i in range(X.shape[0])])
    if total is not None:
        obs[column] = total
    var = pd.DataFrame(index=[f"g{j}" for j in range(X.shape[1])])
    return ad.AnnData(X=csr_matrix(X), obs=obs, var=var)


def _sct_like(X, seed=1):
    """Every cell rescaled to the median depth, then rounded -- the shape of
    SCTransform's corrected counts."""
    rng = np.random.default_rng(seed)
    depth = X.sum(axis=1, keepdims=True)
    scaled = X / depth * np.median(depth)
    return np.round(scaled + rng.uniform(0, 0.4, size=X.shape)).astype(np.float32)


# --- compare_with_recorded_totals ---------------------------------------

def test_no_recorded_total_returns_none():
    a = _adata(_counts())
    assert compare_with_recorded_totals(a.X, a.obs) is None
    assert "no per-cell total" in describe_depth(None)


def test_raw_counts_match_recorded_totals():
    X = _counts()
    a = _adata(X, X.sum(axis=1))
    d = compare_with_recorded_totals(a.X, a.obs)
    assert d["verdict"] == DEPTH_MATCHES
    assert d["column"] == "nCount_RNA"
    assert d["frac_equal"] == 1.0
    assert "equal obs['nCount_RNA']" in describe_depth(d)


def test_gene_subset_is_below_not_rescaled():
    """CELLxGENE drops genes; every row sum is then <= the recorded total."""
    X = _counts()
    a = _adata(X[:, :45], X.sum(axis=1))
    d = compare_with_recorded_totals(a.X, a.obs)
    assert d["verdict"] == DEPTH_BELOW
    assert d["frac_above"] == 0.0
    assert d["ratio_median"] < 1.0
    assert "at or below" in describe_depth(d)


def test_ambient_correction_is_below_not_rescaled():
    """CellBender removes counts from every cell, never adds."""
    rng = np.random.default_rng(3)
    X = _counts()
    corrected = np.floor(X * rng.uniform(0.85, 0.95, size=(N_CELLS, 1)))
    a = _adata(corrected, X.sum(axis=1))
    assert compare_with_recorded_totals(a.X, a.obs)["verdict"] == DEPTH_BELOW


def test_sct_corrected_counts_are_rescaled():
    X = _counts()
    a = _adata(_sct_like(X), X.sum(axis=1))
    d = compare_with_recorded_totals(a.X, a.obs)
    assert d["verdict"] == DEPTH_RESCALED
    assert d["frac_above"] > 0.1 and d["frac_equal"] < 0.05
    assert d["ratio_q01"] < 0.9 < 1.1 < d["ratio_q99"]
    assert "rescaled" in describe_depth(d)


def test_other_total_columns_are_recognised():
    X = _counts()
    for column in ("n_counts", "total_counts", "nUMI", "cellbender_ncount"):
        a = _adata(X, X.sum(axis=1), column=column)
        assert compare_with_recorded_totals(a.X, a.obs)["column"] == column


def test_explicit_column_wins():
    X = _counts()
    a = _adata(X, X.sum(axis=1))
    a.obs["n_counts"] = 1.0
    assert compare_with_recorded_totals(a.X, a.obs, column="nCount_RNA")["verdict"] == DEPTH_MATCHES


# --- describe_count_matrices demotes rescaled integer matrices ------------

def test_rescaled_integer_x_is_not_reported_as_raw_counts():
    X = _counts()
    a = _adata(_sct_like(X), X.sum(axis=1))
    described = describe_count_matrices(a)
    assert described[0]["is_integer"] is True
    assert described[0]["verdict"] == UNCLEAR
    assert described[0]["depth"]["verdict"] == DEPTH_RESCALED
    assert best_counts_slot(a) is None


def test_raw_counts_elsewhere_are_still_found():
    """SCT in X, raw RNA counts in a layer: the layer should win."""
    X = _counts()
    a = _adata(_sct_like(X), X.sum(axis=1))
    a.layers["counts"] = csr_matrix(X)
    by_slot = {d["slot"]: d for d in describe_count_matrices(a)}
    assert by_slot["X"]["verdict"] == UNCLEAR
    assert by_slot["layers:counts"]["verdict"] == RAW_COUNTS
    assert by_slot["layers:counts"]["depth"]["verdict"] == DEPTH_MATCHES
    assert best_counts_slot(a) == "layers:counts"


def test_without_recorded_totals_integer_still_counts_as_raw():
    a = _adata(_sct_like(_counts()))
    described = describe_count_matrices(a)
    assert described[0]["verdict"] == RAW_COUNTS
    assert described[0]["depth"] is None


# --- the import worker refuses rescaled counts -----------------------------

def test_import_worker_refuses_rescaled_counts():
    pytest.importorskip("PyQt6")
    from kosmic.gui.intake.workers import _check_recorded_totals

    X = _counts()
    log = []
    _check_recorded_totals(_adata(X, X.sum(axis=1)), log.append)
    assert log and log[0].startswith("Depth check:")

    bad = _adata(_sct_like(X), X.sum(axis=1))
    with pytest.raises(RuntimeError, match="rescaled counts"):
        _check_recorded_totals(bad, log.append)

    # A matrix the user chose by hand is warned about, not refused.
    log.clear()
    _check_recorded_totals(bad, log.append, user_chose=True)
    assert log and log[0].startswith("WARNING")
