"""Sentinels for candidate-matrix discovery (kosmic.scrna.load.matrices).

Covers the three real depositor layouts this exists for: counts in X alone
(Seurat exports), counts in X *and* a second layer (Broad SCP cardiac), and
normalised X with counts in .raw (CellxGene). Getting the last one wrong hands
log-normalised values to a count model silently, so the verdicts are the point.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from kosmic.scrna.load.matrices import (
    LOG_NORMALISED,
    RAW_COUNTS,
    UNCLEAR,
    best_counts_slot,
    describe_count_matrices,
    promote_matrix,
)

N_CELLS, N_GENES = 40, 12


def _counts(seed=0):
    rng = np.random.default_rng(seed)
    return csr_matrix(rng.integers(0, 30, size=(N_CELLS, N_GENES)).astype(np.float32))


def _lognorm(counts):
    dense = counts.toarray()
    dense = dense / np.maximum(dense.sum(axis=1, keepdims=True), 1) * 1e4
    return csr_matrix(np.log1p(dense).astype(np.float32))


def _obs():
    return pd.DataFrame(index=[f"cell{i}" for i in range(N_CELLS)])


def _var():
    return pd.DataFrame(index=[f"gene{j}" for j in range(N_GENES)])


# --- discovery -------------------------------------------------------------

def test_counts_in_x_only():
    a = ad.AnnData(X=_counts(), obs=_obs(), var=_var())
    described = describe_count_matrices(a)
    assert [d["slot"] for d in described] == ["X"]
    assert described[0]["verdict"] == RAW_COUNTS
    assert described[0]["is_integer"] is True
    assert best_counts_slot(a) == "X"


def test_second_counts_layer_is_listed():
    """Broad SCP layout: CellBender counts in X, CellRanger counts in a layer."""
    a = ad.AnnData(X=_counts(0), obs=_obs(), var=_var())
    a.layers["cellranger_raw"] = _counts(1)
    described = describe_count_matrices(a)
    assert {d["slot"] for d in described} == {"X", "layers:cellranger_raw"}
    assert all(d["verdict"] == RAW_COUNTS for d in described)
    # X already holds counts, so it stays the recommendation
    assert best_counts_slot(a) == "X"


def test_normalised_x_with_raw_counts():
    """CellxGene layout: X normalised, counts in .raw -- the dangerous one."""
    counts = _counts()
    a = ad.AnnData(X=_lognorm(counts), obs=_obs(), var=_var())
    a.raw = ad.AnnData(X=counts, obs=_obs(), var=_var())
    by_slot = {d["slot"]: d for d in describe_count_matrices(a)}
    assert by_slot["X"]["verdict"] == LOG_NORMALISED
    assert by_slot["raw"]["verdict"] == RAW_COUNTS
    assert best_counts_slot(a) == "raw"


def test_log1p_in_uns_forces_log_verdict():
    """Scaled-but-integer-looking data is still log if uns says so."""
    a = ad.AnnData(X=_lognorm(_counts()), obs=_obs(), var=_var())
    a.uns["log1p"] = {"base": None}
    assert describe_count_matrices(a)[0]["verdict"] == LOG_NORMALISED


def test_large_non_integer_is_unclear():
    """Big fractional values are neither counts nor log -- do not guess."""
    rng = np.random.default_rng(3)
    X = csr_matrix((rng.random((N_CELLS, N_GENES)) * 5000).astype(np.float32))
    a = ad.AnnData(X=X, obs=_obs(), var=_var())
    assert describe_count_matrices(a)[0]["verdict"] == UNCLEAR
    assert best_counts_slot(a) is None


# --- promotion -------------------------------------------------------------

def test_promote_raw_puts_counts_in_x():
    counts = _counts()
    a = ad.AnnData(X=_lognorm(counts), obs=_obs(), var=_var())
    a.raw = ad.AnnData(X=counts, obs=_obs(), var=_var())
    a.obs["keepme"] = list(range(N_CELLS))

    promoted = promote_matrix(a, "raw")
    assert np.allclose(promoted.X.toarray(), counts.toarray())
    assert "keepme" in promoted.obs.columns      # obs carried across
    assert describe_count_matrices(promoted)[0]["verdict"] == RAW_COUNTS


def test_promote_layer_replaces_x_and_consumes_layer():
    a = ad.AnnData(X=_counts(0), obs=_obs(), var=_var())
    a.layers["cellranger_raw"] = _counts(1)
    expected = a.layers["cellranger_raw"].toarray()

    promoted = promote_matrix(a, "layers:cellranger_raw")
    assert np.allclose(promoted.X.toarray(), expected)
    assert "cellranger_raw" not in promoted.layers
    assert "cellranger_raw" in a.layers          # original untouched (copy=True)


def test_promote_x_is_a_noop():
    a = ad.AnnData(X=_counts(), obs=_obs(), var=_var())
    assert np.allclose(promote_matrix(a, "X").X.toarray(), a.X.toarray())


def test_promote_clears_log1p_record():
    """After promoting counts, uns must not claim the data is log-transformed."""
    counts = _counts()
    a = ad.AnnData(X=_lognorm(counts), obs=_obs(), var=_var())
    a.raw = ad.AnnData(X=counts, obs=_obs(), var=_var())
    a.uns["log1p"] = {"base": None}
    assert "log1p" not in promote_matrix(a, "raw").uns


@pytest.mark.parametrize("slot", ["raw", "layers:missing", "nonsense"])
def test_promote_rejects_absent_or_bad_slot(slot):
    a = ad.AnnData(X=_counts(), obs=_obs(), var=_var())
    with pytest.raises(ValueError):
        promote_matrix(a, slot)


def test_a_dense_matrix_is_described_not_crashed_on():
    """A dense ndarray has a '.data' attribute too -- a memoryview, not the
    stored values -- so a hasattr check sent it down the sparse branch."""
    import numpy as np

    adata = ad.AnnData(
        X=np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32),
        var=pd.DataFrame(index=["A", "B"]),
        obs=pd.DataFrame(index=["c1", "c2"]))
    described = describe_count_matrices(adata)
    assert [d["slot"] for d in described] == ["X"]
    assert described[0]["verdict"] == RAW_COUNTS
    assert described[0]["max"] == 4.0
