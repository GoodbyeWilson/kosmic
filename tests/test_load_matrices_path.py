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


def test_median_counts_per_cell_tell_corrected_from_uncorrected(tmp_path):
    """A Broad portal file holds CellBender counts in X and Cell Ranger
    counts in a layer; the corrected one has fewer counts per cell."""
    raw = _counts(0)
    corrected = raw.copy()
    corrected.data = np.floor(corrected.data * 0.9)        # background removed
    corrected.eliminate_zeros()
    a = ad.AnnData(X=corrected, obs=_frame(N_CELLS, "cell"), var=_frame(N_GENES, "gene"))
    a.layers["cellranger_raw"] = raw
    path = tmp_path / "portal.h5ad"
    a.write_h5ad(path)
    got = {d["slot"]: d["median_total"] for d in describe_count_matrices_path(path)}
    assert got["layers:cellranger_raw"] == float(np.median(np.asarray(raw.sum(axis=1)).ravel()))
    assert got["X"] == float(np.median(np.asarray(corrected.sum(axis=1)).ravel()))
    assert got["X"] < got["layers:cellranger_raw"]


def test_median_counts_count_empty_cells_as_zero(tmp_path):
    counts = _counts(0).toarray()
    counts[:20] = 0                                         # most cells empty
    a = ad.AnnData(X=csr_matrix(counts), obs=_frame(N_CELLS, "cell"),
                   var=_frame(N_GENES, "gene"))
    path = tmp_path / "sparse_rows.h5ad"
    a.write_h5ad(path)
    assert describe_count_matrices_path(path)[0]["median_total"] == 0.0


def test_the_import_records_the_file_and_matrix_it_used(tmp_path):
    from kosmic.gui.intake.workers import H5adImportWorker

    a = ad.AnnData(X=_counts(0), obs=_frame(N_CELLS, "cell"), var=_frame(N_GENES, "gene"))
    a.layers["cellranger_raw"] = _counts(1)
    src = tmp_path / "portal.h5ad"
    a.write_h5ad(src)
    out_dir = tmp_path / "processed_data"
    out, _ = H5adImportWorker(str(src), str(out_dir), slot="layers:cellranger_raw")._run()
    got = ad.read_h5ad(out)
    assert (got.X.toarray() == a.layers["cellranger_raw"].toarray()).all()
    assert dict(got.uns["kosmic_import"]) == {
        "source_file": "portal.h5ad", "counts_from": "layers['cellranger_raw']"}
    # ...and the import is the study's 'load' stage, with the cells as imported.
    from kosmic import provenance
    load = provenance.last_stage(out_dir, "load")["params"]
    assert load["counts_from"] == "layers['cellranger_raw']"
    assert load["source_file"] == "portal.h5ad" and load["n_cells"] == N_CELLS


def _portal_file(tmp_path):
    a = ad.AnnData(X=_counts(0), obs=_frame(N_CELLS, "cell"), var=_frame(N_GENES, "gene"))
    a.layers["cellranger_raw"] = _counts(1)
    path = tmp_path / "portal.h5ad"
    a.write_h5ad(path)
    return path


def test_both_imports_ask_which_counts_when_there_are_two(tmp_path, monkeypatch):
    """The Project page and Load Data share one question (counts_choice)."""
    from PyQt6.QtWidgets import QInputDialog
    from kosmic.gui.intake.counts_choice import choose_counts_matrix

    asked = {}

    def get_item(parent, title, text, items, current, editable):
        asked.update(text=text, items=items)
        return items[1], True

    monkeypatch.setattr(QInputDialog, 'getItem', staticmethod(get_item))
    slot, described = choose_counts_matrix(None, _portal_file(tmp_path))
    assert slot == 'layers:cellranger_raw'
    assert len(asked['items']) == 2 and 'counts per cell' in asked['text']
    assert {d['slot'] for d in described} == {'X', 'layers:cellranger_raw'}

    monkeypatch.setattr(QInputDialog, 'getItem', staticmethod(lambda *a: ('', False)))
    assert choose_counts_matrix(None, _portal_file(tmp_path))[0] is False


def test_no_question_when_one_matrix_holds_counts(cellxgene_like, monkeypatch):
    from PyQt6.QtWidgets import QInputDialog
    from kosmic.gui.intake.counts_choice import choose_counts_matrix

    def get_item(*a):
        raise AssertionError("should not ask")

    monkeypatch.setattr(QInputDialog, 'getItem', staticmethod(get_item))
    path, _ = cellxgene_like
    assert choose_counts_matrix(None, path)[0] is None


def test_the_load_stage_records_the_matrices_considered(tmp_path):
    from kosmic import provenance
    from kosmic.gui.intake.counts_choice import candidates_summary
    from kosmic.gui.intake.workers import H5adImportWorker

    src = _portal_file(tmp_path)
    out_dir = tmp_path / "processed_data"
    candidates = candidates_summary(describe_count_matrices_path(src))
    H5adImportWorker(str(src), str(out_dir), slot="layers:cellranger_raw",
                     candidates=candidates)._run()
    load = provenance.last_stage(out_dir, "load")["params"]
    assert [c["slot"] for c in load["candidates"]] == ["X", "layers:cellranger_raw"]
    assert all("median_counts_per_cell" in c for c in load["candidates"])
    assert not (tmp_path / "provenance.json").exists()   # nothing outside processed_data
