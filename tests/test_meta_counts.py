"""Stage 4 of ADR-015: a meta-analysis pools results made on one counts
source, and its outputs go to meta_analysis/<counts>/<selection>/.
"""
from __future__ import annotations

import json
import os

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from kosmic.meta_analysis.io import (  # noqa: E402
    UNSPECIFIED_COUNTS, MixedCountsError, meta_selection_path,
    mixed_selection_warning, selection_counts, selection_label, split_selection,
)
from kosmic.paths import de_run_dir, meta_output_dir  # noqa: E402

CELL = "Ventricular_Cardiomyocyte"


@pytest.fixture(scope="module")
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


# -- the Qt-free rules ----------------------------------------------------

def test_one_counts_source_names_the_folder():
    rows = [{"name": "A", "counts": "decontX_counts"},
            {"name": "B", "counts": "decontX_counts"}]
    assert selection_counts(rows) == "decontX_counts"


def test_results_on_different_counts_are_refused_naming_them():
    rows = [{"name": "Chaffin_VCM", "counts": "original_counts"},
            {"name": "Guo_VCM", "counts": "decontX_counts"},
            {"name": "Koenig_VCM", "counts": "decontX_counts"}]
    with pytest.raises(MixedCountsError) as err:
        selection_counts(rows)
    text = str(err.value)
    assert "original counts: Chaffin_VCM" in text
    assert "decontX counts: Guo_VCM, Koenig_VCM" in text


@pytest.mark.parametrize("rows", [
    [{"name": "A", "counts": "original_counts"}, {"name": "X", "external": True}],
    [{"name": "A", "counts": "original_counts"}, {"name": "B", "counts": None}],
    [],
])
def test_unknown_counts_go_to_unspecified_counts(rows):
    assert selection_counts(rows) == UNSPECIFIED_COUNTS


def test_an_external_table_does_not_count_as_mixed():
    rows = [{"name": "A", "counts": "decontX_counts"},
            {"name": "X", "external": True, "counts": "original_counts"}]
    assert selection_counts(rows) == UNSPECIFIED_COUNTS


def test_selection_paths_and_labels():
    sel = meta_selection_path("decontX_counts", CELL)
    assert sel == f"decontX_counts/{CELL}"
    assert split_selection(sel) == ("decontX_counts", CELL)
    assert split_selection(CELL) == (None, CELL)          # before ADR-015
    assert selection_label(sel) == "Ventricular Cardiomyocyte, decontX counts"
    assert selection_label(CELL) == "Ventricular Cardiomyocyte"
    assert meta_output_dir("P", sel) == meta_output_dir("P") / "decontX_counts" / CELL


def test_the_mixed_warning_names_the_counts_folder():
    assert "meta_analysis/original_counts/mixed/" in mixed_selection_warning(
        "original_counts/mixed")
    assert "meta_analysis/mixed/" in mixed_selection_warning("mixed")
    assert mixed_selection_warning(f"original_counts/{CELL}") is None


# -- the Select Studies page -------------------------------------------------

def _write_run(project, study, counts_source):
    folder = de_run_dir(project / study, counts_source, CELL)
    folder.mkdir(parents=True, exist_ok=True)
    genes = [f"G{i}" for i in range(20)]
    lfc = 0.5 if counts_source == "raw" else -0.5
    pd.DataFrame({"names": genes, "logfoldchanges": lfc, "se": 0.2,
                  "pvals": 0.01, "pvals_adj": 0.05}).to_csv(
        folder / f"{study}_{CELL}_DE_deseq2.csv", index=False)


@pytest.fixture
def two_counts_project(tmp_path):
    """Three studies with runs on original and DecontX counts; Youness only
    on the original counts."""
    for study in ("Chaffin", "Guo", "Koenig"):
        _write_run(tmp_path, study, "raw")
        _write_run(tmp_path, study, "decontX_counts")
    _write_run(tmp_path, "Youness", "raw")
    return tmp_path


def _scanned(project):
    from kosmic.gui.meta_analysis.workspace import MetaAnalysisWorkspace
    ws = MetaAnalysisWorkspace()
    ws._project_folder = str(project)
    ws._rescan_project()
    return ws


def _choose_cell_type(ws):
    combo = ws._cell_type_combo
    combo.setCurrentIndex(combo.findData(CELL))


def test_the_counts_box_offers_each_source_original_first(app, two_counts_project):
    ws = _scanned(two_counts_project)
    combo = ws._counts_combo
    items = [(combo.itemData(i), combo.itemText(i)) for i in range(combo.count())]
    assert items == [("original_counts", "Original counts  (4 studies)"),
                     ("decontX_counts", "DecontX counts  (3 studies)")]
    assert combo.currentData() == "original_counts"
    assert combo.isEnabled()


def test_the_pool_holds_results_on_the_chosen_counts_only(app, two_counts_project):
    ws = _scanned(two_counts_project)
    _choose_cell_type(ws)
    assert {d["counts"] for d in ws.datasets} == {"original_counts"}
    assert all(d["df"]["logfoldchanges"].eq(0.5).all() for d in ws.datasets)
    assert ws.output_selection == f"original_counts/{CELL}"
    assert f"meta_analysis/original_counts/{CELL}/" in ws._select_summary.text()

    ws._counts_combo.setCurrentIndex(ws._counts_combo.findData("decontX_counts"))
    _choose_cell_type(ws)
    assert sorted(d["study"] for d in ws.datasets) == ["Chaffin", "Guo", "Koenig"]
    assert all("decontX_counts" in d["path"] for d in ws.datasets)
    assert all(d["df"]["logfoldchanges"].eq(-0.5).all() for d in ws.datasets)
    assert ws.output_selection == f"decontX_counts/{CELL}"


def test_a_study_without_the_chosen_counts_is_named(app, two_counts_project):
    ws = _scanned(two_counts_project)
    logged = []
    ws._log = logged.append
    ws._counts_combo.setCurrentIndex(ws._counts_combo.findData("decontX_counts"))
    assert any("No DE results on decontX counts in: Youness" in m for m in logged)


def test_an_imported_table_moves_the_selection_to_unspecified_counts(app, two_counts_project):
    ws = _scanned(two_counts_project)
    _choose_cell_type(ws)
    ws.datasets.append({"name": "Paper_table", "external": True,
                        "df": ws.datasets[0]["df"]})
    assert ws.output_selection == f"unspecified_counts/{CELL}"


def test_the_gene_run_records_its_counts(app, two_counts_project):
    from kosmic.gui.meta_analysis.pages.gene_ma_page import GeneMAPage
    ws = _scanned(two_counts_project)
    ws._counts_combo.setCurrentIndex(ws._counts_combo.findData("decontX_counts"))
    _choose_cell_type(ws)
    page = GeneMAPage()
    page._project_folder = str(two_counts_project)
    page.set_datasets(ws.datasets, ws.dataset_labels, str(two_counts_project),
                      output_selection=ws.output_selection)
    page._run_output_selection = ws.output_selection
    page._method_keys = ["reml"]
    page._meta_df = pd.DataFrame({"names": ["G1"], "fdr": [0.01]})
    page._auto_save()
    folder = meta_output_dir(two_counts_project, f"decontX_counts/{CELL}")
    sidecar = json.loads(next(folder.glob("consensus_*.json")).read_text())
    assert sidecar["parameters"]["counts"] == "decontX_counts"


# -- saved runs (Figure Export) ----------------------------------------------

def _save_meta_run(project, selection, names, counts=...):
    folder = meta_output_dir(project, selection)
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"names": ["G1"], "logfoldchanges": [0.1], "se": [0.1],
                  "fdr": [0.5]}).to_csv(folder / "consensus_reml_analytical.csv",
                                        index=False)
    params = {"de_method": "deseq2", "study_names": names}
    if counts is not ...:
        params["counts"] = counts
    (folder / "consensus_reml_analytical.json").write_text(
        json.dumps({"parameters": params}), encoding="utf-8")


def test_saved_runs_are_listed_in_both_layouts(two_counts_project):
    from kosmic.meta_analysis.saved import list_meta_runs
    names = [f"{s}_{CELL}" for s in ("Chaffin", "Guo", "Koenig")]
    _save_meta_run(two_counts_project, f"original_counts/{CELL}", names, "original_counts")
    _save_meta_run(two_counts_project, f"decontX_counts/{CELL}", names, "decontX_counts")
    _save_meta_run(two_counts_project, "Endothelial_Cell", names)       # before ADR-015
    runs = {r.selection: r for r in list_meta_runs(two_counts_project)}
    assert set(runs) == {f"original_counts/{CELL}", f"decontX_counts/{CELL}",
                         "Endothelial_Cell"}
    assert runs[f"decontX_counts/{CELL}"].label.startswith(
        "Ventricular Cardiomyocyte, decontX counts: reml")
    assert runs["Endothelial_Cell"].label.startswith("Endothelial Cell: reml")


def test_a_saved_run_draws_the_study_tables_on_its_counts(two_counts_project):
    from kosmic.meta_analysis.saved import figure_inputs, list_meta_runs
    names = [f"{s}_{CELL}" for s in ("Chaffin", "Guo", "Koenig")]
    _save_meta_run(two_counts_project, f"decontX_counts/{CELL}", names, "decontX_counts")
    _save_meta_run(two_counts_project, f"original_counts/{CELL}", names, "original_counts")
    for run in list_meta_runs(two_counts_project):
        lfc = -0.5 if run.selection.startswith("decontX") else 0.5
        studies = figure_inputs(run)["gene_studies"]
        assert [n for n, _ in studies] == names
        assert all(df["logfoldchanges"].eq(lfc).all() for _, df in studies), run.selection


def test_the_counts_choice_does_not_carry_into_another_project(app, two_counts_project, tmp_path_factory):
    ws = _scanned(two_counts_project)
    ws._counts_combo.setCurrentIndex(ws._counts_combo.findData("decontX_counts"))
    ws._rescan_project()
    assert ws._counts_combo.currentData() == "decontX_counts"     # same project: kept
    other = tmp_path_factory.mktemp("other")
    for study in ("A", "B", "C"):
        _write_run(other, study, "raw")
        _write_run(other, study, "decontX_counts")
    ws._project_folder = str(other)
    ws._rescan_project()
    assert ws._counts_combo.currentData() == "original_counts"


def test_the_counts_most_results_were_made_on_are_chosen_first(app, tmp_path):
    """An older project: whole-study runs record their counts, the
    per-cell-type runs saved before records do not."""
    for study in ("A", "B", "C"):
        _write_run(tmp_path, study, "raw")                      # one on original counts
        stats = tmp_path / study / "results" / "de_analysis" / "statistics"
        for ct in ("Endothelial_Cell", "Fibroblast"):
            pd.DataFrame({"names": ["G1"], "logfoldchanges": 0.1, "se": 0.1, "pvals": 0.1,
                          "pvals_adj": 0.1}).to_csv(stats / f"{study}_{ct}_DE_deseq2.csv",
                                                    index=False)
    ws = _scanned(tmp_path)
    assert ws._counts_combo.currentData() == ""                 # counts not recorded
    assert len(ws._all_gene_datasets) == 6
