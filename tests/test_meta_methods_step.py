"""Methods is a workflow step, and every stage records what it did.

The meta workspace had a "View Methods" button opening a modal dialog,
where DE ends on a Methods step. Worse, what it printed spanned every
study and every step -- including Enrichment and Validation, which are
now their own steps -- so it did not belong to the pooling page at all.

These tests pin the step's place in both modes, and pin the provenance
round-trip: a stage recorded by any step must render in the document,
because the whole point is that the Methods text reports what actually
ran rather than what the settings currently say.
"""
from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from kosmic import provenance  # noqa: E402
from kosmic.meta_analysis.io import (  # noqa: E402
    build_combined_methods,
    meta_provenance_dir,
    record_meta_stage,
)


@pytest.fixture(scope="module")
def app():
    """Hold the QApplication.

    'QApplication.instance() or QApplication([])' as a bare statement
    builds one with no reference, Python collects it immediately, and
    Qt then dies part-way through constructing the next widget -- a
    hard interpreter crash, not a test failure.
    """
    pytest.importorskip("PyQt6")
    from PyQt6.QtWidgets import QApplication
    yield QApplication.instance() or QApplication([])


@pytest.fixture
def project(tmp_path):
    """A project with one provenance-bearing study."""
    study = tmp_path / "STUDY_A" / "processed_data"
    study.mkdir(parents=True)
    provenance.record_stage(study, "STUDY_A", "gene_de", {"method": "deseq2"})
    return tmp_path


# -- the shared recorder ----------------------------------------------

def test_record_meta_stage_writes(project):
    assert record_meta_stage(project, "meta_gene", {"n_studies": 2}) is True


def test_record_meta_stage_without_a_project_is_a_no_op():
    """Failing to write a methods note must never lose a result."""
    assert record_meta_stage(None, "meta_gene", {}) is False
    assert record_meta_stage("", "meta_gene", {}) is False


def test_record_meta_stage_survives_a_bad_path(tmp_path):
    bad = tmp_path / "no" / "such" / "file.txt"
    bad.parent.mkdir(parents=True)
    bad.write_text("not a directory")
    assert record_meta_stage(bad, "meta_gene", {}) is False


def test_meta_stages_keep_only_the_latest(project):
    """A run overwrites its output file, so the record must not grow.

    Appending one entry per Run click left 31 recorded pooling runs
    against 8 result files -- most describing results that had already
    been overwritten.
    """
    for i in range(5):
        record_meta_stage(project, "meta_gene", {"run": i})
    rec = provenance.load(meta_provenance_dir(project))
    gene = [st for st in rec["stages"] if st["stage"] == "meta_gene"]
    assert len(gene) == 1
    assert gene[0]["params"] == {"run": 4}


def test_replacing_one_stage_leaves_the_others(project):
    record_meta_stage(project, "meta_enrichment", {"library": "KEGG_2026"})
    record_meta_stage(project, "meta_gene", {"run": 1})
    record_meta_stage(project, "meta_gene", {"run": 2})
    rec = provenance.load(meta_provenance_dir(project))
    names = sorted(st["stage"] for st in rec["stages"])
    assert names == ["meta_enrichment", "meta_gene"]


def test_study_lineage_still_appends(tmp_path):
    """Each study stage transforms the h5ad, so that chain is real."""
    study = tmp_path / "S" / "processed_data"
    study.mkdir(parents=True)
    for i in range(3):
        provenance.record_stage(study, "S", "qc", {"i": i})
    assert len(provenance.load(study)["stages"]) == 3


# -- every step's stage reaches the document --------------------------

@pytest.mark.parametrize("stage,title", [
    ("meta_gene", "Gene meta-analysis"),
    ("meta_enrichment", "Pathway enrichment of the pooled genes"),
    ("meta_reproducibility", "Cross-dataset reproducibility"),
    ("meta_loo", "Leave-one-study-out validation"),
])
def test_stage_renders_with_a_real_title(project, stage, title):
    """An unknown stage renders as its raw key, which reads as a bug."""
    record_meta_stage(project, stage, {"n": 1})
    assert title in build_combined_methods(project, ["STUDY_A"])


def test_parameters_survive_the_round_trip(project):
    record_meta_stage(project, "meta_enrichment", {
        "method": "gsea", "library": "KEGG_2026", "rank_by": "z_stat"})
    record_meta_stage(project, "meta_loo", {"avg_replication_pct": 62.5})
    doc = build_combined_methods(project, ["STUDY_A"])
    for probe in ("gsea", "KEGG_2026", "z_stat", "62.5"):
        assert probe in doc


def test_study_stages_appear_alongside_meta_stages(project):
    record_meta_stage(project, "meta_gene", {"n_studies": 1})
    doc = build_combined_methods(project, ["STUDY_A"])
    assert "===== Study: STUDY_A =====" in doc
    assert "===== Meta-analysis =====" in doc
    assert "Gene differential expression" in doc


def test_no_provenance_yields_empty(tmp_path):
    assert build_combined_methods(tmp_path, []) == ""


# -- the document is a methods section, not the audit trail -----------

def test_long_values_are_summarised(project):
    """One study's rename list was 47,746 chars on a single line."""
    renames = [[f"OLD{i}", f"NEW{i}", "previous_symbol"] for i in range(900)]
    provenance.record_stage(project / "STUDY_A" / "processed_data", "STUDY_A",
                            "gene_names", {"renames": renames, "renamed": 900})
    doc = build_combined_methods(project, ["STUDY_A"])
    longest = max(len(ln) for ln in doc.splitlines())
    assert longest < 400, longest
    assert "900 total" in doc
    assert "provenance.json" in doc


def test_full_mode_keeps_everything(project):
    renames = [[f"OLD{i}", f"NEW{i}", "previous_symbol"] for i in range(900)]
    provenance.record_stage(project / "STUDY_A" / "processed_data", "STUDY_A",
                            "gene_names", {"renames": renames})
    doc = build_combined_methods(project, ["STUDY_A"], full=True)
    assert "OLD899" in doc


def test_superseded_stages_collapse_to_the_last(project):
    """Re-running DE overwrites its output; both settings were listed."""
    study = project / "STUDY_A" / "processed_data"
    provenance.record_stage(study, "STUDY_A", "gene_de", {"min_cells": 100})
    provenance.record_stage(study, "STUDY_A", "gene_de", {"min_cells": 10})
    doc = build_combined_methods(project, ["STUDY_A"])
    assert "min_cells: 10" in doc
    assert "min_cells: 100" not in doc


def test_a_rerun_is_declared_not_hidden(project):
    study = project / "STUDY_A" / "processed_data"
    for i in range(3):
        provenance.record_stage(study, "STUDY_A", "cluster", {"resolution": i})
    doc = build_combined_methods(project, ["STUDY_A"])
    assert "[last of 3 runs]" in doc


def test_only_the_analysed_studies_appear(project):
    """A methods section describes the studies that went in."""
    other = project / "STUDY_B" / "processed_data"
    other.mkdir(parents=True)
    provenance.record_stage(other, "STUDY_B", "gene_de", {"method": "deseq2"})
    doc = build_combined_methods(project, ["STUDY_A"])
    assert "STUDY_A" in doc
    assert "STUDY_B" not in doc


def test_no_names_still_gives_the_project_view(project):
    other = project / "STUDY_B" / "processed_data"
    other.mkdir(parents=True)
    provenance.record_stage(other, "STUDY_B", "gene_de", {"method": "deseq2"})
    doc = build_combined_methods(project)
    assert "STUDY_A" in doc and "STUDY_B" in doc


# -- the step's place in the flow -------------------------------------

def test_methods_is_the_last_step_of_both_modes(app):
    from kosmic.gui.meta_analysis.workspace import MetaAnalysisWorkspace

    for mode in ("exploratory", "hypothesis"):
        steps = [t for t, _ in MetaAnalysisWorkspace.steps_for_mode(mode)]
        assert steps[-1] == "Methods", mode


@pytest.mark.parametrize("mode,attr", [
    ("exploratory", "_EXPLORATORY_MAP"),
    ("hypothesis", "_HYPOTHESIS_MAP"),
    ("comparison", "_COMPARISON_MAP"),
])
def test_step_list_and_page_map_stay_in_step(app, mode, attr):
    from kosmic.gui.meta_analysis.workspace import MetaAnalysisWorkspace

    steps = MetaAnalysisWorkspace.steps_for_mode(mode)
    assert len(steps) == len(getattr(MetaAnalysisWorkspace, attr))


def test_the_old_button_and_dialog_are_gone(app):
    """Methods access moved to the step; the modal must not linger."""
    from PyQt6.QtWidgets import QPushButton
    from kosmic.gui.meta_analysis.pages.gene_ma_page import GeneMAPage

    page = GeneMAPage()
    assert not hasattr(page, "_show_combined_methods")
    labels = [b.text() for b in page.findChildren(QPushButton)]
    assert "View Methods" not in labels


def test_methods_page_reports_when_there_is_nothing_yet(app, tmp_path):
    from kosmic.gui.meta_analysis.pages.methods_page import MetaMethodsPage

    page = MetaMethodsPage()
    assert not page._copy_btn.isEnabled()
    page.set_project(tmp_path, [])
    assert not page._copy_btn.isEnabled()
    assert "No recorded methods yet" in page._status.text()


def test_methods_page_renders_and_enables_copy(app, project):
    from kosmic.gui.meta_analysis.pages.methods_page import MetaMethodsPage

    record_meta_stage(project, "meta_gene", {"n_studies": 1})
    page = MetaMethodsPage()
    page.set_project(project, ["STUDY_A"])
    assert page._copy_btn.isEnabled()
    assert "Gene meta-analysis" in page._view.toPlainText()


# -- per-selection folders (ADR-007) -----------------------------------

@pytest.mark.parametrize("cell_types, folder", [
    (["Endothelial_Cell", "Endothelial_Cell"], "Endothelial_Cell"),
    ([None, None], "all_cells"),
    ([], "all_cells"),
    (["Endothelial_Cell", "Fibroblast"], "mixed"),
    (["Endothelial_Cell", None], "mixed"),
])
def test_selection_names_its_folder(cell_types, folder):
    from kosmic.meta_analysis.io import meta_selection_folder
    assert meta_selection_folder(cell_types) == folder


def test_output_dir_nests_the_selection(tmp_path):
    from kosmic.paths import meta_output_dir
    assert meta_output_dir(tmp_path) == tmp_path / "meta_analysis"
    assert (meta_output_dir(tmp_path, "Fibroblast")
            == tmp_path / "meta_analysis" / "Fibroblast")


def test_each_cell_type_keeps_its_own_record(project):
    """The overwrite in #60: a second cell type must not replace the first."""
    record_meta_stage(project, "meta_gene", {"n_studies": 5}, selection="Endothelial_Cell")
    record_meta_stage(project, "meta_gene", {"n_studies": 4}, selection="Fibroblast")
    endo = provenance.load(meta_provenance_dir(project, "Endothelial_Cell"))
    fib = provenance.load(meta_provenance_dir(project, "Fibroblast"))
    assert endo["stages"][-1]["params"]["n_studies"] == 5
    assert fib["stages"][-1]["params"]["n_studies"] == 4
    assert (meta_provenance_dir(project, "Fibroblast") / "methods.md").exists()


def test_methods_render_the_selections_record(project):
    record_meta_stage(project, "meta_gene", {"marker": "endo"}, selection="Endothelial_Cell")
    record_meta_stage(project, "meta_gene", {"marker": "fib"}, selection="Fibroblast")
    doc = build_combined_methods(project, ["STUDY_A"], selection="Fibroblast")
    assert "fib" in doc and "endo" not in doc


def test_methods_fall_back_to_the_top_level_record(project):
    """Projects from before ADR-007 keep their record at the top level."""
    record_meta_stage(project, "meta_gene", {"marker": "old"})
    doc = build_combined_methods(project, ["STUDY_A"], selection="Fibroblast")
    assert "old" in doc


def test_mixed_selection_is_warned_about():
    from kosmic.meta_analysis.io import mixed_selection_warning
    assert "mixed" in mixed_selection_warning("mixed")
    assert mixed_selection_warning("Fibroblast") is None
