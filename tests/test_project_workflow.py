"""Project workflow steps (ADR-002): completion rules and the page that
renders them."""
from __future__ import annotations

import json
import os

import pytest

from kosmic.paths import MASTER_ACCESSION, scaffold_study
from kosmic.project_workflow import (
    STEP_ATLAS, STEP_DATA, STEP_PROJECT, STEP_REVIEW, STEP_STUDIES, STEPS,
    project_progress, study_kind,
)


def _row(kind="study", raw=False, processed=False, de=False, pb=False):
    return {"kind": kind, "status": {"raw_data": raw, "processed": processed,
                                     "de_results": de, "pseudobulk": pb}}


# ---------------------------------------------------------------------------
# Completion rules
# ---------------------------------------------------------------------------

def test_five_steps_declared():
    assert len(STEPS) == 5
    assert STEPS[STEP_ATLAS][0].endswith("(optional)")


def test_no_project_completes_nothing():
    p = project_progress(None)
    assert p.complete == (False,) * 5
    assert p.status[1] == "warning"


def test_open_empty_project_completes_only_step_one():
    p = project_progress({})
    assert p.complete[STEP_PROJECT]
    assert not p.complete[STEP_STUDIES]
    assert not p.complete[STEP_DATA]      # nothing to have data for
    assert "add at least 2" in p.status[0]


def test_two_studies_without_data():
    p = project_progress({"A": _row(), "B": _row()})
    assert p.complete[STEP_STUDIES]
    assert not p.complete[STEP_DATA]
    assert p.missing_data == ["A", "B"]
    assert p.status == ("2 studies without data", "warning")


def test_data_step_counts_real_studies_only():
    # A subset and the atlas always have processed data; they must not
    # mask a study that lacks it, nor count toward the study floor.
    rows = {"A": _row(processed=True), "B": _row(),
            "A_ECs": _row("subset", processed=True),
            MASTER_ACCESSION: _row("atlas", processed=True)}
    p = project_progress(rows)
    assert p.n_studies == 2
    assert p.missing_data == ["B"]
    assert p.complete[STEP_ATLAS]
    assert not p.complete[STEP_DATA]


def test_review_counts_anything_with_de_and_pseudobulk_except_atlas():
    rows = {"A": _row(processed=True, de=True, pb=True),
            "B": _row(processed=True, de=True),        # no pseudobulk
            "A_ECs": _row("subset", processed=True, de=True, pb=True),
            MASTER_ACCESSION: _row("atlas", processed=True, de=True, pb=True)}
    p = project_progress(rows)
    assert p.n_ready == 2
    assert p.complete[STEP_REVIEW]
    assert p.status == ("2 of 4 ready for meta-analysis", "success")


def test_one_ready_is_not_enough_to_pool():
    rows = {"A": _row(processed=True, de=True, pb=True),
            "B": _row(processed=True)}
    p = project_progress(rows)
    assert not p.complete[STEP_REVIEW]
    assert p.status == ("1 of 2 ready for meta-analysis", "info")


# ---------------------------------------------------------------------------
# Folder kinds come from provenance, not names
# ---------------------------------------------------------------------------

def _write_provenance(study_dir, source=None):
    pd = study_dir / "processed_data"
    pd.mkdir(parents=True, exist_ok=True)
    rec = {"study": study_dir.name, "stages": []}
    if source:
        rec["source"] = source
    (pd / "provenance.json").write_text(json.dumps(rec), encoding="utf-8")


def test_study_kind(tmp_path):
    for name in ("GSE1", "GSE1_ECs", "weird_name", MASTER_ACCESSION):
        scaffold_study(tmp_path / name)
    _write_provenance(tmp_path / "GSE1")
    _write_provenance(
        tmp_path / "GSE1_ECs",
        source=str(tmp_path / "GSE1" / "processed_data" / "GSE1.h5ad"))
    # Underscore in the name is not evidence of a subset.
    assert study_kind(tmp_path / "GSE1") == ("study", None)
    assert study_kind(tmp_path / "weird_name") == ("study", None)
    assert study_kind(tmp_path / "GSE1_ECs") == ("subset", "GSE1")
    assert study_kind(tmp_path / MASTER_ACCESSION) == ("atlas", None)


# ---------------------------------------------------------------------------
# The page
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def qapp():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture
def project(tmp_path):
    for name in ("GSE1", "GSE2"):
        scaffold_study(tmp_path / name)
        _write_provenance(tmp_path / name)
    (tmp_path / "GSE1" / "raw_data" / "x.h5ad").write_bytes(b"")
    return tmp_path


def test_workspace_declares_workflow_contract(qapp):
    from kosmic.gui.project.workspace import ProjectWorkspace
    ws = ProjectWorkspace()
    assert ws.WORKFLOW_STEPS == STEPS
    assert list(ws.iter_completed_steps()) == [(i, False) for i in range(5)]
    assert ws.compute_sidebar_status()[1] == "warning"
    assert ws.current_step == STEP_PROJECT
    ws.deleteLater()


def test_show_project_jumps_to_first_incomplete_step(qapp, project):
    from kosmic.gui.project.workspace import ProjectWorkspace
    ws = ProjectWorkspace()
    seen = []
    ws.step_changed.connect(seen.append)
    ws.show_project(project, None)
    # Two studies exist, neither processed -> Data is the first gap.
    assert ws.current_step == STEP_DATA
    assert seen[-1] == STEP_DATA
    assert ws._side_stack.currentIndex() == STEP_DATA
    done = dict(ws.iter_completed_steps())
    assert done[STEP_PROJECT] and done[STEP_STUDIES] and not done[STEP_DATA]
    assert ws._content_stack.currentIndex() == ws._CONTENT_STUDIES
    ws.deleteLater()


def test_refresh_keeps_the_users_step(qapp, project):
    from kosmic.gui.project.workspace import ProjectWorkspace
    ws = ProjectWorkspace()
    ws.show_project(project, None)
    ws.on_sidebar_step(STEP_REVIEW)
    ws.show_project(project, "GSE1")       # same project: a refresh
    assert ws.current_step == STEP_REVIEW
    ws.deleteLater()


def test_step_click_without_project_stays_on_project_step(qapp):
    from kosmic.gui.project.workspace import ProjectWorkspace
    ws = ProjectWorkspace()
    ws.on_sidebar_step(STEP_REVIEW)
    assert ws.current_step == STEP_PROJECT
    assert ws._content_stack.currentIndex() == ws._CONTENT_WELCOME
    ws.deleteLater()


def test_data_panel_follows_the_selected_study(qapp, project):
    from kosmic.gui.project.workspace import ProjectWorkspace
    ws = ProjectWorkspace()
    ws.show_project(project, None)
    # No active study, but the first row is selected: that is enough.
    assert ws._import_btn.isEnabled()
    assert ws._data_target_label.text() == "Selected study: GSE1"
    assert ws._data_missing_label.text() == "Without data: GSE1, GSE2"
    got = []
    ws.import_requested.connect(lambda a, r: got.append((a, r)))
    ws._select_study_row("GSE2")
    assert ws._data_target_label.text() == "Selected study: GSE2"
    ws._emit_import("h5ad")
    assert got == [("GSE2", "h5ad")]
    ws.deleteLater()


def test_show_empty_returns_to_welcome(qapp, project):
    from kosmic.gui.project.workspace import ProjectWorkspace
    ws = ProjectWorkspace()
    ws.show_project(project, None)
    ws.show_empty()
    assert ws._content_stack.currentIndex() == ws._CONTENT_WELCOME
    assert ws.current_step == STEP_PROJECT
    assert all(not d for _, d in ws.iter_completed_steps())
    ws.deleteLater()


def test_recent_projects_render_in_both_lists(qapp, tmp_path):
    from kosmic.gui.project.workspace import ProjectWorkspace
    ws = ProjectWorkspace()
    ws.set_recent_projects([str(tmp_path), str(tmp_path / "missing")])
    assert ws._recents_layout.count() == 1
    assert ws._side_recents_layout.count() == 1
    ws.deleteLater()
