"""Capture the screenshots the help pages embed, from the running app.

Screenshots taken by hand go stale the moment the interface changes.
This script launches KOSMIC on a real display (fonts render properly --
offscreen rendering draws text as boxes), opens a project, walks to each
documented state and grabs the main window. Re-run it after any UI
change and every image is current.

    python dev/scripts/capture_help_screenshots.py "C:\\dcm raw"
    python dev/scripts/capture_help_screenshots.py "C:\\dcm raw" --only project

Images land in kosmic/gui/help/content/<section>/img/<name>.png and are
referenced from the Markdown by those fixed names. The project you pass
should look like a real one -- several studies at different stages --
because the pages describe what the reader will see.

Each capture is a (help section, image name, setup) triple; 'setup' is
given the AppWindow and puts it in the state to photograph. Add a triple
when a page needs a new picture.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))  # repo root: dev/scripts/<file>

from PyQt6.QtCore import QEvent, QTimer  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

CONTENT = Path(__file__).resolve().parents[2] / "kosmic" / "gui" / "help" / "content"
WINDOW_SIZE = (1600, 1000)
SETTLE_MS = 400     # let layouts and background scans finish before grabbing


# ---------------------------------------------------------------------------
# States to photograph
# ---------------------------------------------------------------------------

def _project_step(window, step: int, select: str | None = None):
    window._nav_rail.set_active("project")
    window._on_mode_selected("project")
    ws = window._project_workspace
    if select:
        ws._select_study_row(select)
    ws.on_sidebar_step(step)


def _scrna_tab(window, study: str, tab: int):
    """Open scRNA on 'study' at step 'tab'. The first call on a study
    starts the h5ad load; later tabs wait for it (see CAPTURES)."""
    window.set_active_study(study)
    window._nav_rail.set_active("scrna")
    window._on_mode_selected("scrna")
    window._scrna_workspace.switch_tab(tab)


def _scrna_load(window, study: str):
    _scrna_tab(window, study, 0)


def _de_setup(window, study: str):
    window.set_active_study(study)
    window._nav_rail.set_active("de")
    window._on_mode_selected("de")
    window._de_workspace.switch_tab(0)


CAPTURES = [
    # section, image name, setup(window), extra settle seconds
    ("project", "step1_project",  lambda w: _project_step(w, 0), 0),
    ("project", "step2_studies",  lambda w: _project_step(w, 1), 0),
    ("project", "step3_data",     lambda w: _project_step(w, 2), 0),
    ("project", "step4_atlas",    lambda w: _project_step(w, 3), 2),   # shared-gene count
    ("project", "step5_review",   lambda w: _project_step(w, 4), 0),
    ("project", "details_subset", lambda w: _project_step(w, 4, "_master_ECs"), 0),
    # scRNA: the first entry pays for the h5ad load; the rest are quick.
    ("scrna",   "load_data",    lambda w: _scrna_tab(w, "GSE292067", 0), 20),
    ("scrna",   "gene_names",   lambda w: _scrna_tab(w, "GSE292067", 1), 4),
    ("scrna",   "inspect",      lambda w: _scrna_tab(w, "GSE292067", 2), 2),
    ("scrna",   "qc",           lambda w: _scrna_tab(w, "GSE292067", 3), 4),
    ("scrna",   "cluster",      lambda w: _scrna_tab(w, "GSE292067", 4), 6),
    ("scrna",   "marker_check", lambda w: _scrna_tab(w, "GSE292067", 5), 4),
    ("scrna",   "decontx",      lambda w: _scrna_tab(w, "GSE292067", 6), 1),
    ("scrna",   "subset",       lambda w: _scrna_tab(w, "GSE292067", 7), 2),
    ("de",      "setup",          lambda w: _de_setup(w, "GSE292067"), 15),
]


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("project", help="project folder to open")
    ap.add_argument("--only", help="capture one section only (e.g. project)")
    args = ap.parse_args()
    project = Path(args.project)
    if not project.is_dir():
        ap.error(f"not a folder: {project}")

    import main as kosmic_main  # the app module, imported for AppWindow
    from kosmic.gui.shared.theme import ThemeMode, apply_theme

    app = QApplication(sys.argv[:1])
    app.setStyle("Fusion")
    apply_theme(app, ThemeMode.DARK)
    window = kosmic_main.AppWindow()
    window.resize(*WINDOW_SIZE)
    window.show()
    window._on_project_open_requested(str(project))
    _pump(app, 1.0)

    todo = [c for c in CAPTURES if not args.only or c[0] == args.only]
    for section, name, setup, extra in todo:
        setup(window)
        _pump(app, SETTLE_MS / 1000 + extra)
        out = CONTENT / section / "img"
        out.mkdir(parents=True, exist_ok=True)
        path = out / f"{name}.png"
        window.grab().save(str(path))
        print(f"wrote {path.relative_to(CONTENT.parent)}")

    QTimer.singleShot(0, app.quit)
    app.exec()
    return 0


def _pump(app, seconds: float) -> None:
    """Run the event loop for a while. Also flushes deferred deletes,
    which processEvents() alone leaves for exec() -- without this a
    widget replaced via deleteLater() stays painted in the grab."""
    end = time.monotonic() + seconds
    while time.monotonic() < end:
        app.processEvents()
        app.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        time.sleep(0.02)


if __name__ == "__main__":
    sys.exit(main())
