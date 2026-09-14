"""Drive the running app through every workspace and step on a real project.

Not a test of results -- a test that nothing blows up when a user clicks
through the whole application in order: Project (five steps), scRNA
(eight steps, on a loaded study), DE (both modes, every step), Meta
(discovery mode, every step), Figures. Every Python exception raised
while driving -- including ones Qt would otherwise print and swallow
inside a slot -- is collected and reported at the end, with the step it
happened on.

    python dev/scripts/smoke_drive.py "C:\\dcm raw"

Exit status is the number of steps that raised. Runs on a real display
(offscreen is fine too; nothing here depends on font metrics).
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))  # repo root: dev/scripts/<file>

from PyQt6.QtCore import QEvent, QTimer  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

STUDY = "GSE292067"          # a fully processed study in the project
errors: list[tuple[str, str]] = []
current = "startup"


def _hook(exc_type, exc, tb):
    errors.append((current, "".join(traceback.format_exception(exc_type, exc, tb))))


def _pump(app, seconds):
    end = time.monotonic() + seconds
    while time.monotonic() < end:
        app.processEvents()
        app.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        time.sleep(0.02)


def _wait_for(app, cond, timeout):
    end = time.monotonic() + timeout
    while time.monotonic() < end:
        _pump(app, 0.2)
        if cond():
            return True
    return False


def main() -> int:
    global current
    ap = argparse.ArgumentParser()
    ap.add_argument("project")
    args = ap.parse_args()

    sys.excepthook = _hook
    import main as kosmic_main
    from kosmic.gui.shared.theme import ThemeMode, apply_theme

    app = QApplication(sys.argv[:1])
    app.setStyle("Fusion")
    apply_theme(app, ThemeMode.DARK)
    w = kosmic_main.AppWindow()
    w.resize(1600, 1000)
    w.show()
    w._on_project_open_requested(str(Path(args.project)))
    _pump(app, 1.0)

    def step(name, fn, settle=0.5):
        global current
        current = name
        try:
            fn()
            _pump(app, settle)
        except Exception:
            errors.append((name, traceback.format_exc()))
        print(f"  {'FAIL' if errors and errors[-1][0] == name else 'ok  '} {name}")

    # ---- Project ----------------------------------------------------------
    pw = w._project_workspace
    for i, label in enumerate(("Project", "Studies", "Data", "Shared Atlas", "Review")):
        step(f"project/{label}", lambda i=i: pw.on_sidebar_step(i))
    step("project/select each row", lambda: [pw._select_study_row(a) for a in list(pw._row_kinds)])
    pw._select_study_row(STUDY)

    # ---- scRNA ------------------------------------------------------------
    def open_scrna():
        w._nav_rail.set_active("scrna")
        w._on_mode_selected("scrna")
    step("scrna/open", open_scrna)
    sc = w._scrna_workspace
    step("scrna/wait for load", lambda: _wait_for(app, lambda: sc.current_adata is not None, 90))
    assert sc.current_adata is not None, "study did not load"
    for i, (label, _) in enumerate(sc.WORKFLOW_STEPS):
        step(f"scrna/{label}", lambda i=i: sc.switch_tab(i), settle=1.5)

    # ---- DE ---------------------------------------------------------------
    def open_de():
        w._nav_rail.set_active("de")
        w._on_mode_selected("de")
    step("de/open", open_de)
    de = w._de_workspace
    step("de/wait for load", lambda: _wait_for(app, lambda: de.current_adata is not None, 90))
    for mode in ("scoring", "discovery"):
        step(f"de/mode {mode}", lambda m=mode: de._on_mode_selected(m))
        for i, (label, _) in enumerate(de.steps_for_mode(mode)):
            step(f"de[{mode}]/{label}", lambda i=i: de.switch_tab(i), settle=1.0)

    # ---- Meta -------------------------------------------------------------
    def open_meta():
        w._nav_rail.set_active("meta")
        w._on_mode_selected("meta")
    step("meta/open", open_meta)
    mw = w._meta_workspace
    step("meta/select studies page", lambda: mw.switch_page(0), settle=2.0)
    step("meta/choose mode page", lambda: mw.switch_page(1))
    step("meta/pick discovery", lambda: mw._on_consensus_mode_selected("exploratory"), settle=1.0)
    for i, (label, _) in enumerate(mw.steps_for_mode(mw.consensus_mode)):
        step(f"meta[discovery]/{label}", lambda i=i: mw.switch_page(i), settle=1.0)

    # ---- Figures ----------------------------------------------------------
    def open_figs():
        w._nav_rail.set_active("figures")
        w._on_mode_selected("figures")
    step("figures/open", open_figs, settle=1.5)

    # ---- back to Project, then F1 target resolution on every page ---------
    step("project/return", lambda: w._on_mode_selected("project"))
    step("help/resolve", lambda: w._resolve_help_id())

    QTimer.singleShot(0, app.quit)
    app.exec()

    print()
    if errors:
        print(f"{len(errors)} step(s) raised:")
        for name, tb in errors:
            print(f"--- {name}\n{tb}")
    else:
        print("no exceptions")
    return len(errors)


if __name__ == "__main__":
    sys.exit(main())
