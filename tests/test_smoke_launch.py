"""Launch smoke test.

Asserts that the KOSMIC entry point can be imported and `AppWindow` can be
constructed without raising. A passing test means `python main.py` can
still start the app; a failing test means a recent change broke the
constructor chain before we notice it interactively.

The AppWindow constructor is exercised in a subprocess with a hard timeout.
Running it inside the pytest process works but leaves Qt timers / event
loops alive and hangs teardown, so subprocess isolation matches how the
app actually launches (a fresh Python process) and gives us a clean exit.

Runs headless via ``QT_QPA_PLATFORM=offscreen``. Does not enter the Qt
event loop and does not trigger network calls. All workspaces are
built eagerly during AppWindow construction (so the first nav-rail
click is just a stack-page swap).

Target: < 20 seconds.
"""
from __future__ import annotations

import os
import subprocess
import sys
import textwrap

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_main_module_importable():
    """`import main` must succeed without side effects."""
    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)
    import importlib
    mod = importlib.import_module("main")
    assert hasattr(mod, "AppWindow"), "main.AppWindow is missing"
    assert hasattr(mod, "main"), "main.main() entry point is missing"


_SMOKE_SCRIPT = textwrap.dedent(
    """
    import os, sys
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    sys.path.insert(0, {root!r})
    from PyQt6.QtWidgets import QApplication
    app = QApplication(sys.argv[:1])
    import main

    window = main.AppWindow()
    assert window._nav_rail is not None
    assert window._explorer is not None
    assert window._stack is not None
    assert window._output_panel is not None

    # Every workspace is built eagerly during AppWindow construction so
    # the first nav-rail click is instant. Catches regression to the old
    # lazy-build pattern that caused per-tab hangs.
    assert window._project_workspace is not None, "Project not built"
    assert window._scrna_workspace is not None, "scRNA not built"
    assert window._de_workspace is not None, "DE not built"
    assert window._meta_workspace is not None, "Meta not built"
    assert window._figures_workspace is not None, "Figures not built"

    window.close()
    # Flush explicitly before the forced exit -- stdout is captured by the
    # parent via pipes and os._exit skips the usual atexit flush.
    sys.stdout.write('SMOKE_OK\\n')
    sys.stdout.flush()
    # Force exit: Qt timers from splash/tutorial can keep Python alive
    # past the final line if we rely on normal shutdown.
    os._exit(0)
    """
)


def test_app_window_constructs_in_subprocess():
    """Spawn a fresh Python process that builds AppWindow and exits cleanly."""
    script = _SMOKE_SCRIPT.format(root=_ROOT)
    try:
        result = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True, text=True, timeout=60,
            cwd=_ROOT,
        )
    except subprocess.TimeoutExpired as exc:
        raise AssertionError(
            f"AppWindow construction timed out after 60s.\n"
            f"stdout: {exc.stdout!r}\nstderr: {exc.stderr!r}"
        )
    assert result.returncode == 0, (
        f"Smoke subprocess exited {result.returncode}.\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    assert "SMOKE_OK" in result.stdout, (
        f"Expected SMOKE_OK sentinel, got:\nstdout: {result.stdout!r}\n"
        f"stderr: {result.stderr!r}"
    )
