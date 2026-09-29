"""Nothing Qt builds during the suite is destroyed while the suite runs.

The Linux runner segfaulted part-way through the suite, always inside
pyqtgraph's menu construction -- ``PlotItem``'s plot-options ``QMenu``
or a ``ViewBoxMenu`` -- while a test built a fresh plot widget. It
happened under the offscreen platform and again under xcb on Xvfb, on
Python 3.11 but not 3.12, with the same PyQt6 and pyqtgraph in both.

The cause is lifetime, not the platform. A test builds a parentless
widget, holds it in a local, and drops it when the test ends. Qt
widgets sit in reference cycles (parent/child links, signal
connections), so they are not freed at that point: the cycle collector
frees them whenever it next runs, which is typically in the middle of
some later test's allocation -- that is, in the middle of building
another plot. Qt then deletes a tree of C++ objects underneath
pyqtgraph while pyqtgraph is building one.

So the fixture below keeps a reference to every widget that exists at
the end of each test. Nothing is destroyed until the process exits and
the operating system reclaims it, which for a test run is free.
Retaining the top-level widgets is enough to retain what they own in
C++, including the ``QGraphicsItem`` trees that are not widgets.

This does not stop a test from checking that something *else* is
released: ``release_dataset`` clears the tabs' references to a study,
and the study is still collected with the workspace widget retained.
"""
import pytest

# Module level, so the references live as long as the interpreter.
_LIVE_WIDGETS: list = []


@pytest.fixture(autouse=True)
def keep_qt_widgets_alive():
    """Retain every widget alive at the end of the test (see module docstring)."""
    yield
    try:
        from PyQt6.QtWidgets import QApplication
    except ImportError:          # a run without PyQt6 installed
        return
    if QApplication.instance() is not None:
        _LIVE_WIDGETS.extend(QApplication.allWidgets())
