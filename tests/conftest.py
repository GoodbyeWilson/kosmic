"""Shared test setup.

Widgets a test drops are destroyed by Python's garbage collector at an
arbitrary later moment, sometimes while the next test is still building
its own widgets inside pyqtgraph's C++ layer, which segfaulted the Linux
CI at random tests (interactive_plot.py, PlotWidget.__init__). After
every test the pending deferred deletions are delivered and the garbage
collector is run, so a test's widgets die before the next test starts.
"""
import gc
import os

import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


@pytest.fixture(autouse=True)
def _tear_down_qt_objects():
    yield
    try:
        from PyQt6.QtCore import QEvent
        from PyQt6.QtWidgets import QApplication
    except ImportError:  # pragma: no cover
        return
    app = QApplication.instance()
    gc.collect()
    if app is not None:
        app.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        app.processEvents()
    gc.collect()
