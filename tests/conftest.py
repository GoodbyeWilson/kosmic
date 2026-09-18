"""Shared test setup.

Widgets a test drops are destroyed whenever Python's garbage collector
next runs, which can be in the middle of a later test constructing its
own widgets inside pyqtgraph's C++ layer; the Linux CI segfaulted there
at random tests. Running the collector at the end of every test makes
that destruction happen at a test boundary instead. Nothing is deleted
by force: Qt's own deferred deletion is left to the event loop, so
module-scoped fixtures other tests still hold are untouched.
"""
import gc
import os

import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


@pytest.fixture(autouse=True)
def _collect_between_tests():
    yield
    gc.collect()
