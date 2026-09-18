"""A finished worker must not keep the dataset it was given.

Tabs keep their last worker as an attribute so the QThread stays alive,
and a worker keeps its input AnnData as an attribute. After the result
has been emitted that input is dead weight -- after a subset, the parent
(16 GB on Reichart) stayed in memory because the finished FilterWorker
still pointed at it. BaseWorker releases AnnData inputs once run() has
finished, on success and on failure alike; the emitted result is
untouched.
"""
import gc
import weakref

import anndata as ad
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PyQt6")


@pytest.fixture(scope="module")
def app():
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _adata():
    return ad.AnnData(X=np.zeros((5, 3), dtype=np.float32),
                      obs=pd.DataFrame(index=list("abcde")),
                      var=pd.DataFrame(index=list("xyz")))


def test_inputs_released_after_success(app):
    from kosmic.gui.shared.widgets.base_worker import BaseWorker

    class W(BaseWorker):
        def __init__(self, adata):
            super().__init__()
            self.adata = adata
            self.note = "kept"

        def _run(self):
            return self.adata[:2].copy()

    parent = _adata()
    ref = weakref.ref(parent)
    w = W(parent)
    results = []
    w.finished_ok.connect(results.append)
    w.run()                      # synchronously, no thread needed
    assert results and results[0].n_obs == 2
    assert w.adata is None and w.note == "kept"
    del parent
    gc.collect()
    assert ref() is None, "the finished worker still held the input"


def test_inputs_released_after_failure(app):
    from kosmic.gui.shared.widgets.base_worker import BaseWorker

    class W(BaseWorker):
        def __init__(self, adata):
            super().__init__()
            self.adata = adata

        def _run(self):
            raise ValueError("boom")

    w = W(_adata())
    failures = []
    w.failed.connect(failures.append)
    w.run()
    assert failures and "boom" in failures[0]
    assert w.adata is None
