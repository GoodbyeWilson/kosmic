"""The Gene DE Run buttons follow the dataset, not only page activation:
a study that finishes loading while the page is open enables them."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pytest
from anndata import AnnData


@pytest.fixture(scope="module")
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_run_buttons_enable_when_the_dataset_finishes_loading(app):
    from kosmic.gui.de_analysis.workspace import DEWorkspace
    ws = DEWorkspace()
    page = ws.gene_de_page
    page.on_activated()                      # opened before the data arrived
    assert not page.run_btn.isEnabled()
    ws.current_adata = AnnData(np.ones((4, 3)))
    ws.dataset_loaded.emit("study.h5ad")     # the background load finishing
    assert page.run_btn.isEnabled() and page.batch_btn.isEnabled()
    ws.release_dataset()
    assert not page.run_btn.isEnabled()
