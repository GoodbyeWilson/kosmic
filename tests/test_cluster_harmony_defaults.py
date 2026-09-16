"""Harmony defaults on the Cluster step: on when a batch column exists,
correcting on the sample; Merge across (the second variable) never
pre-selected. Pre-selecting the condition there pulled disease and
control cells into shared clusters without the user asking.
"""
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


def test_merge_across_blank_and_harmony_on(app):
    from kosmic.gui.scrna.tabs.cluster_tab import ClusterTab

    obs = pd.DataFrame({
        'sample': [f's{i % 6}' for i in range(60)],
        'condition': ['DCM' if i % 2 else 'Donor' for i in range(60)],
    }, index=[f'c{i}' for i in range(60)])
    a = ad.AnnData(X=np.zeros((60, 5), dtype=np.float32), obs=obs,
                   var=pd.DataFrame(index=list('abcde')))
    tab = ClusterTab(None)
    tab.adata = a
    tab._populate_harmony_keys()
    assert tab.batch_combo.currentText() == 'sample'
    assert tab.condition_combo.currentText() == ''
    assert tab.harmony_check.isChecked()
