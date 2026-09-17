"""Manual relabels on the Annotate tab take effect and show at once."""
import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from kosmic.scrna.annotate.labels import set_cell_type

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


def test_set_cell_type_adds_a_missing_category():
    obs = pd.DataFrame({'cell_type': pd.Categorical(['CM', 'CM', 'Fib'])})
    n = set_cell_type(obs, [False, True, True], 'Unknown')
    assert n == 2 and list(obs['cell_type']) == ['CM', 'Unknown', 'Unknown']
    assert 'Unknown' in obs['cell_type'].cat.categories and 'CM' in obs['cell_type'].cat.categories


def test_set_cell_type_object_and_absent_columns():
    obs = pd.DataFrame({'cell_type': ['CM', 'Fib']})
    set_cell_type(obs, [True, False], 'Unknown')
    assert list(obs['cell_type']) == ['Unknown', 'Fib']
    obs2 = pd.DataFrame(index=['a', 'b'])
    set_cell_type(obs2, [True, False], 'CM')
    assert list(obs2['cell_type'].astype(str)) == ['CM', 'nan']


@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_right_click_unknown_updates_obs_and_the_table(app):
    from kosmic.gui.scrna.tabs.annotate_tab import AnnotateTab

    class MW:
        current_adata = None
        _adata_version = 0

    t = AnnotateTab(MW())
    n = 120
    a = ad.AnnData(
        X=sp.csr_matrix(np.ones((n, 5), dtype=np.float32)),
        obs=pd.DataFrame({
            'leiden': pd.Categorical(['0'] * 80 + ['1'] * 40),
            'cell_type': pd.Categorical(['Cardiomyocyte'] * 80 + ['Fibroblast'] * 40),
            'cell_type_auto': pd.Categorical(['Cardiomyocyte'] * 80 + ['Fibroblast'] * 40),
        }, index=[f'c{i}' for i in range(n)]))
    t.adata = a
    t._populate_cluster_table()
    tbl = t.cluster_table
    assert tbl.item(1, 3).text() == 'Fibroblast'

    t._apply_cluster_reannotation(1, 'Unknown')       # what "Mark as Unknown" calls
    assert (a.obs.loc[a.obs['leiden'] == '1', 'cell_type'] == 'Unknown').all()
    assert (a.obs.loc[a.obs['leiden'] == '0', 'cell_type'] == 'Cardiomyocyte').all()
    assert tbl.item(1, 3).text() == 'Unknown'          # visible without any other refresh
    assert t.save_annot_btn.isEnabled()
    assert a.uns['cluster_type_source'] == {'1': 'manual'}
