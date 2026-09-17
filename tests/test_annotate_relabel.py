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
    t.main_window.current_adata = a
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


def test_save_refuses_a_stale_dataset_and_writes_the_current_one(app, tmp_path, monkeypatch):
    """After a study switch the tab can hold the previous study while the
    workspace path names the new one. That pair must never be written."""
    from kosmic.gui.shared import dialogs
    from kosmic.gui.scrna.tabs.annotate_tab import AnnotateTab

    def _study(n, tag):
        return ad.AnnData(
            X=sp.csr_matrix(np.ones((n, 3), dtype=np.float32)),
            obs=pd.DataFrame({'cell_type': pd.Categorical([tag] * n)},
                             index=[f'{tag}{i}' for i in range(n)]))

    warnings, infos = [], []
    monkeypatch.setattr(dialogs, 'warning', lambda *a, **k: warnings.append(a[2]))
    monkeypatch.setattr(dialogs, 'info', lambda *a, **k: infos.append(a[2]))

    class MW:
        current_adata = None
        current_h5ad_path = None
        _adata_version = 0

    ws = MW()
    t = AnnotateTab(ws)
    previous = _study(5, 'youness')
    current = _study(8, 'reichart')
    target = tmp_path / 'reichart.h5ad'
    current.write_h5ad(target)

    # the workspace has moved on; the tab has not refreshed
    ws.current_adata, ws.current_h5ad_path = current, str(target)
    t.adata = previous
    t._save_annotations()
    assert warnings and 'wrong study' in warnings[-1]
    assert ad.read_h5ad(target).n_obs == 8            # untouched

    # in sync: written, atomically, to the workspace path
    t.adata = current
    current.obs['cell_type'] = pd.Categorical(['Unknown'] * 8)
    t._save_annotations()
    assert infos and not (tmp_path / 'reichart.h5ad.tmp').exists()
    assert (ad.read_h5ad(target).obs['cell_type'] == 'Unknown').all()


def test_relabel_acts_on_the_workspace_dataset_not_a_stale_one(app):
    """The study was switched; the tab still holds the old object. A
    relabel must land on the study the workspace has open."""
    from kosmic.gui.scrna.tabs.annotate_tab import AnnotateTab

    class MW:
        current_adata = None
        current_h5ad_path = None
        _adata_version = 3

    def _study(tag, n=20):
        return ad.AnnData(
            X=sp.csr_matrix(np.ones((n, 3), dtype=np.float32)),
            obs=pd.DataFrame({'leiden': pd.Categorical(['0'] * n),
                              'cell_type': pd.Categorical([tag] * n)},
                             index=[f'{tag}{i}' for i in range(n)]))

    ws = MW()
    t = AnnotateTab(ws)
    old, new = _study('old'), _study('new')
    ws.current_adata = new
    t.adata = old
    t._populate_cluster_table()
    t._apply_cluster_reannotation(0, 'Unknown')
    assert (old.obs['cell_type'] == 'old').all()          # untouched
    assert (new.obs['cell_type'] == 'Unknown').all()
    assert t.adata is new
