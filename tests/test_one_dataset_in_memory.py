"""One dataset in memory at a time, and DE loads only what it needs.

A study is tens of gigabytes. Two workspaces each holding a different
one, or a tab keeping the previous study alive after a switch, is what
emptied a 64 GB machine when the DCM atlas was opened for its DE.
"""
import gc
import os
import weakref

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from kosmic.scrna.load.h5ad_meta import read_counts_adata

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


def _study(path, n=30, g=12, with_layer=True):
    rng = np.random.default_rng(0)
    counts = sp.csr_matrix(rng.poisson(2, size=(n, g)).astype(np.float32))
    a = ad.AnnData(X=counts.copy(), obs=pd.DataFrame({'sample': ['s1', 's2', 's3'] * (n // 3),
                                                      '_role': pd.Categorical(['disease', 'control'] * (n // 2))},
                                                     index=[f'c{i}' for i in range(n)]),
                   var=pd.DataFrame(index=[f'G{j}' for j in range(g)]))
    if with_layer:
        a.layers['counts'] = counts
        a.X = sp.csr_matrix(np.log1p(counts.toarray()))
        a.obsm['X_umap'] = np.zeros((n, 2), dtype=np.float32)
        a.obsp['connectivities'] = sp.eye(n, format='csr')
    a.uns['role_map'] = {'x': 'disease'}
    a.write_h5ad(path)
    return counts


def test_read_counts_adata_takes_the_layer_and_nothing_else(tmp_path):
    counts = _study(tmp_path / 's.h5ad')
    a = read_counts_adata(tmp_path / 's.h5ad')
    assert (a.X != counts).nnz == 0                       # counts, not the normalised X
    assert not a.layers and not a.obsm and not a.obsp
    assert a.uns['role_map'] == {'x': 'disease'} and a.uns['counts_loaded_from'] == 'layers:counts'
    assert list(a.obs.columns) == ['sample', '_role'] and a.n_vars == 12


def test_read_counts_adata_falls_back_to_x(tmp_path):
    counts = _study(tmp_path / 'raw.h5ad', with_layer=False)
    a = read_counts_adata(tmp_path / 'raw.h5ad')
    assert (a.X != counts).nnz == 0 and a.uns['counts_loaded_from'] == 'X'


@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_release_dataset_frees_the_study_from_every_tab(app, tmp_path):
    from kosmic.gui.scrna.workspace import ScRNAWorkspace
    ws = ScRNAWorkspace(embedded=True)
    _study(tmp_path / 's.h5ad')
    a = ad.read_h5ad(tmp_path / 's.h5ad')
    ref = weakref.ref(a)
    ws.set_adata(a, str(tmp_path / 's.h5ad'))
    # tabs pick the dataset up as they would when activated
    ws.qc_tab.set_data(a, str(tmp_path / 's.h5ad'))
    ws.cluster_tab.set_data(a, str(tmp_path / 's.h5ad'))
    ws.cluster_tab._umap_widget.set_adata(a)
    ws.annotate_tab.set_adata(a)
    del a
    gc.collect()
    assert ref() is not None                                # held by the tabs, as before

    assert ws.release_dataset() is True
    assert ws.current_adata is None and ws.qc_tab.adata is None
    assert ws.cluster_tab.adata is None and ws.annotate_tab.adata is None
    assert ws.cluster_tab._umap_widget._adata is None
    gc.collect()
    assert ref() is None                                    # actually gone


def test_de_workspace_release(app):
    from kosmic.gui.de_analysis.workspace import DEWorkspace
    ws = DEWorkspace()
    ws.current_adata = ad.AnnData(X=np.ones((3, 2), dtype=np.float32))
    ws.h5ad_path = 'x.h5ad'
    ws._manual_load = True
    assert ws.release_dataset() is True
    assert ws.current_adata is None and ws.h5ad_path is None and ws._manual_load is False
    ws.close()
    ws.deleteLater()
    app.processEvents()


def test_de_setup_releases_before_it_loads(app, tmp_path, monkeypatch):
    """The previous dataset (ours and other workspaces') is released before
    the new file is read, so two studies are never resident together."""
    from PyQt6.QtCore import QEvent
    from kosmic.gui.de_analysis import pages
    from kosmic.gui.de_analysis.workspace import DEWorkspace
    _study(tmp_path / 'new.h5ad')
    ws = DEWorkspace()
    ws.current_adata = ad.AnnData(X=np.ones((3, 2), dtype=np.float32))
    ws.h5ad_path = 'old.h5ad'
    ws._manual_load = True
    order = []
    ws.dataset_loading.connect(lambda p: order.append(('loading', ws.current_adata is None)))
    started = {}
    monkeypatch.setattr(pages.setup_page, 'run_worker', lambda w, **kw: started.setdefault('worker', w) or w)
    ws.setup_page._start_load(str(tmp_path / 'new.h5ad'))
    assert order == [('loading', False)]          # signal fires while the old one is still held...
    assert ws.current_adata is None                 # ...then it is released, before the read starts
    assert 'worker' in started
    ws.close()
    ws.deleteLater()
    app.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    app.processEvents()


def test_scrna_auto_load_waits_until_the_workspace_is_on_screen(app, tmp_path, monkeypatch):
    """Activating a study while another workspace is in front must not load
    the full study behind it (the atlas sat at 27 GB behind DE's 14 GB copy).
    The load happens when the scRNA workspace is opened."""
    from PyQt6.QtCore import QEvent
    from kosmic.gui.scrna.workspace import ScRNAWorkspace
    from kosmic.paths import processed_data_dir
    study = tmp_path / 'study'
    processed_data_dir(study).mkdir(parents=True)
    _study(processed_data_dir(study) / 's.h5ad')
    ws = ScRNAWorkspace(embedded=True)
    loads = []
    monkeypatch.setattr(ws.download_tab, '_load_and_set_adata', lambda p: loads.append(p))

    ws.set_project_directory(str(study))                    # hidden: activated from elsewhere
    assert loads == []
    assert 'loads when this workspace is opened' in ws.download_tab._tab_status_label.text()

    ws.show()                                               # what _activate_scrna does
    ws.download_tab.on_tab_activated()
    assert [os.path.basename(p) for p in loads] == ['s.h5ad']
    ws.close()
    ws.deleteLater()
    app.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    app.processEvents()


def _fake_start(started, worker):
    """Record the load without running it; the worker counts as running
    until wait() is called, as a real thread does at completion time."""
    started.append(worker.h5ad_path)
    worker.isRunning = lambda: worker._fake_running
    worker.wait = lambda *a: setattr(worker, '_fake_running', False)
    worker._fake_running = True
    return worker


def test_de_discards_a_load_that_finishes_after_the_study_changed(app, tmp_path, monkeypatch):
    """A second load request is ignored while one runs, so a study switch
    mid-load used to land the previous study's file as the new study's
    dataset (and its DE results in the new study's folder)."""
    from PyQt6.QtCore import QEvent
    from kosmic.gui.de_analysis import pages
    from kosmic.gui.de_analysis.workspace import DEWorkspace
    from kosmic.paths import processed_data_dir
    for name in ('a', 'b'):
        processed_data_dir(tmp_path / name).mkdir(parents=True)
        _study(processed_data_dir(tmp_path / name) / f'{name}.h5ad')
    ws = DEWorkspace()
    ws.set_project_directory(str(tmp_path / 'a'))
    started = []
    monkeypatch.setattr(pages.setup_page, 'run_worker', lambda w, **kw: _fake_start(started, w))
    ws.setup_page._start_load(str(processed_data_dir(tmp_path / 'a') / 'a.h5ad'))
    ws.set_project_directory(str(tmp_path / 'b'))            # switch while 'a' is in flight
    payload = (ad.read_h5ad(processed_data_dir(tmp_path / 'a') / 'a.h5ad'), started[0], 'a')
    ws.setup_page._on_loaded(payload)
    assert ws.current_adata is None                          # 'a' was not adopted for study 'b'
    assert [os.path.basename(p) for p in started] == ['a.h5ad', 'b.h5ad']
    ws.close()
    ws.deleteLater()
    app.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    app.processEvents()


def test_scrna_discards_a_load_that_finishes_after_the_study_changed(app, tmp_path, monkeypatch):
    from PyQt6.QtCore import QEvent
    from kosmic.gui.scrna.tabs import download_tab as dt
    from kosmic.gui.scrna.workspace import ScRNAWorkspace
    from kosmic.paths import processed_data_dir
    for name in ('a', 'b'):
        processed_data_dir(tmp_path / name).mkdir(parents=True)
        _study(processed_data_dir(tmp_path / name) / f'{name}.h5ad')
    ws = ScRNAWorkspace(embedded=True)
    ws.show()
    started = []
    monkeypatch.setattr(dt, 'run_worker', lambda w, **kw: _fake_start(started, w))
    ws.set_project_directory(str(tmp_path / 'a'))            # visible: auto-load starts
    ws.set_project_directory(str(tmp_path / 'b'))            # switch while 'a' is in flight
    payload = (ad.read_h5ad(processed_data_dir(tmp_path / 'a') / 'a.h5ad'), started[0], '30 cells')
    ws.download_tab._on_h5ad_loaded(payload)
    assert ws.current_adata is None
    assert [os.path.basename(p) for p in started] == ['a.h5ad', 'b.h5ad']
    ws.close()
    ws.deleteLater()
    app.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    app.processEvents()
