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
    assert obs2['cell_type'].iloc[0] == 'CM' and pd.isna(obs2['cell_type'].iloc[1])


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


def test_annotation_and_manual_relabels_are_recorded(app, tmp_path, monkeypatch):
    """An automated pass records the method and per-cluster labels; hand
    relabels are recorded on save, and a new automated pass starts afresh."""
    from kosmic.gui.shared import dialogs
    from kosmic.gui.scrna.tabs.annotate_tab import AnnotateTab
    monkeypatch.setattr(dialogs, 'info', lambda *a, **k: None)

    class MW:
        current_adata = None
        current_h5ad_path = None
        _adata_version = 0
        recorded = []

        def record_provenance(self, stage, params, **kw):
            self.recorded.append((stage, params, kw))

    ws = MW()
    t = AnnotateTab(ws)
    n = 60
    a = ad.AnnData(
        X=sp.csr_matrix(np.ones((n, 3), dtype=np.float32)),
        obs=pd.DataFrame({
            'leiden': pd.Categorical(['0'] * 40 + ['1'] * 20),
            'cell_type': pd.Categorical(['CM'] * 40 + ['Fib'] * 20),
        }, index=[f'c{i}' for i in range(n)]))
    target = tmp_path / 's.h5ad'
    ws.current_adata, ws.current_h5ad_path, t.adata = a, str(target), a

    t._record_annotation({'method': 'CellTypist', 'name': 'Healthy_Adult_Heart.pkl'})
    stage, params, _ = ws.recorded[-1]
    assert stage == 'annotate' and params['name'] == 'Healthy_Adult_Heart.pkl'
    assert params['clusters'] == {'0': {'label': 'CM', 'n_cells': 40},
                                  '1': {'label': 'Fib', 'n_cells': 20}}

    t._populate_cluster_table()
    t._apply_cluster_reannotation(1, 'Unknown')
    t._save_annotations()
    stage, params, kw = ws.recorded[-1]
    assert stage == 'annotate_manual' and kw == {'replace': True}
    assert params['edits'] == [{'cluster': '1', 'from': 'Fib', 'to': 'Unknown', 'n_cells': 20}]

    t._record_annotation({'method': 'CellTypist'})
    assert 'manual_relabels' not in a.uns


def test_a_subset_source_is_annotated_in_memory_and_never_written(app, tmp_path, monkeypatch):
    """Focus mode (ADR-008) points the panel at a subset held in memory.
    Edits, annotation passes and stages stay with the subset; the open
    study, its file and its record are untouched."""
    from kosmic.gui.shared import dialogs
    from kosmic.gui.scrna.tabs.annotate_tab import (
        AnnotateTab, ReferenceAnnotationWorker, SubsetSource)
    import kosmic.scrna.annotate.reference as ref_mod
    warnings = []
    monkeypatch.setattr(dialogs, 'warning', lambda *a, **k: warnings.append(a[2]))

    class MW:
        current_adata = None
        current_h5ad_path = None
        _adata_version = 5
        recorded = []

        def record_provenance(self, stage, params, **kw):
            self.recorded.append(stage)

    def _cells(clusters):
        n = len(clusters)
        return ad.AnnData(
            X=sp.csr_matrix(np.ones((n, 3), dtype=np.float32)),
            obs=pd.DataFrame({'leiden': pd.Categorical(clusters),
                              'cell_type': pd.Categorical(['Fib'] * n),
                              'cell_type_auto': pd.Categorical(['Fib'] * n)},
                             index=[f'c{i}' for i in range(n)]))

    ws = MW()
    study, sub = _cells(['0'] * 30), _cells(['0'] * 20 + ['1'] * 10)
    ws.current_adata, ws.current_h5ad_path = study, str(tmp_path / 's.h5ad')
    t = AnnotateTab(ws)
    src = SubsetSource(sub)
    t.set_source(src)
    assert t.adata is sub and t.save_annot_btn.isHidden()
    assert t.cluster_table.rowCount() == 2

    t._apply_cluster_reannotation(1, 'Fib_activated')
    assert (sub.obs['cell_type'] == 'Fib_activated').sum() == 10 and src.version == 1
    assert (study.obs['cell_type'] == 'Fib').all() and ws._adata_version == 5

    t._record_manual_relabels()
    t._record_manual_relabels()                       # replaced, not appended
    t._record_annotation({'method': 'CellTypist'})
    assert [stage for stage, _ in src.records] == ['annotate_manual', 'annotate']
    assert not ws.recorded

    def fake_run(adata, ref, **kw):
        adata.obs['cell_type'] = pd.Categorical(['FB1'] * 20 + ['FB2'] * 10)
        return adata, {'0': 'FB1', '1': 'FB2'}, {}

    monkeypatch.setattr(ref_mod, 'run_reference_annotation', fake_run)
    out, _ = ReferenceAnnotationWorker(sub, 'any', t._resolve_h5ad_path())._run()
    assert out is sub and (sub.obs['cell_type'] == 'FB2').sum() == 10
    t._save_annotations()
    assert not list(tmp_path.iterdir()) and not warnings

    t.set_source(None)
    assert t.adata is study and not t.save_annot_btn.isHidden()


def test_the_table_shows_agreement_with_each_cells_prediction(app):
    """After majority voting every cell of a cluster has the same label,
    so the table reads agreement from the per-cell predictions."""
    from kosmic.gui.scrna.tabs.annotate_tab import AnnotateTab

    class MW:
        current_adata = None
        _adata_version = 0

    t = AnnotateTab(MW())
    n = 200
    a = ad.AnnData(
        X=sp.csr_matrix(np.ones((n, 5), dtype=np.float32)),
        obs=pd.DataFrame({
            'leiden': pd.Categorical(['0'] * 100 + ['1'] * 100),
            'cell_type': pd.Categorical(['Monocyte'] * 100 + ['Fibroblast'] * 100),
            'cell_type_auto': pd.Categorical(['Monocyte'] * 100 + ['Fibroblast'] * 100),
            'predicted_labels': (['CD16+Mo'] * 40 + ['LYVE1+IGF1+MP'] * 35
                                 + ['DC'] * 25 + ['FB1'] * 100),
        }, index=[f'c{i}' for i in range(n)]))
    t.main_window.current_adata = a
    t.adata = a
    t._populate_cluster_table()
    tbl = t.cluster_table
    assert tbl.item(0, 7).text() == 'Below threshold (40%)'
    assert tbl.item(0, 8).text() == 'Macrophage (35%)'
    assert tbl.item(0, 9).text().startswith('Mixed: 40% of cells')
    assert tbl.item(1, 7).text() == 'High (100%)'
    assert tbl.item(1, 9).text() == ''


def test_a_subclustered_study_lists_and_relabels_its_subclusters(app):
    """After subclustering the table works on 'leiden_sub', so a subcluster
    can be relabelled on its own (ADR-008)."""
    from kosmic.gui.scrna.tabs.annotate_tab import AnnotateTab

    class MW:
        current_adata = None
        _adata_version = 0

    t = AnnotateTab(MW())
    n = 120
    sub = ['0'] * 60 + ['1.0'] * 40 + ['1.1'] * 20
    a = ad.AnnData(
        X=sp.csr_matrix(np.ones((n, 5), dtype=np.float32)),
        obs=pd.DataFrame({
            'leiden': pd.Categorical(['0'] * 60 + ['1'] * 60),
            'leiden_sub': pd.Categorical(sub, categories=['0', '1.0', '1.1']),
            'cell_type': pd.Categorical(['CM'] * 60 + ['T Cell'] * 60),
            'cell_type_auto': pd.Categorical(['CM'] * 60 + ['T Cell'] * 60),
        }, index=[f'c{i}' for i in range(n)]))
    t.main_window.current_adata = a
    t.adata = a
    t._populate_cluster_table()
    tbl = t.cluster_table
    assert [tbl.item(r, 0).text() for r in range(tbl.rowCount())] == ['0', '1.0', '1.1']

    t._apply_cluster_reannotation(2, 'NK Cell')
    assert (a.obs.loc[a.obs['leiden_sub'] == '1.1', 'cell_type'] == 'NK Cell').all()
    assert (a.obs.loc[a.obs['leiden_sub'] == '1.0', 'cell_type'] == 'T Cell').all()
    t._record_manual_relabels()     # no workspace recorder: must not fail
    assert a.uns['manual_relabels']['cluster'] == ['1.1']


def test_fine_type_shows_each_rows_own_predictions(app):
    """'cell_type_fine' keeps the whole cluster's vote after subclustering;
    the column shows the row's most common per-cell prediction instead."""
    from kosmic.gui.scrna.tabs.annotate_tab import AnnotateTab

    class MW:
        current_adata = None
        _adata_version = 0

    t = AnnotateTab(MW())
    n = 100
    a = ad.AnnData(
        X=sp.csr_matrix(np.ones((n, 5), dtype=np.float32)),
        obs=pd.DataFrame({
            'leiden': pd.Categorical(['23'] * n),
            'leiden_sub': pd.Categorical(['23.0'] * 50 + ['23.1'] * 50),
            'cell_type': pd.Categorical(['T Cell'] * 50 + ['NK Cell'] * 50),
            'cell_type_auto': pd.Categorical(['CD4+ T Cell'] * n),
            'cell_type_fine': pd.Categorical(['CD4+T_naive'] * n),   # stale
            'predicted_labels': ['CD8+T_em'] * 50 + ['NK_CD16hi'] * 50,
        }, index=[f'c{i}' for i in range(n)]))
    t.main_window.current_adata = a
    t.adata = a
    t._populate_cluster_table()
    assert [t.cluster_table.item(r, 5).text() for r in range(2)] == ['CD8+T_em', 'NK_CD16hi']
