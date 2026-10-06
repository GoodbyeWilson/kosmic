"""Focus mode in the Cluster step (ADR-008), driven as a user would.

Workers run synchronously and dialogs answer from a script; everything
else is the real scRNA workspace in an offscreen window.
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

CENTRES = {'0': (0, 0), '1': (12, 0), '2': (0, 12)}


def _study(path):
    """Three clusters on a UMAP; cluster '1' holds two planted states."""
    rng = np.random.default_rng(1)
    g = 150
    base = rng.gamma(0.6, 1.0, g)
    kinds = [('0', 'A', 150), ('1', 'B1', 120), ('1', 'B2', 120), ('2', 'C', 110)]
    rows, obs, xy = [], [], []
    for cluster, kind, n in kinds:
        p = base.copy()
        p[{'A': slice(0, 20), 'B1': slice(20, 35), 'B2': slice(35, 50),
           'C': slice(50, 70)}[kind]] *= 12
        rows.append(rng.multinomial(1500, p / p.sum(), size=n))
        obs += [(cluster, kind)] * n
        xy.append(rng.normal(CENTRES[cluster], 1.0, size=(n, 2)))
    n = len(obs)
    a = ad.AnnData(
        X=sp.csr_matrix(np.vstack(rows).astype(np.float32)),
        obs=pd.DataFrame({
            'leiden': pd.Categorical([c for c, _ in obs], categories=['0', '1', '2']),
            'kind': [k for _, k in obs],
            'cell_type': pd.Categorical([{'0': 'CM', '1': 'Fib', '2': 'EC'}[c] for c, _ in obs]),
            'cell_type_auto': pd.Categorical([{'0': 'CM', '1': 'Fib', '2': 'EC'}[c] for c, _ in obs]),
        }, index=[f'c{i}' for i in range(n)]),
        var=pd.DataFrame(index=[f'G{j}' for j in range(g)]))
    a.obsm['X_umap'] = np.vstack(xy).astype(np.float32)
    a.write_h5ad(path)
    return a


@pytest.fixture
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.fixture
def scripted(monkeypatch):
    """Synchronous workers; dialogs answered from 'answers'."""
    from kosmic.gui.shared import dialogs
    from kosmic.gui.scrna.tabs import cluster_tab

    def run_now(worker, *, on_finished, on_failed=None, **_kw):
        on_finished(worker._run())
        return worker

    answers = {'choose': 'Discard', 'confirm': True, 'shown': []}
    monkeypatch.setattr(cluster_tab, 'run_worker', run_now)
    for name in ('info', 'warning', 'error'):
        monkeypatch.setattr(dialogs, name, lambda p, t, m, n=name: answers['shown'].append((n, t)))
    monkeypatch.setattr(dialogs, 'confirm', lambda *a, **k: answers['confirm'])
    monkeypatch.setattr(dialogs, 'choose', lambda *a, **k: answers['choose'])
    return answers


def _click(tab, xy, ctrl=False):
    from PyQt6.QtCore import QPointF, Qt
    vb = tab._umap_widget._plot.getPlotItem().vb

    class Click:
        def button(self):
            return Qt.MouseButton.LeftButton

        def modifiers(self):
            return (Qt.KeyboardModifier.ControlModifier if ctrl
                    else Qt.KeyboardModifier.NoModifier)

        def scenePos(self):
            return vb.mapViewToScene(QPointF(*xy))

    tab._umap_widget._on_plot_clicked(Click())


def test_select_subcluster_annotate_apply_and_close(app, scripted, tmp_path):
    from kosmic import provenance
    from kosmic.gui.scrna.tabs.annotate_tab import SubsetSource, WorkspaceSource
    from kosmic.gui.scrna.workspace import ScRNAWorkspace

    path = tmp_path / 'processed_data' / 'S.h5ad'
    path.parent.mkdir()
    a = _study(path)
    ws = ScRNAWorkspace(embedded=True)
    ws.resize(1400, 900)
    ws.show()
    ws.set_adata(a, str(path))
    tab = ws.cluster_tab
    ws.switch_tab(ws.stack.indexOf(tab))
    tab.set_data(a, str(path))
    tab._show_embedding(a, 'leiden')
    tab._right_tabs.setCurrentIndex(1)
    app.processEvents()
    w = tab._umap_widget
    # the study's own cells, to click on
    one = a.obsm['X_umap'][(a.obs['leiden'] == '1').to_numpy()][0]
    two = a.obsm['X_umap'][(a.obs['leiden'] == '2').to_numpy()][0]

    _click(tab, one)
    assert w.selection() == ('leiden', ['1'], 240) and tab._focus_bar.isVisibleTo(tab)
    _click(tab, two, ctrl=True)
    assert w.selection()[1] == ['1', '2']
    _click(tab, two, ctrl=True)
    _click(tab, (30, 30))                                      # empty space
    assert w.selection() is None and not tab._focus_bar.isVisibleTo(tab)
    _click(tab, one)

    tab.focus_resolution_spin.setValue(0.3)
    tab._run_focus()
    focus = tab._focus
    assert focus is not None and focus.groups == ['1'] and focus.adata.n_obs == 240
    assert isinstance(ws.annotate_tab._source, SubsetSource) and w._adata is focus.adata
    assert not tab._pca_card.isEnabled() and not tab._cluster_card.isEnabled()
    assert tab._focus_apply_btn.isVisibleTo(tab) and not tab._focus_go_btn.isVisibleTo(tab)

    ws.annotate_tab._apply_cluster_reannotation(0, 'Fib_state')     # name subcluster 0
    named = focus.adata.obs_names[(focus.adata.obs['cell_type'] == 'Fib_state').to_numpy()]
    tab.focus_resolution_spin.setValue(0.5)                        # resolution only:
    tab._run_focus()                                               # same population
    assert tab._focus is focus and focus.settings['resolution'] == 0.5
    assert (a.obs['cell_type'] == 'Fib').sum() == 240               # study untouched

    assert tab._apply_focus(wait=True)
    assert tab._focus is None and isinstance(ws.annotate_tab._source, WorkspaceSource)
    assert tab._pca_card.isEnabled()
    saved = ad.read_h5ad(path)
    assert set(saved.obs.loc[saved.obs['leiden'] == '1', 'leiden_sub'].astype(str)) == {
        f'1.{c}' for c in focus.adata.obs['leiden'].cat.categories}
    assert (saved.obs.loc[saved.obs['leiden'] == '0', 'leiden_sub'] == '0').all()
    assert (saved.obs.loc[named, 'cell_type'] == 'Fib_state').all()
    stage = provenance.last_stage(path.parent, 'subcluster')['params']
    assert stage['groups'] == ['1'] and stage['annotation'] == {'method': 'manual'}
    assert stage['manual_edits'][0]['to'] == 'Fib_state'
    assert stage['manual_edits'][0]['resolution'] == 0.3            # made before re-running

    # Closing the study with a focus open asks; Discard leaves it as it was.
    n_stages = len(provenance.load(path.parent)['stages'])
    tab._show_embedding(a, 'leiden')
    _click(tab, two)
    tab._run_focus()
    assert tab._focus is not None
    scripted['choose'] = 'Cancel'
    assert tab.settle_focus(allow_cancel=True) is False and tab._focus is not None
    scripted['choose'] = 'Discard'
    ws.release_dataset()
    assert tab._focus is None and ws.current_adata is None
    assert len(provenance.load(path.parent)['stages']) == n_stages
