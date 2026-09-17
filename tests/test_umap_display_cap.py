"""The interactive embedding draws at most UI_EMBEDDING_MAX_POINTS cells.

Legend counts and the hover lookup still see every cell; only the drawn
points are a fixed random subset, so atlas-sized studies stay responsive.
"""
import os
import time

import numpy as np
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _widget(app):
    from kosmic.gui.scrna.tabs.cluster_tab import _UMAPWidget
    return _UMAPWidget()


def test_draw_order_caps_and_is_stable(monkeypatch):
    from kosmic.gui.scrna.tabs import cluster_tab
    monkeypatch.setattr(cluster_tab, 'UI_EMBEDDING_MAX_POINTS', 1000)
    w = cluster_tab._UMAPWidget
    a = w._draw_order(5000)
    b = w._draw_order(5000)
    assert len(a) == 1000 and np.array_equal(a, b) and len(set(a.tolist())) == 1000
    assert len(w._draw_order(800)) == 800
    monkeypatch.setattr(cluster_tab, 'UI_EMBEDDING_MAX_POINTS', 0)
    assert len(w._draw_order(5000)) == 5000


def test_categorical_plot_counts_every_cell_but_draws_the_cap(app, monkeypatch):
    from kosmic.gui.scrna.tabs import cluster_tab
    monkeypatch.setattr(cluster_tab, 'UI_EMBEDDING_MAX_POINTS', 2000)
    w = _widget(app)
    n = 12_000
    rng = np.random.default_rng(1)
    coords = rng.normal(size=(n, 2))
    labels = np.array(['a'] * 9000 + ['b'] * 3000, dtype=object)
    w.set_data_categorical(coords, labels)
    assert len(w._scatter.data) == 2000
    assert w._legend_entries[0][1].endswith('(9,000)')
    assert '2,000 of 12,000 cells drawn' in w._info_label.text()
    # every drawn point carries its own label
    drawn = np.array([d for d in w._scatter.data['data']], dtype=object)
    assert set(drawn.tolist()) == {'a', 'b'}


def test_continuous_plot_draws_the_cap(app, monkeypatch):
    from kosmic.gui.scrna.tabs import cluster_tab
    monkeypatch.setattr(cluster_tab, 'UI_EMBEDDING_MAX_POINTS', 500)
    w = _widget(app)
    coords = np.random.default_rng(2).normal(size=(3000, 2))
    values = np.linspace(0, 5, 3000)
    w.set_data_continuous(coords, values)
    assert len(w._scatter.data) == 500
    assert '500 of 3,000 cells drawn' in w._info_label.text()
    assert 'range 0.00' in w._info_label.text()


def test_building_a_large_plot_is_fast(app, monkeypatch):
    from kosmic.gui.scrna.tabs import cluster_tab
    monkeypatch.setattr(cluster_tab, 'UI_EMBEDDING_MAX_POINTS', 100_000)
    w = _widget(app)
    n = 300_000
    coords = np.random.default_rng(3).normal(size=(n, 2)).astype(np.float32)
    labels = np.random.default_rng(4).integers(0, 12, n).astype(str).astype(object)
    t0 = time.perf_counter()
    w.set_data_categorical(coords, labels)
    assert time.perf_counter() - t0 < 5.0


def test_unclustered_cells_do_not_break_the_plot(app):
    """Excluded donors carry NaN in obs['leiden'] -> label 'nan' beside '0'..'11'."""
    from kosmic.gui.scrna.tabs.cluster_tab import _label_sort_key
    w = _widget(app)
    labels = np.array(['0'] * 5 + ['11'] * 3 + ['2'] * 2 + ['nan'] * 4, dtype=object)
    w.set_data_categorical(np.random.default_rng(0).normal(size=(14, 2)), labels)
    names = [e[1].split('   ')[0] for e in w._legend_entries]
    assert names == ['0', '2', '11', 'nan']
    assert sorted(['Fibroblast', 'nan', '3', 'Unknown', '10'], key=_label_sort_key) == \
        ['3', '10', 'Fibroblast', 'Unknown', 'nan']
