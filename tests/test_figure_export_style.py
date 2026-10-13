"""Figure export (ADR-009): the chosen format is written, text stays
editable in PDF and SVG, the background is white and dense scatters are
rasterised."""
from __future__ import annotations

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest

import kosmic.visualisation  # noqa: F401  (applies the export rcParams)
from kosmic.visualisation.style import (
    RASTERIZE_THRESHOLD, rasterize_dense_collections, save_figure,
)


def _figure(n_points=10):
    fig, ax = plt.subplots()
    fig.set_facecolor('#1e1e1e')   # dark GUI theme
    ax.scatter(np.arange(n_points), np.arange(n_points))
    ax.set_title('Label')
    return fig


@pytest.mark.parametrize('fmt', ['png', 'pdf', 'svg'])
def test_writes_only_the_chosen_format(tmp_path, fmt):
    out = save_figure(_figure(), tmp_path / f'fig.{fmt}', dpi=72)
    assert out == tmp_path / f'fig.{fmt}'
    assert [p.name for p in tmp_path.iterdir()] == [f'fig.{fmt}']


def test_rejects_unknown_format(tmp_path):
    with pytest.raises(ValueError):
        save_figure(_figure(), tmp_path / 'fig.jpg', dpi=72)


def test_pdf_text_is_not_type3(tmp_path):
    data = save_figure(_figure(), tmp_path / 'fig.pdf', dpi=72).read_bytes()
    assert b'/Type3' not in data
    assert b'/FontFile2' in data   # embedded TrueType (Type 42)


def test_svg_text_is_text(tmp_path):
    svg = save_figure(_figure(), tmp_path / 'fig.svg', dpi=72).read_text()
    assert '<text' in svg and 'Label' in svg


def test_background_is_white(tmp_path):
    out = save_figure(_figure(), tmp_path / 'fig.png', dpi=72)
    corner = plt.imread(out)[0, 0, :3]
    assert np.allclose(corner, 1.0)


def test_only_dense_scatters_are_rasterised():
    sparse, dense = _figure(10), _figure(RASTERIZE_THRESHOLD)
    assert rasterize_dense_collections(sparse) == 0
    assert rasterize_dense_collections(dense) == 1
    plt.close('all')


def test_elbow_plot_does_not_change_global_style():
    # plt.style.use() inside a plotting function turned every later figure
    # dark, which the white export background then made unreadable.
    from kosmic.visualisation.scrna.pca_plots import create_elbow_plot
    before = matplotlib.rcParams['axes.facecolor']
    create_elbow_plot(np.linspace(0.2, 0.01, 30), suggested_pcs=5, dark_mode=True)
    assert matplotlib.rcParams['axes.facecolor'] == before
    plt.close('all')


def test_preview_and_export_draw_under_one_lock(tmp_path):
    """Matplotlib is not thread-safe: a preview (worker thread) and an export
    (main thread) drawing at once failed in mathtext. Both hold MPL_LOCK."""
    import os
    from types import SimpleNamespace
    from unittest import mock
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6.QtWidgets import QApplication, QFileDialog
    from kosmic.gui.figure_export.pages._base import MPL_LOCK, FigurePage, RenderWorker

    app = QApplication.instance() or QApplication([])  # noqa: F841 (held for the test)
    held = []

    def func():
        held.append(MPL_LOCK.locked())
        return _figure()

    RenderWorker(func)._run()

    page = FigurePage(SimpleNamespace(data_version=lambda: 0))
    page._make_render_func = lambda: func
    with mock.patch.object(QFileDialog, 'getSaveFileName',
                           return_value=(str(tmp_path / 'f.png'), '')):
        page._export_current()
    assert held == [True, True]
    assert not MPL_LOCK.locked()


class _Settings:
    """In-memory QSettings, so a test never changes the user's saved settings."""
    store: dict = {}

    def __init__(self, *_):
        pass

    def value(self, key, default=None):
        return self.store.get(key, default)

    def setValue(self, key, value):
        self.store[key] = value


def _export(page, tmp_path, chosen=None):
    """Export *page*, the dialog accepting its suggested name and type unless
    *chosen* (a type picked in the dialog) is given. Returns the dialog mock."""
    from unittest import mock
    from PyQt6.QtWidgets import QFileDialog

    def dialog(_parent, _title, path, _filters, initial):
        return path, chosen or initial

    page._default_export_dir = lambda: tmp_path
    with mock.patch.object(QFileDialog, 'getSaveFileName', side_effect=dialog) as m:
        page._export_current()
    return m


def test_format_dropdown_sets_the_export_format(tmp_path):
    """Every figure page has a Format drop-down beside Export. The dialog
    opens on that format, and the choice is shared by every page."""
    import os
    from types import SimpleNamespace
    from unittest import mock
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6.QtWidgets import QApplication
    from kosmic.gui.figure_export.pages import _base

    app = QApplication.instance() or QApplication([])  # noqa: F841 (held for the test)
    _Settings.store = {}
    with mock.patch.object(_base, 'QSettings', _Settings):
        ws = SimpleNamespace(data_version=lambda: 0)
        page, other = _base.FigurePage(ws), _base.FigurePage(ws)
        page._make_render_func = lambda: _figure
        assert page._format.currentData() == 'png'

        page._format.setCurrentIndex(page._format.findData('pdf'))
        dialog = _export(page, tmp_path)
        assert dialog.call_args.args[4] == 'PDF (*.pdf)'
        assert (tmp_path / 'figure.pdf').exists()

        other.on_activated()
        assert other._format.currentData() == 'pdf'


def test_type_changed_in_the_dialog_sets_the_format(tmp_path):
    """The suggested name has no extension, so a type changed in the dialog
    still decides the format. A suggested '.png' used to stay."""
    import os
    from types import SimpleNamespace
    from unittest import mock
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6.QtWidgets import QApplication
    from kosmic.gui.figure_export.pages import _base

    app = QApplication.instance() or QApplication([])  # noqa: F841 (held for the test)
    _Settings.store = {}
    with mock.patch.object(_base, 'QSettings', _Settings):
        page = _base.FigurePage(SimpleNamespace(data_version=lambda: 0))
        page._make_render_func = lambda: _figure
        dialog = _export(page, tmp_path, chosen='SVG (*.svg)')
    assert dialog.call_args.args[2] == str(tmp_path / 'figure')
    assert (tmp_path / 'figure.svg').exists()


def test_every_figure_page_has_the_format_dropdown():
    import os
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6.QtWidgets import QApplication
    from kosmic.gui.figure_export.workspace import FigureExportWorkspace

    app = QApplication.instance() or QApplication([])  # noqa: F841 (held for the test)
    ws = FigureExportWorkspace()
    assert ws._pages
    for page in ws._pages:
        formats = [page._format.itemData(i) for i in range(page._format.count())]
        assert formats == ['png', 'pdf', 'svg'], page.TITLE


def test_existing_file_is_replaced_only_when_confirmed(tmp_path):
    """The dialog sees the name without its extension and cannot warn that
    figure.pdf exists, so the page asks before replacing it."""
    import os
    from types import SimpleNamespace
    from unittest import mock
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6.QtWidgets import QApplication
    from kosmic.gui.figure_export.pages import _base

    app = QApplication.instance() or QApplication([])  # noqa: F841 (held for the test)
    target = tmp_path / 'figure.pdf'
    target.write_bytes(b'old')
    _Settings.store = {}
    with mock.patch.object(_base, 'QSettings', _Settings):
        page = _base.FigurePage(SimpleNamespace(data_version=lambda: 0))
        page._make_render_func = lambda: _figure
        page._format.setCurrentIndex(page._format.findData('pdf'))
        with mock.patch.object(_base.dialogs, 'confirm', return_value=False) as ask:
            _export(page, tmp_path)
        assert ask.called and target.read_bytes() == b'old'
        with mock.patch.object(_base.dialogs, 'confirm', return_value=True):
            _export(page, tmp_path)
    assert target.read_bytes().startswith(b'%PDF')
