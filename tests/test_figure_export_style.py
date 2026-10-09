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
