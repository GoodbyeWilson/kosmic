"""Cell Types figures: marker dot plot, marker score embeddings, composition."""
from __future__ import annotations

import os
from types import SimpleNamespace
from unittest import mock

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from kosmic.visualisation.scrna.composition import composition_table, create_composition_plot
from kosmic.visualisation.scrna.markers import (
    create_marker_dotplot, create_marker_score_grid, order_rows, panel_scores, scale_scores,
)
from kosmic.visualisation.scrna.umap import _soft_palette, category_colours

PANEL = {'Cardiomyocyte': ['TNNT2', 'MYH6'], 'Fibroblast': ['DCN', 'LUM'],
         'Endothelial': ['VWF', 'PECAM1']}


def _adata(seed=0):
    """Three typed cell populations, each expressing its own panel genes,
    plus filler genes; X is log1p of counts per 10,000."""
    rng = np.random.default_rng(seed)
    types = ['Ventricular Cardiomyocyte', 'Fibroblast', 'Endothelial Cell', 'Unknown']
    genes = [g for gs in PANEL.values() for g in gs] + [f'G{i}' for i in range(60)]
    n = 400
    labels = np.repeat(types, n // 4)
    counts = rng.poisson(1.0, (n, len(genes))).astype(float)
    for t, (_, gs) in zip(types, PANEL.items()):
        for g in gs:
            counts[labels == t, genes.index(g)] += 30
    x = np.log1p(counts / counts.sum(1, keepdims=True) * 1e4)
    obs = pd.DataFrame({
        'cell_type': pd.Categorical(labels),
        '_role': np.tile(['control', 'disease'], n // 2),
        'sample': [f's{i % 6}' for i in range(n)],
    }, index=[f'c{i}' for i in range(n)])
    a = ad.AnnData(X=sp.csr_matrix(x.astype(np.float32)), obs=obs,
                   var=pd.DataFrame(index=genes))
    a.obsm['X_umap'] = rng.normal(size=(n, 2))
    return a


def test_bundled_cardiac_panel_loads():
    from kosmic.reference.marker_panels import builtin_panels
    panel = builtin_panels()['Cardiac (human, curated)']
    assert len(panel) == 14 and panel['Cardiomyocyte'][0] == 'TNNT2'


def test_rows_follow_the_panel_by_expression_not_name():
    mean = pd.DataFrame({'LYVE1+ Macrophage': [0, 0, 3], 'Fibro': [0, 2, 0], 'CM': [4, 0, 0]},
                        index=['TNNT2', 'DCN', 'CD163'])
    owner = {'TNNT2': 'Cardiomyocyte', 'DCN': 'Fibroblast', 'CD163': 'Myeloid'}
    assert order_rows(mean, owner) == ['CM', 'Fibro', 'LYVE1+ Macrophage']


def test_scaled_scores_span_zero_to_one():
    v = scale_scores(np.r_[np.zeros(990), np.linspace(1, 5, 10)])
    assert v.min() == 0 and v.max() == 1


def test_panel_scores_peak_in_the_marked_type():
    a = _adata()
    scores = panel_scores(a, PANEL)
    cm = (a.obs['cell_type'] == 'Ventricular Cardiomyocyte').to_numpy()
    v, genes = scores['Cardiomyocyte']
    assert genes == ['TNNT2', 'MYH6']
    assert v[cm].mean() > 5 * v[~cm].mean()


def test_composition_rows_sum_to_100_largest_first():
    pct = composition_table(['A', 'A', 'A', 'B', 'A', 'B'], ['x', 'x', 'x', 'x', 'y', 'y'],
                            ['y', 'x'])
    assert list(pct.index) == ['y', 'x'] and list(pct.columns) == ['A', 'B']
    assert np.allclose(pct.sum(axis=1), 100)


def test_composition_colours_match_the_embedding_plot():
    labels = ['b', 'a', 'c', 'a']
    assert category_colours(labels) == dict(zip(['a', 'b', 'c'], _soft_palette(3)))


def test_plots_render():
    a = _adata()
    from kosmic.scrna.inspect.gene_group import marker_panel, summarise_gene_group
    genes, owner = marker_panel(PANEL, cell_types=list(PANEL))
    s = summarise_gene_group(a, genes, group_col='cell_type')
    for fig in (create_marker_dotplot(s.mean_matrix, s.pct_matrix, owner),
                create_marker_score_grid(a.obsm['X_umap'], panel_scores(a, PANEL)),
                create_composition_plot(composition_table(a.obs['cell_type'], a.obs['_role']),
                                        category_colours(a.obs['cell_type']),
                                        blocks=[('Control', 1), ('Disease', 1)])):
        assert fig is not None
        plt.close(fig)


@pytest.fixture
def qapp():
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _workspace(a):
    return SimpleNamespace(scrna_ws=SimpleNamespace(current_adata=a, current_project_dir=None),
                           de_ws=None, meta_ws=None, data_version=lambda: 0)


def _export(page, path):
    from PyQt6.QtWidgets import QFileDialog
    with mock.patch.object(QFileDialog, 'getSaveFileName', return_value=(str(path), '')):
        page._export_current()
    assert page._status_lbl.text().startswith('Exported'), page._status_lbl.text()


def test_dot_plot_page_exports_numbers_without_unlabelled_cells(qapp, tmp_path):
    from kosmic.gui.figure_export.pages.cell_types import MarkerDotPlotPage, _expression_label
    a = _adata()
    assert 'log1p(CP10K)' in _expression_label(a)
    page = MarkerDotPlotPage(_workspace(a))
    page._controls.blockSignals(True)
    page.on_activated()
    _export(page, tmp_path / 'dot.png')
    table = pd.read_csv(tmp_path / 'dot.csv')
    assert 'Unknown' not in set(table['cell_type'])
    assert {'gene', 'marks', 'mean_expression', 'pct_expressing'} <= set(table.columns)


def test_score_page_renders(qapp, tmp_path):
    from kosmic.gui.figure_export.pages.cell_types import MarkerScoreUMAPPage
    page = MarkerScoreUMAPPage(_workspace(_adata()))
    page._controls.blockSignals(True)
    assert page.dependencies_met()
    _export(page, tmp_path / 'scores.png')


@pytest.mark.parametrize('by_sample', [False, True])
def test_composition_page(qapp, tmp_path, by_sample):
    from kosmic.gui.figure_export.pages.cell_types import BARS_SAMPLE, CompositionPage
    page = CompositionPage(_workspace(_adata()))
    page._controls.blockSignals(True)
    page.on_activated()
    if by_sample:
        page._controls.bars.setCurrentText(BARS_SAMPLE)
    _export(page, tmp_path / 'comp.png')
    table = pd.read_csv(tmp_path / 'comp.csv', index_col='cell_type')
    assert 'Unknown' not in table.index
    assert list(table.columns) == ([f's{i}' for i in (0, 2, 4, 1, 3, 5)] if by_sample
                                   else ['Control', 'Disease'])
    assert np.allclose(table.sum(), 100, atol=0.01)
