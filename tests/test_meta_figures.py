"""Meta-analysis figures (ADR-010): per-study effects, leave-one-out
re-pooling, the reports beside the figures, and the Figures pages reading
the Meta-Analysis workspace's figure_inputs()."""
from __future__ import annotations

import os
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from kosmic.visualisation.meta.study_forest import (
    create_loo_plot, create_study_forest, forest_report, loo_estimates,
    loo_report, study_effects,
)

EFFECTS = {'S1': 0.6, 'S2': 0.7, 'S3': 1.2, 'S4': 0.6}


def _studies(item='PKD2'):
    return [(s, pd.DataFrame({'names': [item, 'OTHER'], 'logfoldchanges': [v, 0.0],
                              'se': [0.15, 0.1], 'pvals_adj': [0.01, 0.9]}))
            for s, v in EFFECTS.items()]


def _pooled():
    return pd.Series({'names': 'PKD2', 'logfoldchanges': 0.75, 'ci_lower': 0.55,
                      'ci_upper': 0.95, 'tau_squared': 0.02, 'heterogeneity_i2': 0.481,
                      'fdr': 1e-10, 'pooling_method': 'consensus'})


def test_study_effects_skip_studies_without_the_item():
    studies = _studies() + [('S5', pd.DataFrame({'names': ['OTHER'],
                                                 'logfoldchanges': [1.0], 'se': [0.1]}))]
    eff = study_effects('PKD2', studies)
    assert list(eff['study']) == list(EFFECTS)
    assert np.allclose(eff['ci_upper'] - eff['ci_lower'], 2 * 1.959964 * 0.15)


def test_loo_repools_each_fold():
    loo = loo_estimates('PKD2', _studies(), _pooled())
    assert list(loo['omitted']) == list(EFFECTS)
    assert (loo['n_studies'] == 3).all()
    # Leaving out the outlying study pulls the estimate down most.
    assert loo.loc[loo['omitted'] == 'S3', 'logfoldchanges'].iloc[0] == loo['logfoldchanges'].min()


def test_reports_state_percent_i2_and_method():
    eff = study_effects('PKD2', _studies())
    text = forest_report('PKD2', eff, _pooled())
    assert 'I2 = 48.1%' in text and 'weight_pct' in text
    loo_text = loo_report('PKD2', loo_estimates('PKD2', _studies(), _pooled()), _pooled())
    assert 'Leave-one-out re-pooled with: reml' in loo_text


def test_plots_render_and_handle_no_data():
    eff = study_effects('PKD2', _studies())
    for fig in (create_study_forest('PKD2', eff, _pooled()),
                create_loo_plot('PKD2', loo_estimates('PKD2', _studies(), _pooled()), _pooled())):
        assert fig is not None
        plt.close(fig)
    assert create_study_forest('X', eff.iloc[0:0]) is None


@pytest.fixture
def qapp():
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.mark.parametrize('level', ['gene', 'pathway'])
def test_pages_read_figure_inputs(qapp, level, tmp_path):
    from kosmic.gui.figure_export.pages import meta_figures as mf
    results = pd.DataFrame([_pooled()]).assign(names='PKD2')
    inputs = {'gene_results': None, 'pathway_results': None, 'gene_studies': [],
              'pathway_studies': [], 'selection': 'Cardiomyocyte', 'project_folder': None}
    inputs[f'{level}_results'] = results
    inputs[f'{level}_studies'] = _studies()
    ws = SimpleNamespace(meta_ws=SimpleNamespace(figure_inputs=lambda: inputs),
                         data_version=lambda: 0)

    volcano = mf.MetaVolcanoPage if level == 'gene' else mf.PathwayMetaVolcanoPage
    forest = mf.StudyForestPage if level == 'gene' else mf.PathwayForestPage
    pages = [volcano(ws), forest(ws)]
    assert all(p.dependencies_met() for p in pages)
    pages[1]._controls.populate_items(['PKD2'])
    for p in pages:
        plt.close(p._make_render_func()())
    written = pages[1]._export_companions(tmp_path / 'forest.png')
    assert written == [tmp_path / 'forest.txt'] and 'PKD2' in written[0].read_text()

    loo = mf.LeaveOneOutPage(ws)
    # Without signals, so no background re-render runs alongside the test's.
    loo._controls.level.blockSignals(True)
    loo._controls.level.setCurrentText(mf.LEVEL_PATHWAY if level == 'pathway' else mf.LEVEL_GENE)
    loo._controls.populate_items(['PKD2'])
    plt.close(loo._make_render_func()())


def test_study_labels_drop_the_runs_cell_type():
    from kosmic.meta_analysis.io import study_labels
    names = ['Chaffin_Ventricular_Cardiomyocyte', 'GSE183852_Ventricular_Cardiomyocyte']
    assert study_labels(names, 'decontX_counts/Ventricular_Cardiomyocyte') == ['Chaffin', 'GSE183852']
    assert study_labels(names, 'Ventricular_Cardiomyocyte') == ['Chaffin', 'GSE183852']
    # Kept whole: a mixed selection, a name without the cell type, a clash.
    assert study_labels(names, 'decontX_counts/mixed') == names
    assert study_labels(names + ['external'], 'Ventricular_Cardiomyocyte') == names + ['external']
    assert study_labels(['A_B_Cell', 'A_B_Cell'], 'B_Cell') == ['A_B_Cell', 'A_B_Cell']


def test_forest_and_loo_label_studies_without_the_cell_type(qapp, tmp_path):
    from kosmic.gui.figure_export.pages import meta_figures as mf
    studies = [(f'{s}_Ventricular_Cardiomyocyte', df) for s, df in _studies()]
    inputs = {'gene_results': pd.DataFrame([_pooled()]), 'gene_studies': studies,
              'pathway_results': None, 'pathway_studies': [], 'project_folder': None,
              'selection': 'decontX_counts/Ventricular_Cardiomyocyte'}
    ws = SimpleNamespace(meta_ws=SimpleNamespace(figure_inputs=lambda: inputs),
                         data_version=lambda: 0)
    forest, loo = mf.StudyForestPage(ws), mf.LeaveOneOutPage(ws)
    loo._controls.level.blockSignals(True)
    for page in (forest, loo):
        page._controls.populate_items(['PKD2'])
        fig = page._make_render_func()()
        labels = [t.get_text() for t in fig.axes[0].get_yticklabels()]
        plt.close(fig)
        assert set(EFFECTS) <= set(labels)
        assert not any('Cardiomyocyte' in t for t in labels)
        report = page._export_companions(tmp_path / f'{page.TITLE}.png')[0].read_text()
        assert 'Ventricular_Cardiomyocyte' not in report


def test_auto_size_boxes_show_the_size_drawn(qapp):
    """With Auto ticked the boxes show the size the preview was drawn at,
    and unticking Auto keeps it, so the figure does not change size."""
    from kosmic.gui.figure_export.pages import meta_figures as mf
    from kosmic.gui.figure_export.pages._base import _fig_to_png_bytes
    inputs = {'gene_results': pd.DataFrame([_pooled()]), 'gene_studies': _studies(),
              'pathway_results': None, 'pathway_studies': [], 'project_folder': None,
              'selection': 'Cardiomyocyte'}
    page = mf.StudyForestPage(SimpleNamespace(
        meta_ws=SimpleNamespace(figure_inputs=lambda: inputs), data_version=lambda: 0))
    page._controls.populate_items(['PKD2'])
    size = page._controls.figsize
    assert size.get_figsize() is None                      # Auto by default
    png = _fig_to_png_bytes(page._make_render_func()())
    drawn = png.size_inches
    page._on_render_finished(png)
    assert (size.w_spin.value(), size.h_spin.value()) == tuple(round(v, 1) for v in drawn)
    size.auto_cb.blockSignals(True)
    size.auto_cb.setChecked(False)
    assert size.get_figsize() == tuple(round(v, 1) for v in drawn)


def test_pages_unavailable_without_meta_workspace(qapp):
    from kosmic.gui.figure_export.pages.meta_figures import MetaVolcanoPage
    page = MetaVolcanoPage(SimpleNamespace(meta_ws=None, data_version=lambda: 0))
    assert not page.dependencies_met()


def test_meta_volcano_gates_on_fdr_unless_lfc_set(qapp):
    from kosmic.gui.figure_export.pages.meta_figures import MetaVolcanoPage
    results = pd.DataFrame({'names': list('ABCD'), 'logfoldchanges': [0.1, -0.1, 1.0, 0.0],
                            'fdr': [0.01, 0.01, 0.01, 0.5]})
    inputs = {'gene_results': results, 'gene_studies': [], 'selection': None,
              'project_folder': None}
    page = MetaVolcanoPage(SimpleNamespace(
        meta_ws=SimpleNamespace(figure_inputs=lambda: inputs), data_version=lambda: 0))

    def n_coloured():
        fig = page._make_render_func()()
        cols = fig.axes[0].collections   # not significant, down, up
        n = len(cols[1].get_offsets()) + len(cols[2].get_offsets())
        plt.close(fig)
        return n

    assert n_coloured() == 3                 # FDR alone
    # Without signals: a change would also start the page's background
    # re-render, which would draw at the same time as this test.
    page._controls.lfc_gate.blockSignals(True)
    page._controls.lfc_gate.setValue(0.25)
    assert n_coloured() == 1                 # FDR and |log2FC| > 0.25


def test_vote_levels_page_hatches_a_chosen_method_and_exports_the_counts(qapp, tmp_path):
    from kosmic.gui.figure_export.pages.meta_figures import VoteLevelsPage
    results = pd.DataFrame({'names': list('ABCD'), 'n_up': [3, 1, 0, 1],
                            'n_down': [0, 0, 0, 1], 'fdr_reml': [0.01, 0.5, 0.5, 0.5],
                            'fdr': [0.01, 0.5, 0.5, 0.5]})
    inputs = {'gene_results': results, 'gene_studies': [('S1', None), ('S2', None),
                                                         ('S3', None)],
              'selection': 'Cardiomyocyte', 'project_folder': None}
    page = VoteLevelsPage(SimpleNamespace(
        meta_ws=SimpleNamespace(figure_inputs=lambda: inputs), data_version=lambda: 0))
    assert page.dependencies_met()                # any run: levels need only Up / Down
    page._controls.populate_methods(page._methods())
    overlay = page._controls.overlay
    assert overlay.currentText() == 'None'        # off until chosen
    assert page._levels()[1] is None

    overlay.blockSignals(True)       # no background re-render alongside the test's
    overlay.setCurrentIndex(overlay.findData('reml'))
    assert page._levels()[1] == 'REML FDR < 0.05'
    plt.close(page._make_render_func()())
    (csv,) = page._export_companions(tmp_path / 'votes.png')
    counts = pd.read_csv(csv, index_col='level')
    assert counts.loc['3/3', ['up', 'total', 'significant']].tolist() == [1, 1, 1]
    assert counts.loc['conflict', 'total'] == 1


def test_a_page_discarded_while_its_preview_renders_does_not_abort(qapp):
    """The render thread was a Qt child of the page, so discarding the page
    mid-render (a test ending, a window closing) destroyed the running
    thread and Qt aborted the process; CI crashed this way."""
    import gc
    import time

    from kosmic.gui.figure_export.pages import _base
    from kosmic.gui.figure_export.pages.meta_figures import VoteLevelsPage
    results = pd.DataFrame({'names': ['A'], 'n_up': [1], 'n_down': [0], 'fdr': [0.01]})
    inputs = {'gene_results': results, 'gene_studies': [('S1', None)],
              'selection': 'x', 'project_folder': None}

    def open_and_discard():
        page = VoteLevelsPage(SimpleNamespace(
            meta_ws=SimpleNamespace(figure_inputs=lambda: inputs), data_version=lambda: 0))
        page._make_render_func = lambda: (lambda: time.sleep(0.5))   # a slow render
        page._controls.changed.emit()
        assert any(w.isRunning() for w in _base._LIVE_RENDERS)

    open_and_discard()
    gc.collect()
    _base.wait_for_renders()
    qapp.processEvents()             # the finished render reports to a page that is gone
    assert not _base._LIVE_RENDERS
