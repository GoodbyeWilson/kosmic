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
