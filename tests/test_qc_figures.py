"""QC figures: the QC limits and filtering steps read from provenance, the
QC violin and cell-filtering plots, and the QC Violins page's options."""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

import matplotlib  # noqa: E402
matplotlib.use('Agg')

from kosmic import provenance  # noqa: E402
from kosmic.visualisation.scrna.qc_plots import (  # noqa: E402
    create_cell_filtering_plot, create_qc_violin_plots,
)


def _violins(ax):
    from matplotlib.collections import PolyCollection
    return [c for c in ax.collections if isinstance(c, PolyCollection)]


def _record(*stages, n_obs=None):
    return {'stages': [{'stage': s, 'params': p} for s, p in stages],
            'fingerprint': {'n_obs': n_obs} if n_obs is not None else None}


QC = {'min_genes': 300, 'max_genes': 5000, 'min_counts': 500, 'max_counts': 0, 'max_mt': 3}


def test_filtering_steps_follow_the_recorded_counts():
    rec = _record(('load', {'n_cells': 1000}),
                  ('filter_dataset', {'n_cells_before': 1000, 'n_cells_after': 800}),
                  ('qc', {**QC, 'n_cells_before': 800, 'n_cells_after': 750}),
                  ('doublet_filter', {'n_cells_before': 750, 'n_cells_after': 740}),
                  ('qc_recheck', {'n_cells_before': 740, 'n_cells_after': 730}))
    steps = provenance.cell_filtering_steps(rec)
    assert [(s['step'], s['n_cells'], s['n_removed']) for s in steps] == [
        ('Loaded', 1000, 0), ('Pre-analysis filter', 800, 200), ('QC limits', 750, 50),
        ('Doublet removal', 740, 10), ('QC limits on decontaminated counts', 730, 10)]


def test_an_unrecorded_doublet_removal_is_inferred_from_the_detection():
    # Before 'doublet_filter' existed only detection was recorded: a drop
    # equal to its n_doublets is the removal.
    rec = _record(('load', {'n_cells': 1000}),
                  ('qc', {**QC, 'n_cells_before': 1000, 'n_cells_after': 900}),
                  ('doublets', {'n_doublets': 30}),
                  ('propagate_labels', {'n_cells': 870}),
                  ('qc_recheck', {'n_cells_before': 870, 'n_cells_after': 860}))
    steps = provenance.cell_filtering_steps(rec)
    assert [(s['step'], s['n_removed']) for s in steps] == [
        ('Loaded', 0), ('QC limits', 100), ('Doublet removal', 30),
        ('QC limits on decontaminated counts', 10)]


def test_a_drop_no_stage_accounts_for_is_marked_not_recorded():
    rec = _record(('load', {'n_cells': 1000}),
                  ('qc', {**QC, 'n_cells_before': 1000, 'n_cells_after': 900}),
                  n_obs=850)
    assert provenance.cell_filtering_steps(rec)[-1] == {
        'step': 'Not recorded', 'n_cells': 850, 'n_removed': 50}


def test_a_qc_rerun_that_removed_nothing_is_not_a_step():
    rec = _record(('load', {'n_cells': 1000}),
                  ('qc', {**QC, 'n_cells_before': 1000, 'n_cells_after': 900}),
                  ('qc', {**QC, 'n_cells_before': 900, 'n_cells_after': 900}))
    assert len(provenance.cell_filtering_steps(rec)) == 2


def test_a_study_loaded_again_starts_its_steps_again():
    rec = _record(('load', {'n_cells': 1000}),
                  ('qc', {**QC, 'n_cells_before': 1000, 'n_cells_after': 600}),
                  ('load', {'n_cells': 1000}),
                  ('qc', {**QC, 'n_cells_before': 1000, 'n_cells_after': 950}))
    assert [(s['step'], s['n_cells']) for s in provenance.cell_filtering_steps(rec)] == [
        ('Loaded', 1000), ('QC limits', 950)]


def test_no_load_count_gives_no_steps():
    assert provenance.cell_filtering_steps(None) == []
    assert provenance.cell_filtering_steps(_record(('qc', QC))) == []


def test_qc_limits_take_the_last_qc_step_and_skip_unset_limits():
    rec = _record(('qc', {**QC, 'max_mt': 20}), ('qc', QC))
    assert provenance.qc_limits(rec) == {
        'n_genes_by_counts': {'min': 300.0, 'max': 5000.0},
        'total_counts': {'min': 500.0},
        'pct_counts_mt': {'max': 3.0}}
    assert provenance.qc_limits(None) == {}


@pytest.fixture
def obs():
    rng = np.random.default_rng(0)
    n = 600
    return pd.DataFrame({
        'n_genes_by_counts': rng.integers(300, 5000, n),
        'total_counts': rng.integers(500, 20000, n),
        'pct_counts_mt': rng.random(n) * 3,
        'decontX_contamination': rng.random(n) * 0.3,
        'sample': np.repeat([f's{i}' for i in range(6)], n // 6),
        '_role': np.repeat(['control', 'disease'], n // 2),
    })


@pytest.mark.parametrize('grouping', ['pooled', 'condition', 'sample'])
def test_the_violins_draw_one_panel_per_metric(obs, grouping):
    groups = None if grouping == 'pooled' else obs['_role' if grouping == 'condition' else 'sample']
    metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'decontX_contamination',
               'doublet_score']   # absent: left out
    fig = create_qc_violin_plots(obs, metrics=metrics, groups=groups, log_counts=True,
                                 thresholds={'pct_counts_mt': {'max': 3}})
    assert len(fig.axes) == 4
    n_violins = {'pooled': 1, 'condition': 2, 'sample': 6}[grouping]
    assert len(_violins(fig.axes[0])) == n_violins
    # Samples are stacked: the metrics share one x axis of sample names.
    if grouping == 'sample':
        assert [t.get_text() for t in fig.axes[-1].get_xticklabels()] == sorted(obs['sample'].unique())


def test_a_y_range_sets_the_axis_in_the_metrics_units(obs):
    fig = create_qc_violin_plots(obs, log_counts=True,
                                 y_limits={'total_counts': (100, 100_000),
                                           'pct_counts_mt': (None, 10)})
    genes, counts, mt = fig.axes
    # Counts are drawn in log10 units: 100 to 100,000 is 2 to 5.
    assert np.allclose(counts.get_ylim(), (2, 5))
    assert mt.get_ylim()[1] == 10 and mt.get_ylim()[0] < 1
    assert genes.get_ylim()[1] < np.log10(5000) + 0.5   # left on Auto


def test_log_ticks_are_evenly_spaced_steps(obs):
    fig = create_qc_violin_plots(obs, log_counts=True,
                                 y_limits={'total_counts': (400, 50_000)})
    labels = [t.get_text() for t in fig.axes[1].get_yticklabels()]
    assert labels == ['1,000', '3,000', '10,000', '30,000']


def test_the_violins_return_none_without_qc_metrics():
    assert create_qc_violin_plots(pd.DataFrame({'x': [1, 2]})) is None


def test_the_filtering_flow_chart_has_a_box_per_count_and_per_step():
    steps = [{'step': 'Loaded', 'n_cells': 1000, 'n_removed': 0},
             {'step': 'QC limits', 'n_cells': 900, 'n_removed': 100},
             {'step': 'Doublet removal', 'n_cells': 890, 'n_removed': 10}]
    fig = create_cell_filtering_plot(steps)
    # Three counts down the left, and a box for each of the two steps.
    assert len(fig.axes[0].patches) == 5
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert 'QC limits\n100 removed' in texts and '890 cells' in texts
    assert create_cell_filtering_plot([]) is None


@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_the_pages_read_the_studys_provenance(app, tmp_path):
    from kosmic.gui.figure_export.pages.qc_violins import (
        GROUP_POOLED, GROUP_SAMPLE, START_FILTERED, CellFilteringPage, QCViolinsPage,
    )
    from kosmic.scrna.load.h5ad_meta import read_figure_meta
    from tests.study_folder import figure_workspace, make_study_folder

    study = make_study_folder(tmp_path)
    adata = read_figure_meta(next((study / 'processed_data').glob('*.h5ad')))
    ws = figure_workspace(scrna_adata=adata, study=study)

    violins = QCViolinsPage(ws)
    violins._controls.blockSignals(True)
    size = violins._controls.figsize
    assert size.get_figsize() == (12, 6)
    violins._controls.group.setCurrentText(GROUP_SAMPLE)
    assert size.get_figsize() == (12, 10)        # stacked rows get more height
    violins._controls.group.setCurrentText(GROUP_POOLED)
    size.h_spin.setValue(8)                      # a height set by hand is kept
    violins._controls.group.setCurrentText(GROUP_SAMPLE)
    assert size.get_figsize() == (12, 8)
    assert violins.dependencies_met()
    fig = violins._make_render_func()()
    assert len(_violins(fig.axes[0])) == adata.obs['sample'].nunique()
    # The fixture's QC step set max_mt 4: one dotted line on the MT panel.
    assert any(line.get_linestyle() == ':' for line in fig.axes[2].get_lines())
    assert violins._controls.y_limits() == {}
    violins._controls.y_range['pct_counts_mt'][1].setText('5')
    assert violins._controls.y_limits() == {'pct_counts_mt': (None, 5.0)}
    assert violins._make_render_func()().axes[2].get_ylim()[1] == 5
    violins._controls.y_range['pct_counts_mt'][0].setText('0')    # 0 is a limit, not Auto
    assert violins._controls.y_limits() == {'pct_counts_mt': (0.0, 5.0)}
    assert violins._make_render_func()().axes[2].get_ylim() == (0, 5)
    violins._controls.y_range['total_counts'][1].setText('50,000')
    violins._controls.y_range['n_genes_by_counts'][0].setText('abc')    # not a number: Auto
    assert violins._controls.y_limits() == {'pct_counts_mt': (0.0, 5.0),
                                            'total_counts': (None, 50000.0)}
    # The DecontX and doublet ranges are used only with their panels.
    assert not violins._controls.y_range['doublet_score'][0].isEnabled()
    table = violins._export_companions(tmp_path / 'qc.png')[0]
    assert pd.read_csv(table)['cells'].sum() == adata.n_obs

    filtering = CellFilteringPage(ws)
    filtering._controls.blockSignals(True)
    assert filtering.dependencies_met()
    assert [s['n_removed'] for s in filtering._steps()] == [0, 20]
    filtering._controls.start.setCurrentText(START_FILTERED)
    assert len(filtering._steps()) == 2     # no pre-analysis filter recorded
