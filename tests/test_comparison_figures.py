"""The Methods Comparison Venn diagram and its Figure Export page (ADR-017 §6)."""
import json
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from kosmic.comparison.agreement import Rule, region_counts, venn_regions  # noqa: E402
from kosmic.comparison.bundle import ComparisonBundle  # noqa: E402
from kosmic.comparison.standard import COLUMNS  # noqa: E402

SETS = {'pb_merged': {'A', 'B', 'C', 'D'}, 'ma_reml': {'A', 'B', 'E'}, 'ma_vote': {'A'}}


def test_region_counts_cover_every_region():
    assert region_counts(SETS) == {(1, 0, 0): 2, (0, 1, 0): 1, (0, 0, 1): 0,
                                   (1, 1, 0): 1, (1, 0, 1): 0, (0, 1, 1): 0, (1, 1, 1): 1}
    assert region_counts({'x': {'A'}, 'y': {'A', 'B'}}) == {(1, 0): 0, (0, 1): 1, (1, 1): 1}


def test_venn_regions_name_the_methods_that_call_each_gene():
    t = venn_regions(SETS).set_index('gene')
    assert t.loc['A', 'region'] == 'pb_merged + ma_reml + ma_vote'
    assert t.loc['E', 'region'] == 'ma_reml' and not t.loc['E', 'pb_merged']
    assert list(venn_regions(SETS, universe={'A', 'E'})['gene']) == ['A', 'E']


FOUR = {'a': set('ABCDEFGH'), 'b': set('ABCDE'), 'c': set('ABXY'), 'd': set('AZ')}


@pytest.mark.parametrize('n', [2, 3, 4])
def test_the_venn_plot_writes_each_non_empty_region_and_a_key(n):
    from kosmic.visualisation.comparison.venn import create_venn_plot
    sets = dict(list(FOUR.items())[:n])
    fig = create_venn_plot(sets, {m: m.upper() for m in sets}, n_compared=20, title='T',
                           rule_text='adjusted p < 0.05')
    texts = sorted(t.get_text() for t in fig.axes[0].texts)
    regions = [v for v in region_counts(sets).values() if v]
    assert len(texts) <= len(regions)            # a region too thin to draw has no label
    assert set(texts) <= {f"{v:,}" for v in regions}
    legend = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert legend[0] == 'A (8 genes)' and len(legend) == n
    note = next(t.get_text() for t in fig.texts if t.get_text().startswith('Called'))
    assert 'Called: adjusted p < 0.05.' in note and '20 genes compared' in note
    assert not any('$' in t for t in texts + legend)


def test_the_shapes_are_drawn_in_proportion_to_the_genes_called():
    from kosmic.visualisation.venn_layout import venn_layout
    geom, _ = venn_layout({'big': set(range(400)), 'small': set(range(100))})
    (_, _, r_big, _, _), (_, _, r_small, _, _) = geom
    assert (r_big / r_small) ** 2 == pytest.approx(4.0)


def test_each_count_is_written_inside_its_own_region():
    # A crescent's mean point lies in the circle it wraps; the label must not.
    from kosmic.visualisation.comparison.venn import _inside, region_positions
    from kosmic.visualisation.venn_layout import venn_layout
    geom, _ = venn_layout({'outer': set(range(1000)), 'inner': set(range(500))})
    x, y = region_positions(geom, 2)[(1, 0)]
    inside = _inside(geom, np.array([x]), np.array([y]))[:, 0]
    assert inside.tolist() == [True, False]


def test_the_venn_plot_needs_two_to_four_sets():
    from kosmic.visualisation.comparison.venn import create_venn_plot
    assert create_venn_plot({'a': {'A'}}) is None
    assert create_venn_plot({k: {'A'} for k in 'abcde'}) is None


def _bundle(folder=None):
    rows = []
    for m, called in SETS.items():
        for i, g in enumerate('ABCDE'):
            rows.append({'gene': g, 'method': m, 'effect': 1.0, 'direction': 1, 'stat': 0.0,
                         'pvalue': 0.01 if g in called else 0.5,
                         'padj': (0.01 if g in called else 0.5) if m != 'ma_vote' else np.nan,
                         'tested': True, 'called': g in called, 'rank': i + 1})
    b = ComparisonBundle.build(pd.DataFrame(rows)[COLUMNS], Rule(),
                               {'cell_type': 'Ventricular Cardiomyocyte'})
    b.folder = folder
    return b


@pytest.fixture(scope='module')
def qapp():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _workspace(bundle=None):
    from kosmic.gui.figure_export.workspace import FigureExportWorkspace
    ws = FigureExportWorkspace()
    if bundle is not None:
        ws.set_comparison_workspace(SimpleNamespace(
            figure_inputs=lambda: {'bundle': bundle, 'project_folder': None}))
    return ws


def test_the_comparison_section_lists_the_venn_greyed_without_a_comparison(qapp):
    from kosmic.gui.figure_export.pages.methods_comparison import MethodsVennPage
    ws = _workspace()
    ws.show_section('compare')
    page = next(p for p in ws._pages if isinstance(p, MethodsVennPage))
    assert not page.dependencies_met()
    assert 'Methods Comparison workspace' in page.unavailable_reason()


def test_the_venn_page_draws_the_chosen_methods_and_exports_its_regions(qapp, tmp_path):
    from kosmic.gui.figure_export.pages import _base
    from kosmic.gui.figure_export.pages.methods_comparison import MethodsVennPage
    ws = _workspace(_bundle(tmp_path))
    ws.show_section('compare')
    page = next(p for p in ws._pages if isinstance(p, MethodsVennPage))
    assert page.dependencies_met()
    page._controls.populate(page._bundle().methods)
    # The defaults present in this comparison (it has no Stouffer).
    assert page._controls.chosen() == ['pb_merged', 'ma_reml', 'ma_vote']
    for box in page._controls.method_boxes[2:]:
        box.blockSignals(True)
        box.setCurrentIndex(box.findData(None))    # two methods
        box.blockSignals(False)
    fig = page._make_render_func()()
    assert len(fig.axes[0].get_legend().get_texts()) == 2
    assert page._default_export_dir() == tmp_path
    csv, params = page._export_companions(tmp_path / 'methods_venn.pdf')
    regions = pd.read_csv(csv)
    assert set(regions['region']) == {'Merged pseudobulk + REML', 'Merged pseudobulk', 'REML'}
    saved = json.loads(params.read_text())
    assert saved['methods_shown'] == ['pb_merged', 'ma_reml']
    assert saved['rule']['rule'] == 'significance'
    _base.wait_for_renders()


def test_figures_reads_a_saved_comparison_and_offers_a_choice(qapp, tmp_path):
    first = _bundle()
    first.save(tmp_path, 'decontX_counts/Ventricular_Cardiomyocyte')
    second = _bundle()
    second.save(tmp_path, 'decontX_counts/Fibroblast')
    ws = _workspace()
    ws.set_project(tmp_path)
    ws.show_section('compare')
    assert ws.comparison_inputs()['bundle'] is not None
    labels = [ws._run_combo.itemText(i) for i in range(ws._run_combo.count())]
    assert set(labels) == {'Ventricular Cardiomyocyte, decontX counts',
                           'Fibroblast, decontX counts'}
    assert ws._run_row.isVisibleTo(ws)
    ws._run_combo.setCurrentIndex(labels.index('Fibroblast, decontX counts'))
    assert ws.comparison_inputs()['bundle'].folder.name == 'Fibroblast'
