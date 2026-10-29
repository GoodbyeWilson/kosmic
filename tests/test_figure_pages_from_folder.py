"""Stage 4 of ADR-014: the scRNA and DE figure pages draw a study from its
folder when no workspace holds it, and draw the same figure either way.
"""
from __future__ import annotations

import os
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

import matplotlib  # noqa: E402
matplotlib.use('Agg')

from kosmic.de import saved  # noqa: E402
from kosmic.paths import de_figures_dir, de_pathway_scoring_dir, scrna_figures_dir  # noqa: E402
from kosmic.scrna.load import h5ad_meta  # noqa: E402
from tests.study_folder import (  # noqa: E402
    figure_workspace, make_study_folder, write_de_outputs,
)

GENE_SETS = {'Immune': ['PTPRC', 'CD3E', 'CD14', 'LYZ'],
             'Vascular': ['PECAM1', 'VWF', 'ACTA2', 'MYH11']}


@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _complete(study, *, counts_source='raw', details=True):
    """Add what the DE steps of this change save: the first run's record
    (with its gene sets), the ORA table and the GSEA curves."""
    run = saved.list_gene_de_runs(study)[0]
    saved.save_run_record(run.path, {
        'method': run.method, 'cell_type': run.cell_type, 'counts_source': counts_source,
        'sample_column': 'sample', 'condition_column': 'condition',
        'disease_label': 'Disease', 'control_label': 'Control', 'min_cells': 3,
        'pathway_gene_sets': GENE_SETS, 'geneset_label': '2 gene sets'})
    scoring = de_pathway_scoring_dir(study)
    saved.save_enrichment(scoring, pd.DataFrame({
        'Term': ['immune response', 'angiogenesis'], 'FDR': [0.001, 0.2],
        'Fold_Enrichment': [3.0, 1.2], 'Gene_Count': [3, 2],
        'Genes': ['PTPRC, CD3E, CD14', 'VWF, PECAM1']}))
    if details:
        saved.save_gsea_details(saved.gsea_library_path(study), {
            t: {'es': e, 'nes': n, 'RES': np.linspace(0, e, 50).tolist(), 'hits': [3, 10, 20]}
            for t, e, n in (('HALLMARK_A', 0.5, 1.8), ('HALLMARK_B', -0.3, -1.2))},
            'Hallmark')
    return run


@pytest.fixture
def study(tmp_path):
    s = write_de_outputs(make_study_folder(tmp_path, form='log'))
    _complete(s)
    return s


def _workspace(app, study):
    """Figure Export with no dataset in scRNA or DE, on 'study', once its
    file's metadata has been read."""
    from kosmic.gui.figure_export.workspace import FigureExportWorkspace
    ws = FigureExportWorkspace()
    ws.set_study(study)
    ws.source.scrna_adata()
    ws.source.wait()
    app.processEvents()
    return ws


def _draw(app, page):
    from kosmic.gui.figure_export.pages._base import wait_for_renders
    page.on_activated()
    wait_for_renders()
    app.processEvents()
    return page._preview._pixmap is not None


def _export(page, path):
    from PyQt6.QtWidgets import QFileDialog
    with mock.patch.object(QFileDialog, 'getSaveFileName', return_value=(str(path), '')):
        page._export_current()
    assert page._status_lbl.text().startswith('Exported'), page._status_lbl.text()


# -- Every figure from the folder ---------------------------------------------------

def test_every_scrna_and_de_figure_draws_from_the_folder(app, study):
    ws = _workspace(app, study)
    assert ws.scrna_ws is None and ws.de_ws is None
    drawn, unavailable = [], []
    for section, _group, page in ws._entries:
        if section == 'meta':
            continue
        if not ws._available(page):
            unavailable.append(page.TITLE)
            continue
        assert _draw(app, page), (page.TITLE, page._status_lbl.text())
        drawn.append(page.TITLE)
    assert unavailable == [], unavailable
    assert len(drawn) == 19


def test_a_figure_whose_input_was_not_saved_says_what_to_run(app, tmp_path):
    study = write_de_outputs(make_study_folder(tmp_path))
    _complete(study, details=False)                       # GSEA run before curves were kept
    ws = _workspace(app, study)
    from kosmic.gui.figure_export.pages.library_gsea_mountain import LibraryGseaMountainPage
    page = next(p for p in ws._pages if isinstance(p, LibraryGseaMountainPage))
    assert not ws._available(page)
    assert 'Run GSEA again' in page.unavailable_reason()


def test_the_scrna_figures_say_why_while_the_file_is_not_read(app, study, tmp_path,
                                                              monkeypatch):
    import threading
    from kosmic.gui.figure_export.workspace import FigureExportWorkspace
    from kosmic.gui.figure_export.pages.umap import UMAPPage
    release = threading.Event()
    real = h5ad_meta.read_figure_meta
    monkeypatch.setattr(h5ad_meta, 'read_figure_meta',
                        lambda p: (release.wait(10), real(p))[1])
    ws = FigureExportWorkspace()
    page = next(p for p in ws._pages if isinstance(p, UMAPPage))
    ws.set_study(study)
    assert not ws._available(page)
    assert page.unavailable_reason().startswith('Reading StudyA.h5ad')
    release.set()
    ws.source.wait()
    app.processEvents()
    assert ws._available(page)
    ws.set_study(tmp_path / 'Renamed')
    assert page.unavailable_reason() == f"The study folder {tmp_path / 'Renamed'} cannot be found."


# -- The same figure either way -------------------------------------------------------

@pytest.mark.parametrize('form', ['log', 'counts'])
def test_the_marker_dot_plot_numbers_are_the_same_held_or_from_the_file(app, tmp_path, form):
    import anndata as ad
    from kosmic.gui.figure_export.pages.cell_types import MarkerDotPlotPage
    study = make_study_folder(tmp_path, form=form)
    held = ad.read_h5ad(h5ad_meta.registered_h5ad(study))
    out = {}
    for name, ws in (('held', figure_workspace(scrna_adata=held, study=study)),
                     ('file', figure_workspace(study=study))):
        ws.source.scrna_adata()
        ws.source.wait()
        page = MarkerDotPlotPage(ws)
        page._controls.source.blockSignals(True)
        page._controls.source.setCurrentText('General (built-in)')
        page._controls.source.blockSignals(False)
        assert _draw(app, page), page._status_lbl.text()
        _export(page, tmp_path / f'{name}.png')
        out[name] = pd.read_csv(tmp_path / f'{name}.csv')
        out[name + '_label'] = page._cache[next(iter(page._cache))][3]
    pd.testing.assert_frame_equal(out['held'], out['file'])
    assert out['held_label'] == out['file_label'] == 'Mean expression\nlog1p(CP10K)'


def _heatmap_page(ws):
    from kosmic.gui.figure_export.pages.gene_heatmap import GeneHeatmapPage
    return GeneHeatmapPage(ws)


def test_the_heatmap_table_is_the_same_held_or_from_the_file(app, study):
    run = saved.list_gene_de_runs(study)[0]
    h5ad = h5ad_meta.registered_h5ad(study)
    de = SimpleNamespace(
        current_adata=h5ad_meta.read_counts_adata(h5ad), de_results=saved.read_gene_de(run),
        pathway_gene_sets=GENE_SETS, sample_col='sample', condition_col='condition',
        min_cells=3, disease_label='Disease', control_label='Control', h5ad_path=h5ad)
    held = _heatmap_page(figure_workspace(de=de, study=study))
    from_file = _heatmap_page(figure_workspace(study=study))
    with mock.patch.object(h5ad_meta, 'read_genes', wraps=h5ad_meta.read_genes) as reads:
        sd_held, avail_held = held._prepare(held._inputs())
        assert reads.call_count == 0                     # the DE workspace's counts
        sd_file, avail_file = from_file._prepare(from_file._inputs())
        assert reads.call_count == 1                     # only the gene sets' genes
    assert set(reads.call_args.args[1]) <= {g for gs in GENE_SETS.values() for g in gs}
    pd.testing.assert_frame_equal(sd_held, sd_file)
    assert avail_held == avail_file and avail_file


def test_the_heatmap_from_the_file_uses_the_runs_counts(app, tmp_path):
    study = write_de_outputs(make_study_folder(tmp_path))
    _complete(study, counts_source='decontX_counts')
    page = _heatmap_page(figure_workspace(study=study))
    sd, _avail = page._prepare(page._inputs())
    from kosmic.visualisation.de.gene_heatmap import prepare_heatmap_data
    want, _ = prepare_heatmap_data(
        h5ad_meta.read_counts_adata(h5ad_meta.registered_h5ad(study),
                                    counts_layer='decontX_counts'),
        GENE_SETS, 'sample', 'condition', min_cells=3,
        allowed_genes=set(saved.read_gene_de(saved.list_gene_de_runs(study)[0])['names']))
    pd.testing.assert_frame_equal(sd, want)


# -- The workspace ------------------------------------------------------------------

def test_the_run_picker_chooses_the_de_table(app, study):
    from kosmic.gui.figure_export.pages.volcano import VolcanoPage
    ws = _workspace(app, study)
    ws.show_section('de')
    runs = saved.list_gene_de_runs(study)
    combo = ws._run_combo
    assert [combo.itemText(i) for i in range(combo.count())] == [r.label for r in runs]
    assert not ws._run_row.isHidden()
    page = next(p for p in ws._pages if isinstance(p, VolcanoPage))
    pd.testing.assert_frame_equal(page._de().de_results, saved.read_gene_de(runs[0]))
    before = ws.data_version()
    combo.setCurrentIndex(2)
    assert ws.source.de_run() == runs[2] and ws.data_version() > before
    pd.testing.assert_frame_equal(page._de().de_results, saved.read_gene_de(runs[2]))
    ws.show_section('scrna')
    assert ws._run_row.isHidden()


def test_the_run_picker_offers_the_de_workspaces_run_first(app, study):
    from kosmic.gui.figure_export.workspace import FigureExportWorkspace
    live = pd.DataFrame({'names': ['VWF'], 'logfoldchanges': [1.0], 'pvals_adj': [0.01]})
    ws = FigureExportWorkspace(de_workspace=SimpleNamespace(project_dir=study, de_results=live))
    ws.set_study(study)
    ws.show_section('de')
    assert ws._run_combo.itemText(0) == 'Run in the DE workspace'
    assert ws._run_combo.count() == 1 + len(saved.list_gene_de_runs(study))


def test_filling_the_list_reads_no_table(app, study, monkeypatch):
    reads = []
    for name in ('read_gene_de', 'read_pathway_de', 'read_per_donor_scores', 'read_fgsea',
                 'read_gsea_library', 'read_gsea_details', 'read_enrichment'):
        monkeypatch.setattr(saved, name, lambda *a, n=name, **k: reads.append(n))
    for name in ('read_genes', 'read_cells', 'read_var', 'read_counts_adata'):
        monkeypatch.setattr(h5ad_meta, name, lambda *a, n=name, **k: reads.append(n))
    ws = _workspace(app, study)
    # Opening the figure a section selects reads its data; only the list
    # (availability and tooltips for every figure) is measured here.
    ws._on_selection_changed = lambda: None
    for section in ('scrna', 'de', 'meta'):
        ws.show_section(section)
        ws._rebuild_list()
    assert reads == []


def test_results_older_than_the_dataset_carry_a_note(app, study):
    from kosmic.gui.figure_export.pages.volcano import VolcanoPage
    ws = _workspace(app, study)
    page = next(p for p in ws._pages if isinstance(p, VolcanoPage))
    run = saved.list_gene_de_runs(study)[0]
    h5ad = h5ad_meta.registered_h5ad(study)
    t = h5ad.stat().st_mtime
    os.utime(run.path, (t + 60, t + 60))
    page._update_staleness()
    assert page._staleness_banner.text() == ''
    os.utime(run.path, (t - 60, t - 60))
    page._update_staleness()
    assert page._staleness_banner.text().startswith(
        "These results were made before the study's file was last saved")


def test_figures_export_into_the_studys_folder_whoever_holds_another(app, study, tmp_path):
    import anndata as ad
    from kosmic.gui.figure_export.workspace import FigureExportWorkspace
    from kosmic.gui.figure_export.pages.umap import UMAPPage
    from kosmic.gui.figure_export.pages.volcano import VolcanoPage
    other = SimpleNamespace(current_adata=ad.AnnData(np.zeros((2, 2))),
                            current_project_dir=str(tmp_path / 'Other'))
    ws = FigureExportWorkspace(scrna_workspace=other,
                               de_workspace=SimpleNamespace(project_dir=tmp_path / 'Other'))
    ws.set_study(study)
    umap = next(p for p in ws._pages if isinstance(p, UMAPPage))
    volcano = next(p for p in ws._pages if isinstance(p, VolcanoPage))
    assert umap._default_export_dir() == scrna_figures_dir(study)
    assert volcano._default_export_dir() == de_figures_dir(study)
    ws.source.wait()
