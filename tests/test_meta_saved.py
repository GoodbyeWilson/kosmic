"""Stage 5b of ADR-014: the Meta figures draw a meta-analysis run saved in
the project folder when the Meta-Analysis workspace holds none, and a
copied project (without raw_data/) exports every figure.
"""
from __future__ import annotations

import json
import os
import shutil
import time

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

import matplotlib  # noqa: E402
matplotlib.use('Agg')

from kosmic.de import saved as de_saved  # noqa: E402
from kosmic.meta_analysis import saved  # noqa: E402
from kosmic.meta_analysis.io import discover_de_results, load_study_de  # noqa: E402
from kosmic.paths import meta_output_dir  # noqa: E402
from tests.study_folder import (  # noqa: E402
    complete_de_outputs, make_study_folder, write_de_outputs,
)

STUDIES = ('StudyA', 'StudyB', 'StudyC')
SELECTION = 'Endothelial'


@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _pooled(genes, *, seed=0):
    """A pooled table with the columns a gene-level run saves."""
    rng = np.random.default_rng(seed)
    lfc = rng.normal(size=len(genes))
    se = rng.random(len(genes)) * 0.2 + 0.1
    p = rng.random(len(genes)) ** 4
    return pd.DataFrame({
        'names': genes, 'logfoldchanges': lfc, 'se': se,
        'ci_lower': lfc - 1.96 * se, 'ci_upper': lfc + 1.96 * se,
        'pvals_pooled': p, 'pvals_reml': p, 'fdr': np.minimum(p * 3, 1.0),
        'fdr_reml': np.minimum(p * 3, 1.0), 'tau_squared': 0.02,
        'heterogeneity_i2': 0.3, 'n_studies': 3, 'pooling_method': 'consensus',
        'n_up': rng.integers(0, 4, len(genes)), 'n_down': 0})


def _write_run(folder, stem, table, params):
    folder.mkdir(parents=True, exist_ok=True)
    table.to_csv(folder / f'{stem}.csv', index=False)
    (folder / f'{stem}.json').write_text(json.dumps({'parameters': params}), encoding='utf-8')


def make_project(root, *, gsea_details=True):
    """Three studies with DE outputs and the meta-analysis runs a session
    saves: two gene-level runs and a pathway run of one selection, the
    Enrichment step's ORA and library GSEA, and an older run of 'mixed'."""
    project = root / 'Project'
    for name in STUDIES:
        complete_de_outputs(write_de_outputs(make_study_folder(project, name=name)))
    found = discover_de_results(project)
    gene_names = [f'{s}_Endothelial' for s in STUDIES]
    pathway_names = sorted({e['accession'] for e in found if e['file_type'] == 'pathway_de'})
    genes = list(pd.read_csv(found[0]['path'])['names'])
    sel = meta_output_dir(project, SELECTION)
    _write_run(meta_output_dir(project, 'mixed'), 'consensus_reml_analytical',
               _pooled(genes, seed=5), {'de_method': 'deseq2', 'study_names': gene_names})
    time.sleep(0.02)
    _write_run(sel, 'consensus_reml_cc', _pooled(genes, seed=2),
               {'de_method': 'deseq2', 'study_names': gene_names})
    time.sleep(0.02)
    _write_run(sel, 'consensus_reml_fisher_vote_count_analytical', _pooled(genes),
               {'de_method': 'deseq2', 'study_names': gene_names})
    _write_run(sel, 'pathway_ma_reml', _pooled(['Glycolysis', 'OXPHOS'], seed=3),
               {'study_names': pathway_names})
    de_saved.save_enrichment(sel, pd.DataFrame({
        'Term': ['angiogenesis', 'other'], 'FDR': [0.001, 0.2], 'Fold_Enrichment': [3.0, 1.2],
        'Gene_Count': [2, 1], 'Genes': ['VWF, PECAM1', 'LYZ']}))
    gsea = sel / 'gsea_library_results_Hallmark.csv'
    pd.DataFrame({'names': ['HALLMARK_A', 'HALLMARK_B'], 'nes': [1.8, -1.2], 'es': [0.5, -0.3],
                  'pvals': [0.001, 0.04], 'pvals_adj': [0.01, 0.2],
                  'leading_edge': ['VWF;PECAM1', 'LYZ']}).to_csv(gsea, index=False)
    if gsea_details:
        de_saved.save_gsea_details(gsea, {
            t: {'es': e, 'nes': n, 'RES': np.linspace(0, e, 50).tolist(), 'hits': [3, 10, 20]}
            for t, e, n in (('HALLMARK_A', 0.5, 1.8), ('HALLMARK_B', -0.3, -1.2))}, 'Hallmark')
    return project


@pytest.fixture
def project(tmp_path):
    return make_project(tmp_path)


# -- Listing and reading runs ----------------------------------------------------------

def test_runs_are_listed_by_selection_newest_first(project):
    runs = saved.list_meta_runs(project)
    assert [(r.selection, r.gene_path.name) for r in runs] == [
        (SELECTION, 'consensus_reml_fisher_vote_count_analytical.csv'),
        (SELECTION, 'consensus_reml_cc.csv'),
        ('mixed', 'consensus_reml_analytical.csv')]
    assert runs[0].pathway_path.name == 'pathway_ma_reml.csv'
    assert runs[2].pathway_path is None
    assert runs[0].label == 'Endothelial: reml + fisher + vote count (analytical)'
    assert runs[1].label == 'Endothelial: reml (calibrated)'


def test_a_saved_run_gives_the_figure_inputs_keys(project, app):
    from kosmic.gui.meta_analysis.workspace import MetaAnalysisWorkspace
    run = saved.list_meta_runs(project)[0]
    inputs = saved.figure_inputs(run)
    held = MetaAnalysisWorkspace().figure_inputs()
    assert set(held) <= set(inputs)
    pd.testing.assert_frame_equal(inputs['gene_results'], pd.read_csv(run.gene_path))
    assert inputs['selection'] == SELECTION and inputs['project_folder'] == str(project)
    assert [n for n, _ in inputs['gene_studies']] == [f'{s}_Endothelial' for s in STUDIES]
    assert len(inputs['pathway_studies']) == 3 and inputs['missing_gene_studies'] == []
    assert list(inputs['enrichment_results']['Term']) == ['angiogenesis', 'other']
    assert inputs['gsea_library_name'] == 'Hallmark'
    assert set(inputs['gsea_library_details']) == {'HALLMARK_A', 'HALLMARK_B'}
    assert inputs['enrichment_ranking'] is inputs['gene_results']


def test_availability_reads_no_table(project, monkeypatch):
    reads = []
    monkeypatch.setattr(pd, 'read_csv', lambda *a, **k: reads.append(a))
    inputs = saved.figure_inputs(saved.list_meta_runs(project)[0])
    assert all(inputs.available(k) for k in (
        'gene_results', 'pathway_results', 'gene_studies', 'enrichment_results',
        'gsea_library_results', 'gsea_library_details'))
    assert inputs['missing_gene_studies'] == []       # from the file listing
    assert reads == []
    assert not saved.figure_inputs(saved.list_meta_runs(project)[2]).available('pathway_results')


def test_study_tables_are_the_meta_workspaces(project, app):
    """The forest plot reads the same numbers from a saved run as from the
    workspace: both load through load_study_de."""
    from kosmic.gui.meta_analysis.workspace import MetaAnalysisWorkspace
    ws = MetaAnalysisWorkspace()
    ws._project_folder = str(project)       # not set_project_directory_external: no QSettings
    ws._rescan_project()
    held = {d['name']: d['df'] for d in ws._all_gene_datasets}
    names = [f'{s}_Endothelial' for s in STUDIES]
    tables, missing = saved.study_tables(project, names, de_method='deseq2')
    assert missing == []
    for name, df in tables:
        pd.testing.assert_frame_equal(df, held[name])


def test_the_moved_loader_uses_the_moderated_se(tmp_path):
    df = pd.DataFrame({'names': ['A', 'B', 'C'], 'logfoldchanges': [1.0, -1.0, 0.5],
                       'pvals': 0.01, 'pvals_adj': 0.02, 'se': [0.5, 0.5, 0.5],
                       'se_moderated': [0.2, np.nan, 0.0]})
    path = tmp_path / 'S_DE_welch_cpm_eb.csv'
    df.to_csv(path, index=False)
    assert list(load_study_de(path, 'welch_cpm_eb')['se']) == [0.2, 0.5, 0.5]
    assert list(load_study_de(path, 'deseq2')['se']) == [0.5, 0.5, 0.5]


def test_a_study_without_a_table_is_reported(project):
    shutil.rmtree(project / 'StudyB' / 'results')
    tables, missing = saved.study_tables(
        project, [f'{s}_Endothelial' for s in STUDIES], de_method='deseq2')
    assert missing == ['StudyB_Endothelial'] and len(tables) == 2


# -- Figure Export -----------------------------------------------------------------

def _workspace(app, project, study='StudyA'):
    from kosmic.gui.figure_export.workspace import FigureExportWorkspace
    ws = FigureExportWorkspace()
    ws.set_project(project)
    ws.set_study(project / study)
    ws.source.scrna_adata()
    ws.source.wait()
    app.processEvents()
    return ws


def _draw(app, page):
    from kosmic.gui.figure_export.pages._base import wait_for_renders
    if hasattr(page._controls, 'populate_items'):
        df = page._results()
        page._controls.populate_items(list(df.sort_values('fdr')['names'].astype(str)))
    page.on_activated()
    wait_for_renders()
    app.processEvents()
    return page._preview._pixmap is not None


def test_every_meta_figure_draws_a_saved_run(app, project):
    ws = _workspace(app, project)
    assert ws.meta_ws is None
    ws.show_section('meta')
    for section, _group, page in ws._entries:
        if section != 'meta':
            continue
        assert ws._available(page), (page.TITLE, page.unavailable_reason())
        assert _draw(app, page), (page.TITLE, page._status_lbl.text())
    # The enrichment figures are listed there too, drawing the saved run.
    from kosmic.gui.figure_export.pages.library_gsea_mountain import LibraryGseaMountainPage
    assert any(isinstance(p, LibraryGseaMountainPage)
               for section, _g, p in ws._entries if section == 'meta')


def test_the_run_picker_chooses_the_meta_run(app, project):
    from kosmic.gui.figure_export.pages.meta_figures import MetaVolcanoPage
    ws = _workspace(app, project)
    ws.show_section('meta')
    runs = saved.list_meta_runs(project)
    combo = ws._run_combo
    assert [combo.itemText(i) for i in range(combo.count())] == [r.label for r in runs]
    assert not ws._run_row.isHidden()
    page = next(p for p in ws._pages if isinstance(p, MetaVolcanoPage))
    pd.testing.assert_frame_equal(page._results(), pd.read_csv(runs[0].gene_path))
    before = ws.data_version()
    combo.setCurrentIndex(2)
    assert ws.data_version() > before
    pd.testing.assert_frame_equal(page._results(), pd.read_csv(runs[2].gene_path))
    assert page._default_export_dir() == meta_output_dir(project, 'mixed')


def test_the_meta_workspaces_run_comes_first(app, project):
    from types import SimpleNamespace
    held = {'gene_results': _pooled(['A', 'B']),
            'pathway_results': None, 'gene_studies': [], 'pathway_studies': [],
            'selection': 'x', 'project_folder': str(project)}
    ws = _workspace(app, project)
    ws.set_meta_workspace(SimpleNamespace(figure_inputs=lambda: held))
    ws.show_section('meta')
    assert ws._run_combo.itemText(0) == 'Run in the Meta-Analysis workspace'
    assert ws.meta_inputs() is held
    ws._run_combo.setCurrentIndex(1)
    assert ws.meta_inputs() is not held


def test_a_missing_study_greys_only_the_forest_figures(app, project):
    from kosmic.gui.figure_export.pages.meta_figures import (
        LeaveOneOutPage, MetaVolcanoPage, StudyForestPage,
    )
    shutil.rmtree(project / 'StudyB' / 'results' / 'de_analysis' / 'statistics')
    ws = _workspace(app, project)
    pages = {type(p): p for p in ws._pages}
    assert ws._available(pages[MetaVolcanoPage])
    for cls in (StudyForestPage, LeaveOneOutPage):
        assert not ws._available(pages[cls])
        assert 'StudyB_Endothelial' in pages[cls].unavailable_reason()


def test_gsea_curves_saved_before_they_were_kept_grey_the_mountain(app, tmp_path):
    project = make_project(tmp_path, gsea_details=False)
    inputs = saved.figure_inputs(saved.list_meta_runs(project)[0])
    assert inputs.available('gsea_library_results')
    assert not inputs.available('gsea_library_details')
    assert inputs['gsea_library_details'] is None and inputs['gsea_library_name'] == 'Hallmark'


def test_a_new_study_redraws_although_an_id_is_reused(app, project):
    """Previews are cached against data_version: a study change bumps it,
    so a preview of the last study is not reused when a new object happens
    to get the old one's id."""
    from kosmic.gui.figure_export.pages.umap import UMAPPage
    ws = _workspace(app, project)
    page = next(p for p in ws._pages if isinstance(p, UMAPPage))
    assert _draw(app, page)
    version = ws.data_version()
    assert page._cached_data_version == version
    ws.set_study(project / 'StudyB')
    assert ws.data_version() > version and page._cached_pixmap is None
    ws.source.wait()


def test_a_study_folder_renamed_while_shown_greys_its_figures(app, project):
    ws = _workspace(app, project)
    (project / 'StudyA').rename(project / 'StudyA_old')
    for section, _group, page in ws._entries:
        if section == 'meta':
            continue
        assert not ws._available(page), page.TITLE
        assert page.unavailable_reason() == (
            f"The study folder {project / 'StudyA'} cannot be found."), page.TITLE


def test_a_copied_project_without_raw_data_exports_every_figure(app, project, tmp_path):
    """Lab workstation to laptop: the project folder is copied without
    raw_data/, and every figure is available from the copy and exports into
    the copy's folders."""
    from unittest import mock
    from PyQt6.QtWidgets import QFileDialog
    copy = tmp_path / 'Laptop' / 'Project'
    shutil.copytree(project, copy, ignore=shutil.ignore_patterns('raw_data'))
    shutil.rmtree(project)                       # nothing may come from the original
    assert not any(p.name == 'raw_data' for p in copy.rglob('*'))
    ws = _workspace(app, copy)
    for section in ('scrna', 'de', 'meta'):
        ws.show_section(section)
        for sec, _group, page in ws._entries:
            if sec != section:
                continue
            assert ws._available(page), (page.TITLE, page.unavailable_reason())
            assert _draw(app, page), (page.TITLE, page._status_lbl.text())
            target = page._default_export_dir() / f'{page._default_export_basename()}.png'
            with mock.patch.object(QFileDialog, 'getSaveFileName', return_value=(str(target), '')):
                page._export_current()
            assert page._status_lbl.text().startswith('Exported'), (page.TITLE,
                                                                   page._status_lbl.text())
            assert target.is_file() and copy in target.parents
