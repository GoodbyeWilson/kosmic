"""Stage 1 of ADR-016: the outputs of the analyses made on one gene DE run
(enrichment, GSEA, pathway DE) are read from a folder for that run's
counts and cell type, and from the flat 'pathway_scoring/' folder of
earlier versions for a whole-study run."""
from __future__ import annotations

import os

import pandas as pd
import pytest

from kosmic.de import saved
from kosmic.paths import de_downstream_dir, de_pathway_scoring_dir, de_run_dir
from tests.study_folder import MARKER_GENES, de_table, make_study_folder, write_de_outputs

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

ORA = pd.DataFrame({'Term': ['immune response'], 'FDR': [0.001], 'Fold_Enrichment': [3.0],
                    'Gene_Count': [3], 'Genes': ['PTPRC, CD3E, CD14']})
PATHWAYS = pd.DataFrame({'names': ['Glycolysis'], 'logfoldchanges': [0.5], 'pvals': [0.01],
                         'pvals_adj': [0.02], 'se': [0.2]})


def _ora(term):
    return ORA.assign(Term=[term])


def _new_layout(study, counts='decontX_counts', cell_types=('Endothelial', 'CD8 T Cell')):
    """Gene DE runs and their downstream outputs in the ADR-015/016 folders:
    a whole-study run and one run per cell type, each with its own ORA
    table and pathway DE table."""
    from kosmic.de.batch import slugify_cell_type
    name = study.name
    whole = de_run_dir(study, counts)
    whole.mkdir(parents=True, exist_ok=True)
    de_table(MARKER_GENES).to_csv(whole / f'{name}_DE_deseq2.csv', index=False)
    saved.save_enrichment(de_downstream_dir(study, counts), _ora('whole study'))
    PATHWAYS.to_csv(de_downstream_dir(study, counts) / 'pathway_de_results.csv', index=False)
    for ct in cell_types:
        slug = slugify_cell_type(ct)
        run = de_run_dir(study, counts, ct)
        run.mkdir(parents=True, exist_ok=True)
        de_table(MARKER_GENES).to_csv(run / f'{name}_{slug}_DE_deseq2.csv', index=False)
        out = de_downstream_dir(study, counts, ct)
        saved.save_enrichment(out, _ora(ct))
        PATHWAYS.assign(names=[f'{ct} pathway']).to_csv(
            out / 'pathway_de_results.csv', index=False)
    return study


# -- Paths -------------------------------------------------------------------------

def test_the_folder_of_a_run_is_named_after_its_counts_and_cell_type(tmp_path):
    scoring = de_pathway_scoring_dir(tmp_path)
    assert de_downstream_dir(tmp_path) == scoring / 'original_counts'
    assert de_downstream_dir(tmp_path, 'raw') == scoring / 'original_counts'
    assert de_downstream_dir(tmp_path, 'original_counts') == scoring / 'original_counts'
    assert (de_downstream_dir(tmp_path, 'decontX_counts', 'CD8 T Cell')
            == scoring / 'decontX_counts' / 'CD8_T_Cell')
    # The same cell-type folder name as the run's statistics/ folder.
    assert (de_downstream_dir(tmp_path, 'decontX_counts', 'CD8 T Cell').name
            == de_run_dir(tmp_path, 'decontX_counts', 'CD8 T Cell').name)


# -- Readers -----------------------------------------------------------------------

def test_two_cell_types_are_read_separately(tmp_path):
    study = _new_layout(make_study_folder(tmp_path, form='counts'))
    for ct in ('Endothelial', 'CD8 T Cell'):
        folder = saved.downstream_dir(study, 'decontX_counts', ct)
        assert list(saved.read_enrichment(folder)['Term']) == [ct]
        assert list(saved.read_pathway_de(study, counts='decontX_counts',
                                          cell_type=ct)['names']) == [f'{ct} pathway']
    whole = saved.downstream_dir(study, 'decontX_counts')
    assert list(saved.read_enrichment(whole)['Term']) == ['whole study']
    assert saved.pathway_levels(study, 'decontX_counts', 'Endothelial') == ['']


def test_a_run_whose_counts_are_not_known_is_found_in_any_counts_folder(tmp_path):
    study = _new_layout(make_study_folder(tmp_path, form='counts'))
    folder = saved.downstream_dir(study, None, 'Endothelial')
    assert folder == de_downstream_dir(study, 'decontX_counts', 'Endothelial')
    assert saved.downstream_dir(study) == de_downstream_dir(study, 'decontX_counts')


def test_other_counts_are_not_read_for_a_run(tmp_path):
    study = _new_layout(make_study_folder(tmp_path, form='counts'))
    assert saved.read_pathway_de(study, counts='original_counts',
                                 cell_type='Endothelial') is None
    assert saved.read_enrichment(saved.downstream_dir(study, 'raw', 'Endothelial')) is None


def test_the_flat_folder_is_read_for_a_whole_study_run_only(tmp_path):
    study = write_de_outputs(make_study_folder(tmp_path, form='counts'))  # flat layout
    flat = de_pathway_scoring_dir(study)
    assert saved.downstream_dir(study) == flat
    assert saved.downstream_dir(study, 'original_counts') == flat
    assert list(saved.read_pathway_de(study, counts='decontX_counts')['names']) == [
        'Glycolysis', 'OXPHOS']
    # A cell type's outputs are never taken from the flat folder: which
    # cell type they were made on was not recorded.
    assert saved.read_pathway_de(study, cell_type='Endothelial') is None
    assert saved.read_fgsea(study, 'original_counts', 'Endothelial') is None


def test_a_whole_study_runs_own_folder_is_preferred_to_the_flat_one(tmp_path):
    study = write_de_outputs(make_study_folder(tmp_path, form='counts'))
    own = de_downstream_dir(study, 'original_counts')
    saved.save_enrichment(own, _ora('own folder'))
    assert saved.downstream_dir(study, 'original_counts') == own
    assert list(saved.read_enrichment(saved.downstream_dir(study))['Term']) == ['own folder']


# -- Pathway-level discovery -------------------------------------------------------

def test_pathway_discovery_names_each_run_with_its_counts(tmp_path):
    from kosmic.meta_analysis.io import discover_de_results

    study = _new_layout(make_study_folder(tmp_path, 'StudyB', form='counts'))
    entries = [e for e in discover_de_results(tmp_path) if e['file_type'] == 'pathway_de']
    found = sorted((e['accession'], e['counts']) for e in entries)
    assert found == [('StudyB', 'decontX_counts'),
                     ('StudyB_CD8_T_Cell', 'decontX_counts'),
                     ('StudyB_Endothelial', 'decontX_counts')]
    by_acc = {e['accession']: e for e in entries}
    assert by_acc['StudyB_Endothelial']['path'].parent.name == 'Endothelial'
    assert all(e['analysis_folder'] == study.name for e in entries)


def test_pathway_discovery_reads_the_flat_folder_without_counts(tmp_path):
    from kosmic.meta_analysis.io import discover_de_results

    write_de_outputs(make_study_folder(tmp_path, form='counts'))
    entries = [e for e in discover_de_results(tmp_path) if e['file_type'] == 'pathway_de']
    assert sorted(e['path'].name for e in entries) == [
        'parent_pathway_de_results.csv', 'pathway_de_results.csv']
    # The whole-study run's accession, not one of the cell types'.
    assert {(e['accession'], e['counts']) for e in entries} == {('StudyA', None)}


# -- Figure Export -----------------------------------------------------------------

@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_figure_export_reads_the_selected_runs_folder(app, tmp_path):
    from kosmic.gui.figure_export.shared.study_source import StudySource

    study = _new_layout(make_study_folder(tmp_path, form='counts'))
    src = StudySource(lambda: None, lambda: None)
    src.set_study(study)
    runs = {r.cell_type: r for r in saved.list_gene_de_runs(study)}
    for ct in (None, 'Endothelial', 'CD8 T Cell'):
        src.set_de_run(runs[ct])
        st = src.de_state()
        assert st.available('enrichment_results')
        assert list(st.enrichment_results['Term']) == [ct or 'whole study']
        assert st.available('pathway_de_results')
        assert not st.available('fgsea_results')
        assert src.de_state().saved_path('enrichment_results').parent == \
            de_downstream_dir(study, 'decontX_counts', ct)


def test_the_meta_pathway_pool_holds_each_run_on_the_chosen_counts(app, tmp_path):
    from kosmic.gui.meta_analysis.workspace import MetaAnalysisWorkspace

    cell = 'Ventricular_Cardiomyocyte'
    for study in ('Chaffin', 'Guo'):
        for counts, lfc in (('raw', 0.5), ('decontX_counts', -0.5)):
            run = de_run_dir(tmp_path / study, counts, cell)
            run.mkdir(parents=True)
            de_table(MARKER_GENES).to_csv(run / f'{study}_{cell}_DE_deseq2.csv', index=False)
            out = de_downstream_dir(tmp_path / study, counts, cell)
            out.mkdir(parents=True)
            PATHWAYS.assign(logfoldchanges=[lfc]).to_csv(
                out / 'pathway_de_results.csv', index=False)
    ws = MetaAnalysisWorkspace()
    ws._project_folder = str(tmp_path)
    ws._rescan_project()
    pool = {d['name']: d for d in ws._all_pathway_datasets}
    assert sorted(pool) == [f'Chaffin_{cell}', f'Guo_{cell}']
    assert all(d['df']['logfoldchanges'].eq(0.5).all() for d in pool.values())

    ws._counts_combo.setCurrentIndex(ws._counts_combo.findData('decontX_counts'))
    pool = {d['name']: d for d in ws._all_pathway_datasets}
    assert sorted(pool) == [f'Chaffin_{cell}', f'Guo_{cell}']
    assert all(d['df']['logfoldchanges'].eq(-0.5).all() for d in pool.values())
