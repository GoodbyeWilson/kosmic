"""Stage 2 of ADR-014: reading what the DE steps saved to a study folder,
and saving what they did not (run records, over-representation tables,
GSEA running-enrichment curves)."""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
import pytest

from kosmic.de import saved
from tests.study_folder import make_study_folder, write_de_outputs

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


@pytest.fixture
def study(tmp_path):
    return write_de_outputs(make_study_folder(tmp_path, form='counts'))


def test_runs_are_listed_with_cell_types_and_methods(study):
    runs = saved.list_gene_de_runs(study)
    assert [(r.cell_type, r.method) for r in runs] == [
        (None, 'welch_cpm_eb'),                   # the whole study first
        ('Endothelial', 'deseq2'),                # from the batch summary
        ('Fibroblast', 'deseq2'),
        ('Myeloid NN', 'deseq2'),                 # lower-case file name, no summary row
    ]
    assert runs[0].label == 'StudyA (welch_cpm_eb)'
    # Tables only: the significant subsets and pseudobulk files are not runs.
    assert all('_DE_' in r.path.name for r in runs)
    # The summary's de_path points to another computer; paths are rebuilt
    # from the study folder.
    assert all(r.path.parent == study / 'results' / 'de_analysis' / 'statistics'
               for r in runs)
    assert all(r.path.is_file() for r in runs)


def test_a_run_table_reads_back(study):
    run = saved.list_gene_de_runs(study)[1]
    df = saved.read_gene_de(run)
    assert 'names' in df.columns and len(df) == 60
    run.path.write_text('not,a\ntable')
    assert saved.read_gene_de(run) is None


def test_a_run_record_round_trips_and_names_the_count_source(study):
    from kosmic.scrna.load.h5ad_meta import read_counts_adata, registered_h5ad

    h5ad = registered_h5ad(make_study_folder(study.parent, 'Other', form='log'))
    adata = read_counts_adata(h5ad, counts_layer='decontX_counts')
    record = saved.build_run_record(
        adata, method='deseq2', sample_col='sample', condition_col='condition',
        disease_label='Disease', control_label='Control', min_cells=10, min_counts=0,
        pathway_gene_sets={'Glycolysis': ['PKM', 'HK2']}, geneset_label='1 pathway',
        h5ad_path=h5ad)
    assert record['counts_source'] == 'decontX_counts'
    assert record['dataset_file'] == h5ad.name          # the name only
    run = saved.list_gene_de_runs(study)[1]
    assert saved.read_run_record(run.path) is None      # saved before records
    saved.save_run_record(run.path, {**record, 'cell_type': 'Endothelial'})
    back = saved.read_run_record(run.path)
    assert back['pathway_gene_sets'] == {'Glycolysis': ['PKM', 'HK2']}
    assert back['counts_source'] == 'decontX_counts'
    json.loads(run.record_path.read_text())             # plain JSON
    # The record's cell type is used for the run's label.
    saved.save_run_record(run.path, {**record, 'cell_type': 'Endothelial cell'})
    assert saved.list_gene_de_runs(study)[1].cell_type == 'Endothelial cell'
    # A record is not taken for a run of its own.
    assert len(saved.list_gene_de_runs(study)) == 4


def test_pathway_outputs_read_back_at_both_levels(study):
    assert saved.pathway_levels(study) == ['', 'parent_']
    assert list(saved.read_pathway_de(study)['names']) == ['Glycolysis', 'OXPHOS']
    assert list(saved.read_pathway_de(study, 'parent_')['names']) == ['Glycolysis']
    scores = saved.read_per_donor_scores(study)
    assert scores['names'] == ['Glycolysis', 'OXPHOS']
    assert saved.read_per_donor_scores(study, 'parent_') is None
    assert list(saved.read_fgsea(study)['names']) == ['HALLMARK_A', 'HALLMARK_B']
    assert list(saved.read_gsea_library(study)['nes']) == [1.8, -1.2]


def test_gsea_details_round_trip(study):
    path = saved.gsea_library_path(study)
    assert saved.read_gsea_details(path) == (None, '')  # a run saved before
    details = {'HALLMARK_A': {'RES': np.linspace(0, 0.5, 60), 'hits': [3, 7, 40],
                              'es': 0.5, 'nes': 1.8},
               'odd name / with, punctuation': {'RES': np.zeros(60), 'hits': [],
                                                'es': -0.1, 'nes': -0.4}}
    saved.save_gsea_details(path, details, 'MSigDB Hallmark')
    back, library = saved.read_gsea_details(path)
    assert library == 'MSigDB Hallmark'
    assert list(back) == list(details)
    for term, d in details.items():
        np.testing.assert_array_equal(back[term]['RES'], d['RES'])
        assert back[term]['hits'] == d['hits']
        assert back[term]['es'] == d['es'] and back[term]['nes'] == d['nes']


def test_enrichment_table_round_trips(tmp_path):
    df = pd.DataFrame({'Direction': ['2/3 up'], 'Term': ['GO:1 glycolysis'],
                       'P_value': [0.001], 'FDR': [0.01], 'Fold_Enrichment': [3.2],
                       'Gene_Count': [3], 'Term_Size': [40], 'Genes': ['PKM, HK2, LDHA']})
    assert saved.read_enrichment(tmp_path) is None
    saved.save_enrichment(tmp_path / 'x', df)
    pd.testing.assert_frame_equal(saved.read_enrichment(tmp_path / 'x'), df)


def test_batch_de_saves_a_record_beside_each_table(tmp_path):
    from kosmic.de import batch as B
    from tests.test_de_batch import _adata

    record = {'counts_source': 'raw', 'geneset_label': 'x', 'pathway_gene_sets': {}}
    runs = B.run_de_by_cell_type(
        _adata(), 'cell_type', 'donor', 'condition', tmp_path / 'GSE1', 'GSE1',
        de_method='ttest', min_cells=10, detection_min_pct=0.0, run_record=record)
    for r in runs:
        rec = saved.read_run_record(r.de_path)
        assert rec['cell_type'] == r.cell_type
        assert rec['cell_type_column'] == 'cell_type'
        assert rec['method'] == 'welch_cpm' and rec['counts_source'] == 'raw'
    listed = saved.list_gene_de_runs(tmp_path / 'GSE1')
    assert {r.cell_type for r in listed} == {r.cell_type for r in runs}


def test_the_whole_study_export_saves_a_record(app, tmp_path):
    from kosmic.gui.de_analysis.workers import MetaExportWorker
    from tests.study_folder import de_table

    w = MetaExportWorker(None, 'sample', 'condition', 10, 'ttest', True,
                         de_table(['A', 'B']), tmp_path / 'S', 'S',
                         run_record={'counts_source': 'raw', 'cell_type': None})
    w._run()                         # the pseudobulk step fails without data; reported
    (run,) = saved.list_gene_de_runs(tmp_path / 'S')
    assert run.method == 'welch_cpm_eb' and run.cell_type is None
    assert saved.read_run_record(run.path)['method'] == 'welch_cpm_eb'


@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _enrichment_page(study):
    from kosmic.gui.de_analysis.pages.enrichment_page import EnrichmentPage
    from kosmic.gui.meta_analysis.pages.enrichment_step import _MetaEnrichmentContext

    ws = _MetaEnrichmentContext()
    ws.set_meta_results(pd.DataFrame({'names': ['PKM', 'HK2'], 'logfoldchanges': [1.0, -1.0],
                                      'fdr': [0.01, 0.2], 'pvals_pooled': [0.001, 0.1]}))
    ws.project_dir = study
    return EnrichmentPage(ws), ws


def test_the_enrichment_page_saves_an_ora_run(app, study):
    page, _ws = _enrichment_page(study)
    page._annotate_direction_and_display(pd.DataFrame({
        'Term': ['T1', 'T2'], 'P_value': [0.001, 0.2], 'FDR': [0.01, 0.4],
        'Fold_Enrichment': [3.0, 1.1], 'Gene_Count': [2, 1], 'Term_Size': [20, 30],
        'Genes': ['PKM, HK2', 'HK2']}))
    from kosmic.paths import de_pathway_scoring_dir
    back = saved.read_enrichment(de_pathway_scoring_dir(study))
    assert list(back['Term']) == ['T1', 'T2'] and 'Direction' in back.columns


def test_the_meta_enrichment_step_saves_an_ora_run_in_its_selection(app, tmp_path):
    from kosmic.gui.meta_analysis.pages.enrichment_step import MetaEnrichmentPage
    from kosmic.paths import meta_output_dir

    p = MetaEnrichmentPage()
    p.set_project_folder(tmp_path, selection='Cardiomyocyte')
    p.set_meta_results(pd.DataFrame({'names': ['G1', 'G2'], 'logfoldchanges': [1.0, -1.0],
                                     'fdr': [0.01, 0.2], 'pvals_pooled': [0.001, 0.1]}))
    p._on_enrichr_finished((pd.DataFrame({
        'Term': ['a'], 'P_value': [1e-5], 'FDR': [3e-5], 'Fold_Enrichment': [3.0],
        'Gene_Count': [2], 'Genes': [['G1', 'G2']]}), 'done'))
    back = saved.read_enrichment(meta_output_dir(tmp_path, 'Cardiomyocyte'))
    assert list(back['Term']) == ['a']


def test_the_enrichment_page_saves_gsea_curves(app, study):
    page, ws = _enrichment_page(study)
    ws.gsea_library_name = 'MSigDB Hallmark'
    ws.gsea_library_details = {'HALLMARK_A': {'RES': np.linspace(0, 1, 5), 'hits': [1],
                                              'es': 0.9, 'nes': 2.0}}
    page._show_gsea_results(pd.DataFrame({
        'names': ['HALLMARK_A'], 'nes': [2.0], 'es': [0.9], 'pvals': [0.001],
        'pvals_adj': [0.01], 'leading_edge': ['PKM']}))
    details, library = saved.read_gsea_details(saved.gsea_library_path(study))
    assert library == 'MSigDB Hallmark' and list(details) == ['HALLMARK_A']
