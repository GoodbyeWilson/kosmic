"""Stage 3 of ADR-014: where Figure Export's study figures get their data.

The source serves a workspace's data when that workspace holds the study
and the study's folder otherwise; it caches what it read against the file,
reads the file's metadata off the GUI thread, and answers 'is this
available?' without reading a table.
"""
from __future__ import annotations

import os
import threading
import time
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from kosmic.de import saved  # noqa: E402
from kosmic.scrna.load import h5ad_meta  # noqa: E402
from tests.study_folder import make_study_folder, write_de_outputs  # noqa: E402


@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.fixture
def study(tmp_path):
    return write_de_outputs(make_study_folder(tmp_path, form='log'))


def _source(scrna=None, de=None):
    from kosmic.gui.figure_export.shared.study_source import StudySource
    return StudySource(lambda: scrna, lambda: de)


def _settle(app, src):
    src.wait()
    app.processEvents()


# -- scRNA ----------------------------------------------------------------------

def test_the_scrna_workspaces_dataset_is_used_when_it_holds_the_study(app, study):
    held = ad.AnnData(np.zeros((2, 2)))
    src = _source(scrna=SimpleNamespace(current_adata=held, current_project_dir=str(study)))
    src.set_study(study)
    assert src.scrna_adata() is held
    assert src.holds_scrna() and src.scrna_status() is None
    assert src.scrna_expression(['VWF']) is held


def test_the_file_is_read_in_the_background_when_no_workspace_holds_it(app, study, tmp_path):
    other = SimpleNamespace(current_adata=ad.AnnData(np.zeros((2, 2))),
                            current_project_dir=str(tmp_path / 'Other'))
    src = _source(scrna=other)                      # scRNA holds a different study
    changes = []
    src.changed.connect(lambda: changes.append(1))
    src.set_study(study)
    assert src.scrna_adata() is None                # not read yet...
    assert src.scrna_status().startswith('Reading StudyA.h5ad')
    _settle(app, src)
    meta = src.scrna_adata()                        # ...then the file's
    assert meta is not None and meta.X is None and 'X_umap' in meta.obsm
    assert src.scrna_status() is None
    assert len(changes) == 2                        # the study set, the read landed


def test_reads_are_cached_against_the_file(app, study, monkeypatch):
    calls = []
    real = h5ad_meta.read_figure_meta
    monkeypatch.setattr(h5ad_meta, 'read_figure_meta',
                        lambda p: calls.append(p) or real(p))
    src = _source()
    src.set_study(study)
    src.scrna_adata()
    _settle(app, src)
    first = src.scrna_adata()
    assert src.scrna_adata() is first and len(calls) == 1
    # Saved again: read again.
    path = src.h5ad_path()
    time.sleep(0.01)
    h5ad_meta.save_h5ad(ad.read_h5ad(path), path)
    assert src.scrna_adata() is None
    _settle(app, src)
    assert src.scrna_adata() is not first and len(calls) == 2
    # Another study: nothing kept.
    src.set_study(None)
    assert src.scrna_adata() is None and src.scrna_status() is None


def test_expression_from_the_file_is_read_once_and_retried_after_a_save(app, study,
                                                                         monkeypatch):
    src = _source()
    src.set_study(study)
    calls = []
    real = h5ad_meta.read_genes

    def _flaky(*a, **kw):
        calls.append(1)
        if len(calls) == 1:
            raise h5ad_meta.StudyFileChanged('saved meanwhile')
        return real(*a, **kw)

    monkeypatch.setattr(h5ad_meta, 'read_genes', _flaky)
    got = src.scrna_expression(['VWF', 'PTPRC'])
    assert list(got.var_names) == ['VWF', 'PTPRC'] and len(calls) == 2
    assert src.scrna_expression(['VWF', 'PTPRC']) is got and len(calls) == 2


def test_an_unreadable_file_is_reported_not_raised(app, study):
    src = _source()
    src.set_study(study)
    path = src.h5ad_path()
    path.write_bytes(path.read_bytes()[:500])
    src.scrna_adata()
    _settle(app, src)
    assert src.scrna_adata() is None
    assert src.scrna_status().startswith('Could not read StudyA.h5ad')


def test_a_missing_study_folder_is_reported(app, tmp_path):
    src = _source()
    src.set_study(tmp_path / 'Renamed')
    assert src.scrna_adata() is None
    assert 'cannot be found' in src.scrna_status()


def test_a_source_discarded_while_reading_does_not_abort(app, study, monkeypatch):
    from PyQt6.QtCore import QEvent
    release = threading.Event()
    real = h5ad_meta.read_figure_meta
    monkeypatch.setattr(h5ad_meta, 'read_figure_meta',
                        lambda p: (release.wait(10), real(p))[1])
    src = _source()
    src.set_study(study)
    src.scrna_adata()
    thread = src._meta_thread
    src.deleteLater()
    app.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    del src
    release.set()
    thread.join(10)
    app.processEvents()
    assert not thread.is_alive()


# -- DE -------------------------------------------------------------------------

def test_de_attributes_come_from_the_folder_when_de_is_elsewhere(app, study, tmp_path):
    de = SimpleNamespace(project_dir=tmp_path / 'Other', de_results=pd.DataFrame({'x': [1]}),
                         current_adata=None)
    src = _source(de=de)
    src.set_study(study)
    st = src.de_state()
    first = saved.list_gene_de_runs(study)[0]
    pd.testing.assert_frame_equal(st.de_results, saved.read_gene_de(first))
    assert st.run_label == first.label
    # Labels and columns from the manifest when the run has no record.
    assert (st.disease_label, st.control_label) == ('Disease', 'Control')
    assert (st.sample_col, st.condition_col) == ('sample', 'condition')
    assert st.counts_source == 'raw'
    assert list(st.pathway_de_results['names']) == ['Glycolysis', 'OXPHOS']
    assert st.pathway_score_data['names'] == ['Glycolysis', 'OXPHOS']
    assert list(st.gsea_library_results['names']) == ['HALLMARK_A', 'HALLMARK_B']
    assert st.enrichment_results is None and not st.available('enrichment_results')
    assert st.project_dir == study and st.h5ad_path.name == 'StudyA.h5ad'
    assert st.current_adata is None


def test_a_run_record_supplies_the_runs_settings(app, study):
    src = _source()
    src.set_study(study)
    run = saved.list_gene_de_runs(study)[1]
    saved.save_run_record(run.path, {
        'cell_type': 'Endothelial', 'disease_label': 'DCM', 'control_label': 'NF',
        'counts_source': 'decontX_counts', 'geneset_label': '1 pathway',
        'pathway_gene_sets': {'Glycolysis': ['VWF', 'PTPRC', 'NOPE']}})
    src.set_de_run(run)
    st = src.de_state()
    assert (st.disease_label, st.control_label) == ('DCM', 'NF')
    assert st.counts_source == 'decontX_counts' and st.geneset_label == '1 pathway'
    assert st.pathway_gene_sets == {'Glycolysis': ['VWF', 'PTPRC', 'NOPE']}
    assert st.available('pathway_coverage')
    assert st.pathway_coverage['Glycolysis']['available_genes'] == 2


def test_the_de_workspaces_run_is_used_when_it_is_on_the_study(app, study):
    live = pd.DataFrame({'names': ['A'], 'logfoldchanges': [1.0]})
    de = SimpleNamespace(project_dir=study, de_results=live, disease_label='DCM',
                         control_label='NF', pathway_de_results=None, current_adata=None)
    src = _source(de=de)
    src.set_study(study)
    st = src.de_state()
    assert st.de_results is live and st.disease_label == 'DCM'
    assert st.run_label == 'Run in the DE workspace'
    # What the workspace has not got this session comes from the folder.
    assert list(st.pathway_de_results['names']) == ['Glycolysis', 'OXPHOS']
    # A saved run chosen instead serves its own table.
    run = saved.list_gene_de_runs(study)[2]
    src.set_de_run(run)
    assert src.de_state().de_results is not live


def test_availability_reads_no_table(app, study, monkeypatch):
    reads = []
    for name in ('read_gene_de', 'read_pathway_de', 'read_per_donor_scores', 'read_fgsea',
                 'read_gsea_library', 'read_gsea_details', 'read_enrichment'):
        monkeypatch.setattr(saved, name, lambda *a, n=name, **k: reads.append(n))
    src = _source()
    src.set_study(study)
    st = src.de_state()
    for name in ('de_results', 'pathway_de_results', 'pathway_score_data', 'fgsea_results',
                 'gsea_library_results', 'gsea_library_details', 'enrichment_results',
                 'disease_label', 'sample_col', 'pathway_gene_sets'):
        st.available(name)
    assert reads == []
    assert st.available('de_results') and st.available('gsea_library_results')
    assert not st.available('gsea_library_details')       # saved before curves were kept


def test_de_counts_come_from_the_runs_count_source(app, study):
    src = _source()
    src.set_study(study)
    got = src.de_counts(['CD14', 'G1'], 'decontX_counts')
    want = h5ad_meta.read_genes(src.h5ad_path(), ['CD14', 'G1'], source='counts',
                                counts_layer='decontX_counts')
    np.testing.assert_array_equal(got.X.toarray(), want.X.toarray())
    assert got.uns['counts_loaded_from'] == 'layers:decontX_counts'


# -- The workspace ----------------------------------------------------------------

def test_changing_study_invalidates_previews_and_refills_the_list(app, study, tmp_path):
    from kosmic.gui.figure_export.workspace import FigureExportWorkspace
    ws = FigureExportWorkspace()
    before = ws.data_version()
    rebuilt = []
    ws._rebuild_list = lambda: rebuilt.append(1)
    ws.set_study(study)
    assert ws.data_version() > before and rebuilt
    v = ws.data_version()
    ws.set_study(study)                                   # the same study: nothing
    assert ws.data_version() == v
    ws.source.wait()


def test_labels_come_from_the_dataset_when_no_record_or_manifest_has_them(app, study):
    """An atlas has no manifest and its runs may predate records: the
    labels come from the dataset's role map, else its '_role' column."""
    import anndata as ad
    from kosmic.manifest import manifest_path
    manifest_path(study).unlink()
    path = h5ad_meta.registered_h5ad(study)
    a = ad.read_h5ad(path)
    a.uns['role_map'] = {'NF': 'control', 'DCM': 'disease'}
    h5ad_meta.save_h5ad(a, path)
    src = _source()
    src.set_study(study)
    st = src.de_state()
    assert st.disease_label == ''                    # metadata not read yet: nothing waits
    _settle(app, src)
    assert (st.disease_label, st.control_label) == ('DCM', 'NF')
    del a.uns['role_map']
    h5ad_meta.save_h5ad(a, path)
    _settle(app, src)
    src.scrna_adata()
    _settle(app, src)
    st = src.de_state()
    assert (st.disease_label, st.control_label) == ('Disease', 'Control')
