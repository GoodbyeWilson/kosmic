"""Stage 3 of ADR-016: Enrichment and Pathway DE work on the gene DE run
chosen in the run bar. They save in that run's folder, clear and reload
when another run is chosen, score a cell type's cells only, and record
the run's cell type and counts in the provenance."""
from __future__ import annotations

import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from kosmic import provenance  # noqa: E402
from kosmic.de import saved  # noqa: E402
from kosmic.paths import de_downstream_dir, de_pathway_scoring_dir, de_run_dir  # noqa: E402

GENES = [f'G{i}' for i in range(6)]
ORA = pd.DataFrame({'Term': ['T1', 'T2'], 'P_value': [0.001, 0.2], 'FDR': [0.01, 0.4],
                    'Fold_Enrichment': [3.0, 1.1], 'Gene_Count': [2, 1], 'Term_Size': [20, 30],
                    'Genes': ['G0, G1', 'G1']})


@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _table():
    return pd.DataFrame({
        'names': GENES, 'logfoldchanges': [1.0, -1.0, 0.3, -0.3, 0.1, 0.0],
        'pvals': [0.001] * 6, 'pvals_adj': [0.01, 0.01, 0.5, 0.5, 0.5, 0.5],
        'disease_samples': [3] * 6, 'control_samples': [3] * 6})


def _run(study, cell_type=None, **record):
    folder = de_run_dir(study, 'raw', cell_type)
    folder.mkdir(parents=True, exist_ok=True)
    name = 'S1' + (f"_{cell_type.replace(' ', '_')}" if cell_type else '')
    path = folder / f'{name}_DE_deseq2.csv'
    _table().to_csv(path, index=False)
    saved.save_run_record(path, {'method': 'deseq2', 'cell_type': cell_type,
                                 'counts_source': 'raw', 'fdr': 0.05, 'lfc_threshold': 0.0,
                                 **record})


def _adata():
    """Six donors, two cell types; the cell types are in 'annotation'."""
    rng = np.random.default_rng(0)
    obs = pd.DataFrame({
        'sample': np.repeat([f'd{i}' for i in range(6)], 20),
        'condition': np.repeat(['Control'] * 3 + ['Disease'] * 3, 20),
        'annotation': np.tile(['Endothelial Cell'] * 12 + ['Fibroblast'] * 8, 6),
    })
    obs.index = [f'c{i}' for i in range(len(obs))]
    return ad.AnnData(X=rng.poisson(5, size=(len(obs), len(GENES))).astype(np.float32),
                      obs=obs, var=pd.DataFrame(index=GENES))


@pytest.fixture
def ws(app, tmp_path):
    from kosmic.gui.de_analysis.workspace import DEWorkspace
    study = tmp_path / 'S1'
    _run(study)
    for ct, floor in (('Endothelial Cell', 3), ('Fibroblast', 7)):
        _run(study, ct, cell_type_column='annotation', min_cells=floor, min_counts=0)
    w = DEWorkspace()
    w.project_dir = study
    w.h5ad_path = study / 'processed_data' / 'S1.h5ad'
    w.analysis_mode = 'discovery'
    w.current_adata = _adata()
    w.sample_col, w.condition_col = 'sample', 'condition'
    w.control_label, w.disease_label = 'Control', 'Disease'
    w.gene_de_page._try_load_previous_results()
    return w


def _choose(w, name):
    combo = w.run_bar.combo
    combo.setCurrentIndex(next(i for i in range(combo.count())
                               if combo.itemText(i).startswith(name)))


# -- Enrichment ----------------------------------------------------------------------

def test_ora_runs_of_two_cell_types_are_saved_separately(ws):
    page = ws.enrichment_page
    for ct, term in (('Endothelial Cell', 'endo'), ('Fibroblast', 'fibro')):
        _choose(ws, ct)
        page._annotate_direction_and_display(ORA.assign(Term=[term, 'T2']))
    study = ws.project_dir
    for ct, term in (('Endothelial Cell', 'endo'), ('Fibroblast', 'fibro')):
        back = saved.read_enrichment(de_downstream_dir(study, 'original_counts', ct))
        assert list(back['Term']) == [term, 'T2']
    # Nothing is written to the flat folder of earlier versions.
    assert not (de_pathway_scoring_dir(study) / saved.ENRICHMENT_RESULTS_FILE).exists()


def test_choosing_another_run_clears_and_shows_its_saved_results(ws):
    page = ws.enrichment_page
    _choose(ws, 'Endothelial Cell')
    page._annotate_direction_and_display(ORA)
    assert page._table.rowCount() == 2

    _choose(ws, 'Fibroblast')                       # nothing saved for it
    assert page._table.rowCount() == 0
    assert ws.enrichment_results is None
    assert 'No saved run' in page._status_label.text()

    _choose(ws, 'Endothelial Cell')                 # its saved run comes back
    assert page._table.rowCount() == 2
    assert list(ws.enrichment_results['Term']) == ['T1', 'T2']
    assert 'Endothelial Cell' in page._status_label.text()


def test_enrichment_provenance_keeps_one_entry_per_cell_type(ws):
    page = ws.enrichment_page
    for ct in ('Endothelial Cell', 'Fibroblast', 'Endothelial Cell'):
        _choose(ws, ct)
        page._annotate_direction_and_display(ORA)
        page._record_study_enrichment('ora')
    stages = [st for st in provenance.load(ws.h5ad_path.parent)['stages']
              if st['stage'] == 'enrichment']
    # The rerun of Endothelial Cell replaced its entry; Fibroblast's stays.
    assert sorted(st['params']['cell_type'] for st in stages) == [
        'Endothelial Cell', 'Fibroblast']
    assert {st['params']['counts'] for st in stages} == {'original_counts'}


# -- Pathway DE ------------------------------------------------------------------------

def test_a_cell_types_cells_are_selected_by_its_records_column(ws):
    _choose(ws, 'Fibroblast')
    cells = ws.run_adata()
    assert cells.n_obs == 48
    assert set(cells.obs['annotation']) == {'Fibroblast'}
    assert ws.run_donor_floors() == (7, 0)
    _choose(ws, 'Whole study')
    assert ws.run_adata() is ws.current_adata


def test_pathway_scoring_of_a_cell_type_uses_its_cells_and_donors(ws, monkeypatch):
    from kosmic.gui.de_analysis.pages import pathway_de_page
    started = []
    monkeypatch.setattr(pathway_de_page, 'run_worker', lambda w, **k: started.append(w))
    ws.pathway_gene_sets = {'P': GENES[:4]}
    _choose(ws, 'Fibroblast')
    ws.pathway_de_page._start_scoring()
    worker = started[0]
    assert set(worker.adata.obs['annotation']) == {'Fibroblast'}
    assert (worker.min_cells, worker.min_counts) == (7, 0)
    folder = de_downstream_dir(ws.project_dir, 'original_counts', 'Fibroblast')
    assert worker.scoring_dir == folder

    # The worker saves there.
    worker._save_and_emit(pd.DataFrame({'names': ['P'], 'logfoldchanges': [0.4],
                                        'pvals': [0.01], 'pvals_adj': [0.01]}), ['P'])
    assert (folder / 'pathway_de_results.csv').is_file()


def test_pathway_results_follow_the_chosen_run(ws):
    folder = de_downstream_dir(ws.project_dir, 'original_counts', 'Fibroblast')
    folder.mkdir(parents=True)
    pd.DataFrame({'names': ['P'], 'logfoldchanges': [0.4], 'pvals': [0.01],
                  'pvals_adj': [0.01]}).to_csv(folder / 'pathway_de_results.csv', index=False)
    _choose(ws, 'Fibroblast')
    assert list(ws.pathway_de_results['names']) == ['P']
    _choose(ws, 'Endothelial Cell')
    assert ws.pathway_de_results is None
    assert 'No saved pathway results for Endothelial Cell' in ws.pathway_de_page.status_label.text()
