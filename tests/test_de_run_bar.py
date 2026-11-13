"""Stage 2 of ADR-016: the DE workspace's run bar lists the study's gene DE
runs, and the run chosen in it is the one every page after gene DE uses."""
from __future__ import annotations

import os

import pandas as pd
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from kosmic.de import saved  # noqa: E402
from kosmic.de.batch import CellTypeRun, summary_path  # noqa: E402
from kosmic.paths import de_run_dir  # noqa: E402


@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _table(n_sig, lfc=1.0):
    genes = [f'G{i}' for i in range(6)]
    return pd.DataFrame({
        'names': genes, 'logfoldchanges': [lfc, -lfc, 0.3, -0.3, 0.1, 0.0],
        'pvals': [0.001] * 6, 'pvals_adj': [0.01 if i < n_sig else 0.5 for i in range(6)],
        'disease_samples': [4] * 6, 'control_samples': [5] * 6})


def _save(folder, name, table, rule, cell_type=None):
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f'{name}_DE_deseq2.csv'
    table.to_csv(path, index=False)
    fdr, lfc = rule
    saved.save_run_record(path, {'method': 'deseq2', 'cell_type': cell_type,
                                 'counts_source': 'raw', 'fdr': fdr, 'lfc_threshold': lfc})
    return path


@pytest.fixture
def study(tmp_path):
    """S1: a whole-study run, two cell types with their own rules, and a
    cell type the batch run skipped."""
    s = tmp_path / 'S1'
    _save(de_run_dir(s, 'raw'), 'S1', _table(4), (0.05, 0.0))
    _save(de_run_dir(s, 'raw', 'Endothelial Cell'), 'S1_Endothelial_Cell', _table(3),
          (0.05, 0.5), 'Endothelial Cell')
    _save(de_run_dir(s, 'raw', 'Fibroblast'), 'S1_Fibroblast', _table(2), (0.1, 0.0),
          'Fibroblast')
    pd.DataFrame({'cell_type': ['Endothelial Cell', 'Fibroblast', 'Myeloid'],
                  'status': ['ok', 'ok', 'skipped'],
                  'message': ['', '', 'too few donors: 1 control'],
                  'method': ['deseq2'] * 3}).to_csv(
        summary_path(s, 'S1', 'raw'), index=False)
    return s


@pytest.fixture
def ws(app, study):
    from kosmic.gui.de_analysis.workspace import DEWorkspace
    w = DEWorkspace()
    w.project_dir = study
    w.analysis_mode = 'discovery'
    w.gene_de_page._try_load_previous_results()
    return w


def _entries(w):
    combo = w.run_bar.combo
    return [(combo.itemText(i).split('  ')[0], combo.model().item(i).isEnabled())
            for i in range(combo.count())]


def test_the_runs_are_listed_whole_study_first(ws):
    assert _entries(ws) == [('Whole study', True), ('Endothelial Cell', True),
                            ('Fibroblast', True), ('Myeloid', False)]
    assert 'too few donors' in ws.run_bar.combo.itemText(3)
    # The whole-study run is shown first, with its own rule.
    assert ws.run_label() == 'Whole study'
    assert len(ws.de_results) == 6 and len(ws.significant_genes) == 4
    assert ws.significance == (0.05, 0.0, False)


def test_choosing_a_cell_type_changes_the_results_and_rule(ws, study):
    ws.run_bar.combo.setCurrentIndex(2)                 # Fibroblast
    assert ws.run_label() == 'Fibroblast'
    expected = pd.read_csv(de_run_dir(study, 'raw', 'Fibroblast') / 'S1_Fibroblast_DE_deseq2.csv')
    pd.testing.assert_frame_equal(ws.de_results, expected)
    assert ws.significance == (0.1, 0.0, False)
    assert len(ws.significant_genes) == 2
    detail = ws.run_bar.detail.text()
    assert 'DESeq2' in detail and 'Original counts' in detail
    assert '5 control, 4 disease donors' in detail
    assert '2 significant (FDR < 0.1)' in detail

    ws.run_bar.combo.setCurrentIndex(1)                 # Endothelial Cell
    assert ws.significance == (0.05, 0.5, False)
    assert '|log2FC| > 0.5' in ws.run_bar.detail.text()


def test_a_cell_type_that_did_not_run_cannot_be_chosen(ws):
    ws.select_run(3)
    assert ws.run_label() == 'Whole study'
    assert ws.run_bar.combo.currentIndex() == 0


def test_the_bar_is_shown_on_the_pages_after_gene_de(ws):
    ws.switch_tab(0)                                    # Load Data
    assert ws.run_bar.isHidden()
    for step in (2, 3, 4):                              # Gene DE, Enrichment, Methods
        ws.switch_tab(step)
        assert not ws.run_bar.isHidden()
    ws.clear_runs()
    assert ws.run_bar.isHidden()


def test_the_status_line_names_the_run(ws):
    ws.switch_tab(2)
    ws.run_bar.combo.setCurrentIndex(2)
    text, _state = ws.compute_sidebar_status()
    assert text.startswith('Fibroblast: DE complete')


def test_a_batch_run_lists_its_cell_types_and_shows_the_first(ws, monkeypatch):
    from kosmic.gui.shared import dialogs
    monkeypatch.setattr(dialogs, 'info', lambda *a, **k: None)
    runs = [CellTypeRun(cell_type='T cell', slug='T_cell', accession='S1_T_cell', status='ok',
                        n_samples=9, n_significant=1, de_results=_table(1)),
            CellTypeRun(cell_type='Atrial', slug='Atrial', accession='S1_Atrial',
                        status='failed', message='no donors'),
            CellTypeRun(cell_type='B cell', slug='B_cell', accession='S1_B_cell', status='ok',
                        n_samples=9, n_significant=2, de_results=_table(2))]
    ws.gene_de_page._on_batch_finished((runs, None, 'Batch DE complete'))
    assert _entries(ws) == [('Whole study', True), ('Atrial', False),
                            ('B cell', True), ('T cell', True)]
    assert ws.run_label() == 'B cell'
    assert ws.de_results is runs[2].de_results
