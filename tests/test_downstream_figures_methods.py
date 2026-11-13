"""Stage 4 of ADR-016: Figure Export, the provenance's methods and the DE
Methods page follow the gene DE run, per cell type."""
from __future__ import annotations

import os

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from kosmic import provenance  # noqa: E402
from kosmic.de import saved  # noqa: E402
from kosmic.paths import de_downstream_dir  # noqa: E402
from tests import test_downstream_per_run as per_run  # noqa: E402

ORA, _choose = per_run.ORA, per_run._choose
app, ws = per_run.app, per_run.ws          # the stage 3 fixtures


def _source(w):
    from kosmic.gui.figure_export.shared.study_source import StudySource
    src = StudySource(lambda: None, lambda: w)
    src.set_study(w.project_dir)
    return src


# -- Figure Export -----------------------------------------------------------------

def test_figure_export_names_and_reads_the_run_chosen_in_the_bar(ws):
    _choose(ws, 'Fibroblast')
    ws.enrichment_page._annotate_direction_and_display(ORA.assign(Term=['fibro', 'T2']))
    st = _source(ws).de_state()
    assert st.run_label == 'Fibroblast (DE workspace)'
    assert list(st.enrichment_results['Term']) == ['fibro', 'T2']
    assert st.cell_type_filter == ('annotation', 'Fibroblast')


def test_a_saved_run_picked_in_figure_export_is_not_given_the_workspaces_results(ws):
    # Endothelial Cell has a saved ORA run; the DE workspace shows Fibroblast's.
    _choose(ws, 'Endothelial Cell')
    ws.enrichment_page._annotate_direction_and_display(ORA.assign(Term=['endo', 'T2']))
    _choose(ws, 'Fibroblast')
    ws.enrichment_page._annotate_direction_and_display(ORA.assign(Term=['fibro', 'T2']))
    src = _source(ws)
    endo = next(r for r in saved.list_gene_de_runs(ws.project_dir)
                if r.cell_type == 'Endothelial Cell')
    src.set_de_run(endo)
    st = src.de_state()
    assert list(st.enrichment_results['Term']) == ['endo', 'T2']
    assert st.cell_type_filter == ('annotation', 'Endothelial Cell')
    assert st.saved_path('enrichment_results').parent == de_downstream_dir(
        ws.project_dir, 'original_counts', 'Endothelial Cell')


def test_a_live_cell_type_run_without_results_in_memory_reads_its_own_folder(ws):
    _choose(ws, 'Fibroblast')
    folder = de_downstream_dir(ws.project_dir, 'original_counts', 'Fibroblast')
    saved.save_enrichment(folder, ORA.assign(Term=['saved', 'T2']))
    ws.enrichment_results = None            # nothing held in memory
    st = _source(ws).de_state()
    assert list(st.enrichment_results['Term']) == ['saved', 'T2']


# -- Provenance's methods ------------------------------------------------------------

def _stage(stage, **params):
    return {'stage': stage, 'timestamp': '', 'params': params}


def test_a_whole_study_de_rerun_leaves_the_cell_types_enrichments_current():
    stages = [_stage('enrichment', method='ora', library='L', cell_type='Fibroblast',
                     counts='original_counts'),
              _stage('enrichment', method='ora', library='L', cell_type=None,
                     counts='original_counts'),
              _stage('gene_de', count_source='raw')]
    final = provenance._final_stages(stages)
    kept = [st['params'].get('cell_type') for st in final if st['stage'] == 'enrichment']
    assert kept == ['Fibroblast']           # the whole-study enrichment is superseded
    gene_de = next(st for st in final if st['stage'] == 'gene_de')
    assert gene_de['_n_superseded'] == {'enrichment': 1}
    # A rerun on other counts does not supersede it either.
    stages[2] = _stage('gene_de', count_source='decontX_counts')
    final = provenance._final_stages(stages)
    assert sum(st['stage'] == 'enrichment' for st in final) == 2


def test_the_methods_keep_the_last_pathway_de_of_each_run():
    stages = [_stage('pathway_de', scoring_method='a', cell_type='Fibroblast', counts='c'),
              _stage('pathway_de', scoring_method='b', cell_type='Endothelial Cell',
                     counts='c'),
              _stage('pathway_de', scoring_method='c', cell_type='Fibroblast', counts='c')]
    final = provenance._final_stages(stages)
    assert sorted((st['params']['cell_type'], st['params']['scoring_method'], st['_n_runs'])
                  for st in final) == [('Endothelial Cell', 'b', 1), ('Fibroblast', 'c', 2)]


# -- DE Methods page --------------------------------------------------------------------

def test_the_methods_sheet_describes_the_chosen_run_and_names_the_others(ws):
    page = ws.enrichment_page
    for ct in ('Endothelial Cell', 'Fibroblast'):
        _choose(ws, ct)
        page._annotate_direction_and_display(ORA)
        page._record_study_enrichment('ora')
    ws.switch_tab(4)                                    # Methods (discovery)
    text = ws.methods_page._text.toPlainText()
    assert 'Results for:' in text and 'Fibroblast' in text
    assert 'Cell-type column:' in text and 'annotation' in text
    assert 'Over-representation:' in text and '2 genes in' in text
    assert 'Endothelial Cell (original counts): enrichment' in text

    _choose(ws, 'Endothelial Cell')                     # refreshed for the new run
    text = ws.methods_page._text.toPlainText()
    assert 'Fibroblast (original counts): enrichment' in text
