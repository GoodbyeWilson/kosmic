"""ADR-015 stage 1: DE outputs are found in both layouts -- flat in
'statistics/' (earlier versions) and in 'statistics/<counts>/<cell type>/'
-- and each run knows the counts it used."""
from __future__ import annotations

import json

import pandas as pd
import pytest

from kosmic.de.saved import find_de_result, list_gene_de_runs, run_counts
from kosmic.meta_analysis.io import discover_de_results
from kosmic.paths import (
    counts_folder_name, counts_folder_of, de_output_dirs, de_run_dir, de_stats_dir,
)

CELL_TYPES = ['Fibroblast', 'Ventricular Cardiomyocyte']


def _table(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({'names': ['A', 'B'], 'logfoldchanges': [1.0, -1.0],
                  'pvals_adj': [0.01, 0.5], 'se': [0.1, 0.1]}).to_csv(path, index=False)


def _record(path, **fields):
    path.with_suffix('.json').write_text(json.dumps(fields), encoding='utf-8')


def _study(root, layout, counts_source='decontX_counts', record=True):
    """A study 'S1' with a DE run per cell type, written in one layout."""
    study = root / 'S1'
    for ct in CELL_TYPES:
        slug = ct.replace(' ', '_')
        folder = (de_run_dir(study, counts_source, ct) if layout == 'new'
                  else de_stats_dir(study))
        de = folder / f'S1_{slug}_DE_deseq2.csv'
        _table(de)
        _table(folder / f'S1_{slug}_pseudobulk.csv')
        if record:
            _record(de, cell_type=ct, counts_source=counts_source)
    return study


def test_counts_folder_names():
    assert counts_folder_name(None) == counts_folder_name('raw') == 'original_counts'
    assert counts_folder_name('decontX_counts') == 'decontX_counts'


def test_run_folder_and_its_counts(tmp_path):
    d = de_run_dir(tmp_path / 'S1', 'raw', 'Ventricular Cardiomyocyte')
    assert d.parts[-3:] == ('statistics', 'original_counts', 'Ventricular_Cardiomyocyte')
    assert counts_folder_of(d / 'x.csv') == 'original_counts'
    assert counts_folder_of(de_stats_dir(tmp_path / 'S1') / 'x.csv') is None
    assert de_run_dir(tmp_path / 'S1', 'decontX_counts').name == 'decontX_counts'


def test_output_folders_cover_both_layouts(tmp_path):
    study = _study(tmp_path, 'new')
    _table(de_stats_dir(study) / 'S1_DE_deseq2.csv')          # an old flat file too
    names = [d.name for d in de_output_dirs(de_stats_dir(study))]
    assert names == ['statistics', 'decontX_counts', 'Fibroblast', 'Ventricular_Cardiomyocyte']


@pytest.mark.parametrize('layout', ['old', 'new'])
def test_both_layouts_list_the_same_runs(tmp_path, layout):
    runs = list_gene_de_runs(_study(tmp_path, layout))
    assert [(r.name, r.cell_type, r.method) for r in runs] == [
        ('S1_Fibroblast', 'Fibroblast', 'deseq2'),
        ('S1_Ventricular_Cardiomyocyte', 'Ventricular Cardiomyocyte', 'deseq2')]
    assert {r.counts for r in runs} == {'decontX_counts'}     # folder, or record


def test_counts_are_unknown_without_folder_or_record(tmp_path):
    runs = list_gene_de_runs(_study(tmp_path, 'old', record=False))
    assert {r.counts for r in runs} == {None}
    assert all('counts' not in r.label for r in runs)


def test_two_count_sources_side_by_side(tmp_path):
    study = _study(tmp_path, 'new', 'raw')
    _study(tmp_path, 'new', 'decontX_counts')
    runs = list_gene_de_runs(study)
    assert len(runs) == 4
    assert sorted({r.counts for r in runs}) == ['decontX_counts', 'original_counts']
    only = list_gene_de_runs(study, counts='decontX_counts')
    assert {r.counts for r in only} == {'decontX_counts'} and len(only) == 2
    assert 'decontX counts' in only[0].label


def test_counts_come_from_the_record_in_the_flat_layout(tmp_path):
    study = _study(tmp_path, 'old', counts_source='raw')
    de = next(de_stats_dir(study).glob('*_DE_deseq2.csv'))
    assert run_counts(de) == 'original_counts'


def test_whole_study_run_is_found_in_either_layout(tmp_path):
    study = tmp_path / 'S1'
    _table(de_stats_dir(study) / 'S1_DE_deseq2.csv')
    assert find_de_result(study, 'S1', 'deseq2') == de_stats_dir(study) / 'S1_DE_deseq2.csv'
    new = de_run_dir(study, 'decontX_counts') / 'S1_DE_deseq2.csv'
    _table(new)
    assert find_de_result(study, 'S1', 'deseq2', counts='decontX_counts') == new
    assert find_de_result(study, 'S2', 'deseq2') is None


@pytest.mark.parametrize('layout', ['old', 'new'])
def test_meta_discovery_finds_runs_and_their_counts(tmp_path, layout):
    _study(tmp_path, layout)
    found = [e for e in discover_de_results(tmp_path) if e['file_type'] == 'de_results']
    assert sorted(e['accession'] for e in found) == ['S1_Fibroblast',
                                                     'S1_Ventricular_Cardiomyocyte']
    assert {e['counts'] for e in found} == {'decontX_counts'}


# -- Stage 2: writers and summary ------------------------------------------

def _adata(decontx=False):
    """Two cell types, three donors per arm, and one donor too small to keep."""
    import anndata
    import numpy as np

    rng = np.random.default_rng(0)
    rows = []
    for role, donors in (('disease', ['D1', 'D2', 'D3']), ('control', ['C1', 'C2', 'C3'])):
        for d in donors:
            for ct in ('Endothelial', 'Fibroblast'):
                rows += [(d, role, ct)] * 30
    rows += [('C4', 'control', 'Endothelial')] * 5        # below min_cells: excluded
    obs = pd.DataFrame(rows, columns=['donor', '_role', 'cell_type'])
    obs['condition'] = np.where(obs['_role'] == 'disease', 'DCM', 'Donor')
    obs.index = [f'c{i}' for i in range(len(obs))]
    x = rng.poisson(6, size=(len(obs), 40)).astype(np.float32)
    a = anndata.AnnData(X=x, obs=obs, var=pd.DataFrame(index=[f'G{i}' for i in range(40)]))
    if decontx:
        a.layers['decontX_counts'] = np.floor(x * 0.9)
    return a


def _batch(tmp_path, **kw):
    from kosmic.de import batch as B
    return B.run_de_by_cell_type(
        _adata(decontx='counts_layer' in kw), 'cell_type', 'donor', 'condition',
        tmp_path / 'GSE1', 'GSE1', de_method='ttest', min_cells=10,
        detection_min_pct=0.0, run_record={'counts_source': 'raw'}, **kw)


def test_batch_writes_each_cell_type_into_its_counts_and_cell_type_folder(tmp_path):
    runs = _batch(tmp_path)
    for r in runs:
        folder = de_run_dir(tmp_path / 'GSE1', 'raw', r.cell_type)
        written = [r.de_path, r.de_path.with_suffix('.json'), r.pseudobulk_path,
                   r.significant_path]
        assert all(p.parent == folder and p.is_file() for p in written)
        assert [p.name for p in written[2:]] == [f'GSE1_{r.slug}_pseudobulk.csv',
                                                 f'GSE1_{r.slug}_significant_genes.csv']
        record = json.loads(r.de_path.with_suffix('.json').read_text())
        assert record['counts_source'] == 'raw' and record['saved'] == r.saved
    # Nothing is left in the flat layout.
    assert not list(de_stats_dir(tmp_path / 'GSE1').glob('*.csv'))


def test_a_corrected_layer_gets_its_own_counts_folder(tmp_path):
    runs = _batch(tmp_path, counts_layer='decontX_counts')
    assert {r.de_path.parent.parent.name for r in runs} == {'decontX_counts'}
    assert {r.counts for r in runs} == {'decontX_counts'}


def test_the_summary_describes_each_run(tmp_path):
    from kosmic.de import batch as B
    runs = _batch(tmp_path)
    s = B.batch_summary(runs).set_index('cell_type')
    endo = s.loc['Endothelial']
    assert endo['counts'] == 'original_counts' and endo['folder'] == 'Endothelial'
    assert (endo['n_donors_control'], endo['n_donors_disease'],
            endo['n_donors_excluded']) == (3, 3, 1)            # C4 has 5 cells
    assert endo['median_cells_per_donor'] == 30
    assert endo['n_up'] + endo['n_down'] == endo['n_significant']
    assert endo['min_cells'] == 10 and endo['fdr'] == 0.05
    assert not any(':' in str(v) and '\\' in str(v) for v in s['folder'])   # no absolute paths
    out = B.summary_path(tmp_path / 'GSE1', 'GSE1', 'raw')
    assert out.parent == de_run_dir(tmp_path / 'GSE1', 'raw')


def test_written_runs_are_read_back_after_the_project_is_moved(tmp_path):
    import shutil
    _batch(tmp_path)
    moved = tmp_path / 'elsewhere' / 'GSE1'
    shutil.move(str(tmp_path / 'GSE1'), str(moved))
    runs = list_gene_de_runs(moved)
    assert [r.cell_type for r in runs] == ['Endothelial', 'Fibroblast']
    assert {r.counts for r in runs} == {'original_counts'}
    found = [e for e in discover_de_results(moved.parent) if e['file_type'] == 'de_results']
    assert len(found) == 2


# -- Stage 3: significance is a setting of each run ------------------------

def test_the_gate_is_off_by_default_and_not_in_the_config():
    import tomllib
    from pathlib import Path

    import kosmic
    assert kosmic.DEFAULT_LFC_THRESHOLD == 0.0
    cfg = tomllib.loads((Path(kosmic.__file__).parent.parent / 'config.toml').read_text())
    assert 'lfc_threshold' not in cfg['significance']


def _n_fdr(r, lfc=0.0):
    t = r.de_results
    return int(((t['pvals_adj'] < 0.05) & (t['logfoldchanges'].abs() > lfc)).sum())


def test_a_run_with_the_gate_off_counts_on_fdr_alone(tmp_path):
    from kosmic.de import batch as B
    runs = _batch(tmp_path)
    s = B.batch_summary(runs).set_index('cell_type')
    for r in runs:
        assert r.n_significant == _n_fdr(r)
        assert len(pd.read_csv(r.significant_path)) == r.n_significant
        record = json.loads(r.de_path.with_suffix('.json').read_text())
        assert (record['fdr'], record['lfc_threshold']) == (0.05, 0.0)
        assert s.loc[r.cell_type, 'lfc_threshold'] == 0.0


def test_a_chosen_gate_is_recorded_and_used(tmp_path):
    runs = _batch(tmp_path, lfc=0.5)
    for r in runs:
        assert r.n_significant == _n_fdr(r, 0.5)
        assert len(pd.read_csv(r.significant_path)) == r.n_significant
        assert json.loads(r.de_path.with_suffix('.json').read_text())['lfc_threshold'] == 0.5


def test_a_run_without_a_recorded_rule_is_read_as_assumed():
    from kosmic.de.saved import run_significance
    assert run_significance({'counts_source': 'raw'}) == (0.05, 0.25, True)
    assert run_significance(None) == (0.05, 0.25, True)
    assert run_significance({'fdr': 0.1, 'lfc_threshold': 0.0}) == (0.1, 0.0, False)


@pytest.fixture(scope='module')
def app():
    import os
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _spin(value):
    from PyQt6.QtWidgets import QDoubleSpinBox
    s = QDoubleSpinBox()
    s.setRange(0, 10)
    s.setDecimals(4)
    s.setValue(value)
    return s


def test_reopening_a_saved_run_sets_the_de_page_to_its_rule(app, tmp_path):
    from types import SimpleNamespace

    from kosmic.gui.de_analysis.pages.gene_de_page import GeneDEPage
    de = tmp_path / 'S1_DE_deseq2.csv'
    _table(de)
    _record(de, fdr=0.05, lfc_threshold=0.5)
    page = SimpleNamespace(pval_filter=_spin(0.05), fc_filter=_spin(0.0),
                           ws=SimpleNamespace(log_message=SimpleNamespace(emit=print)))
    GeneDEPage._apply_run_significance(page, de)
    assert page.fc_filter.value() == 0.5 and page.ws.significance == (0.05, 0.5, False)
    _record(de, counts_source='raw')                      # an older run: no rule
    GeneDEPage._apply_run_significance(page, de)
    assert page.fc_filter.value() == 0.25 and page.ws.significance[2] is True


def test_enrichment_starts_from_the_runs_rule_once_per_run(app):
    from types import SimpleNamespace

    from kosmic.gui.de_analysis.pages.enrichment_page import EnrichmentPage
    page = SimpleNamespace(_fdr_spin=_spin(0.05), _lfc_spin=_spin(0.25),
                           ws=SimpleNamespace(significance=(0.05, 0.0, False)))
    EnrichmentPage._sync_significance(page)
    assert page._lfc_spin.value() == 0.0
    page._lfc_spin.setValue(1.0)                          # typed by the user
    EnrichmentPage._sync_significance(page)
    assert page._lfc_spin.value() == 1.0                  # kept: same run
    page.ws.significance = (0.05, 0.5, False)             # another run shown
    EnrichmentPage._sync_significance(page)
    assert page._lfc_spin.value() == 0.5


def test_figure_gate_starts_from_the_runs_rule(app, tmp_path):
    from kosmic.de import saved
    from kosmic.gui.figure_export.pages.volcano import VolcanoPage
    from tests.study_folder import figure_workspace, make_study_folder, write_de_outputs

    study = write_de_outputs(make_study_folder(tmp_path))
    ws = figure_workspace(study=study)
    run = saved.list_gene_de_runs(study)[0]
    ws.source.set_de_run(run)
    _record(run.path, lfc_threshold=0.5, fdr=0.05)
    page = VolcanoPage(ws)
    page._controls.blockSignals(True)
    page._sync_gate_from_run()
    assert page._controls.lfc_gate.value() == 0.5
    page._controls.lfc_gate.setValue(0.0)                 # set here: kept for this run
    page._sync_gate_from_run()
    assert page._controls.lfc_gate.value() == 0.0


def test_pseudobulk_is_the_one_beside_each_run(tmp_path, monkeypatch):
    from kosmic.gui.meta_analysis import pseudobulk_discovery as pd_disc
    monkeypatch.setattr(pd_disc, '_load_persisted_cache', lambda: {})
    monkeypatch.setattr(pd_disc, '_save_persisted_cache', lambda cache: None)
    pd_disc.clear_memory_cache()
    _study(tmp_path, 'new')
    found = [e for e in discover_de_results(tmp_path) if e['file_type'] == 'de_results']
    datasets = [{'path': str(e['path']), 'name': e['accession']} for e in found]
    paths = pd_disc.find_pseudobulk_paths(datasets, str(tmp_path))
    assert [p.split('\\')[-1].split('/')[-1] for p in paths] == [
        f"{e['accession']}_pseudobulk.csv" for e in found]
