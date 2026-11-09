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
