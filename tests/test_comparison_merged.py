"""Merged pseudobulk: finding the atlas's saved run and fitting it (ADR-017 §3)."""
import json

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")

from kosmic.comparison.design import Design  # noqa: E402
from kosmic.comparison.merged import find_saved_merged, run_merged  # noqa: E402

CT = 'Ventricular Cardiomyocyte'


def _design(**kw):
    base = dict(covariates=['sex'], mapping={'sex_inferred': 'sex'}, sample_column='sample',
                cell_type_column='cell_type', cell_type=CT, counts='original_counts',
                min_cells=5, min_counts=0)
    return Design(**{**base, **kw})


def _save_run(project, counts, covariates, method='deseq2'):
    folder = project / '_master/results/de_analysis/statistics' / counts / 'Ventricular_Cardiomyocyte'
    folder.mkdir(parents=True)
    path = folder / f'_master_Ventricular_Cardiomyocyte_DE_{method}.csv'
    pd.DataFrame({'names': ['A'], 'logfoldchanges': [1.0], 'pvals': [0.01],
                  'pvals_adj': [0.02], 'se': [0.3]}).to_csv(path, index=False)
    path.with_suffix('.json').write_text(json.dumps({
        'cell_type': CT, 'covariates': covariates, 'sample_column': 'sample',
        'min_cells': 5, 'min_counts': 0, 'counts_source': 'raw'}), encoding='utf-8')
    return path


def test_a_matching_saved_run_is_found(tmp_path):
    _save_run(tmp_path, 'original_counts', ['sex', 'study'])
    saved = find_saved_merged(tmp_path, _design())
    assert saved is not None and saved.matches
    assert saved.read()['names'].tolist() == ['A']


def test_a_saved_run_without_study_is_found_but_does_not_match(tmp_path):
    _save_run(tmp_path, 'original_counts', ['sex'])
    saved = find_saved_merged(tmp_path, _design())
    assert saved is not None and not saved.matches
    assert any('study' in d for d in saved.differences)


def test_a_run_on_other_counts_is_not_the_selections(tmp_path):
    _save_run(tmp_path, 'decontX_counts', ['sex', 'study'])
    assert find_saved_merged(tmp_path, _design()) is None


def test_no_atlas_folder_means_no_saved_run(tmp_path):
    assert find_saved_merged(tmp_path, _design()) is None


def _write_atlas(project, seed=0):
    """Two studies, 3 disease and 3 control donors each, both sexes in each arm."""
    rng = np.random.default_rng(seed)
    rows = []
    for study in ('S1', 'S2'):
        for role in ('disease', 'control'):
            for i in range(3):
                for _ in range(20):
                    rows.append((study, f'{study}_{role}{i}', role, ('female', 'male')[i % 2], CT))
    obs = pd.DataFrame(rows, columns=['study', 'sample', '_role', 'sex', 'cell_type'])
    obs.index = [f'c{i}' for i in range(len(obs))]
    X = rng.poisson(8, size=(len(obs), 60)).astype(np.float32)
    X[obs['_role'].values == 'disease', 0] *= 4   # one clearly changed gene
    a = anndata.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f'G{i}' for i in range(60)]))
    path = project / '_master/processed_data/master.h5ad'
    path.parent.mkdir(parents=True)
    a.write_h5ad(path)


def test_run_merged_fits_the_design_and_saves_it_for_the_atlas(tmp_path):
    _write_atlas(tmp_path)
    saved = run_merged(tmp_path, _design())
    assert saved.matches
    assert saved.record['design'] == '~sex + study + condition'
    table = saved.read()
    assert 'G0' in set(table['names']) and len(table) > 10
    assert saved.path.parent.name == 'Ventricular_Cardiomyocyte'


def test_run_merged_refuses_a_covariate_the_atlas_lacks(tmp_path):
    _write_atlas(tmp_path)
    with pytest.raises(ValueError, match='batch'):
        run_merged(tmp_path, _design(unresolved=['batch']))


def test_run_merged_refuses_counts_the_atlas_does_not_hold(tmp_path):
    _write_atlas(tmp_path)
    with pytest.raises(ValueError, match='decontX counts'):
        run_merged(tmp_path, _design(counts='decontX_counts'))
