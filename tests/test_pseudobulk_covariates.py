"""Covariates in the pseudobulk CSV (ADR-018, stage 1).

The exact case-control permutation has to fit the design the DE run used,
so the run's sample-level covariates travel in the pseudobulk file as
'cov:<name>' columns between 'role' and the genes. The tests cover the
schema the writer produces, the loader that must not read a covariate as a
gene, and the files written before the change.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")

from kosmic.de import batch as B  # noqa: E402
from kosmic.meta_analysis.cc_permutation import (  # noqa: E402
    load_pseudobulk, load_pseudobulk_covariates)
from kosmic.meta_analysis.io import load_pseudobulk_set  # noqa: E402
from kosmic.paths import PSEUDOBULK_COVARIATE_PREFIX  # noqa: E402


def _sample_df(with_sex=True):
    df = pd.DataFrame({
        'sample': ['D1', 'D2', 'C1', 'C2'],
        'condition': ['DCM', 'DCM', 'Donor', 'Donor'],
        'n_cells': [10, 12, 15, 9],
        'total_counts': [50_000, 60_000, 70_000, 40_000],
        'role': ['disease', 'disease', 'control', 'control'],
    })
    if with_sex:
        df['sex_inferred'] = ['male', 'female', 'male', 'female']
    return df


def _write(tmp_path, with_sex=True, name='X_pseudobulk.csv'):
    path = tmp_path / name
    matrix = np.arange(8, dtype=float).reshape(4, 2)
    B.pseudobulk_frame(matrix, _sample_df(with_sex), ['GENE_A', 'GENE_B']).to_csv(path)
    return path


def test_covariates_are_written_after_role_with_the_prefix():
    df = B.pseudobulk_frame(np.zeros((4, 2)), _sample_df(), ['A', 'B'])
    assert list(df.columns) == [
        'condition', 'n_cells', 'total_counts', 'role',
        f'{PSEUDOBULK_COVARIATE_PREFIX}sex_inferred', 'A', 'B']
    assert df['cov:sex_inferred'].tolist() == ['male', 'female', 'male', 'female']


def test_the_prefix_cannot_be_a_gene_symbol_prefix():
    # A colon never occurs in the symbols KOSMIC harmonises to (HGNC, MGI).
    assert PSEUDOBULK_COVARIATE_PREFIX.endswith(':')


def test_load_pseudobulk_does_not_read_a_covariate_as_a_gene(tmp_path):
    path = _write(tmp_path)
    expr, conds, genes, samples, n_cells, roles = load_pseudobulk(path)
    assert genes == ['GENE_A', 'GENE_B']
    assert expr.shape == (4, 2)
    assert roles.tolist() == ['disease', 'disease', 'control', 'control']
    assert samples == ['D1', 'D2', 'C1', 'C2']


def test_covariates_round_trip(tmp_path):
    cov = load_pseudobulk_covariates(_write(tmp_path))
    assert list(cov.columns) == ['sex_inferred']
    assert list(cov.index) == ['D1', 'D2', 'C1', 'C2']
    assert cov['sex_inferred'].tolist() == ['male', 'female', 'male', 'female']


def test_a_file_without_covariates_gives_none(tmp_path):
    path = _write(tmp_path, with_sex=False)
    assert load_pseudobulk_covariates(path) is None
    # and still loads exactly as before
    assert load_pseudobulk(path)[2] == ['GENE_A', 'GENE_B']


def test_load_pseudobulk_set_carries_samples_and_covariates(tmp_path):
    with_cov = _write(tmp_path, True, 'A_pseudobulk.csv')
    without = _write(tmp_path, False, 'B_pseudobulk.csv')
    a, b = load_pseudobulk_set([with_cov, without], labels=['A', 'B'])
    assert a['samples'] == ['D1', 'D2', 'C1', 'C2']
    assert a['covariates']['sex_inferred'].tolist() == [
        'male', 'female', 'male', 'female']
    assert b['covariates'] is None
    assert a['gene_names'] == b['gene_names'] == ['GENE_A', 'GENE_B']


def _adata_with_sex():
    rng = np.random.default_rng(0)
    rows = []
    sex = {'D1': 'male', 'D2': 'female', 'D3': 'male',
           'C1': 'female', 'C2': 'male', 'C3': 'female'}
    for donor, role in (('D1', 'disease'), ('D2', 'disease'), ('D3', 'disease'),
                        ('C1', 'control'), ('C2', 'control'), ('C3', 'control')):
        for _ in range(30):
            rows.append((donor, role, 'Endothelial', sex[donor]))
    obs = pd.DataFrame(rows, columns=['donor', '_role', 'cell_type', 'sex_inferred'])
    obs['condition'] = obs['_role'].map({'disease': 'DCM', 'control': 'Donor'})
    obs.index = [f'cell{i}' for i in range(len(obs))]
    X = rng.poisson(6, size=(len(obs), 30)).astype(np.float32)
    return anndata.AnnData(
        X=X, obs=obs, var=pd.DataFrame(index=[f'GENE{i}' for i in range(30)]))


def test_batch_run_writes_the_runs_covariates_into_the_pseudobulk(tmp_path):
    runs = B.run_de_by_cell_type(
        _adata_with_sex(), 'cell_type', 'donor', 'condition', tmp_path / 'GSE1',
        'GSE1', cell_types=['Endothelial'], de_method='ttest', min_cells=10,
        detection_min_pct=0.0, covariates=['sex_inferred'])
    cov = load_pseudobulk_covariates(runs[0].pseudobulk_path)
    assert cov is not None
    assert set(cov.index) == {'D1', 'D2', 'D3', 'C1', 'C2', 'C3'}
    assert cov.loc['D1', 'sex_inferred'] == 'male'
    # and the genes are still only the genes
    assert load_pseudobulk(runs[0].pseudobulk_path)[2] == [
        f'GENE{i}' for i in range(30)]


def test_batch_run_without_covariates_writes_none(tmp_path):
    runs = B.run_de_by_cell_type(
        _adata_with_sex(), 'cell_type', 'donor', 'condition', tmp_path / 'GSE1',
        'GSE1', cell_types=['Endothelial'], de_method='ttest', min_cells=10,
        detection_min_pct=0.0)
    assert load_pseudobulk_covariates(runs[0].pseudobulk_path) is None


# --- covariates a study's DE run used but its file lacks (ADR-018) -------------

def _with_record(tmp_path, design, with_covariates):
    import json
    path = _write(tmp_path, with_sex=with_covariates, name='S_pseudobulk.csv')
    (tmp_path / 'S_DE_deseq2.json').write_text(
        json.dumps({'method': 'deseq2', 'design': design}), encoding='utf-8')
    return path


def test_a_file_without_the_covariates_its_run_used_is_marked(tmp_path):
    path = _with_record(tmp_path, '~sex_inferred + condition', with_covariates=False)
    (pb,) = load_pseudobulk_set([path], labels=['S'])
    assert pb['covariates'] is None
    assert pb['missing_covariates'] == ['sex_inferred']


def test_a_file_that_carries_them_is_not_marked(tmp_path):
    path = _with_record(tmp_path, '~sex_inferred + condition', with_covariates=True)
    (pb,) = load_pseudobulk_set([path], labels=['S'])
    assert pb['missing_covariates'] == []


def test_a_condition_only_design_needs_no_covariates(tmp_path):
    path = _with_record(tmp_path, '~condition', with_covariates=False)
    (pb,) = load_pseudobulk_set([path], labels=['S'])
    assert pb['missing_covariates'] == []


def test_without_a_run_record_nothing_is_claimed_missing(tmp_path):
    path = _write(tmp_path, with_sex=False, name='S_pseudobulk.csv')
    (pb,) = load_pseudobulk_set([path], labels=['S'])
    assert pb['missing_covariates'] == []


def test_the_exact_engine_refuses_such_a_study(tmp_path):
    from kosmic.meta_analysis.cc_permutation import exact_study_state
    path = _with_record(tmp_path, '~sex_inferred + condition', with_covariates=False)
    (pb,) = load_pseudobulk_set([path], labels=['S'])
    with pytest.raises(ValueError, match="study 'S'.*sex_inferred.*DE run"):
        exact_study_state(pb, seed=1, study_index=0)
