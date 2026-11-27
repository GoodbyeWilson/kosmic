"""The vectorised case-control permutation engine (ADR-019, stage 2).

'deseq2_vec' takes the shuffle, covariates, saved permutations and refusals of the
exact engine and differs only in how one study is fitted on one shuffled labelling.
These tests hold it to the exact engine: the same permuted fit, the same null, the
same saved-run behaviour, and a clear refusal when scipy's private optimiser interface
is not usable.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("pydeseq2")

from kosmic.de.de_analysis import _run_deseq2  # noqa: E402
from kosmic.meta_analysis import cc_permutation as cc  # noqa: E402
from kosmic.meta_analysis import vec_lbfgsb  # noqa: E402
from kosmic.meta_analysis.pooling.reml import reml_fast  # noqa: E402


def _study(name, seed, n_dis=5, n_ctl=5, n_excl=1, n_genes=60):
    rng = np.random.default_rng(seed)
    n = n_dis + n_ctl + n_excl
    role = np.array(['disease'] * n_dis + ['control'] * n_ctl + ['exclude'] * n_excl)
    samples = [f'{name}_{i}' for i in range(n)]
    mean = rng.gamma(2.0, 60.0, n_genes) + 25
    mu = mean[None, :] * rng.uniform(0.7, 1.4, n)[:, None]
    mu[:n_dis, :8] *= 2.2
    counts = rng.negative_binomial(8, 8 / (8 + mu)).astype(float)
    sex = ['male', 'female'] * (n // 2) + ['male'] * (n % 2)
    cov = pd.DataFrame({'sex_inferred': sex}, index=samples)
    return {'dataset_name': name, 'expression': counts,
            'gene_names': [f'G{i}' for i in range(n_genes)], 'role': role,
            'samples': samples, 'covariates': cov, 'conditions': role}


@pytest.fixture(scope='module')
def studies():
    return [_study('A', 1), _study('B', 2), _study('C', 3)]


@pytest.fixture(scope='module')
def observed(studies):
    tables = []
    for pb in studies:
        sdf = pd.DataFrame({'sample': pb['samples'], 'role': pb['role'],
                            'sex_inferred': pb['covariates']['sex_inferred'].values})
        t = _run_deseq2(pb['expression'], sdf, pb['gene_names'], covariates=['sex_inferred'])
        t['dataset'] = pb['dataset_name']
        tables.append(t)
    return {'reml': reml_fast(tables, min_studies=2)}


def _run(studies, observed, engine, n_perms, resume_dir=None, **kw):
    msgs = []
    out = cc.consensus_cc_permutation(
        studies, observed, n_perms=n_perms, seed=11, de_method=engine,
        min_studies=2, n_workers=1, resume_dir=resume_dir, status_cb=msgs.append, **kw)
    return out, msgs


def test_the_engine_has_its_own_name():
    assert cc.normalise_de_method('deseq2_vec') == 'deseq2_vec'
    assert cc.CC_ENGINE_VEC in cc._REFIT_ENGINES and cc.CC_ENGINE_EXACT in cc._REFIT_ENGINES
    assert cc.CC_ENGINE_FAST not in cc._REFIT_ENGINES


def test_a_vectorised_permuted_fit_equals_the_exact_permuted_fit(studies):
    pb = studies[1]
    exact = cc.exact_study_state(pb, seed=5, study_index=1)
    vec = dict(cc.exact_study_state(pb, seed=5, study_index=1), engine=cc.CC_ENGINE_VEC)
    for perm_i in (0, 3):
        lfc_e, p_e, se_e, ok_e = cc._exact_de(exact, perm_i)
        lfc_v, p_v, se_v, ok_v = cc._exact_de(vec, perm_i)
        assert ok_e and ok_v
        np.testing.assert_array_equal(np.isnan(p_v), np.isnan(p_e))
        np.testing.assert_allclose(lfc_v, lfc_e, rtol=1e-5, atol=1e-8, equal_nan=True)
        np.testing.assert_allclose(se_v, se_e, rtol=1e-5, atol=1e-8, equal_nan=True)
        np.testing.assert_allclose(p_v, p_e, rtol=1e-4, equal_nan=True)


def test_a_vectorised_run_gives_the_exact_runs_null(studies, observed, tmp_path):
    _run(studies, observed, 'deseq2_exact', 3, resume_dir=tmp_path / 'exact')
    out, msgs = _run(studies, observed, 'deseq2_vec', 3, resume_dir=tmp_path / 'vec')
    info = out['reml'].attrs['cc_permutation']
    assert info['engine'] == 'deseq2_vec' and info['stratified'] is True
    assert info['approximation'] is False and info['failed_study_fits'] == 0
    assert info['n_permutations'] == 3
    assert any('vectorised' in m for m in msgs)
    for i in range(3):
        with np.load(tmp_path / 'exact' / f'perm_{i:05d}.npz') as a, \
                np.load(tmp_path / 'vec' / f'perm_{i:05d}.npz') as b:
            np.testing.assert_allclose(b['reml'], a['reml'], rtol=1e-4, atol=1e-6)


def test_a_resumed_vectorised_run_equals_an_uninterrupted_one(studies, observed, tmp_path):
    whole, _ = _run(studies, observed, 'deseq2_vec', 3, resume_dir=tmp_path / 'whole')
    _run(studies, observed, 'deseq2_vec', 2, resume_dir=tmp_path / 'split')
    resumed, msgs = _run(studies, observed, 'deseq2_vec', 3, resume_dir=tmp_path / 'split')
    assert any('2 of 3' in m for m in msgs)
    assert resumed['reml'].attrs['cc_permutation']['resumed_permutations'] == 2
    a = whole['reml'].set_index('names')['pvals_pooled'].sort_index()
    b = resumed['reml'].set_index('names')['pvals_pooled'].sort_index()
    np.testing.assert_allclose(a.to_numpy(), b.to_numpy())


def test_saved_permutations_of_the_other_engine_are_not_reused(studies, observed, tmp_path):
    folder = tmp_path / 'cc'
    _run(studies, observed, 'deseq2_exact', 2, resume_dir=folder)
    _out, msgs = _run(studies, observed, 'deseq2_vec', 2, resume_dir=folder)
    assert any('other settings' in m for m in msgs)


def test_a_study_without_the_covariates_its_run_used_is_refused(studies, observed):
    bad = [dict(studies[0], covariates=None, missing_covariates=['sex_inferred'])] + studies[1:]
    with pytest.raises(ValueError, match="study 'A'.*sex_inferred.*DE run"):
        _run(bad, observed, 'deseq2_vec', 2)


def test_the_engine_refuses_to_run_when_the_optimiser_interface_is_unusable(
        studies, observed, monkeypatch):
    monkeypatch.setattr(vec_lbfgsb, '_SETULB', None)
    monkeypatch.setattr(vec_lbfgsb, '_AVAILABLE', None)
    with pytest.raises(RuntimeError, match="private.*L-BFGS-B.*exact engine"):
        _run(studies, observed, 'deseq2_vec', 2)
