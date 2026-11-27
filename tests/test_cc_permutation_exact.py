"""The exact case-control permutation (ADR-018, stage 4).

Covers what makes the engine exact and what makes a run usable: the shuffle
(roles move only among the donors in the contrast, within covariate groups),
that a permuted fit is the study's own DE on the shuffled labels, the
permutation p-value, resume and cancel, and the refusals that stop a run
which could not reproduce the study's model.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("pydeseq2")

from kosmic.de.de_analysis import _run_deseq2, fit_deseq2_wald  # noqa: E402
from kosmic.meta_analysis import cc_permutation as cc  # noqa: E402
from kosmic.meta_analysis.pooling.reml import reml_fast  # noqa: E402


# --- the shuffle -----------------------------------------------------------

ROLES = np.array(['disease'] * 4 + ['control'] * 4 + ['exclude'] * 2)


def test_only_donors_in_the_contrast_change_role():
    rng = np.random.default_rng(0)
    out = cc.shuffle_roles(ROLES, None, rng)
    assert (out[ROLES == 'exclude'] == 'exclude').all()
    assert sorted(out[:8]) == sorted(ROLES[:8])        # same arms, same sizes


def test_shuffle_within_strata_keeps_each_stratum_s_table():
    strata = np.array([0, 0, 1, 1, 0, 0, 1, 1, 0, 1])
    for seed in range(20):
        out = cc.shuffle_roles(ROLES, strata, np.random.default_rng(seed))
        for s in (0, 1):
            sel = (strata == s) & (ROLES != 'exclude')
            assert sorted(out[sel]) == sorted(ROLES[sel])


def test_the_shuffle_is_reproducible_and_does_move_roles():
    a = cc.shuffle_roles(ROLES, None, np.random.default_rng([7, 3, 1]))
    b = cc.shuffle_roles(ROLES, None, np.random.default_rng([7, 3, 1]))
    assert (a == b).all()
    moved = [(cc.shuffle_roles(ROLES, None, np.random.default_rng(s)) != ROLES).any()
             for s in range(10)]
    assert any(moved)


# --- the permutation p-value ------------------------------------------------

def test_permutation_p_counts_the_observed_statistic_as_one_draw():
    null = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    p = cc.permutation_p(null, [0.0, 3.0, 5.0, 6.0])
    # (1 + #null >= obs) / (1 + 5)
    np.testing.assert_allclose(p, [6 / 6, 4 / 6, 2 / 6, 1 / 6])


def test_permutation_p_is_never_zero():
    assert cc.permutation_p(np.arange(100.0), [1e9])[0] == pytest.approx(1 / 101)


# --- engine names ------------------------------------------------------------

def test_engine_names():
    assert cc.normalise_de_method('deseq2') == 'deseq2_fast'
    assert cc.normalise_de_method('deseq2_fast') == 'deseq2_fast'
    assert cc.normalise_de_method('deseq2_exact') == 'deseq2_exact'
    assert cc.normalise_de_method('welch_cpm') == 'welch_cpm'


# --- fixtures: three small studies ---------------------------------------------

def _study(name, seed, n_dis=4, n_ctl=4, n_excl=1, n_genes=40, sex=None):
    rng = np.random.default_rng(seed)
    n = n_dis + n_ctl + n_excl
    role = np.array(['disease'] * n_dis + ['control'] * n_ctl + ['exclude'] * n_excl)
    samples = [f'{name}_{i}' for i in range(n)]
    mean = rng.gamma(2.0, 60.0, n_genes) + 25
    mu = mean[None, :] * rng.uniform(0.7, 1.4, n)[:, None]
    mu[:n_dis, :8] *= 2.2
    counts = rng.negative_binomial(8, 8 / (8 + mu)).astype(float)
    sex = sex or ['male', 'female'] * (n // 2) + ['male'] * (n % 2)
    cov = pd.DataFrame({'sex_inferred': sex}, index=samples)
    return {'dataset_name': name, 'expression': counts,
            'gene_names': [f'G{i}' for i in range(n_genes)], 'role': role,
            'samples': samples, 'covariates': cov, 'conditions': role}


@pytest.fixture(scope='module')
def studies():
    return [_study('A', 1), _study('B', 2), _study('C', 3)]


def _observed(studies):
    tables = []
    for pb in studies:
        sdf = pd.DataFrame({'sample': pb['samples'], 'role': pb['role'],
                            'sex_inferred': pb['covariates']['sex_inferred'].values})
        t = _run_deseq2(pb['expression'], sdf, pb['gene_names'],
                        covariates=['sex_inferred'])
        t['dataset'] = pb['dataset_name']
        tables.append(t)
    return {'reml': reml_fast(tables, min_studies=2)}


@pytest.fixture(scope='module')
def observed(studies):
    return _observed(studies)


# --- the state and its refusals ----------------------------------------------

def test_state_records_the_covariates_the_design_will_use(studies):
    st = cc.exact_study_state(studies[0], seed=1, study_index=0)
    assert st['covariates'] == ['sex_inferred']
    assert st['usable'] == ['sex_inferred']
    assert st['strata'] is not None
    assert st['n_labelings'] >= 2


def test_state_without_covariates_shuffles_among_all_donors_in_the_contrast(studies):
    pb = dict(studies[0], covariates=None)
    st = cc.exact_study_state(pb, seed=1, study_index=0)
    assert st['covariates'] == [] and st['strata'] is None
    assert st['n_labelings'] == 70                      # C(8, 4)


def test_a_run_that_used_covariates_the_file_lacks_is_refused(studies):
    pb = dict(studies[0], covariates=None, missing_covariates=['sex_inferred'])
    with pytest.raises(ValueError, match="study 'A'.*sex_inferred.*DE run"):
        cc.exact_study_state(pb, seed=1, study_index=0)


def test_a_covariate_collinear_with_condition_is_dropped_as_in_the_study_s_own_fit():
    # every disease donor male, every control donor female: the DE drops sex
    # from the design, so the permutation shuffles among all donors
    pb = _study('Z', 9, n_excl=0, sex=['male'] * 4 + ['female'] * 4)
    st = cc.exact_study_state(pb, seed=1, study_index=0)
    assert st['covariates'] == ['sex_inferred'] and st['usable'] == []
    assert st['strata'] is None and st['n_labelings'] == 70


def test_no_movable_donor_is_refused():
    # Each covariate has both roles at every level, but no sex-by-batch
    # combination does: within a group nobody can change role.
    pb = _study('Z', 9, n_excl=0)
    sex = ['male', 'male', 'female', 'female', 'male', 'male', 'female', 'female']
    batch = ['a', 'a', 'b', 'b', 'b', 'b', 'a', 'a']
    pb['covariates'] = pd.DataFrame({'sex_inferred': sex, 'batch': batch},
                                    index=pb['samples'])
    with pytest.raises(ValueError, match="nothing to shuffle"):
        cc.exact_study_state(pb, seed=1, study_index=0)


def test_a_study_needs_both_arms():
    pb = _study('Z', 9, n_ctl=0)
    with pytest.raises(ValueError, match="at least one disease and one control"):
        cc.exact_study_state(pb, seed=1, study_index=0)


# --- a permuted fit is the study's own DE on the shuffled labels ---------------

def test_a_permuted_fit_is_the_study_s_own_fit_on_the_shuffled_roles(studies):
    pb = studies[1]
    st = cc.exact_study_state(pb, seed=5, study_index=1)
    lfc, p, se, ok = cc._exact_de(st, perm_i=2)
    assert ok

    rng = np.random.default_rng([5, 2, 1])
    roles = cc.shuffle_roles(pb['role'], st['strata'], rng)
    sdf = pd.DataFrame({'sample': pb['samples'], 'role': roles,
                        'sex_inferred': pb['covariates']['sex_inferred'].values})
    ref = fit_deseq2_wald(pb['expression'], sdf, pb['gene_names'],
                          covariates=['sex_inferred'], quiet=True).ds.results_df
    pos = [pb['gene_names'].index(g) for g in ref.index]
    np.testing.assert_allclose(lfc[pos], ref['log2FoldChange'], rtol=1e-8)
    np.testing.assert_allclose(se[pos], ref['lfcSE'], rtol=1e-8)
    np.testing.assert_allclose(p[pos], ref['pvalue'], rtol=1e-8, equal_nan=True)


# --- a run ----------------------------------------------------------------------

def _run(studies, observed, n_perms, resume_dir=None, **kw):
    msgs = []
    out = cc.consensus_cc_permutation(
        studies, observed, n_perms=n_perms, seed=11, de_method='deseq2_exact',
        min_studies=2, n_workers=1, resume_dir=resume_dir,
        status_cb=msgs.append, **kw)
    return out, msgs


def test_an_exact_run_calibrates_and_records_how(studies, observed, tmp_path):
    out, msgs = _run(studies, observed, 3, resume_dir=tmp_path / 'cc')
    df = out['reml']
    info = df.attrs['cc_permutation']
    assert info['engine'] == 'deseq2_exact' and info['stratified'] is True
    assert info['n_permutations'] == 3 and info['approximation'] is False
    assert info['covariates'] == {'A': ['sex_inferred'], 'B': ['sex_inferred'],
                                  'C': ['sex_inferred']}
    assert info['failed_study_fits'] == 0
    p = df['pvals_pooled'].to_numpy()
    assert np.all(p >= info['smallest_calibrated_p'] - 1e-15) and np.all(p <= 1)
    assert 'pvals_analytical' in df.columns
    assert any('exact' in m for m in msgs)
    assert len(list((tmp_path / 'cc').glob('perm_*.npz'))) == 3


def test_a_resumed_run_equals_an_uninterrupted_one(studies, observed, tmp_path):
    whole, _ = _run(studies, observed, 3, resume_dir=tmp_path / 'whole')
    _run(studies, observed, 2, resume_dir=tmp_path / 'split')
    resumed, msgs = _run(studies, observed, 3, resume_dir=tmp_path / 'split')
    assert any('2 of 3' in m for m in msgs)
    assert resumed['reml'].attrs['cc_permutation']['resumed_permutations'] == 2
    a = whole['reml'].set_index('names')['pvals_pooled']
    b = resumed['reml'].set_index('names')['pvals_pooled']
    np.testing.assert_allclose(a.sort_index().to_numpy(), b.sort_index().to_numpy())


def test_a_cancelled_run_keeps_what_it_finished(studies, observed, tmp_path):
    folder = tmp_path / 'cc'
    with pytest.raises(cc.CCPermutationCancelled):
        _run(studies, observed, 3, resume_dir=folder, cancelled=lambda: True)
    # the check runs before each permutation, so nothing finished yet ...
    assert len(list(folder.glob('perm_*.npz'))) == 0
    # ... and the same folder serves a later complete run
    out, _ = _run(studies, observed, 2, resume_dir=folder)
    assert out['reml'].attrs['cc_permutation']['n_permutations'] == 2


def test_saved_permutations_of_other_settings_are_not_reused(studies, observed, tmp_path):
    folder = tmp_path / 'cc'
    _run(studies, observed, 2, resume_dir=folder)
    _out, msgs = _run(studies, observed, 2, resume_dir=folder, stratify=False)
    assert any('other settings' in m for m in msgs)


# --- the fast engine's name -------------------------------------------------------

def test_the_fast_engine_is_selected_by_its_own_name(studies, observed):
    msgs = []
    out = cc.consensus_cc_permutation(
        studies, observed, n_perms=2, seed=3, de_method='deseq2_fast',
        min_studies=2, n_workers=1, status_cb=msgs.append)
    assert any('fast, an approximation' in m for m in msgs)
    assert out['reml'].attrs['cc_permutation']['approximation'] is True


def test_worker_processes_give_the_same_result_as_one_process(studies, observed, tmp_path):
    one, _ = _run(studies, observed, 3, resume_dir=tmp_path / 'one')
    msgs = []
    many = cc.consensus_cc_permutation(
        studies, observed, n_perms=3, seed=11, de_method='deseq2_exact',
        min_studies=2, n_workers=2, resume_dir=tmp_path / 'many',
        status_cb=msgs.append)
    assert any('2 processes' in m for m in msgs)
    a = one['reml'].set_index('names')['pvals_pooled'].sort_index()
    b = many['reml'].set_index('names')['pvals_pooled'].sort_index()
    np.testing.assert_allclose(a.to_numpy(), b.to_numpy())
