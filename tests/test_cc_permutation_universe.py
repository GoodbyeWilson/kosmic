"""Calibrating against a gene universe, and permutation fits without a second filter.

'gene_universe': every study is still fitted on all its genes, but only the
genes the observed results test are pooled into the null and calibrated.
'refilter=False': permutation fits skip the DE expression filter, for
pseudobulk whose genes already passed it. Both are off by default, and runs
with the defaults keep the signature their saved permutations were written
under.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("pydeseq2")

from kosmic.de.de_analysis import _run_deseq2  # noqa: E402
from kosmic.meta_analysis import cc_permutation as cc  # noqa: E402
from kosmic.meta_analysis.pooling.reml import reml_fast  # noqa: E402
from kosmic.meta_analysis.pooling.sumrank import sumrank_fast  # noqa: E402

N_GENES = 60
LOW = 'G_low'


def _study(name, seed, n_dis=5, n_ctl=5, n_excl=1, n_genes=N_GENES):
    rng = np.random.default_rng(seed)
    n = n_dis + n_ctl + n_excl
    role = np.array(['disease'] * n_dis + ['control'] * n_ctl + ['exclude'] * n_excl)
    samples = [f'{name}_{i}' for i in range(n)]
    mean = rng.gamma(2.0, 60.0, n_genes) + 25
    mu = mean[None, :] * rng.uniform(0.7, 1.4, n)[:, None]
    mu[:n_dis, :8] *= 2.2
    counts = rng.negative_binomial(8, 8 / (8 + mu)).astype(float)
    # one gene below the filter's CPM cutoff (~10 counts at the median library
    # size) but above its summed-count rule (15), so only the CPM rule drops it
    low = rng.poisson(2.0, (n, 1)).astype(float) + 1
    counts = np.hstack([counts, low])
    sex = ['male', 'female'] * (n // 2) + ['male'] * (n % 2)
    cov = pd.DataFrame({'sex_inferred': sex}, index=samples)
    return {'dataset_name': name, 'expression': counts,
            'gene_names': [f'G{i}' for i in range(n_genes)] + [LOW], 'role': role,
            'samples': samples, 'covariates': cov, 'conditions': role}


@pytest.fixture(scope='module')
def studies():
    return [_study('A', 1), _study('B', 2), _study('C', 3)]


@pytest.fixture(scope='module')
def tables(studies):
    out = []
    for pb in studies:
        sdf = pd.DataFrame({'sample': pb['samples'], 'role': pb['role'],
                            'sex_inferred': pb['covariates']['sex_inferred'].values})
        t = _run_deseq2(pb['expression'], sdf, pb['gene_names'], covariates=['sex_inferred'])
        t['dataset'] = pb['dataset_name']
        out.append(t)
    return out


def _observed(tables, genes=None):
    if genes is not None:
        tables = [t[t['names'].isin(genes)] for t in tables]
    return {'reml': reml_fast(tables, min_studies=2),
            'sumrank': sumrank_fast(tables, min_studies=2)}


def _run(studies, observed, n_perms=3, resume_dir=None, engine='deseq2_vec', **kw):
    msgs = []
    out = cc.consensus_cc_permutation(
        studies, observed, n_perms=n_perms, seed=11, de_method=engine,
        min_studies=2, n_workers=1, resume_dir=resume_dir, status_cb=msgs.append, **kw)
    return out, msgs


def _p(frame):
    return frame.set_index('names')['pvals_pooled'].sort_index()


def test_the_defaults_leave_the_saved_signature_unchanged(studies, tables, tmp_path):
    out, _ = _run(studies, _observed(tables), resume_dir=tmp_path / 'cc')
    sig = json.loads((tmp_path / 'cc' / 'manifest.json').read_text())
    assert 'gene_universe' not in sig and 'refilter' not in sig
    info = out['reml'].attrs['cc_permutation']
    assert info['gene_universe'] is None and info['refilter'] is True
    assert info['rows_outside_universe_dropped'] == 0


def test_a_universe_of_every_pooled_gene_changes_nothing(studies, tables):
    obs = _observed(tables)
    base, _ = _run(studies, obs)
    every = set(obs['reml']['names']) | set(obs['sumrank']['names'])
    same, _ = _run(studies, obs, gene_universe=every)
    for k in obs:
        np.testing.assert_allclose(_p(same[k]).to_numpy(), _p(base[k]).to_numpy())
        assert (same[k].attrs['cc_permutation']['null_values']
                == base[k].attrs['cc_permutation']['null_values'])


def test_only_universe_genes_enter_the_null_and_are_calibrated(studies, tables):
    universe = [f'G{i}' for i in range(30)]
    obs = _observed(tables)                       # observed tables still hold every gene
    out, _ = _run(studies, obs, n_perms=4, gene_universe=universe)
    for k, frame in out.items():
        assert set(frame['names']) <= set(universe)
        info = frame.attrs['cc_permutation']
        assert info['gene_universe'] == 30
        assert info['rows_outside_universe_dropped'] == len(obs[k]) - len(frame)
        assert info['null_values'] <= 4 * 30
        assert info['smallest_calibrated_p'] >= 1 / (4 * 30 + 1)


def test_rank_methods_rank_among_universe_genes_only(studies, tables):
    """With a universe, permuted genes are ranked among the universe genes (as
    the observed SumRank was), not among every fitted gene, so the null and
    the calibrated p-values differ from a run without one."""
    universe = [f'G{i}' for i in range(20)]
    inside, _ = _run(studies, _observed(tables, universe), gene_universe=universe)
    unrestricted, _ = _run(studies, _observed(tables, universe))
    # same observed statistics, different nulls: the unrestricted run ranks each
    # permuted gene among all 61, the universe run among the 20
    a = _p(inside['sumrank'])
    b = _p(unrestricted['sumrank']).reindex(a.index)
    assert not np.allclose(a.to_numpy(), b.to_numpy())


def test_without_refiltering_a_sparse_gene_is_still_fitted(studies):
    pb = studies[0]
    state = dict(cc.exact_study_state(pb, seed=5, study_index=0), engine=cc.CC_ENGINE_VEC)
    i = list(state['genes']).index(LOW)
    _, p_default, _, ok = cc._exact_de(state, 0)
    assert ok and not np.isfinite(p_default[i])
    _, p_all, _, ok = cc._exact_de(dict(state, filter_min_count=0), 0)
    assert ok and np.isfinite(p_all[i])


def test_refilter_false_is_recorded_and_changes_the_signature(studies, tables, tmp_path):
    obs = _observed(tables)
    _run(studies, obs, n_perms=2, resume_dir=tmp_path / 'cc')
    out, msgs = _run(studies, obs, n_perms=2, resume_dir=tmp_path / 'cc', refilter=False)
    assert any('other settings' in m for m in msgs)
    sig = json.loads((tmp_path / 'cc' / 'manifest.json').read_text())
    assert sig['refilter'] is False
    assert out['reml'].attrs['cc_permutation']['refilter'] is False


def test_a_universe_changes_the_signature(studies, tables, tmp_path):
    obs = _observed(tables)
    _run(studies, obs, n_perms=2, resume_dir=tmp_path / 'cc')
    _out, msgs = _run(studies, obs, n_perms=2, resume_dir=tmp_path / 'cc',
                      gene_universe=[f'G{i}' for i in range(30)])
    assert any('other settings' in m for m in msgs)


def test_refilter_false_needs_a_refitting_engine(studies, tables):
    with pytest.raises(ValueError, match='refilter'):
        _run(studies, _observed(tables), engine='deseq2_fast', refilter=False)


def test_a_universe_with_no_pseudobulk_gene_is_refused(studies, tables):
    with pytest.raises(ValueError, match='gene_universe'):
        _run(studies, _observed(tables), gene_universe=['NOT_A_GENE'])
