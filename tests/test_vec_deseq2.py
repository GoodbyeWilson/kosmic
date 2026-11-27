"""The vectorised per-study fit (ADR-019, stage 1).

It must return what the study's own DE ('fit_deseq2_wald', PyDESeq2) returns on the
same labels: the unshrunk log2 fold change, its standard error and the Wald p-value,
with the same genes left without a p-value by the Cook's filter. The lockstep
L-BFGS-B it relies on must equal scipy's own, and must report itself unavailable
when scipy's private interface is missing.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("pydeseq2")

from kosmic.de.de_analysis import fit_deseq2_wald  # noqa: E402
from kosmic.meta_analysis import vec_deseq2 as vd  # noqa: E402
from kosmic.meta_analysis import vec_lbfgsb  # noqa: E402


def simulate(n_dis, n_ctl, n_genes, seed, outliers=15):
    """NB pseudobulk counts with a depth effect, a sex effect, some DE genes and outliers."""
    rng = np.random.default_rng(seed)
    n = n_dis + n_ctl
    role = np.array(['disease'] * n_dis + ['control'] * n_ctl)
    sex = np.where(rng.random(n) < 0.5, 'F', 'M')
    depth = rng.uniform(0.6, 1.6, n)
    base = rng.lognormal(3.0, 1.6, n_genes)
    disp = 0.05 + 0.4 / np.sqrt(base)
    lfc = np.where(rng.random(n_genes) < 0.15, rng.normal(0, 1.0, n_genes), 0.0)
    sex_effect = rng.normal(0, 0.3, n_genes)
    mu = (base[None] * depth[:, None]
          * np.exp((role == 'disease')[:, None] * lfc[None]
                   + (sex == 'M')[:, None] * sex_effect[None]))
    r = 1 / disp
    y = rng.negative_binomial(r[None], r[None] / (r[None] + mu)).astype(float)
    for g in rng.choice(n_genes, outliers, replace=False):
        y[rng.integers(n), g] *= 25
    sample_df = pd.DataFrame({'sample': [f'D{i}' for i in range(n)], 'role': role, 'sex': sex})
    return y.round(), sample_df, [f'G{i}' for i in range(n_genes)]


# --- the lockstep optimiser ------------------------------------------------

def test_lockstep_lbfgsb_equals_scipy_minimize():
    from scipy.optimize import minimize
    rng = np.random.default_rng(5)
    a = rng.uniform(0.5, 3.0, size=(8, 2))
    c = rng.uniform(-2.0, 2.0, size=(8, 2))

    def evaluate(idx, x):
        d = x - c[idx]
        return (a[idx] * d ** 2).sum(axis=1) + 0.2 * (x ** 4).sum(axis=1), \
            2 * a[idx] * d + 0.8 * x ** 3

    x0 = rng.uniform(-1, 1, size=(8, 2))
    x, ok = vec_lbfgsb.lbfgsb_lockstep(evaluate, x0, -1.5, 1.5)
    assert ok.all()
    for k in range(8):
        ref = minimize(lambda z, k=k: tuple(v[0] for v in evaluate(np.array([k]), z[None])),
                       x0[k], jac=True, method='L-BFGS-B', bounds=[(-1.5, 1.5)] * 2)
        np.testing.assert_array_equal(x[k], ref.x)


def test_available_when_the_private_interface_works():
    assert vec_lbfgsb.available() is True


def test_unavailable_when_the_private_interface_is_missing(monkeypatch):
    monkeypatch.setattr(vec_lbfgsb, '_SETULB', None)
    monkeypatch.setattr(vec_lbfgsb, '_AVAILABLE', None)
    assert vec_lbfgsb.available() is False


def test_unavailable_when_the_private_interface_has_changed(monkeypatch):
    def changed(*args, **kwargs):
        raise TypeError("setulb() takes a different set of arguments")
    monkeypatch.setattr(vec_lbfgsb, '_SETULB', changed)
    monkeypatch.setattr(vec_lbfgsb, '_AVAILABLE', None)
    assert vec_lbfgsb.available() is False


# --- the fit ------------------------------------------------------------------

@pytest.mark.parametrize("n_dis, n_ctl, seed, covariates", [
    (8, 8, 1, ['sex']),
    (12, 10, 2, ['sex']),
    (6, 9, 3, []),
])
def test_fit_study_equals_fit_deseq2_wald(n_dis, n_ctl, seed, covariates):
    y, sample_df, genes = simulate(n_dis, n_ctl, 300, seed)
    fit = fit_deseq2_wald(y, sample_df, genes, covariates=covariates, n_cpus=1, quiet=True)
    res = fit.ds.results_df
    lfc, se, pval = vd.fit_study(y, sample_df, genes, covariates=covariates)

    pos = np.array([genes.index(g) for g in res.index])
    # the same genes were fitted, and the same genes lost their p-value to the Cook's filter
    assert np.isfinite(lfc).sum() == len(res)
    np.testing.assert_array_equal(np.isnan(pval[pos]), res['pvalue'].isna().to_numpy())
    np.testing.assert_allclose(lfc[pos], res['log2FoldChange'], rtol=1e-5, atol=1e-8)
    np.testing.assert_allclose(se[pos], res['lfcSE'], rtol=1e-5, atol=1e-8)
    ok = res['pvalue'].notna().to_numpy()
    np.testing.assert_allclose(pval[pos][ok], res['pvalue'].to_numpy()[ok], rtol=1e-4)


def test_fit_study_drops_donors_marked_exclude_like_fit_deseq2_wald():
    y, sample_df, genes = simulate(7, 7, 120, 4)
    sample_df = pd.concat([sample_df, pd.DataFrame(
        {'sample': ['X1', 'X2'], 'role': ['exclude', 'exclude'], 'sex': ['F', 'M']})],
        ignore_index=True)
    y = np.vstack([y, y[:2] * 3])
    fit = fit_deseq2_wald(y, sample_df, genes, covariates=['sex'], n_cpus=1, quiet=True)
    lfc, _se, _p = vd.fit_study(y, sample_df, genes, covariates=['sex'])
    pos = np.array([genes.index(g) for g in fit.ds.results_df.index])
    np.testing.assert_allclose(lfc[pos], fit.ds.results_df['log2FoldChange'],
                               rtol=1e-5, atol=1e-8)


def test_fit_study_gives_none_when_no_donor_pair_remains():
    y, sample_df, genes = simulate(1, 1, 30, 6, outliers=0)
    sample_df['role'] = ['exclude', 'disease']
    assert vd.fit_study(y, sample_df, genes) is None
