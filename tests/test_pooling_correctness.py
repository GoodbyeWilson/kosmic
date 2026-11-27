"""Correctness sentinels for the pooling methods exposed in the GUI.

KOSMIC offers six pooling methods across three families:

- Effect-size: DL, REML
- Rank: SumRank, gwOP
- P-value: Fisher, Stouffer

DL has its own file (``test_dl_correctness.py``) covering the τ²
truncation behaviour and HKSJ variant. This file adds one canonical
correctness sentinel per remaining method.

The rank and p-value methods are tested in a realistic gene-set context
(SumRank needs a real
rank space). All such tests use a 30-gene fixture.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats as scipy_stats

from kosmic.meta_analysis.pooling.fisher import fisher_fast
from kosmic.meta_analysis.pooling.gwop import gwop_fast
from kosmic.meta_analysis.pooling.reml import reml_fast
from kosmic.meta_analysis.pooling.stouffer import stouffer_fast
from kosmic.meta_analysis.pooling.sumrank import sumrank_fast


# ----------------------------------------------------------------------
# REML: single-gene cross-check against statsmodels
# ----------------------------------------------------------------------
def test_reml_matches_statsmodels():
    """REML pooled effect must agree with statsmodels' iterated REML
    estimator on a heterogeneous dataset where τ² is well-defined.
    Tolerance is rtol=2e-3 because the two implementations use slightly
    different convergence criteria."""
    from statsmodels.stats.meta_analysis import combine_effects

    effects = np.array([0.2, 1.5, 0.3, 2.0, -0.4, 1.8, 0.1, 1.2])
    variances = np.array([0.05, 0.08, 0.03, 0.07, 0.04, 0.06, 0.05, 0.09])
    ses = np.sqrt(variances)
    datasets = [
        pd.DataFrame({
            'names': ['G1'], 'logfoldchanges': [effects[i]],
            'se': [ses[i]], 'pvals': [0.01],
            'dataset': [f'Study_{i+1}'],
        })
        for i in range(len(effects))
    ]

    sm = combine_effects(effects, variances, method_re='iterated')
    out = reml_fast(datasets, min_studies=2)
    row = out.iloc[0]

    np.testing.assert_allclose(row['logfoldchanges'], sm.mean_effect_re,
                                rtol=2e-3)
    np.testing.assert_allclose(row['se'], sm.sd_eff_w_re, rtol=2e-3)


# ----------------------------------------------------------------------
# Prediction interval: t(k-2), matching metafor predict(pi.type="Riley")
# ----------------------------------------------------------------------
def _single_gene_datasets(effects, variances):
    return [
        pd.DataFrame({'names': ['G1'], 'logfoldchanges': [e], 'se': [np.sqrt(v)],
                      'pvals': [0.01], 'dataset': [f'Study_{i + 1}']})
        for i, (e, v) in enumerate(zip(effects, variances))
    ]


# Reference values from metafor 4.8-0 (rma(method=...), predict(pi.type="Riley")).
_PI_EFFECTS = [-0.82, -1.14, -0.59, -0.75, -0.89, 0.10]
_PI_VARIANCES = [0.043, 0.058, 0.028, 0.018, 0.058, 0.05]


@pytest.mark.parametrize("k, method, tau2, lo, hi", [
    (5, "reml", 0.0, -1.0446099202, -0.5231923222),
    (6, "reml", 0.1170077380, -1.7293649236, 0.3730126766),
    (3, "dl", 0.0321662791, -3.8348272653, 2.2081755650),
])
def test_prediction_interval_matches_metafor(k, method, tau2, lo, hi):
    """The prediction interval uses t with k-2 degrees of freedom, as in
    metafor's pi.type="Riley". A normal (1.96) interval is narrower, most
    of all for few studies."""
    from kosmic.meta_analysis.pooling.dl import dl_fast

    fn = reml_fast if method == "reml" else dl_fast
    row = fn(_single_gene_datasets(_PI_EFFECTS[:k], _PI_VARIANCES[:k]),
             min_studies=2).iloc[0]
    np.testing.assert_allclose(row['tau_squared'], tau2, atol=2e-5)
    np.testing.assert_allclose(row['pi_lower'], lo, atol=2e-4)
    np.testing.assert_allclose(row['pi_upper'], hi, atol=2e-4)
    normal_half = 1.96 * np.sqrt(row['se'] ** 2 + row['tau_squared'])
    assert (row['pi_upper'] - row['pi_lower']) / 2 > normal_half


def test_prediction_interval_degrees_of_freedom_floor():
    """Two studies give df = 0; the degrees of freedom are floored at 1
    (as in metafor), so the interval is finite and very wide."""
    row = reml_fast(_single_gene_datasets([0.2, 0.5], [0.04, 0.05]),
                    min_studies=2).iloc[0]
    half = (row['pi_upper'] - row['pi_lower']) / 2
    np.testing.assert_allclose(
        half, scipy_stats.t.ppf(0.975, 1) * np.sqrt(row['se'] ** 2 + row['tau_squared']))


# ----------------------------------------------------------------------
# Multi-gene fixture: 30 genes per study, 3 studies, with a planted SIG
# gene (p≪0.01) and a planted NULL gene (p≈0.7). Other 28 genes are
# random.
# ----------------------------------------------------------------------
@pytest.fixture
def three_studies_30_genes():
    rng = np.random.default_rng(42)
    n_genes = 30
    studies = []
    for k in range(3):
        names = [f'G{i:02d}' for i in range(n_genes)]
        # Planted significant gene at index 0
        # Planted null gene at index 1
        # Random for the rest
        pvals = np.concatenate([
            [0.0001 + 0.0001 * k],   # SIG (index 0): consistently significant
            [0.6 + 0.05 * k],        # NULL (index 1): consistently non-significant
            rng.uniform(0.01, 0.99, n_genes - 2),
        ])
        effects = np.concatenate([
            [1.5 + 0.1 * k],         # SIG: large positive
            [0.05],                  # NULL: tiny
            rng.normal(0, 0.3, n_genes - 2),
        ])
        ses = rng.uniform(0.15, 0.35, n_genes)
        studies.append(pd.DataFrame({
            'names': names,
            'logfoldchanges': effects,
            'se': ses,
            'pvals': pvals,
            'dataset': [f'S{k+1}'] * n_genes,
        }))
    return studies


def _row(df, gene_name):
    return df[df['names'] == gene_name].iloc[0]


# ----------------------------------------------------------------------
# Fisher: matches scipy on the planted SIG gene
# ----------------------------------------------------------------------
def test_fisher_pooled_pvalue_orders_sig_before_null(three_studies_30_genes):
    """Fisher: planted SIG gene must have a smaller pooled p-value than
    the planted NULL gene. Sanity check that χ² combined-p is doing
    something."""
    out = fisher_fast(three_studies_30_genes, min_studies=2)
    sig = _row(out, 'G00')
    null = _row(out, 'G01')
    assert sig['pvals_pooled'] < null['pvals_pooled']


def test_fisher_matches_scipy_on_the_raw_pvalues(three_studies_30_genes):
    """KOSMIC's Fisher is the textbook test: each study's p used as it is,
    so the combined p equals scipy.stats.combine_pvalues."""
    sig_raw = np.array([
        ds[ds['names'] == 'G00'].iloc[0]['pvals']
        for ds in three_studies_30_genes
    ])
    expected = scipy_stats.combine_pvalues(sig_raw, method='fisher').pvalue

    out = fisher_fast(three_studies_30_genes, min_studies=2)
    sig = _row(out, 'G00')
    np.testing.assert_allclose(sig['pvals_pooled'], expected, rtol=1e-10)


# ----------------------------------------------------------------------
# Stouffer: planted SIG gene must beat planted NULL gene
# ----------------------------------------------------------------------
def test_stouffer_pooled_pvalue_orders_sig_before_null(three_studies_30_genes):
    """Stouffer's z-method: SIG must outrank NULL."""
    out = stouffer_fast(three_studies_30_genes, min_studies=2)
    sig = _row(out, 'G00')
    null = _row(out, 'G01')
    assert sig['pvals_pooled'] < null['pvals_pooled']


# ----------------------------------------------------------------------
# Rank methods: SIG outranks NULL
# ----------------------------------------------------------------------
def test_sumrank_orders_sig_before_null(three_studies_30_genes):
    out = sumrank_fast(three_studies_30_genes, min_studies=2)
    sig = _row(out, 'G00')
    null = _row(out, 'G01')
    assert sig['pvals_pooled'] < null['pvals_pooled']


def test_gwop_orders_sig_before_null(three_studies_30_genes):
    out = gwop_fast(three_studies_30_genes, min_studies=2)
    sig = _row(out, 'G00')
    null = _row(out, 'G01')
    assert sig['pvals_pooled'] < null['pvals_pooled']
