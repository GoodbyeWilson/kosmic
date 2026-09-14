"""Correctness of DL (DerSimonian-Laird) and HKSJ pooling vs statsmodels.

Statsmodels is the reference; KOSMIC's vectorised dl_fast must agree at
machine precision. Plus the τ²-truncation behaviour that statsmodels
doesn't share (we truncate Q < df to FE; statsmodels returns NaN).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

from kosmic.meta_analysis.pooling.dl import dl_fast


def _datasets(effects, ses, gene='GENE1'):
    return [
        pd.DataFrame({
            'names': [gene], 'logfoldchanges': [effects[i]],
            'se': [ses[i]], 'pvals': [0.01],
            'dataset': [f'Study_{i + 1}'],
        })
        for i in range(len(effects))
    ]


def test_dl_matches_statsmodels_at_machine_precision():
    """Pooled effect, SE, I², p-value all match statsmodels under
    heterogeneous data (Q >> df, τ² well-defined)."""
    from statsmodels.stats.meta_analysis import combine_effects

    effects = np.array([0.2, 1.5, 0.3, 2.0, -0.4, 1.8, 0.1, 1.2])
    variances = np.array([0.05, 0.08, 0.03, 0.07, 0.04, 0.06, 0.05, 0.09])
    ses = np.sqrt(variances)

    sm = combine_effects(effects, variances, method_re='dl')
    out = dl_fast(_datasets(effects, ses), min_studies=2)
    row = out.iloc[0]

    np.testing.assert_allclose(row['logfoldchanges'], sm.mean_effect_re,
                                rtol=1e-12)
    np.testing.assert_allclose(row['se'], sm.sd_eff_w_re, rtol=1e-12)
    np.testing.assert_allclose(row['heterogeneity_i2'], sm.i2, rtol=1e-12)

    z_sm = sm.mean_effect_re / sm.sd_eff_w_re
    p_sm = 2 * scipy_stats.norm.sf(abs(z_sm))
    np.testing.assert_allclose(row['pvals_pooled'], p_sm, rtol=1e-12)


def test_dl_hksj_matches_statsmodels():
    """HKSJ-corrected SE matches statsmodels' sd_eff_w_re_hksj. p-value
    comes from t-distribution with k-1 df."""
    from statsmodels.stats.meta_analysis import combine_effects

    effects = np.array([0.2, 1.5, 0.3, 2.0, -0.4, 1.8, 0.1, 1.2])
    variances = np.array([0.05, 0.08, 0.03, 0.07, 0.04, 0.06, 0.05, 0.09])
    ses = np.sqrt(variances)
    k = len(effects)

    sm = combine_effects(effects, variances, method_re='dl')
    out = dl_fast(_datasets(effects, ses), min_studies=2, hksj=True)
    row = out.iloc[0]

    np.testing.assert_allclose(row['se'], sm.sd_eff_w_re_hksj, rtol=1e-12)
    t_stat = sm.mean_effect_re / sm.sd_eff_w_re_hksj
    p_sm = 2 * scipy_stats.t.sf(abs(t_stat), k - 1)
    np.testing.assert_allclose(row['pvals_pooled'], p_sm, rtol=1e-12)


def test_dl_truncates_negative_tau2_to_fixed_effects():
    """When Q < df the method-of-moments τ² is negative; KOSMIC truncates
    to 0 (Borenstein 2009, ch. 16) -- under truncation RE collapses to FE
    with I²=0. Statsmodels returns NaN here, so we cross-check internally."""
    effects = np.array([0.42, 0.55, 0.28, 0.61, 0.38])
    ses = np.array([0.2, 0.3, 0.17, 0.22, 0.28])
    datasets = _datasets(effects, ses)

    re_out = dl_fast(datasets, min_studies=2)
    fe_out = dl_fast(datasets, min_studies=2, method='fixed')

    np.testing.assert_allclose(
        re_out.iloc[0]['logfoldchanges'],
        fe_out.iloc[0]['logfoldchanges'],
        rtol=1e-12)
    assert re_out.iloc[0]['heterogeneity_i2'] == 0.0
