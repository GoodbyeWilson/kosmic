"""Textbook checks for REML's I², Stouffer's direction in the comparison,
BH with a p of 0, and the GSEA gene-set size defaults."""
import numpy as np
import pandas as pd
import pytest

# dat.bcg (Colditz et al. 1994), as shipped with metafor: tpos, tneg, cpos, cneg.
BCG = [(4, 119, 11, 128), (6, 300, 29, 274), (3, 228, 11, 209),
       (62, 13536, 248, 12619), (33, 5036, 47, 5761), (180, 1361, 372, 1079),
       (8, 2537, 10, 619), (505, 87886, 499, 87892), (29, 7470, 45, 7232),
       (17, 1699, 65, 1600), (186, 50448, 141, 27197), (5, 2493, 3, 2338),
       (27, 16886, 29, 17825)]


def _bcg_datasets():
    out = []
    for i, (a, b, c, d) in enumerate(BCG):
        yi = np.log((a / (a + b)) / (c / (c + d)))
        vi = 1 / a - 1 / (a + b) + 1 / c - 1 / (c + d)
        out.append(pd.DataFrame({'names': ['g'], 'logfoldchanges': [yi],
                                 'se': [np.sqrt(vi)], 'pvals': [0.5]}))
    return out


def test_reml_i2_matches_metafor_dat_bcg():
    """metafor::rma(yi, vi, data=dat.bcg, method='REML'): tau² 0.3132,
    I² 92.22%, estimate -0.7145."""
    from kosmic.meta_analysis.pooling.reml import reml_fast
    r = reml_fast(_bcg_datasets(), min_studies=2).iloc[0]
    assert r['tau_squared'] == pytest.approx(0.3132, abs=5e-5)
    assert r['logfoldchanges'] == pytest.approx(-0.7145, abs=5e-5)
    assert r['heterogeneity_i2'] == pytest.approx(0.9222, abs=5e-5)


def test_dl_i2_matches_metafor_dat_bcg():
    """method='DL': tau² 0.3088, I² 92.12%, which for DL is (Q - df) / Q."""
    from kosmic.meta_analysis.pooling.dl import dl_fast
    r = dl_fast(_bcg_datasets(), min_studies=2).iloc[0]
    assert r['tau_squared'] == pytest.approx(0.3088, abs=5e-5)
    assert r['heterogeneity_i2'] == pytest.approx(0.9212, abs=5e-5)


def test_reml_i2_zero_when_no_heterogeneity():
    from kosmic.meta_analysis.pooling.reml import reml_fast
    ds = [pd.DataFrame({'names': ['g'], 'logfoldchanges': [0.5], 'se': [0.2],
                        'pvals': [0.01]}) for _ in range(4)]
    assert reml_fast(ds, min_studies=2).iloc[0]['heterogeneity_i2'] == 0.0


def test_stouffer_direction_is_the_sign_of_z():
    """Effect and Z can disagree in sign; the comparison takes Stouffer's own."""
    from kosmic.comparison.standard import stouffer_to_standard
    df = pd.DataFrame({'names': ['a', 'b', 'c'], 'logfoldchanges': [0.3, -0.2, 0.1],
                       'z_stat': [-2.5, 3.0, 0.0], 'pvals_pooled': [0.01, 0.003, 1.0]})
    out = stouffer_to_standard(df).set_index('gene')
    assert out.loc['a', 'direction'] == -1
    assert out.loc['b', 'direction'] == 1
    assert out.loc['c', 'direction'] == 1  # Z = 0: falls back to the effect


def test_bh_counts_a_p_of_zero():
    """p = 0 is the smallest p, not a missing one: it stays in the BH
    denominator and gets q = 0, as statsmodels gives."""
    from statsmodels.stats.multitest import multipletests
    from kosmic.numerical import bh_fdr
    p = np.array([0.0, 0.01, 0.02, 0.5, np.nan])
    q = bh_fdr(p)
    expect = multipletests(p[:4], method='fdr_bh')[1]
    np.testing.assert_allclose(q[:4], expect)
    assert q[0] == 0.0 and q[4] == 1.0


def test_gsea_size_defaults_are_gseas():
    import inspect
    from kosmic.de import fgsea
    assert (fgsea.GSEA_MIN_SIZE, fgsea.GSEA_MAX_SIZE) == (15, 500)
    sig = inspect.signature(fgsea.run_fgsea)
    assert sig.parameters['min_size'].default == 15
    assert sig.parameters['max_size'].default == 500


def test_sumrank_top_r_monte_carlo_p_is_never_zero():
    """A rank sum below every simulated one gets (0 + 1) / (n + 1), not 0."""
    from kosmic.meta_analysis.pooling.sumrank import _top_r_sum_cdf
    p = _top_r_sum_cdf(1e-9, 5, 3)
    assert 0 < p < 1e-6
