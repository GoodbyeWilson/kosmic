"""Stouffer, Fisher and gwOP are the textbook methods: each study's p-value is
used as it is (no floor at 1/n_genes), a study with p = 1 counts (z = 0), and a
p of 0 (an underflow) counts as the most significant value rather than being
dropped. The permutation engine pools the same way as the analytical method,
so a calibration's null matches the observed statistic. SumRank's z takes the
sign of the direction its rank test found."""
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from kosmic.meta_analysis.cc_permutation import _pool_fisher_pvals, _pool_stouffer_pvals
from kosmic.meta_analysis.pooling.fisher import fisher_fast
from kosmic.meta_analysis.pooling.stouffer import stouffer_fast
from kosmic.meta_analysis.pooling.sumrank import sumrank_fast

# Gene G0: very strong in every study (p far below 1/n_genes); G1: p = 1 in one
# study; G2: p = 0 (underflow) in one study; the rest moderate.
P = np.array([
    [1e-30, 1e-25, 1e-40],
    [0.01, 1.0, 0.02],
    [0.0, 0.03, 0.04],
    [0.2, 0.5, 0.7],
    [0.04, 0.3, 0.9],
])
LFC = np.array([
    [1.0, 1.2, 0.9],
    [0.5, 0.0, 0.4],
    [0.8, 0.6, 0.5],
    [-0.2, 0.1, -0.3],
    [0.3, -0.1, 0.05],
])
NAMES = [f'G{i}' for i in range(len(P))]


def _studies():
    return [pd.DataFrame({'names': NAMES, 'logfoldchanges': LFC[:, j], 'se': 0.2,
                          'pvals': P[:, j], 'dataset': f'S{j}'}) for j in range(P.shape[1])]


def _textbook_stouffer(p, lfc):
    z = np.sign(lfc) * stats.norm.isf(np.clip(p, 1e-300, 1.0) / 2)
    return z.sum() / np.sqrt(len(p))


def _row(df, g):
    return df.set_index('names').loc[g]


def test_stouffer_uses_each_p_as_it_is():
    out = stouffer_fast(_studies(), min_studies=2)
    for i, g in enumerate(NAMES):
        np.testing.assert_allclose(_row(out, g)['z_stat'], _textbook_stouffer(P[i], LFC[i]),
                                   rtol=1e-12)
    # No cap: G0's z is far above the 1/n_genes ceiling of the old floor
    # (isf(1/(2*5)) * sqrt(3) ~ 2.8).
    assert _row(out, 'G0')['z_stat'] > 15


def test_a_study_with_p_1_counts():
    out = stouffer_fast(_studies(), min_studies=2)
    assert _row(out, 'G1')['n_studies'] == 3
    # z = (z1 + 0 + z3) / sqrt(3), not / sqrt(2)
    np.testing.assert_allclose(_row(out, 'G1')['z_stat'], _textbook_stouffer(P[1], LFC[1]))


def test_fisher_uses_each_p_as_it_is():
    out = fisher_fast(_studies(), min_studies=2)
    for i, g in enumerate(NAMES):
        expected = stats.combine_pvalues(np.clip(P[i], 1e-300, 1.0), method='fisher').pvalue
        np.testing.assert_allclose(_row(out, g)['pvals_pooled'], max(expected, 1e-300), rtol=1e-9)


@pytest.mark.parametrize('pooler,analytical', [(_pool_stouffer_pvals, stouffer_fast),
                                               (_pool_fisher_pvals, fisher_fast)])
def test_the_permutation_engine_pools_as_the_analytical_method(pooler, analytical):
    valid = np.ones_like(P, dtype=bool)
    engine = pooler(LFC, np.full_like(P, 0.04), valid, pval_matrix=P)
    out = analytical(_studies(), min_studies=2).set_index('names').loc[NAMES]
    np.testing.assert_allclose(engine, out['pvals_pooled'].to_numpy(float), rtol=1e-9)


def test_sumrank_z_takes_the_direction_of_its_test():
    # G0 ranks top in every study (up) but its mean log2FC is made negative by
    # one study's large negative fold change.
    lfc = LFC.copy()
    lfc[0] = [1.0, 1.2, -9.0]
    studies = [pd.DataFrame({'names': NAMES, 'logfoldchanges': lfc[:, j], 'se': 0.2,
                             'pvals': P[:, j], 'dataset': f'S{j}'}) for j in range(3)]
    out = sumrank_fast(studies, min_studies=2).set_index('names')
    assert out.at['G0', 'logfoldchanges'] < 0
    assert (out['z_stat'] > 0).equals(out['direction'] == 'up')
