"""The one per-study model function (ADR-018, stage 2).

'fit_deseq2_wald' defines the per-study DESeq2 model. The DE workspace
('_run_deseq2') and the exact case-control permutation both call it, so the
tests pin down that the DE results are the function's results, and the
rules that decide the model: which donors, which covariates, which design.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("pydeseq2")

from kosmic.de.de_analysis import (  # noqa: E402
    Deseq2Fit, _run_deseq2, fit_deseq2_wald)


def _data(seed=0, n_genes=50, sex=None):
    """Four disease, four control and one excluded donor; NB counts."""
    rng = np.random.default_rng(seed)
    samples = ['D1', 'D2', 'D3', 'D4', 'C1', 'C2', 'C3', 'C4', 'X1']
    role = ['disease'] * 4 + ['control'] * 4 + ['exclude']
    mean = rng.gamma(2.0, 60.0, n_genes) + 20
    depth = rng.uniform(0.7, 1.4, len(samples))
    disease_up = np.where(np.arange(n_genes) < 10, 2.0, 1.0)
    mu = mean[None, :] * depth[:, None]
    mu[:4] *= disease_up[None, :]
    counts = rng.negative_binomial(8, 8 / (8 + mu)).astype(float)
    sample_df = pd.DataFrame({
        'sample': samples, 'condition': role, 'n_cells': 100, 'role': role,
        'sex_inferred': sex or ['male', 'female', 'male', 'female',
                                'female', 'male', 'female', 'male', 'male'],
    })
    return counts, sample_df, [f'G{i}' for i in range(n_genes)]


def test_the_fit_has_the_design_and_drops_the_excluded_donor():
    counts, sdf, genes = _data()
    fit = fit_deseq2_wald(counts, sdf, genes, covariates=['sex_inferred'],
                          quiet=True)
    assert isinstance(fit, Deseq2Fit)
    assert fit.design == '~sex_inferred + condition'
    assert fit.used_covariates == ['sex_inferred']
    assert len(fit.sample_df) == 8 and 'X1' not in set(fit.sample_df['sample'])
    assert 'log2FoldChange' in fit.ds.results_df.columns


def test_a_covariate_collinear_with_condition_is_left_out_of_the_design():
    sex = ['male'] * 4 + ['female'] * 4 + ['male']
    counts, sdf, genes = _data(sex=sex)
    fit = fit_deseq2_wald(counts, sdf, genes, covariates=['sex_inferred'],
                          quiet=True)
    assert fit.design == '~condition' and fit.used_covariates == []


def test_run_deseq2_reports_the_fits_unshrunk_estimates():
    counts, sdf, genes = _data()
    fit = fit_deseq2_wald(counts, sdf, genes, covariates=['sex_inferred'],
                          quiet=True)
    out = _run_deseq2(counts, sdf, genes, covariates=['sex_inferred']
                      ).set_index('names')
    ref = fit.ds.results_df
    common = ref.index.intersection(out.index)
    assert len(common) > 40
    np.testing.assert_allclose(out.loc[common, 'logfoldchanges'],
                               ref.loc[common, 'log2FoldChange'], rtol=1e-9)
    np.testing.assert_allclose(out.loc[common, 'se'],
                               ref.loc[common, 'lfcSE'], rtol=1e-9)
    np.testing.assert_allclose(out.loc[common, 'pvals'],
                               ref.loc[common, 'pvalue'], rtol=1e-9)
    assert out.attrs['design'] == '~sex_inferred + condition'


def test_one_core_and_quiet_do_not_change_the_fit():
    counts, sdf, genes = _data()
    a = fit_deseq2_wald(counts, sdf, genes, covariates=['sex_inferred'],
                        quiet=True).ds.results_df
    b = fit_deseq2_wald(counts, sdf, genes, covariates=['sex_inferred'],
                        n_cpus=1, quiet=True).ds.results_df
    np.testing.assert_allclose(a['log2FoldChange'], b['log2FoldChange'])
    np.testing.assert_allclose(a['pvalue'], b['pvalue'], equal_nan=True)


def test_too_few_donors_gives_none():
    counts, sdf, genes = _data()
    sdf = sdf.copy()
    sdf['role'] = ['disease'] + ['exclude'] * 8
    assert fit_deseq2_wald(counts, sdf, genes, quiet=True) is None
    assert _run_deseq2(counts, sdf, genes).empty


def test_missing_role_column_is_an_error():
    counts, sdf, genes = _data()
    with pytest.raises(ValueError, match="role"):
        fit_deseq2_wald(counts, sdf.drop(columns='role'), genes, quiet=True)
