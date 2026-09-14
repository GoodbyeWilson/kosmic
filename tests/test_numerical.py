"""Tests for kosmic.numerical: bh_fdr and neg_log10."""
from __future__ import annotations

import numpy as np
import pytest
from statsmodels.stats.multitest import multipletests

from kosmic.numerical import bh_fdr, neg_log10


def test_bh_fdr_matches_statsmodels():
    """Proof that we do Benjamini-Hochberg correctly."""
    rng = np.random.default_rng(42)
    pvals = rng.uniform(0.0001, 0.99, size=50)
    expected = multipletests(pvals, method="fdr_bh")[1]
    np.testing.assert_allclose(bh_fdr(pvals), expected)


def test_bh_fdr_excludes_invalid_from_denominator():
    """NaN/inf must not inflate the BH denominator -- that would silently
    deflate every q-value in the array."""
    valid_only = bh_fdr([0.01, 0.05, 0.1])
    mixed = bh_fdr([0.01, np.nan, 0.05, 0.1])
    assert mixed[0] == pytest.approx(valid_only[0])
    assert mixed[2] == pytest.approx(valid_only[1])
    assert mixed[3] == pytest.approx(valid_only[2])


def test_neg_log10_of_zero_is_finite():
    """The whole point of this function: matplotlib / pyqtgraph can't plot inf,
    so neg_log10 floors p at 1e-300 by default."""
    out = neg_log10([0.0, -0.5, 1.5, 0.01])
    assert np.all(np.isfinite(out))
    assert out[3] == pytest.approx(2.0)
