"""Tests for kosmic.scrna.qc.decontx: supervised per-sample decontamination.

Sentinels: the three prerequisite gates (annotation, >=2 cell types,
sample column) and that a supervised run reduces cross-cell-type ambient
contamination while writing the decontaminated layer non-destructively.
"""
from __future__ import annotations

import importlib.util

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy.sparse import csr_matrix

from kosmic.scrna.qc.decontx import (
    DECONTX_CONTAMINATION,
    DECONTX_LAYER,
    check_decontx_prerequisites,
    run_decontx,
)

# DecontX is optional (not on PyPI); the app reports its absence and
# carries on, so the tests do the same rather than fail an install
# that is otherwise complete. The prerequisite-gate tests below do not
# need the package, so only the ones that run it are skipped.
requires_decontx = pytest.mark.skipif(
    importlib.util.find_spec("decontx") is None,
    reason="decontx not installed (optional; not on PyPI)")


def _contaminated_adata(n_genes=120, contam=0.12, seed=0):
    """Two cell types (CM genes 0..59, EC genes 60..119), two samples.

    EC cells carry ambient CM counts so DecontX has something to remove.
    """
    rng = np.random.default_rng(seed)
    half = n_genes // 2
    rows, cts, samples = [], [], []

    def _add(x, ct, s):
        rows.append(x)
        cts.append(ct)
        samples.append(s)

    for s in ('S1', 'S2'):
        for _ in range(120):  # CM cells
            x = np.zeros(n_genes)
            x[:half] += rng.poisson(8, half)
            x[half:] += rng.poisson(0.2, half)
            _add(x, 'CM', s)
        for _ in range(120):  # EC cells with CM ambient leak
            x = np.zeros(n_genes)
            x[half:] += rng.poisson(8, half)
            x[:half] += rng.poisson(8 * contam, half)
            _add(x, 'EC', s)
    X = csr_matrix(np.vstack(rows).astype(np.float32))
    obs = pd.DataFrame({'cell_type': cts, 'sample': samples})
    obs.index = [f'c{i}' for i in range(X.shape[0])]
    var = pd.DataFrame(index=[f'g{i}' for i in range(n_genes)])
    return AnnData(X=X, obs=obs, var=var)


def test_prerequisite_gates():
    a = _contaminated_adata()

    # All three satisfied.
    ok, reason = check_decontx_prerequisites(a, sample_col='sample')
    assert ok and reason == ''

    # Gate 1: no annotation.
    noann = a.copy()
    del noann.obs['cell_type']
    ok, reason = check_decontx_prerequisites(noann, sample_col='sample')
    assert not ok and 'nnotate' in reason

    # Gate 2: single cell type (already subset).
    sub = a[a.obs.cell_type == 'EC'].copy()
    ok, reason = check_decontx_prerequisites(sub, sample_col='sample')
    assert not ok and 'before Subset' in reason

    # Gate 3: no sample column.
    ok, reason = check_decontx_prerequisites(a, sample_col=None)
    assert not ok and 'per sample' in reason
    ok, reason = check_decontx_prerequisites(a, sample_col='nope')
    assert not ok


@requires_decontx
def test_run_reduces_contamination_nondestructively():
    a = _contaminated_adata()
    half = a.n_vars // 2

    ec = a[a.obs.cell_type == 'EC'].X.toarray()
    frac_before = ec[:, :half].sum() / ec.sum()

    raw_X = a.X.copy()
    _, info = run_decontx(a, sample_col='sample')

    # Decontaminated layer written; raw X untouched.
    assert DECONTX_LAYER in a.layers
    assert DECONTX_CONTAMINATION in a.obs.columns
    assert (a.X != raw_X).nnz == 0

    dec = a.layers[DECONTX_LAYER]
    dec = dec.toarray() if hasattr(dec, 'toarray') else np.asarray(dec)
    ec_dec = dec[(a.obs.cell_type == 'EC').values]
    frac_after = ec_dec[:, :half].sum() / ec_dec.sum()

    # Cross-type ambient markedly reduced.
    assert frac_after < frac_before * 0.5
    assert 0.0 <= info['mean_contamination'] <= 1.0


def test_run_refuses_when_gated():
    sub = _contaminated_adata()
    sub = sub[sub.obs.cell_type == 'EC'].copy()
    with pytest.raises(ValueError, match='before Subset'):
        run_decontx(sub, sample_col='sample')
