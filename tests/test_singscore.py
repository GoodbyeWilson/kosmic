"""Tests for the singscore rank-based pathway score (Foroutan et al. 2018).

Pins the [-0.5, +0.5] bounds and the extreme behaviour: a set of the
highest-expressed genes scores near +0.5, the lowest-expressed near -0.5.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad

from kosmic.de.pathway_scoring import score_pathways


def _adata_with_high_and_low_sets(seed=0):
    rng = np.random.default_rng(seed)
    N, cells = 40, 120
    genes = [f"GENE{i}" for i in range(N)]
    high, low = genes[:5], genes[-5:]
    X = np.rint(rng.uniform(1, 3, (cells, N)))          # mid genes: ties
    X[:, :5] += 50                                       # high set: top ranks
    X[:, -5:] = 0                                         # low set: bottom ranks
    a = ad.AnnData(X, var=pd.DataFrame(index=genes))
    a.obs["sample"] = ["s0"] * cells
    a.obs["condition"] = ["A"] * cells
    a.raw = a.copy()                                     # raw = counts
    return a, high, low


def test_singscore_bounds_and_extremes():
    a, high, low = _adata_with_high_and_low_sets()
    work, avail = score_pathways(a, {"HI": high, "LO": low}, method="singscore")
    assert "HI" in avail and "LO" in avail

    hi = work.obs["HI_score"].to_numpy()
    lo = work.obs["LO_score"].to_numpy()

    # Bounded to [-0.5, +0.5].
    assert np.all(hi <= 0.5 + 1e-9) and np.all(hi >= -0.5 - 1e-9)
    assert np.all(lo <= 0.5 + 1e-9) and np.all(lo >= -0.5 - 1e-9)

    # Highest-expressed set near +0.5, lowest near -0.5.
    assert hi.mean() > 0.3
    assert lo.mean() < -0.3


def test_singscore_rank_based_invariant_to_normalisation():
    """Ranks are monotonic under log1p, so the score is identical whether the
    matrix is raw counts or log-normalised."""
    a, high, _ = _adata_with_high_and_low_sets()
    work, _ = score_pathways(a, {"HI": high}, method="singscore")
    from_counts = work.obs["HI_score"].to_numpy()

    # Pre-normalise so _ensure_lognorm sees floats and skips its own step.
    import scanpy as sc
    b = a.copy()
    b.raw = None
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    b.raw = b.copy()
    work2, _ = score_pathways(b, {"HI": high}, method="singscore")
    assert np.allclose(from_counts, work2.obs["HI_score"].to_numpy())
