"""Tests for per-cell pathway scoring: log-normalisation of the raw-counts
copy, and faithful delegation of `score_genes` to scanpy.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

from kosmic.de.pathway_scoring import score_pathways, _ensure_lognorm


def _counts_adata(seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"GENE{i}" for i in range(60)]
    pw = genes[:8]
    depth = rng.uniform(0.5, 2.0, (240, 1))
    counts = np.rint(rng.poisson(3.0, size=(240, len(genes))) * depth).astype(float)
    a = ad.AnnData(counts.copy(), var=pd.DataFrame(index=genes))
    a.obs["sample"] = [f"s{i % 8}" for i in range(240)]
    a.obs["condition"] = ["A" if i % 8 < 4 else "B" for i in range(240)]
    a.raw = a.copy()                       # raw = counts (standard pipeline)
    return a, pw


def test_ensure_lognorm_normalises_counts():
    a, _ = _counts_adata()
    work = a.raw.to_adata()
    assert np.allclose(work.X, np.rint(work.X))     # counts before
    _ensure_lognorm(work)
    assert not np.allclose(work.X, np.rint(work.X))  # log-norm floats after
    assert "log1p" in work.uns


def test_ensure_lognorm_skips_already_normalised():
    a, _ = _counts_adata()
    work = a.raw.to_adata()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    before = work.X.copy()
    _ensure_lognorm(work)                            # must be a no-op
    assert np.allclose(work.X, before)


def test_score_genes_matches_scanpy_reference():
    """score_pathways('score_genes') == direct sc.tl.score_genes on the
    log-normalised full-gene matrix (same ctrl_size, same seed)."""
    a, pw = _counts_adata()

    work_k, avail = score_pathways(a, {"PW": pw}, method="score_genes")
    assert "PW" in avail

    ref = a.raw.to_adata()
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)
    sc.tl.score_genes(ref, gene_list=pw, ctrl_size=50,
                      score_name="PW_score", random_state=0, use_raw=False)

    assert np.allclose(work_k.obs["PW_score"].to_numpy(),
                       ref.obs["PW_score"].to_numpy())
