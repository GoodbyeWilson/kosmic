"""Tests for `per_donor_pathway_score` — the single per-donor score used by
the Pathway Score plot and the spillover diagnostic.

Pins that it (a) returns one value per donor for each scoring method by
reusing that method's own computation, (b) the `deseq2` method produces a
DESeq2 size-factor-log2 per-donor score distinct from CPM, and (c) reports
coverage failure as empty.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
import pytest

from kosmic.de.pathway_scoring import per_donor_pathway_score


def _make_adata(seed=0):
    rng = np.random.default_rng(seed)
    n_samples, cells = 8, 60
    pw = [f"P{i}" for i in range(6)]
    genes = pw + [f"G{i}" for i in range(20)]
    blocks, samples, conds = [], [], []
    for s in range(n_samples):
        depth = 1.0 if s % 2 == 0 else 3.0          # depth differences
        X = rng.poisson(4.0 * depth, size=(cells, len(genes))).astype(float)
        blocks.append(X)
        samples += [f"s{s}"] * cells
        conds += (["A"] if s < 4 else ["B"]) * cells
    counts = np.vstack(blocks)
    adata = ad.AnnData(X=counts.copy(), var=pd.DataFrame(index=genes))
    adata.obs["sample"] = samples
    adata.obs["condition"] = conds
    adata.layers["counts"] = counts.copy()
    adata.raw = adata.copy()
    lib = np.clip(counts.sum(1, keepdims=True), 1, None)
    adata.X = np.log1p(counts / lib * 1e4)          # processed state
    return adata, pw


@pytest.mark.parametrize("method",
                         ["pseudobulk_cpm", "deseq2",
                          "mean", "zscore", "singscore"])
def test_returns_one_score_per_donor(method):
    adata, pw = _make_adata()
    scores, sdf = per_donor_pathway_score(
        adata, "PW", pw, "sample", "condition", method)
    assert scores.shape == (8,)
    assert list(sdf.columns[:2]) == ["sample", "condition"]
    assert len(sdf) == 8
    assert np.all(np.isfinite(scores))


def test_deseq2_differs_from_cpm():
    """deseq2 produces a DESeq2 size-factor-log2 per-donor score, not CPM."""
    adata, pw = _make_adata()
    cpm, _ = per_donor_pathway_score(adata, "PW", pw, "sample", "condition",
                                     "pseudobulk_cpm")
    dsq, _ = per_donor_pathway_score(adata, "PW", pw, "sample", "condition",
                                     "deseq2")
    assert not np.allclose(cpm, dsq)


def test_coverage_failure_returns_empty():
    adata, _ = _make_adata()
    scores, sdf = per_donor_pathway_score(
        adata, "PW", ["NOT_A_GENE"], "sample", "condition", "pseudobulk_cpm")
    assert scores.size == 0 and sdf.empty
