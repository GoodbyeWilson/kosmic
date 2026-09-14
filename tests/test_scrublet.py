"""Doublet detection finds planted cross-type doublets and nothing else.

The wrapper moved from the standalone 'scrublet' package to scanpy's
port of the same algorithm, because the standalone package hard-depends
on 'annoy', which has no Windows wheel and so demanded a C++ compiler
at install time. These tests pin that the swap kept the behaviour:
doublets made by summing cells of *different* types are the case the
method exists to catch, and it should catch nearly all of them while
calling almost no singlets.
"""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")
pytest.importorskip("scanpy")

from kosmic.scrna.qc.scrublet import run_scrublet  # noqa: E402


def _planted(seed=11, n_per=400, n_genes=600, n_types=3, n_dbl=90):
    """Three populations with distinct marker programmes, plus doublets
    formed only across populations."""
    rng = np.random.default_rng(seed)
    cells, labels = [], []
    for t in range(n_types):
        mu = np.full(n_genes, 0.3)
        mu[t * 120:(t + 1) * 120] = 4.0
        cells.append(rng.poisson(mu, size=(n_per, n_genes)))
        labels += [t] * n_per
    base = np.vstack(cells).astype(np.float32)
    labels = np.array(labels)
    i = rng.integers(0, len(base), n_dbl)
    j = rng.integers(0, len(base), n_dbl)
    keep = labels[i] != labels[j]
    i, j = i[keep], j[keep]
    X = np.vstack([base, base[i] + base[j]])
    truth = np.r_[np.zeros(len(base), bool), np.ones(len(i), bool)]
    return X, truth


def test_finds_planted_doublets_and_spares_singlets():
    X, truth = _planted()
    a = anndata.AnnData(X)
    a, res = run_scrublet(a, expected_doublet_rate=0.06)
    called = a.obs["predicted_doublet"].to_numpy().astype(bool)
    assert called[truth].mean() > 0.9          # recall
    assert called[~truth].mean() < 0.02        # false positives
    assert res["n_doublets"] == int(called.sum())
    assert res["n_total"] == len(truth)
    assert not res["used_manual_threshold"]


def test_scores_rank_doublets_above_singlets():
    X, truth = _planted(seed=3)
    a = anndata.AnnData(X)
    a, _ = run_scrublet(a)
    s = a.obs["doublet_score"].to_numpy()
    # Every planted doublet outscores the median singlet.
    assert (s[truth] > np.median(s[~truth])).all()


def test_uses_raw_counts_when_present():
    """Scrublet wants counts; normalised .X must not be what it sees."""
    X, truth = _planted(seed=5)
    a = anndata.AnnData(X)
    a.raw = a
    a.X = np.log1p(a.X / a.X.sum(axis=1, keepdims=True) * 1e4)  # normalised
    a, _ = run_scrublet(a)
    called = a.obs["predicted_doublet"].to_numpy().astype(bool)
    assert called[truth].mean() > 0.9


def test_does_not_mutate_caller_matrix():
    X, _ = _planted(seed=9)
    a = anndata.AnnData(X)
    before = np.asarray(a.X).copy()
    run_scrublet(a)
    assert np.array_equal(np.asarray(a.X), before)


def test_low_count_cells_are_not_scored():
    X, _ = _planted(seed=2)
    X[:5] = 0                                   # five empty cells
    a = anndata.AnnData(X)
    a, _ = run_scrublet(a, min_counts=2)
    assert (a.obs["doublet_score"].to_numpy()[:5] == 0).all()
    assert not a.obs["predicted_doublet"].to_numpy()[:5].any()


def test_runs_without_annoy(monkeypatch):
    """The whole point of the swap: no annoy, no C++ compiler."""
    import sys
    monkeypatch.setitem(sys.modules, "annoy", None)   # `import annoy` raises
    X, truth = _planted(seed=1)
    a = anndata.AnnData(X)
    a, _ = run_scrublet(a)
    assert a.obs["predicted_doublet"].to_numpy()[truth].mean() > 0.9
