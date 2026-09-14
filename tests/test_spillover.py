"""Tests for the ambient-spillover diagnostic.

Pins that `spillover_attribution`:
  (a) reports ~0% attributable when contamination is unrelated to the score;
  (b) reports a high % attributable when the disease effect is *entirely*
      driven by contamination;
  (c) reports a partial fraction and keeps a real effect significant when
      the effect is part-real, part-contamination (the KOSMIC EC OXPHOS case);
  (d) validates inputs.
"""
from __future__ import annotations

import numpy as np
import pytest

import anndata as ad
import pandas as pd

from kosmic.de.spillover import spillover_attribution, pathway_spillover


def _labels(n_disease, n_control):
    return np.array([1] * n_disease + [0] * n_control)


def test_pure_contamination_high_attribution():
    """Score is a deterministic function of contamination only; disease acts
    solely through contamination -> most of the effect is attributable."""
    rng = np.random.default_rng(0)
    nd, nc = 15, 20
    d = _labels(nd, nc)
    # Disease has higher contamination; score tracks contamination exactly.
    contam = np.where(d == 1, 0.30, 0.15) + rng.normal(0, 0.01, nd + nc)
    score = 5.0 * contam + rng.normal(0, 1e-3, nd + nc)

    r = spillover_attribution(score, contam, d)
    assert r["pct_attributable"] > 80
    assert abs(r["adjusted_effect"]) < abs(r["raw_effect"])


def test_no_contamination_coupling_low_attribution():
    """Contamination is random noise unrelated to disease or score -> the
    disease effect barely moves under adjustment."""
    rng = np.random.default_rng(1)
    nd, nc = 15, 20
    d = _labels(nd, nc)
    contam = rng.uniform(0.1, 0.3, nd + nc)          # unrelated to disease
    score = np.where(d == 1, -0.6, 0.0) + rng.normal(0, 0.1, nd + nc)

    r = spillover_attribution(score, contam, d)
    assert r["raw_effect"] < 0                        # disease lowers score
    assert r["pct_attributable"] < 20


def test_partial_real_effect_survives():
    """Part contamination, part real (the EC OXPHOS case): the effect
    attenuates but a real component remains and stays significant."""
    rng = np.random.default_rng(2)
    nd, nc = 13, 25
    d = _labels(nd, nc)
    contam = np.where(d == 1, 0.14, 0.18) + rng.normal(0, 0.02, nd + nc)
    # Score = genuine disease drop  +  contamination coupling  + noise.
    score = (-0.5 * d) + (2.0 * contam) + rng.normal(0, 0.15, nd + nc)

    r = spillover_attribution(score, contam, d)
    assert r["raw_effect"] < 0
    assert 0 < r["pct_attributable"] < 100
    # A genuine component survives adjustment and stays significant.
    assert r["adjusted_effect"] < 0
    assert r["adjusted_p"] < 0.05


def test_input_validation():
    with pytest.raises(ValueError):
        spillover_attribution([1, 2, 3], [1, 2, 3], [0, 1, 0])   # n < 4
    with pytest.raises(ValueError):
        spillover_attribution([1, 2, 3, 4], [1, 2, 3], [0, 1, 0, 1])  # mismatch


def _make_pathway_adata(seed=0):
    """Synthetic adata with a pathway, a contamination metric, and a
    disease effect that is partly contamination-driven."""
    rng = np.random.default_rng(seed)
    n_samples, cells = 20, 60
    pw_genes = [f"OX{i}" for i in range(8)]
    other = [f"G{i}" for i in range(12)]
    genes = pw_genes + other

    blocks, sample_ids, conditions = [], [], []
    for s in range(n_samples):
        disease = s < 10
        cell_contam = rng.uniform(0.10, 0.20) if disease else rng.uniform(0.15, 0.25)
        X = rng.poisson(5.0, size=(cells, len(genes))).astype(float)
        # Pathway lower in disease (real) + tracks contamination (soup).
        X[:, :len(pw_genes)] += (0 if disease else 4) + cell_contam * 10
        blocks.append(X)
        sample_ids += [f"s{s}"] * cells
        conditions += (["DCM"] if disease else ["Donor"]) * cells

    counts = np.vstack(blocks)
    adata = ad.AnnData(X=counts.copy(), var=pd.DataFrame(index=genes))
    adata.obs["sample"] = sample_ids
    adata.obs["condition"] = conditions
    adata.layers["counts"] = counts.copy()
    adata.raw = adata.copy()
    # A contamination metric that differs by condition, with per-cell
    # noise so per-sample means aren't perfectly collinear with condition.
    base = np.array([0.15 if c == "DCM" else 0.20 for c in conditions])
    adata.obs["pct_counts_cardiomyocyte"] = base + rng.normal(0, 0.02, len(base))
    return adata, {"OXPHOS": pw_genes}


def test_pathway_spillover_end_to_end():
    adata, gene_sets = _make_pathway_adata()
    r = pathway_spillover(
        adata, "OXPHOS", gene_sets["OXPHOS"], "sample", "condition",
        disease_label="DCM", control_label="Donor",
        contam_key="pct_counts_cardiomyocyte")
    # 20 samples, 10 per arm.
    assert r["n_disease"] == 10 and r["n_control"] == 10
    assert len(r["score"]) == 20 and len(r["contamination"]) == 20
    for k in ("raw_effect", "adjusted_effect", "pct_attributable", "adjusted_p"):
        assert k in r
    assert r["pathway"] == "OXPHOS"


def test_pathway_spillover_missing_metric_raises():
    adata, gene_sets = _make_pathway_adata()
    del adata.obs["pct_counts_cardiomyocyte"]
    with pytest.raises(ValueError, match="score a contamination panel"):
        pathway_spillover(adata, "OXPHOS", gene_sets["OXPHOS"], "sample",
                          "condition", "DCM", "Donor",
                          "pct_counts_cardiomyocyte")


def test_result_schema():
    rng = np.random.default_rng(3)
    d = _labels(10, 10)
    contam = rng.uniform(0.1, 0.3, 20)
    score = rng.normal(0, 1, 20)
    r = spillover_attribution(score, contam, d)
    for k in ("raw_effect", "adjusted_effect", "pct_attributable",
              "raw_effect_sd", "adjusted_effect_sd", "adjusted_p",
              "spearman_score_contam", "n_disease", "n_control"):
        assert k in r
    assert r["n_disease"] == 10 and r["n_control"] == 10
