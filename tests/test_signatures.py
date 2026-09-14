"""Tests for ambient-contamination signature scoring.

Pins that:
  (a) `compute_signature_pct` returns the correct per-cell percent and
      reads from `adata.raw` so it works after HVG subsetting;
  (b) `run_qc_pipeline` stores `pct_counts_<key>` and the metric survives
      cell + gene filtering and a later cell-type subset;
  (c) the bundled cardiomyocyte panel loads.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad

from kosmic.scrna.qc.signatures import (
    compute_signature_pct, flag_signature_genes, add_signature_metrics)
from kosmic.scrna.qc.filter import run_qc_pipeline
from kosmic.reference.contamination import builtin_panels, get_panel


def _make_adata(n_cells=200, seed=0):
    """Counts matrix with a known 'contaminant' panel fraction per cell.

    Panel genes contribute exactly `panel_counts` per cell; background
    genes contribute `bg_counts`, so the true panel fraction is known.
    """
    rng = np.random.default_rng(seed)
    panel = ["RYR2", "MYBPC3", "TECRL"]          # contaminant panel
    background = [f"G{i:02d}" for i in range(20)]
    genes = panel + background

    bg = rng.poisson(5.0, size=(n_cells, len(background))).astype(float)
    pan = rng.poisson(1.0, size=(n_cells, len(panel))).astype(float)
    counts = np.hstack([pan, bg])

    adata = ad.AnnData(X=counts.copy(), var=pd.DataFrame(index=genes))
    adata.obs_names = [f"c{i}" for i in range(n_cells)]
    return adata, panel


def test_compute_signature_pct_matches_manual():
    adata, panel = _make_adata()
    pct = compute_signature_pct(adata, panel, key="cm", store=True)

    panel_idx = [adata.var_names.get_loc(g) for g in panel]
    manual = (adata.X[:, panel_idx].sum(axis=1)
              / np.clip(adata.X.sum(axis=1), 1, None) * 100.0)

    assert np.allclose(pct, manual)
    assert "pct_counts_cm" in adata.obs
    assert np.allclose(adata.obs["pct_counts_cm"].to_numpy(), manual)
    # A small-but-nonzero, bounded fraction.
    assert (pct >= 0).all() and (pct <= 100).all()
    assert 0 < pct.mean() < 50


def test_flag_signature_genes():
    adata, panel = _make_adata()
    matched = flag_signature_genes(adata, panel + ["NOT_A_GENE"], key="cm")
    assert set(matched) == set(panel)
    assert adata.var["cm"].sum() == len(panel)


def test_missing_panel_genes_give_zero():
    adata, _ = _make_adata()
    pct = compute_signature_pct(adata, ["ABSENT1", "ABSENT2"], key="x")
    assert np.allclose(pct, 0.0)


def test_reads_from_raw_after_hvg_subset():
    """The whole point: panel genes gone from X but present in raw."""
    adata, panel = _make_adata()
    adata.raw = adata.copy()
    # Simulate HVG subsetting: drop the panel genes from the working matrix.
    keep = [g for g in adata.var_names if g not in panel]
    sub = adata[:, keep].copy()
    assert not any(g in sub.var_names for g in panel)   # gone from X

    pct = compute_signature_pct(sub, panel, key="cm", use_raw=True)
    # Recover the true fraction from raw despite X lacking the genes.
    raw = adata.raw
    panel_idx = [list(raw.var_names).index(g) for g in panel]
    manual = (raw.X[:, panel_idx].sum(axis=1)
              / np.clip(raw.X.sum(axis=1), 1, None) * 100.0)
    assert np.allclose(pct, np.asarray(manual).ravel())


def test_qc_pipeline_stores_and_survives_subset():
    adata, panel = _make_adata(n_cells=300)
    params = {"min_genes": 0, "min_cells": 1,
              "signature_panels": {"cm": panel}}
    adata, stats = run_qc_pipeline(adata, params)

    assert "pct_counts_cm" in adata.obs
    assert stats["signature_matched"]["cm"] == len(panel)

    # Subset to an arbitrary half -- the metric must ride along unchanged.
    half = adata[adata.obs_names[:100]].copy()
    assert "pct_counts_cm" in half.obs
    assert half.obs["pct_counts_cm"].notna().all()


def test_add_signature_metrics_returns_match_counts():
    adata, panel = _make_adata()
    matched = add_signature_metrics(adata, {"cm": panel, "empty": ["ZZZ"]},
                                    use_raw=False)
    assert matched["cm"] == len(panel)
    assert matched["empty"] == 0
    assert "pct_counts_cm" in adata.obs
    assert "pct_counts_empty" in adata.obs


def test_bundled_cardiomyocyte_panel_loads():
    panels = builtin_panels()
    assert "cardiomyocyte" in panels
    cm = get_panel("cardiomyocyte")
    assert cm["dominant_cell_type"] == "Cardiomyocyte"
    assert len(cm["genes"]) == 34
    assert len(set(cm["genes"])) == 34            # no duplicates
    assert "RYR2" in cm["genes"] and "MYBPC3" in cm["genes"]
    # No leaky broadly-expressed genes that would over-attribute the effect.
    for leaky in ("MYL2", "MYL3", "MYH7", "TTN", "MYH6"):
        assert leaky not in cm["genes"]
