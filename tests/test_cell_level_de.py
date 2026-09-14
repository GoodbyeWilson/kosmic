"""Sentinel tests for cell-level (ignore-donor) gene DE.

Cell mode runs scanpy rank_genes_groups (Wilcoxon) on individual cells
split by disease/control, instead of per-donor pseudobulk. These tests
pin that it (a) returns the standard DE schema, (b) recovers a planted
signal, and (c) runs even when there are too few donors for pseudobulk.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad

from kosmic.de.de_analysis import run_cell_level_de, run_de_pipeline


def _make_adata(n_disease_donors, n_control_donors, cells_per=150,
                n_genes=40, n_planted=5, seed=0):
    """Synthetic adata with `n_planted` genes up in disease.

    Mirrors KOSMIC's processed state: adata.raw + layers['counts'] hold raw
    counts, adata.X holds log1p(CPM). run_cell_level_de must normalise the
    counts itself rather than testing the raw counts directly.
    """
    rng = np.random.default_rng(seed)
    genes = [f"G{i:02d}" for i in range(n_genes)]
    blocks, roles, donors = [], [], []

    for d in range(n_disease_donors):
        X = rng.poisson(2.0, size=(cells_per, n_genes)).astype(float)
        X[:, :n_planted] += rng.poisson(7.0, size=(cells_per, n_planted))
        blocks.append(X)
        roles += ["disease"] * cells_per
        donors += [f"d{d}"] * cells_per
    for c in range(n_control_donors):
        X = rng.poisson(2.0, size=(cells_per, n_genes)).astype(float)
        blocks.append(X)
        roles += ["control"] * cells_per
        donors += [f"c{c}"] * cells_per

    counts = np.vstack(blocks)
    adata = ad.AnnData(X=counts.copy(), var=pd.DataFrame(index=genes))
    adata.obs["_role"] = pd.Categorical(
        roles, categories=["control", "disease", "exclude"])
    adata.obs["donor"] = donors
    adata.layers["counts"] = counts.copy()
    adata.raw = adata.copy()                       # raw holds COUNTS

    # adata.X -> log1p(CPM): the processed state the DE workspace operates on.
    lib = np.clip(counts.sum(axis=1, keepdims=True), 1, None)
    adata.X = np.log1p(counts / lib * 1e4)
    return adata, genes[:n_planted]


def test_cell_level_de_schema_and_signal():
    adata, planted = _make_adata(2, 2)
    res = run_cell_level_de(adata)

    # Standard DE schema the rest of the pipeline expects.
    for col in ("names", "logfoldchanges", "pvals", "pvals_adj"):
        assert col in res.columns
    assert len(res) == adata.n_vars

    # Per-arm counts are cell counts, not donor counts.
    assert res["disease_samples"].iloc[0] == int((adata.obs["_role"] == "disease").sum())
    assert res["control_samples"].iloc[0] == int((adata.obs["_role"] == "control").sum())

    planted_rows = res[res["names"].isin(planted)]
    assert (planted_rows["pvals_adj"] < 0.05).all()
    assert (planted_rows["logfoldchanges"] > 0).all()

    # The planted genes are cleanly separated as the most up-regulated.
    # (Cell mode flags many genes overall — CPM normalization plus huge
    # cell n makes even small shifts significant — so we check effect-size
    # separation rather than null specificity.)
    null_rows = res[~res["names"].isin(planted)]
    assert planted_rows["logfoldchanges"].min() > null_rows["logfoldchanges"].max()


def test_cell_level_runs_with_one_donor_per_arm():
    """Pseudobulk needs >=2 samples/arm; cell mode does not care about donors."""
    adata, planted = _make_adata(1, 1)
    de_results, sig, coverage, sample_df = run_de_pipeline(
        adata, sample_col="donor", condition_col="_role",
        pathway_gene_sets={}, full_genome=True, unit="cell",
    )
    assert not de_results.empty
    assert set(planted).issubset(set(sig["names"]))
    # sample_df reports per-role cell counts, not pseudobulk samples.
    assert set(sample_df["role"]) == {"disease", "control"}


def test_cell_level_pathway_filter():
    """Without full_genome, results are filtered to the pathway gene set."""
    adata, planted = _make_adata(2, 2)
    pathway = {"planted": list(planted) + ["G10", "G11"]}
    de_results, sig, coverage, _ = run_de_pipeline(
        adata, sample_col="donor", condition_col="_role",
        pathway_gene_sets=pathway, full_genome=False, unit="cell",
    )
    assert not de_results.empty
    assert set(de_results["names"]).issubset(set(pathway["planted"]))


def test_hypothesis_mode_restricts_fdr_to_committed_genes():
    """Discovery corrects FDR genome-wide; hypothesis corrects over the
    committed gene list only. The fit is identical (raw p-values match), but
    the smaller BH denominator lifts power for the committed set."""
    adata, planted = _make_adata(3, 3)
    committed = set(planted) | {"G10", "G11", "G12"}

    common = dict(sample_col="donor", condition_col="_role",
                  pathway_gene_sets={}, full_genome=True, unit="sample",
                  detection_min_pct=0.0)

    disc, _, _, _ = run_de_pipeline(adata, fdr_genes=None, **common)
    hyp, _, _, _ = run_de_pipeline(adata, fdr_genes=committed, **common)

    # Hypothesis results are restricted to the committed set.
    assert not hyp.empty
    assert set(hyp["names"]).issubset(committed)

    # Same underlying fit: raw p-values match for the shared genes.
    d = disc.set_index("names")
    h = hyp.set_index("names")
    shared = list(h.index)
    assert np.allclose(d.loc[shared, "pvals"].values, h.loc[shared, "pvals"].values)

    # Smaller multiple-testing burden -> committed FDR no larger than genome-wide,
    # and strictly smaller for the strongest gene.
    assert (h["pvals_adj"].values <= d.loc[shared, "pvals_adj"].values + 1e-9).all()
    assert h["pvals_adj"].min() < d.loc[shared, "pvals_adj"].min()

    # The genome-wide ranking is preserved for consumers whose null is the whole
    # transcriptome (fgsea): it spans more genes than the restricted table.
    gw = hyp.attrs.get("genome_wide_ranking")
    assert gw is not None
    assert set(hyp["names"]).issubset(set(gw["names"]))
    assert len(gw) > len(hyp)


def test_hypothesis_mode_species_converts_committed_genes():
    """Committed pathway genes are usually human symbols (HK1); a mouse
    dataset names the same genes Hk1. Hypothesis mode must map the committed
    list into the dataset's species namespace before restricting -- otherwise
    every gene drops out and the run returns empty ("No genes passed
    filtering"). Regression for that bug.
    """
    rng = np.random.default_rng(3)
    mouse_genes = ["Hk1", "Pkm", "Ldha", "Gapdh", "Eno1", "Pgk1", "Aldoa",
                   "Tpi1", "Pfkl", "Gpi1", "Slc2a1", "Cs", "Idh1", "Mdh1",
                   "Sdha", "Fh1", "Aco2", "Ogdh", "Sucla2", "Pkm2"]
    n_planted = 4
    blocks, roles, donors = [], [], []
    for d in range(3):
        X = rng.poisson(3.0, size=(120, len(mouse_genes))).astype(float)
        X[:, :n_planted] += rng.poisson(9.0, size=(120, n_planted))
        blocks.append(X)
        roles += ["disease"] * 120
        donors += [f"d{d}"] * 120
    for c in range(3):
        blocks.append(rng.poisson(3.0, size=(120, len(mouse_genes))).astype(float))
        roles += ["control"] * 120
        donors += [f"c{c}"] * 120

    counts = np.vstack(blocks)
    adata = ad.AnnData(X=counts.copy(), var=pd.DataFrame(index=mouse_genes))
    adata.obs["_role"] = pd.Categorical(
        roles, categories=["control", "disease", "exclude"])
    adata.obs["donor"] = donors
    adata.layers["counts"] = counts.copy()
    adata.raw = adata.copy()
    lib = np.clip(counts.sum(axis=1, keepdims=True), 1, None)
    adata.X = np.log1p(counts / lib * 1e4)

    # Committed list in HUMAN symbol case, as pathway .gmt files ship them.
    committed_human = {g.upper() for g in mouse_genes}

    hyp, sig, _, _ = run_de_pipeline(
        adata, sample_col="donor", condition_col="_role",
        pathway_gene_sets={}, full_genome=True, unit="sample",
        detection_min_pct=0.0, fdr_genes=committed_human,
    )
    assert not hyp.empty, "committed human symbols were not mapped to mouse"
    # Results are the mouse-cased gene names, restricted to the committed set.
    assert set(hyp["names"]).issubset(set(mouse_genes))
    assert set(mouse_genes[:n_planted]).issubset(set(sig["names"]))


def test_deseq2_surfaces_independent_and_cooks_filtering():
    """DESeq2's independent + Cook's filtering is exposed, not hidden: results
    carry a 'filter_status' column and a filter-count summary, and turning
    independent filtering off removes the 'independent_filter' verdicts."""
    adata, _ = _make_adata(4, 4)
    common = dict(sample_col="donor", condition_col="_role",
                  pathway_gene_sets={}, full_genome=True, unit="sample",
                  de_method="deseq2", detection_min_pct=0.0)

    on, _, _, _ = run_de_pipeline(adata, deseq2_independent_filter=True, **common)
    off, _, _, _ = run_de_pipeline(adata, deseq2_independent_filter=False, **common)

    # Surfaced (column + counts), and the counts partition every gene.
    assert "filter_status" in on.columns
    counts = on.attrs.get("deseq2_filter_counts")
    assert counts is not None and sum(counts.values()) == len(on)

    # Off -> no independent-filter verdicts; results are not silently all-1.0.
    assert (off["filter_status"] == "independent_filter").sum() == 0
    assert on["pvals_adj"].notna().any()


def test_cell_mode_inflates_false_positives_vs_pseudobulk():
    """Pure-null data with per-donor batch effects: pseudobulk treats the
    donor as the unit and controls false positives, while cell mode ignores
    donor grouping and mistakes between-donor variance for signal. This is
    the pseudoreplication failure that makes cell mode opt-in."""
    rng = np.random.default_rng(7)
    n_genes = 50
    cells_per = 120
    blocks, roles, donors = [], [], []
    # No true disease effect; each donor gets its own multiplicative batch
    # offset so between-donor variance dominates.
    for grp, tag in (("disease", "d"), ("control", "c")):
        for k in range(4):
            scale = rng.lognormal(0.0, 0.4, size=n_genes)
            X = rng.poisson(2.0 * scale, size=(cells_per, n_genes)).astype(float)
            blocks.append(X)
            roles += [grp] * cells_per
            donors += [f"{tag}{k}"] * cells_per

    counts = np.vstack(blocks)
    genes = [f"G{i:02d}" for i in range(n_genes)]
    adata = ad.AnnData(X=counts.copy(), var=pd.DataFrame(index=genes))
    adata.obs["_role"] = pd.Categorical(
        roles, categories=["control", "disease", "exclude"])
    adata.obs["donor"] = donors
    adata.layers["counts"] = counts.copy()
    adata.raw = adata.copy()
    lib = np.clip(counts.sum(axis=1, keepdims=True), 1, None)
    adata.X = np.log1p(counts / lib * 1e4)

    _, sig_cell, _, _ = run_de_pipeline(
        adata, sample_col="donor", condition_col="_role",
        pathway_gene_sets={}, full_genome=True, unit="cell",
    )
    _, sig_pb, _, _ = run_de_pipeline(
        adata, sample_col="donor", condition_col="_role",
        pathway_gene_sets={}, full_genome=True, unit="sample",
        de_method="ttest", min_cells=10,
    )
    # Pseudobulk should find essentially nothing; cell mode over-calls.
    assert len(sig_cell) > len(sig_pb)
