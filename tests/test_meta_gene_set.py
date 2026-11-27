"""The genes a meta-analysis tests, and the rules that choose them
(kosmic.meta_analysis.gene_set), plus the merged arm's per-cell-type FDR list
(kosmic.de.batch.run_de_by_cell_type with a mapping).

Every rule is off by default, so a run that sets none pools what it always
pooled; each rule's effect is recorded.
"""
import os
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from kosmic.meta_analysis.gene_set import (
    build_gene_set, current_symbols, gene_set_path, inconsistent_merge_genes,
    meta_gene_sets, newest_gene_set, read_gene_set, save_gene_set,
    study_folders, genes_tested)
from kosmic.meta_analysis.saved import _GENE_RUN
from kosmic.provenance import record_stage


def _study(names, pct):
    return {"name": "s", "df": pd.DataFrame({
        "names": names, "pct_disease": pct, "pct_control": pct})}


DATASETS = [
    _study(["TTN", "NPPA", "HNRNPU", "ENSG000001", "C11orf58"], [.5, .5, .5, .5, .01]),
    _study(["TTN", "NPPA", "HNRNPU", "ENSG000001", "C11orf58"], [.5, .5, .5, .5, .01]),
]


# --- building the set --------------------------------------------------------

def test_no_rule_set_means_no_gene_set():
    genes, record = build_gene_set(DATASETS, min_studies=2)
    assert genes is None
    assert record == {"n_genes_reported": 5, "rules": [], "n_genes": 5}


def test_each_rule_removes_its_genes_and_is_recorded_in_order():
    genes, record = build_gene_set(
        DATASETS, min_studies=2, filter_threshold=0.05,
        filter_stat="min_studies_above_pct",
        symbols={"TTN", "NPPA", "HNRNPU", "C11orf58"},
        exclude={"HNRNPU"}, exclude_reason="merged inconsistently across studies")
    assert genes == {"TTN", "NPPA"}
    assert [(s["rule"], s["n_before"], s["n_removed"]) for s in record["rules"]] == [
        ("expression filter", 5, 1), ("current symbols only", 4, 1),
        ("excluded genes", 3, 1)]
    assert record["rules"][2]["reason"] == "merged inconsistently across studies"
    assert record["n_genes"] == 2


def test_current_human_symbols_keep_their_lower_case_orf():
    symbols = current_symbols("human")
    assert {"TTN", "C11orf58", "MT-ND1"} <= symbols
    assert "ENSG000001" not in symbols and "C11ORF58" not in symbols


def test_tested_genes_are_those_every_method_returned():
    dfs = {"reml": pd.DataFrame({"names": ["A", "B", "C"]}),
           "stouffer": pd.DataFrame({"names": ["A", "B"]}),
           "vote_count": pd.DataFrame({"names": ["A", "B", "C", "D"]})}
    assert genes_tested(dfs, ["reml", "stouffer", "vote_count"]) == ["A", "B"]


# --- saving and finding it ---------------------------------------------------

def test_the_gene_set_sits_beside_the_results_and_is_not_a_run(tmp_path):
    csv = tmp_path / "consensus_reml_stouffer_analytical.csv"
    path = save_gene_set(gene_set_path(csv), {"B", "A", "C11orf58"})
    assert path.name == "consensus_reml_stouffer_analytical_genes.csv"
    assert read_gene_set(path) == {"A", "B", "C11orf58"}
    assert not _GENE_RUN.match(path.name)


def test_the_newest_saved_gene_set_is_found(tmp_path):
    old = save_gene_set(tmp_path / "consensus_reml_analytical_genes.csv", {"A"})
    new = save_gene_set(tmp_path / "consensus_reml_cc_genes.csv", {"B"})
    os.utime(old, (time.time() - 100, time.time() - 100))
    assert newest_gene_set(tmp_path) == new
    assert newest_gene_set(tmp_path / "missing") is None


def test_each_cell_type_finds_its_own_meta_analysis(tmp_path):
    folder = tmp_path / "meta_analysis" / "decontX_counts" / "Endothelial_Cell"
    path = save_gene_set(folder / "consensus_reml_analytical_genes.csv", {"A"})
    found = meta_gene_sets(tmp_path, "decontX_counts", ["Endothelial Cell", "Fibroblast"])
    assert found == {"Endothelial Cell": path, "Fibroblast": None}
    # The original counts are their own folder.
    assert meta_gene_sets(tmp_path, "raw", ["Endothelial Cell"]) == {"Endothelial Cell": None}


# --- inconsistent merges -----------------------------------------------------

def _write(path, genes):
    path.parent.mkdir(parents=True, exist_ok=True)
    ad.AnnData(X=csr_matrix(np.ones((3, len(genes)), dtype=np.float32)),
               var=pd.DataFrame(index=list(genes)),
               obs=pd.DataFrame(index=["c0", "c1", "c2"])).write_h5ad(path)


def test_the_studies_of_a_meta_analysis_are_checked_for_inconsistent_merges(tmp_path):
    for accession, merged in (("S1", {"HNRNPU": ["HNRNPU", "HNRNPU-AS1"]}),
                              ("S2", {}), ("S3", None)):
        processed = tmp_path / accession / "processed_data"
        _write(processed / f"{accession}.h5ad", ["TTN", "HNRNPU"])
        if merged is not None:
            record_stage(processed, accession, "gene_names", params={"merged": merged})
    genes, unknown = inconsistent_merge_genes(tmp_path, ["S1", "S2", "S3"])
    assert genes == {"HNRNPU"}
    assert unknown == ["S3"]


def test_dataset_names_map_to_their_study_folders(tmp_path):
    datasets = [
        {"name": "S1_Endothelial_Cell",
         "path": str(tmp_path / "S1" / "statistics" / "de.csv")},
        {"name": "imported", "path": str(tmp_path / "x.csv"), "external": True},
        {"name": "elsewhere", "path": str(tmp_path.parent / "other.csv")},
    ]
    assert study_folders(datasets, tmp_path) == {
        "S1_Endothelial_Cell": "S1", "imported": None, "elsewhere": None}


# --- the merged arm: one FDR list per cell type ------------------------------

def _adata(n_genes=20, seed=0):
    rng = np.random.default_rng(seed)
    rows = [(d, role, ct) for role, donors in (("disease", ["D1", "D2", "D3"]),
                                               ("control", ["C1", "C2", "C3"]))
            for d in donors for ct in ("Endothelial", "Fibroblast") for _ in range(20)]
    obs = pd.DataFrame(rows, columns=["donor", "_role", "cell_type"])
    obs["condition"] = obs["_role"].map({"disease": "DCM", "control": "Donor"})
    obs.index = [f"cell{i}" for i in range(len(obs))]
    X = rng.poisson(6, size=(len(obs), n_genes)).astype(np.float32)
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes)]))


def test_a_mapping_gives_each_cell_type_its_own_fdr_list(tmp_path):
    from kosmic.de.batch import run_de_by_cell_type
    lists = {"Endothelial": {f"GENE{i}" for i in range(5)},
             "Fibroblast": {f"GENE{i}" for i in range(5, 12)}}
    runs = run_de_by_cell_type(
        _adata(), "cell_type", "donor", "condition", tmp_path / "GSE1", "GSE1",
        de_method="ttest", min_cells=10, detection_min_pct=0.0,
        fdr_genes=lists, write=False)
    by = {r.cell_type: r for r in runs}
    assert set(by["Endothelial"].de_results["names"]) <= lists["Endothelial"]
    assert set(by["Fibroblast"].de_results["names"]) <= lists["Fibroblast"]
    assert by["Fibroblast"].n_genes_tested > by["Endothelial"].n_genes_tested


def test_a_cell_type_missing_from_the_mapping_is_skipped_not_run_genome_wide(tmp_path):
    from kosmic.de.batch import run_de_by_cell_type
    runs = run_de_by_cell_type(
        _adata(), "cell_type", "donor", "condition", tmp_path / "GSE1", "GSE1",
        de_method="ttest", min_cells=10, detection_min_pct=0.0,
        fdr_genes={"Endothelial": {"GENE0", "GENE1"}}, write=False)
    by = {r.cell_type: r for r in runs}
    assert by["Endothelial"].status == "ok"
    assert by["Fibroblast"].status == "skipped"
    assert by["Fibroblast"].message == "no FDR gene list for this cell type"
