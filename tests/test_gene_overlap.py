"""Cross-study gene overlap.

Combining studies is an inner join on gene names, so the master dataset only
ever carries genes present in every study. These cover the ways that number
goes wrong quietly: a study annotated against an older build, a file that
cannot be read, and a project where nothing has been processed yet.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from kosmic.scrna.load.gene_overlap import (
    gene_overlap,
    project_gene_overlap,
    read_var_names,
)

CORE = ["TTN", "NPPA", "MYH7"]


def _write(path, genes):
    a = ad.AnnData(
        X=csr_matrix(np.eye(len(genes), dtype=np.float32)),
        var=pd.DataFrame(index=list(genes)))
    a.write_h5ad(path)
    return path


def test_read_var_names_returns_the_index(tmp_path):
    p = _write(tmp_path / "a.h5ad", CORE)
    assert read_var_names(p) == CORE


def test_shared_is_the_intersection_not_the_smallest_study(tmp_path):
    a = _write(tmp_path / "a.h5ad", CORE + ["GENE_A"])
    b = _write(tmp_path / "b.h5ad", CORE + ["GENE_B"])
    res = gene_overlap({"a": a, "b": b})
    assert res["shared"] == 3
    assert res["union"] == 5
    assert res["studies"]["a"]["n_genes"] == 4


def test_protein_coding_subset_is_reported(tmp_path):
    p = _write(tmp_path / "a.h5ad", CORE + ["MALAT1"])
    res = gene_overlap({"a": p}, protein_coding=CORE)
    assert res["studies"]["a"]["n_protein_coding"] == 3
    assert res["shared_protein_coding"] == 3


def test_harmonisation_recovers_a_retired_symbol(tmp_path):
    """The whole point of the panel: an old build costs you shared genes."""
    a = _write(tmp_path / "a.h5ad", CORE)
    b = _write(tmp_path / "b.h5ad", ["TTN", "NPPA", "MYHCB"])  # retired MYH7
    res = gene_overlap({"a": a, "b": b}, lookup={"MYHCB": "MYH7"})
    assert res["shared"] == 2
    assert res["harmonised"]["shared"] == 3


def test_an_unreadable_study_is_reported_not_fatal(tmp_path):
    good = _write(tmp_path / "a.h5ad", CORE)
    bad = tmp_path / "b.h5ad"
    bad.write_bytes(b"not an h5ad")
    res = gene_overlap({"a": good, "b": bad})
    assert res["studies"]["b"]["error"]
    assert res["studies"]["b"]["n_genes"] is None
    assert res["shared"] == 3          # the readable study still counts


def test_no_readable_studies_gives_zero_not_a_crash(tmp_path):
    bad = tmp_path / "b.h5ad"
    bad.write_bytes(b"not an h5ad")
    res = gene_overlap({"b": bad})
    assert res["shared"] == 0
    assert res["harmonised"] is None


def test_project_walk_uses_the_newest_processed_file(tmp_path):
    for accession, genes in (("S1", CORE + ["A"]), ("S2", CORE + ["B"])):
        proc = tmp_path / accession / "processed_data"
        proc.mkdir(parents=True)
        _write(proc / f"{accession}.h5ad", genes)
    res = project_gene_overlap(tmp_path)
    assert set(res["studies"]) == {"S1", "S2"}
    assert res["shared"] == 3


def test_project_skips_studies_with_nothing_processed(tmp_path):
    proc = tmp_path / "S1" / "processed_data"
    proc.mkdir(parents=True)
    _write(proc / "S1.h5ad", CORE)
    (tmp_path / "S2" / "raw_data").mkdir(parents=True)
    res = project_gene_overlap(tmp_path)
    assert set(res["studies"]) == {"S1"}


@pytest.mark.parametrize("accessions", [["S1"], []])
def test_project_honours_an_explicit_study_list(tmp_path, accessions):
    for accession in ("S1", "S2"):
        proc = tmp_path / accession / "processed_data"
        proc.mkdir(parents=True)
        _write(proc / f"{accession}.h5ad", CORE)
    res = project_gene_overlap(tmp_path, accessions=accessions)
    assert set(res["studies"]) == set(accessions)
