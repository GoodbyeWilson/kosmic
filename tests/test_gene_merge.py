"""Genes that harmonisation merged inconsistently across studies.

After harmonisation both studies carry a column called 'HNRNPU'; only the
provenance record says whether it is one feature or two summed together.
These cover the reading of that record and the comparison across studies.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from kosmic.combine import concat_studies
from kosmic.combine.gene_merge import (
    find_merge_asymmetry,
    merge_map_from_provenance,
    project_merge_asymmetry,
)
from kosmic.provenance import record_stage

GENES = ["TTN", "NPPA", "HNRNPU"]


def _write(path, genes, n_obs=4):
    path.parent.mkdir(parents=True, exist_ok=True)
    a = ad.AnnData(
        X=csr_matrix(np.ones((n_obs, len(genes)), dtype=np.float32)),
        var=pd.DataFrame(index=list(genes)),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n_obs)]))
    a.write_h5ad(path)
    return path


# --- reading the record -----------------------------------------------------

def test_no_record_at_all_reads_as_unknown(tmp_path):
    assert merge_map_from_provenance(tmp_path) is None


def test_harmonised_with_no_merges_reads_as_empty_not_unknown(tmp_path):
    record_stage(tmp_path, "S1", "gene_names",
                 params={"renamed": 3, "merged": {}})
    assert merge_map_from_provenance(tmp_path) == {}


def test_merge_map_round_trips(tmp_path):
    record_stage(tmp_path, "S1", "gene_names",
                 params={"merged": {"HNRNPU": ["HNRNPU", "HNRNPU-AS1"]}})
    assert merge_map_from_provenance(tmp_path) == {
        "HNRNPU": ["HNRNPU", "HNRNPU-AS1"]}


def test_the_latest_harmonisation_wins(tmp_path):
    record_stage(tmp_path, "S1", "gene_names", params={"merged": {"A": ["A", "B"]}})
    record_stage(tmp_path, "S1", "gene_names", params={"merged": {"C": ["C", "D"]}})
    assert merge_map_from_provenance(tmp_path) == {"C": ["C", "D"]}


# --- the comparison ---------------------------------------------------------

def test_two_features_in_one_study_one_in_another_is_flagged():
    result = find_merge_asymmetry({
        "a": (GENES, {"HNRNPU": ["HNRNPU", "HNRNPU-AS1"]}),
        "b": (GENES, {}),
    })
    assert result["genes"] == {"HNRNPU": {"a": 2, "b": 1}}


def test_the_same_merge_everywhere_is_not_an_asymmetry():
    merges = {"HNRNPU": ["HNRNPU", "HNRNPU-AS1"]}
    result = find_merge_asymmetry({"a": (GENES, merges), "b": (GENES, merges)})
    assert result["genes"] == {}


def test_a_plain_rename_is_not_an_asymmetry():
    """AARS -> AARS1 is one feature in both studies. That is the point of
    harmonising, not a problem with it."""
    result = find_merge_asymmetry({
        "a": (["AARS1"], {"AARS1": ["AARS"]}),
        "b": (["AARS1"], {}),
    })
    assert result["genes"] == {}


def test_a_gene_missing_from_one_study_is_not_flagged():
    """The inner join drops it anyway, so it is not an asymmetry to fix."""
    result = find_merge_asymmetry({
        "a": (GENES, {"HNRNPU": ["HNRNPU", "HNRNPU-AS1"]}),
        "b": (["TTN", "NPPA"], {}),
    })
    assert result["genes"] == {}


def test_studies_without_a_record_are_named_not_guessed_at():
    result = find_merge_asymmetry({
        "a": (GENES, {"HNRNPU": ["HNRNPU", "HNRNPU-AS1"]}),
        "b": (GENES, None),
    })
    assert result["unknown"] == ["b"]


def test_a_single_study_has_nothing_to_compare_against():
    assert find_merge_asymmetry({"a": (GENES, {})})["genes"] == {}


# --- end to end over a project ---------------------------------------------

def test_project_walk_finds_the_asymmetry(tmp_path):
    for accession, merged in (("S1", {"HNRNPU": ["HNRNPU", "HNRNPU-AS1"]}),
                              ("S2", {})):
        sdir = tmp_path / accession
        _write(sdir / "processed_data" / f"{accession}.h5ad", GENES)
        record_stage(sdir / "processed_data", accession, "gene_names",
                     params={"merged": merged})
    result = project_merge_asymmetry(tmp_path)
    assert set(result["genes"]) == {"HNRNPU"}
    assert result["unknown"] == []


def test_project_walk_reports_unharmonised_studies(tmp_path):
    for accession in ("S1", "S2"):
        _write(tmp_path / accession / "processed_data" / f"{accession}.h5ad",
               GENES)
    result = project_merge_asymmetry(tmp_path)
    assert result["unknown"] == ["S1", "S2"]


# --- the drop itself --------------------------------------------------------

def test_concat_drops_the_named_genes(tmp_path):
    paths = [_write(tmp_path / s / "processed_data" / f"{s}.h5ad", GENES)
             for s in ("S1", "S2")]
    out = concat_studies(paths, tmp_path / "master.h5ad",
                         drop_genes=["HNRNPU"])
    master = ad.read_h5ad(out)
    assert list(master.var_names) == ["TTN", "NPPA"]
    assert master.n_obs == 8


def test_concat_ignores_a_drop_name_no_study_has(tmp_path):
    paths = [_write(tmp_path / s / "processed_data" / f"{s}.h5ad", GENES)
             for s in ("S1", "S2")]
    out = concat_studies(paths, tmp_path / "master.h5ad",
                         drop_genes=["NOT_A_GENE"])
    assert ad.read_h5ad(out).n_vars == 3


@pytest.mark.parametrize("drop", [None, []])
def test_concat_without_a_drop_list_is_unchanged(tmp_path, drop):
    paths = [_write(tmp_path / s / "processed_data" / f"{s}.h5ad", GENES)
             for s in ("S1", "S2")]
    out = concat_studies(paths, tmp_path / "master.h5ad", drop_genes=drop)
    assert ad.read_h5ad(out).n_vars == 3


# --- leaving out excluded cells --------------------------------------------

def _write_roles(path, genes, roles):
    path.parent.mkdir(parents=True, exist_ok=True)
    n = len(roles)
    a = ad.AnnData(
        X=csr_matrix(np.ones((n, len(genes)), dtype=np.float32)),
        var=pd.DataFrame(index=list(genes)),
        obs=pd.DataFrame({"_role": roles},
                         index=[f"c{i}" for i in range(n)]))
    a.write_h5ad(path)
    return path


def test_excluded_cells_are_left_out_when_asked(tmp_path):
    paths = [
        _write_roles(tmp_path / "S1" / "processed_data" / "S1.h5ad", GENES,
                     ["disease", "control", "exclude", "exclude"]),
        _write_roles(tmp_path / "S2" / "processed_data" / "S2.h5ad", GENES,
                     ["disease", "control"]),
    ]
    out = concat_studies(paths, tmp_path / "master.h5ad", drop_excluded=True)
    master = ad.read_h5ad(out)
    assert master.n_obs == 4                       # 2 kept + 2
    assert "exclude" not in set(master.obs["_role"].astype(str))


def test_excluded_cells_are_kept_by_default(tmp_path):
    """Unchanged behaviour unless the caller asks for the drop."""
    paths = [
        _write_roles(tmp_path / "S1" / "processed_data" / "S1.h5ad", GENES,
                     ["disease", "exclude"]),
        _write_roles(tmp_path / "S2" / "processed_data" / "S2.h5ad", GENES,
                     ["control", "exclude"]),
    ]
    out = concat_studies(paths, tmp_path / "master.h5ad")
    assert ad.read_h5ad(out).n_obs == 4


def test_a_study_with_no_role_column_is_unaffected(tmp_path):
    s1 = tmp_path / "S1" / "processed_data" / "S1.h5ad"
    s1.parent.mkdir(parents=True)
    _write(s1, GENES, n_obs=3)                     # no _role at all
    s2 = _write_roles(tmp_path / "S2" / "processed_data" / "S2.h5ad", GENES,
                      ["disease", "exclude"])
    out = concat_studies([s1, s2], tmp_path / "master.h5ad",
                         drop_excluded=True)
    assert ad.read_h5ad(out).n_obs == 4            # 3 + 1
