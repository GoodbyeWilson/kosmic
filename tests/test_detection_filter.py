"""Tests for the detection pre-filter's units.

The pooled rule asks what fraction of an arm's *cells* carry a gene,
with every donor thrown in together. That lets a gene pass on one
donor's enthusiasm, and on a combined object it blends cohorts so a gene
visible in one study and absent from another passes on the average.
These pin the donor- and study-stratified rules that fix both, and the
ordering rule that keeps thin donors out of the calculation.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")

from kosmic.de.de_analysis import (  # noqa: E402
    donors_meeting_min_cells, genes_passing_detection,
)


def build(spec, cells_per_donor=100, donors_per_arm=4, studies=None):
    """Build an AnnData where each gene's per-donor detection is exact.

    'spec' maps gene name -> {donor: fraction of that donor's cells with
    a non-zero count}. Donors are 'D0..' (disease) and 'C0..' (control).
    """
    donors, roles = [], []
    for i in range(donors_per_arm):
        donors += [f"D{i}"] * cells_per_donor
        roles += ["disease"] * cells_per_donor
    for i in range(donors_per_arm):
        donors += [f"C{i}"] * cells_per_donor
        roles += ["control"] * cells_per_donor
    donors = np.array(donors)
    roles = np.array(roles)

    names = list(spec)
    X = np.zeros((len(donors), len(names)), dtype=np.float32)
    for j, gene in enumerate(names):
        for donor, frac in spec[gene].items():
            idx = np.where(donors == donor)[0]
            X[idx[:int(round(frac * len(idx)))], j] = 7.0

    obs = pd.DataFrame({"sample": donors, "_role": roles})
    obs["condition"] = np.where(roles == "disease", "DCM", "Donor")
    if studies is not None:
        obs["study"] = [studies[d] for d in donors]
    obs.index = [f"c{i}" for i in range(len(donors))]
    return anndata.AnnData(X=X, obs=obs, var=pd.DataFrame(index=names))


ALL_D = {f"D{i}": 0.06 for i in range(4)}
ONE_D = {"D0": 0.60}


# --- the pooled rule (unchanged default) --------------------------------

def test_pooled_rule_is_the_default():
    a = build({"EVEN": ALL_D, "ONE_DONOR": ONE_D})
    kept = genes_passing_detection(a, ["EVEN", "ONE_DONOR"], 0.05)
    # 0.60 in one of four donors pools to 15% -- passes on one donor.
    assert set(kept) == {"EVEN", "ONE_DONOR"}


def test_disabled_returns_input_unchanged():
    a = build({"EVEN": ALL_D})
    assert genes_passing_detection(a, ["EVEN", "ABSENT"], 0) == ["EVEN", "ABSENT"]


def test_missing_roles_returns_input_unchanged():
    a = build({"EVEN": ALL_D})
    del a.obs["_role"]
    assert genes_passing_detection(a, ["EVEN"], 0.05) == ["EVEN"]


def test_genes_absent_from_the_matrix_are_kept():
    a = build({"EVEN": ALL_D})
    assert "NOT_MEASURED" in genes_passing_detection(
        a, ["EVEN", "NOT_MEASURED"], 0.05)


# --- donor stratification ------------------------------------------------

def test_donor_rule_drops_the_one_donor_gene():
    a = build({"EVEN": ALL_D, "ONE_DONOR": ONE_D})
    kept = genes_passing_detection(
        a, ["EVEN", "ONE_DONOR"], 0.05,
        donor_col="sample", min_donor_frac=0.5)
    assert kept == ["EVEN"]


def test_donor_rule_respects_the_fraction():
    # Two of four disease donors clear 5%.
    half = {"D0": 0.20, "D1": 0.20}
    a = build({"HALF": half})
    for frac, expected in ((0.25, ["HALF"]), (0.5, ["HALF"]),
                           (0.75, []), (1.0, [])):
        assert genes_passing_detection(
            a, ["HALF"], 0.05, donor_col="sample",
            min_donor_frac=frac) == expected, frac


def test_donor_rule_needs_only_one_arm():
    """A gene on in every disease donor and no control donor survives."""
    a = build({"DISEASE_ONLY": ALL_D})
    kept = genes_passing_detection(
        a, ["DISEASE_ONLY"], 0.05, donor_col="sample", min_donor_frac=1.0)
    assert kept == ["DISEASE_ONLY"]


def test_zero_fraction_falls_back_to_pooled():
    a = build({"ONE_DONOR": ONE_D})
    assert genes_passing_detection(
        a, ["ONE_DONOR"], 0.05, donor_col="sample",
        min_donor_frac=0.0) == ["ONE_DONOR"]


def test_unknown_donor_column_falls_back_to_pooled():
    a = build({"ONE_DONOR": ONE_D})
    assert genes_passing_detection(
        a, ["ONE_DONOR"], 0.05, donor_col="nope",
        min_donor_frac=0.5) == ["ONE_DONOR"]


# --- study stratification ------------------------------------------------

STUDIES = {"D0": "S1", "D1": "S1", "D2": "S2", "D3": "S2",
           "C0": "S1", "C1": "S1", "C2": "S2", "C3": "S2"}


def test_study_rule_drops_a_gene_only_one_study_can_see():
    # 20% in both S1 disease donors, absent from S2 -> pools to 10%.
    one_study = {"D0": 0.20, "D1": 0.20}
    a = build({"ONE_STUDY": one_study, "BOTH": ALL_D}, studies=STUDIES)
    pooled = genes_passing_detection(a, ["ONE_STUDY", "BOTH"], 0.05)
    assert set(pooled) == {"ONE_STUDY", "BOTH"}

    stratified = genes_passing_detection(
        a, ["ONE_STUDY", "BOTH"], 0.05, study_col="study")
    assert stratified == ["BOTH"]


def test_study_rule_combines_with_the_donor_rule():
    # Passes per-study, but inside S1 only one of two donors carries it.
    lopsided = {"D0": 0.40, "D2": 0.20, "D3": 0.20}
    a = build({"LOPSIDED": lopsided}, studies=STUDIES)
    assert genes_passing_detection(
        a, ["LOPSIDED"], 0.05, study_col="study") == ["LOPSIDED"]
    assert genes_passing_detection(
        a, ["LOPSIDED"], 0.05, study_col="study",
        donor_col="sample", min_donor_frac=1.0) == []


def test_a_study_with_no_cells_does_not_veto():
    studies = dict(STUDIES)
    a = build({"EVEN": ALL_D}, studies=studies)
    # Add a third study level that has no cells in either arm.
    a.obs["study"] = a.obs["study"].astype(str)
    a.obs.loc[a.obs.index[:0], "study"] = "S3"
    assert genes_passing_detection(
        a, ["EVEN"], 0.05, study_col="study") == ["EVEN"]


def test_unknown_study_column_is_ignored():
    a = build({"EVEN": ALL_D})
    assert genes_passing_detection(
        a, ["EVEN"], 0.05, study_col="nope") == ["EVEN"]


# --- cell mask / ordering ------------------------------------------------

def test_cell_mask_excludes_dropped_donors():
    """A thin donor's cells must not decide the gene list."""
    a = build({"THIN_DONOR_ONLY": {"D0": 0.80}})
    assert genes_passing_detection(a, ["THIN_DONOR_ONLY"], 0.05) == [
        "THIN_DONOR_ONLY"]
    mask = (a.obs["sample"] != "D0").to_numpy()
    assert genes_passing_detection(
        a, ["THIN_DONOR_ONLY"], 0.05, cell_mask=mask) == []


def test_donors_meeting_min_cells():
    a = build({"EVEN": ALL_D}, cells_per_donor=100)
    assert donors_meeting_min_cells(a, "sample", 50) is None   # all pass
    assert donors_meeting_min_cells(a, "sample", 0) is None    # disabled
    assert donors_meeting_min_cells(a, "nope", 50) is None     # no column

    obs = a.obs.copy()
    keep = np.ones(a.n_obs, bool)
    keep[np.where((obs["sample"] == "D0").to_numpy())[0][5:]] = False
    thin = a[keep].copy()
    mask = donors_meeting_min_cells(thin, "sample", 50)
    assert mask is not None
    assert set(thin.obs.loc[~mask, "sample"]) == {"D0"}


def test_mask_and_donor_rule_together():
    spec = {"CARRIED_BY_THIN": {"D0": 0.90, "D1": 0.06, "D2": 0.06, "D3": 0.06}}
    a = build(spec)
    mask = (a.obs["sample"] != "D0").to_numpy()
    # Without the mask D0 lifts the whole arm; with it, three donors at 6%
    # still clear a 50% donor requirement.
    assert genes_passing_detection(
        a, list(spec), 0.05, donor_col="sample", min_donor_frac=0.5,
        cell_mask=mask) == ["CARRIED_BY_THIN"]
    assert genes_passing_detection(
        a, list(spec), 0.07, donor_col="sample", min_donor_frac=0.5,
        cell_mask=mask) == []
