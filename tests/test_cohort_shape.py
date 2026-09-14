"""Is the loaded object one study with several donors, or several
studies pooled?

That single fact decides whether 'study' belongs in the DE design and
whether the detection filter has to run per study, so it has to be
right, and it has to be visible rather than inferred from a column list.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")

from kosmic.scrna.inspect.batch import (  # noqa: E402
    cohort_headline, cohort_shape, detect_study_column,
)


def build(donors, study_of=None, cells=10, sample_col="sample"):
    """donors: {donor_id: role}."""
    rows = []
    for donor, role in donors.items():
        rows += [(donor, role)] * cells
    obs = pd.DataFrame(rows, columns=[sample_col, "_role"])
    if study_of is not None:
        obs["study"] = [study_of[d] for d in obs[sample_col]]
    obs.index = [f"c{i}" for i in range(len(obs))]
    X = np.ones((len(obs), 3), dtype=np.float32)
    return anndata.AnnData(X=X, obs=obs, var=pd.DataFrame(index=list("ABC")))


SINGLE = {"D1": "disease", "D2": "disease", "C1": "control", "C2": "control"}
ATLAS_STUDIES = {"D1": "S1", "D2": "S2", "C1": "S1", "C2": "S2"}


# --- single study --------------------------------------------------------

def test_single_study_is_not_an_atlas():
    shape = cohort_shape(build(SINGLE))
    assert shape["is_atlas"] is False
    assert shape["n_studies"] == 1
    assert shape["study_col"] is None
    assert shape["n_donors"] == 4
    assert (shape["n_disease"], shape["n_control"]) == (2, 2)


def test_single_study_headline():
    line1, line2 = cohort_headline(cohort_shape(build(SINGLE)))
    assert line1 == "Single study"
    assert line2 == "4 donors: 2 disease / 2 control"


def test_constant_study_column_is_still_a_single_study():
    """A combined object carrying one study label is not an atlas."""
    adata = build(SINGLE, study_of={d: "S1" for d in SINGLE})
    shape = cohort_shape(adata)
    assert shape["is_atlas"] is False
    assert shape["study_col"] is None
    assert detect_study_column(adata) is None


# --- atlas ---------------------------------------------------------------

def test_atlas_is_detected_and_broken_down():
    shape = cohort_shape(build(SINGLE, study_of=ATLAS_STUDIES))
    assert shape["is_atlas"] is True
    assert shape["n_studies"] == 2
    assert shape["study_col"] == "study"
    assert [s["study"] for s in shape["studies"]] == ["S1", "S2"]
    assert all(s["n_donors"] == 2 for s in shape["studies"])
    assert all(s["n_disease"] == 1 and s["n_control"] == 1
               for s in shape["studies"])


def test_atlas_headline_names_a_short_list():
    """Two or three studies are named; the point is recognising the cohort."""
    line1, line2 = cohort_headline(
        cohort_shape(build(SINGLE, study_of=ATLAS_STUDIES)))
    assert line1 == "Atlas — S1 + S2"
    assert line2 == "4 donors: 2 disease / 2 control"


def test_atlas_headline_counts_a_long_list():
    donors = {f"D{i}": "disease" for i in range(4)}
    donors.update({f"C{i}": "control" for i in range(4)})
    studies = {d: f"S{i}" for i, d in enumerate(donors)}
    line1, _ = cohort_headline(cohort_shape(build(donors, study_of=studies)))
    assert line1 == "Atlas — 8 studies"


def test_per_study_cell_counts():
    adata = build(SINGLE, study_of=ATLAS_STUDIES, cells=7)
    shape = cohort_shape(adata)
    assert sum(s["n_cells"] for s in shape["studies"]) == adata.n_obs
    assert all(s["n_cells"] == 14 for s in shape["studies"])


def test_unbalanced_studies():
    donors = {"D1": "disease", "D2": "disease", "D3": "disease",
              "C1": "control"}
    studies = {"D1": "S1", "D2": "S1", "D3": "S2", "C1": "S2"}
    shape = cohort_shape(build(donors, study_of=studies))
    by_study = {s["study"]: s for s in shape["studies"]}
    assert by_study["S1"]["n_disease"] == 2 and by_study["S1"]["n_control"] == 0
    assert by_study["S2"]["n_disease"] == 1 and by_study["S2"]["n_control"] == 1


# --- roles ---------------------------------------------------------------

def test_excluded_donors_are_counted_separately():
    donors = dict(SINGLE, X1="exclude", X2="exclude")
    shape = cohort_shape(build(donors))
    assert shape["n_excluded"] == 2
    assert shape["n_donors"] == 6
    assert (shape["n_disease"], shape["n_control"]) == (2, 2)
    assert "(2 excluded)" in cohort_headline(shape)[1]


def test_missing_roles_degrade_gracefully():
    adata = build(SINGLE)
    del adata.obs["_role"]
    shape = cohort_shape(adata)
    assert shape["n_donors"] == 4
    assert (shape["n_disease"], shape["n_control"]) == (0, 0)
    assert "no roles set" in cohort_headline(shape)[1]


# --- degenerate inputs ---------------------------------------------------

def test_no_adata():
    shape = cohort_shape(None)
    assert shape["n_donors"] == 0
    assert cohort_headline(shape)[0] == "No donors resolved"


def test_no_sample_column():
    adata = build(SINGLE, sample_col="weird_name")
    shape = cohort_shape(adata)
    assert shape["n_donors"] == 0
    assert shape["sample_col"] is None
    # Naming it explicitly recovers everything.
    named = cohort_shape(adata, sample_col="weird_name")
    assert named["n_donors"] == 4


def test_explicit_study_column_overrides_detection():
    adata = build(SINGLE, study_of=ATLAS_STUDIES)
    adata.obs["site"] = np.where(
        adata.obs["sample"].isin(["D1", "C1"]), "London", "Glasgow")
    shape = cohort_shape(adata, study_col="site")
    assert shape["study_col"] == "site"
    assert [s["study"] for s in shape["studies"]] == ["Glasgow", "London"]


def test_unknown_study_column_falls_back_to_single_study():
    shape = cohort_shape(build(SINGLE), study_col="nope")
    assert shape["is_atlas"] is False
    assert shape["study_col"] is None
