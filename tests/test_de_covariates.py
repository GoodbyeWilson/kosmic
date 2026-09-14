"""Adjusting pooled DE for a study covariate.

The mega arm fits one model over several cohorts. Without the cohort in
the design its offset lands in the residual, inflating dispersion and
shrinking every effect toward zero. These cover carrying the covariate
onto the pseudobulk table, and refusing the designs that cannot be fitted.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from kosmic.de.de_analysis import _usable_covariates, create_pseudobulk

GENES = ["G1", "G2", "G3"]


def _adata(samples, conditions, studies=None, cells_per=5):
    rows, obs = [], []
    rng = np.random.default_rng(0)
    for i, (s, c) in enumerate(zip(samples, conditions)):
        for _ in range(cells_per):
            rows.append(rng.poisson(5, len(GENES)))
            rec = {"sample": s, "condition": c, "_role": c}
            if studies is not None:
                rec["study"] = studies[i]
            obs.append(rec)
    return ad.AnnData(
        X=csr_matrix(np.array(rows, dtype=np.float32)),
        var=pd.DataFrame(index=GENES),
        obs=pd.DataFrame(obs, index=[f"c{i}" for i in range(len(obs))]))


# --- carrying the covariate through -----------------------------------------

def test_a_sample_level_covariate_reaches_sample_df():
    a = _adata(["S1", "S2"], ["disease", "control"], ["A", "B"])
    _, sample_df, _ = create_pseudobulk(
        a, GENES, "sample", "condition", min_cells=1, aggregate="sum",
        covariates=["study"])
    assert set(sample_df["study"]) == {"A", "B"}


def test_no_covariates_leaves_the_table_as_it_was():
    a = _adata(["S1", "S2"], ["disease", "control"], ["A", "B"])
    _, sample_df, _ = create_pseudobulk(
        a, GENES, "sample", "condition", min_cells=1, aggregate="sum")
    assert "study" not in sample_df.columns


def test_a_covariate_that_varies_within_a_sample_is_skipped():
    """Not a sample-level property; taking the first cell's value would
    invent one."""
    a = _adata(["S1"], ["disease"], ["A"], cells_per=4)
    a.obs["messy"] = ["x", "x", "y", "y"]
    _, sample_df, _ = create_pseudobulk(
        a, GENES, "sample", "condition", min_cells=1, aggregate="sum",
        covariates=["messy"])
    assert "messy" not in sample_df.columns


def test_a_missing_column_is_ignored():
    a = _adata(["S1", "S2"], ["disease", "control"])
    _, sample_df, _ = create_pseudobulk(
        a, GENES, "sample", "condition", min_cells=1, aggregate="sum",
        covariates=["not_a_column"])
    assert "not_a_column" not in sample_df.columns


# --- which covariates can actually be fitted --------------------------------

def _df(**cols):
    return pd.DataFrame(cols)


def test_a_balanced_covariate_is_usable():
    df = _df(study=["A", "A", "B", "B"])
    conds = ["disease", "control", "disease", "control"]
    assert _usable_covariates(df, ["study"], conds) == ["study"]


def test_a_constant_covariate_is_dropped():
    """No information, and it makes the model matrix singular."""
    df = _df(study=["A", "A", "A", "A"])
    conds = ["disease", "control", "disease", "control"]
    assert _usable_covariates(df, ["study"], conds) == []


def test_a_covariate_unique_per_sample_is_dropped():
    """Fits the data perfectly, leaving nothing for condition."""
    df = _df(donor=["d1", "d2", "d3", "d4"])
    conds = ["disease", "control", "disease", "control"]
    assert _usable_covariates(df, ["donor"], conds) == []


def test_a_covariate_collinear_with_condition_is_dropped():
    """One cohort all disease, the other all control: study and
    condition explain the same split, so the contrast is not estimable."""
    df = _df(study=["A", "A", "B", "B"])
    conds = ["disease", "disease", "control", "control"]
    assert _usable_covariates(df, ["study"], conds) == []


def test_order_is_preserved():
    # Both must cross condition, or they are dropped as collinear -- so
    # neither can simply alternate with it.
    df = _df(study=["A", "A", "A", "B", "B", "B"],
             sex=["m", "m", "f", "f", "m", "f"])
    conds = ["disease", "control", "disease",
             "control", "disease", "control"]
    assert _usable_covariates(df, ["study", "sex"], conds) == ["study", "sex"]
    assert _usable_covariates(df, ["sex", "study"], conds) == ["sex", "study"]


@pytest.mark.parametrize("covariates", [None, []])
def test_nothing_requested_gives_nothing(covariates):
    df = _df(study=["A", "B"])
    assert _usable_covariates(df, covariates, ["disease", "control"]) == []


# --- the fitted design ------------------------------------------------------

def _pooled(n_genes=60, seed=0):
    """Two cohorts, both with disease and control, plus a cohort offset."""
    rng = np.random.default_rng(seed)
    study = np.array(["A"] * 8 + ["B"] * 6)
    role = np.array(["disease"] * 4 + ["control"] * 4
                    + ["disease"] * 3 + ["control"] * 3)
    base = rng.lognormal(4, 0.8, n_genes)
    counts = rng.poisson(base[None, :] * np.ones((len(role), 1)))
    counts = (counts * np.where(study[:, None] == "B", 3.0, 1.0)).astype(float)
    df = pd.DataFrame({"sample": [f"D{i}" for i in range(len(role))],
                       "role": role, "study": study, "n_cells": 500})
    return counts, df, [f"G{i}" for i in range(n_genes)]


def test_the_covariate_enters_the_design():
    from kosmic.de.de_analysis import _run_deseq2
    counts, df, genes = _pooled()
    res = _run_deseq2(counts, df, genes, covariates=["study"])
    assert res.attrs["design"] == "~study + condition"
    assert res.attrs["covariates_used"] == ["study"]


def test_without_covariates_the_design_is_condition_only():
    from kosmic.de.de_analysis import _run_deseq2
    counts, df, genes = _pooled()
    res = _run_deseq2(counts, df, genes)
    assert res.attrs["design"] == "~condition"
    assert res.attrs["covariates_used"] == []


def test_shrinkage_is_applied_to_condition_not_the_covariate():
    """With a covariate present, taking the first non-intercept
    coefficient picks 'study' -- which would leave the reported fold
    changes unshrunk while claiming otherwise."""
    from kosmic.de.de_analysis import _run_deseq2
    counts, df, genes = _pooled()
    res = _run_deseq2(counts, df, genes, covariates=["study"])
    assert "condition" in res.attrs["shrinkage_status"]
    assert "study" not in res.attrs["shrinkage_status"]


def test_public_and_private_names_are_the_same_rule():
    """The DE page filters its 'Adjust for' list with this same function.

    If the alias ever diverged, the GUI would offer covariates the model
    then silently drops -- the user would read a design that was never
    fitted.
    """
    from kosmic.de.de_analysis import _usable_covariates, usable_covariates

    assert _usable_covariates is usable_covariates
