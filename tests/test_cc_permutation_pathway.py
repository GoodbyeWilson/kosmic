"""Pathway-level case/control permutation calibration.

Not wired into the GUI yet, so this is the only thing that exercises it.
Planted data: one pathway shifted up in disease in every study, one null
pathway; the calibrated p-value must separate them, and a permutation
null built from label shuffles must leave the null pathway near 1.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from kosmic.meta_analysis.cc_permutation import (
    fast_pathway_de, pathway_cc_permutation,
)

# Scores are CPM-based, so any planted increase dilutes every other
# gene and a "null" pathway picks up a small, consistent negative shift
# that pooling across studies then calls significant. The plant is
# therefore balanced: 8 genes doubled, 16 halved, library unchanged in
# expectation; the NULL pathway is drawn from the untouched genes.
GENES = [f"G{i}" for i in range(400)]
PATHWAYS = {"UP": GENES[:8], "NULL": GENES[200:208]}


def _study(rng, name, n_per_group=6, shift=2.0):
    n = 2 * n_per_group
    counts = rng.poisson(50, size=(n, len(GENES))).astype(float)
    is_disease = np.r_[np.ones(n_per_group, bool), np.zeros(n_per_group, bool)]
    counts[np.ix_(is_disease, np.arange(8))] *= shift        # UP pathway genes
    counts[np.ix_(is_disease, np.arange(8, 24))] /= shift   # balance the library
    return {
        "expression": counts,
        "conditions": np.where(is_disease, "disease", "control"),
        "gene_names": GENES,
        "role": np.where(is_disease, "disease", "control"),
        "dataset_name": name,
    }


def _stouffer(per_study_dfs):
    """Pool per-study pathway p-values (signed by logFC) with Stouffer."""
    frames = [df.set_index("names") for df in per_study_dfs]
    names = sorted(set.intersection(*(set(f.index) for f in frames)))
    z = np.zeros(len(names))
    for f in frames:
        p = np.clip(f.loc[names, "pvals"].to_numpy(), 1e-300, 1 - 1e-16)
        sign = np.sign(f.loc[names, "logfoldchanges"].to_numpy())
        z += sign * stats.norm.isf(p / 2)
    z /= np.sqrt(len(frames))
    return pd.Series(2 * stats.norm.sf(np.abs(z)), index=names)


@pytest.fixture(scope="module")
def planted():
    rng = np.random.default_rng(0)
    return [_study(rng, f"S{i}") for i in range(3)]


def _observed(pb_data):
    dfs = []
    for s in pb_data:
        is_d = s["role"] == "disease"
        names, lfc, p, se = fast_pathway_de(
            s["expression"], is_d, np.ones(len(is_d), bool), s["gene_names"], PATHWAYS)
        dfs.append(pd.DataFrame({"names": names, "logfoldchanges": lfc,
                                 "se": se, "pvals": p, "pvals_adj": p,
                                 "dataset": s["dataset_name"]}))
    return _stouffer(dfs)


def test_fast_pathway_de_finds_the_planted_pathway(planted):
    obs = _observed(planted)
    assert set(obs.index) == {"UP", "NULL"}
    assert obs["UP"] < 1e-4
    assert obs["NULL"] > 0.05


def test_calibration_separates_planted_from_null(planted):
    obs = _observed(planted)
    seen = []
    cal = pathway_cc_permutation(
        planted, {"stouffer": obs}, PATHWAYS, {"stouffer": _stouffer},
        n_perms=60, seed=1, progress_cb=lambda i, n: seen.append((i, n)))
    assert set(cal) == {"stouffer"}
    s = cal["stouffer"]
    assert list(s.index) == list(obs.index)
    assert ((s > 0) & (s <= 1)).all()
    # Planted effect: nothing in 60 label shuffles matches it -> 1/61.
    assert s["UP"] == pytest.approx(1 / 61)
    # Null pathway: a fair share of shuffles beat the observed statistic.
    assert s["NULL"] > 0.1
    assert seen and seen[-1] == (60, 60)


def test_calibration_is_deterministic_for_a_seed(planted):
    obs = _observed(planted)
    a = pathway_cc_permutation(planted, {"m": obs}, PATHWAYS, {"m": _stouffer},
                               n_perms=20, seed=7)["m"]
    b = pathway_cc_permutation(planted, {"m": obs}, PATHWAYS, {"m": _stouffer},
                               n_perms=20, seed=7)["m"]
    pd.testing.assert_series_equal(a, b)
