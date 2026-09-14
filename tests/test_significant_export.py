"""The significant-gene table is written automatically, not on a button.

'significant_subset' is the single definition of "significant" shared by
the in-app counts, the auto-exported CSV and the batch runner, so these
pin its edge cases -- particularly that a blank FDR (independent-filtered
or a Cook's outlier) is not silently counted as a hit.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from kosmic.de.de_analysis import significant_subset
from kosmic.paths import de_stats_dir, significant_path


def frame():
    return pd.DataFrame({
        "names": ["HIT", "SMALL_FC", "NOT_SIG", "NAN_FDR", "NEG_HIT"],
        "pvals_adj": [0.001, 0.001, 0.9, np.nan, 0.01],
        "logfoldchanges": [1.5, 0.1, 2.0, 3.0, -1.5],
        "abs_logfoldchange": [1.5, 0.1, 2.0, 3.0, 1.5],
    })


def test_applies_both_thresholds():
    out = significant_subset(frame(), fdr=0.05, lfc=0.25)
    assert set(out["names"]) == {"HIT", "NEG_HIT"}


def test_negative_fold_changes_count():
    out = significant_subset(frame(), fdr=0.05, lfc=0.25)
    assert "NEG_HIT" in set(out["names"])


def test_blank_fdr_is_not_significant():
    """An independent-filtered or Cook's-outlier gene has no p-value."""
    out = significant_subset(frame(), fdr=0.05, lfc=0.25)
    assert "NAN_FDR" not in set(out["names"])


def test_thresholds_are_honoured():
    assert len(significant_subset(frame(), fdr=0.05, lfc=3.0)) == 0
    assert len(significant_subset(frame(), fdr=1.0, lfc=0.0)) == 4  # NaN still out


def test_falls_back_to_logfoldchanges_when_abs_missing():
    df = frame().drop(columns=["abs_logfoldchange"])
    out = significant_subset(df, fdr=0.05, lfc=0.25)
    assert set(out["names"]) == {"HIT", "NEG_HIT"}


def test_empty_and_none_are_passed_through():
    assert significant_subset(None) is None
    empty = pd.DataFrame()
    assert significant_subset(empty).empty


def test_frame_without_fdr_column_yields_nothing():
    df = pd.DataFrame({"names": ["A"], "logfoldchanges": [2.0]})
    assert significant_subset(df).empty


def test_result_is_a_copy():
    df = frame()
    out = significant_subset(df, fdr=0.05, lfc=0.25)
    out.loc[out.index[0], "names"] = "CHANGED"
    assert "CHANGED" not in set(df["names"])


# --- the batch runner writes one per cell type --------------------------

def test_batch_writes_a_significant_file_per_cell_type(tmp_path):
    anndata = pytest.importorskip("anndata")
    from kosmic.de import batch as B

    rng = np.random.default_rng(0)
    rows = []
    for role, donors in (("disease", ["D1", "D2", "D3"]),
                         ("control", ["C1", "C2", "C3"])):
        for d in donors:
            for ct in ("Endothelial", "Fibroblast"):
                rows += [(d, role, ct)] * 30
    obs = pd.DataFrame(rows, columns=["donor", "_role", "cell_type"])
    obs["condition"] = np.where(obs["_role"] == "disease", "DCM", "Donor")
    obs.index = [f"c{i}" for i in range(len(obs))]
    adata = anndata.AnnData(
        X=rng.poisson(6, size=(len(obs), 40)).astype(np.float32),
        obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(40)]))

    runs = B.run_de_by_cell_type(
        adata, "cell_type", "donor", "condition", tmp_path / "GSE1", "GSE1",
        de_method="ttest", min_cells=10, detection_min_pct=0.0)

    for r in runs:
        assert r.significant_path is not None
        assert r.significant_path.exists()
        sig = pd.read_csv(r.significant_path)
        # The file must agree with the count the run reported.
        assert len(sig) == r.n_significant

    names = {p.name for p in de_stats_dir(tmp_path / "GSE1").iterdir()}
    assert "GSE1_Endothelial_significant_genes.csv" in names


def test_significant_filename_is_not_mistaken_for_a_de_input():
    """Meta-analysis discovery keys on '_DE_<method>.csv'."""
    from kosmic.meta_analysis.io import _DE_RESULTS_PATTERN

    name = significant_path("/study", "GSE1_Endothelial").name
    assert name == "GSE1_Endothelial_significant_genes.csv"
    assert _DE_RESULTS_PATTERN.match(name) is None
