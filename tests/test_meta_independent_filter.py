"""Sentinel tests for the cross-study independent filter on the
meta-analysis input universe.

The filter drops genes that are dropouts across studies (low
expression everywhere), so BH FDR after pooling sees a smaller
universe. These tests check the core contract: drop the obviously
dropout-y genes, keep the obviously expressed ones, respect the
chosen statistic.
"""
from __future__ import annotations

import pandas as pd
import pytest

from kosmic.meta_analysis.independent_filter import (
    apply_filter, filter_datasets_to_universe,
)


def _make_study(name: str, gene_pcts: dict[str, tuple[float, float]]) -> dict:
    """Build a {name, df} dict with pct_disease / pct_control per gene."""
    rows = [
        {"names": gene, "pct_disease": d, "pct_control": c}
        for gene, (d, c) in gene_pcts.items()
    ]
    return {"name": name, "df": pd.DataFrame(rows)}


def test_mean_max_pct_drops_dropout_genes_and_keeps_expressed():
    """Genes expressed in only one study fall below the mean threshold."""
    expressed_everywhere = {"name": "EXPR", "pcts": (0.3, 0.4)}
    dropout_one_study = {"name": "DROP", "pcts": (0.5, 0.5)}

    studies = [
        _make_study("S1", {
            expressed_everywhere["name"]: expressed_everywhere["pcts"],
            dropout_one_study["name"]: dropout_one_study["pcts"],  # only here
        }),
        _make_study("S2", {
            expressed_everywhere["name"]: (0.35, 0.42),
            "DROP": (0.0, 0.0),
        }),
        _make_study("S3", {
            expressed_everywhere["name"]: (0.4, 0.45),
            "DROP": (0.0, 0.0),
        }),
    ]

    # DROP averages max_pct = 0.5 / 3 = 0.167 across studies; EXPR
    # averages ~0.42. A 0.20 threshold separates them cleanly:
    kept = apply_filter(studies, threshold=0.20, stat="mean_max_pct")
    assert "EXPR" in kept
    assert "DROP" not in kept

    # And a permissive 0.10 threshold keeps both:
    kept_loose = apply_filter(studies, threshold=0.10, stat="mean_max_pct")
    assert "EXPR" in kept_loose
    assert "DROP" in kept_loose


def test_min_studies_above_pct_requires_replication():
    """min_studies_above_pct demands per-study expression in K studies."""
    studies = [
        _make_study("S1", {"REPL": (0.3, 0.3), "ONESHOT": (0.9, 0.9)}),
        _make_study("S2", {"REPL": (0.3, 0.3), "ONESHOT": (0.0, 0.0)}),
        _make_study("S3", {"REPL": (0.3, 0.3), "ONESHOT": (0.0, 0.0)}),
    ]

    kept = apply_filter(
        studies, threshold=0.10, stat="min_studies_above_pct",
        min_studies=2,
    )
    assert "REPL" in kept           # >= 0.10 in 3 studies
    assert "ONESHOT" not in kept    # >= 0.10 in only 1 study


def test_filter_handles_studies_missing_pct_columns():
    """Studies without pct columns are skipped, not fatal."""
    df_with = pd.DataFrame([
        {"names": "G", "pct_disease": 0.5, "pct_control": 0.5},
    ])
    df_without = pd.DataFrame([{"names": "G", "logfoldchanges": 0.1}])

    studies = [{"name": "S1", "df": df_with}, {"name": "S2", "df": df_without}]
    kept = apply_filter(studies, threshold=0.10, stat="mean_max_pct")
    assert kept == {"G"}


def test_filter_returns_empty_when_no_data():
    """Empty input -> empty output, not an error."""
    assert apply_filter([], threshold=0.05) == set()
    assert apply_filter(
        [{"name": "S", "df": pd.DataFrame()}], threshold=0.05
    ) == set()


def test_filter_rejects_unknown_stat():
    studies = [_make_study("S", {"G": (0.5, 0.5)})]
    with pytest.raises(ValueError, match="Unknown stat"):
        apply_filter(studies, threshold=0.1, stat="bogus")


def test_filter_datasets_to_universe_restricts_each_df():
    """The dataset filter keeps only rows whose names are in the universe."""
    studies = [
        _make_study("S1", {"A": (0.5, 0.5), "B": (0.1, 0.1), "C": (0.0, 0.0)}),
        _make_study("S2", {"A": (0.4, 0.4), "B": (0.2, 0.2)}),
    ]
    out = filter_datasets_to_universe(studies, gene_universe={"A"})
    assert [list(s["df"]["names"]) for s in out] == [["A"], ["A"]]


def test_filter_datasets_to_universe_none_passthrough():
    """gene_universe=None returns the input unchanged."""
    studies = [_make_study("S", {"G": (0.5, 0.5)})]
    assert filter_datasets_to_universe(studies, None) is studies
