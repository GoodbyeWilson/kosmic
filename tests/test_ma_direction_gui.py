"""The meta-analysis results tables show direction conflicts (issue #62).

The gene Results tab shows Up / Down / Conflict columns, reports the
number of significant genes with opposite-direction effects, and can be
limited to those genes. The pathway Results tab shows the same columns
when the pooled table has them and omits them otherwise.
"""
from __future__ import annotations

import os

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from kosmic.gui.meta_analysis.pages.gene_ma.results import ResultsTab  # noqa: E402
from kosmic.gui.meta_analysis.pages.pathway_ma_page import PathwayMAPage  # noqa: E402


@pytest.fixture(scope="module")
def app():
    inst = QApplication.instance() or QApplication([])
    yield inst


@pytest.fixture
def results_df():
    return pd.DataFrame({
        "names": ["A", "B", "C", "D"],
        "logfoldchanges": [0.5, 1.2, 0.9, 3.1],
        "fdr": [0.3, 0.004, 0.02, 1e-6],
        "n_studies": [3, 3, 3, 3],
        "n_up": [1, 1, 2, 0],
        "n_down": [1, 1, 0, 2],
        "direction_conflict": [True, True, False, False],
    })


def _headers(table):
    return [table.horizontalHeaderItem(i).text()
            for i in range(table.columnCount())]


def _names(table):
    return sorted(table.item(r, 0).text() for r in range(table.rowCount()))


def test_gene_results_columns_label_and_filter(app, results_df):
    tab = ResultsTab()
    tab.populate(results_df, ["fisher"])

    assert {"Up", "Down", "Conflict"} <= set(_headers(tab.table))
    # A is a conflict but not significant; B is both.
    assert "3 consensus genes" in tab.sig_label.text()
    assert "1 with opposite-direction effects" in tab.sig_label.text()
    assert tab.conflict_only_cb.isEnabled()

    tab.conflict_only_cb.setChecked(True)
    assert _names(tab.table) == ["A", "B"]
    tab.conflict_only_cb.setChecked(False)
    assert _names(tab.table) == ["A", "B", "C", "D"]


def test_gene_results_without_direction_columns(app, results_df):
    """A result saved before the columns existed still displays."""
    tab = ResultsTab()
    tab.conflict_only_cb.setChecked(True)
    old = results_df.drop(columns=["n_up", "n_down", "direction_conflict"])
    tab.populate(old, ["fisher"])

    assert "Conflict" not in _headers(tab.table)
    assert not tab.conflict_only_cb.isEnabled()
    assert not tab.conflict_only_cb.isChecked()
    assert tab.table.rowCount() == 4
    assert "opposite" not in tab.sig_label.text()


def test_pathway_results_columns_follow_the_table(app, results_df):
    page = PathwayMAPage()
    page._populate_results_table(results_df)
    assert {"Up", "Down", "Conflict"} <= set(_headers(page._results_table))

    page._populate_results_table(
        results_df.drop(columns=["n_up", "n_down", "direction_conflict"]))
    assert "Conflict" not in _headers(page._results_table)


def test_deseq2_vif_summaries_carry_adjusted_pvalues(app):
    """The default pathway scoring must yield countable per-study tables."""
    from kosmic.gui.meta_analysis.pages.pathway_ma_page import PathwayMAWorker

    genes = [f'G{i}' for i in range(6)]
    up = pd.DataFrame({'names': genes, 'se': 0.1,
                       'logfoldchanges': [2.0, 2.2, 1.8, 0.1, -0.1, 0.0]})
    down = up.assign(logfoldchanges=-up['logfoldchanges'])
    worker = PathwayMAWorker(
        [], [up, down], None, ['reml'], 2, scoring_method='deseq2_vif',
        pathway_gene_sets={'PW_A': genes[:3], 'PW_B': genes[3:]})
    per_study = worker._build_deseq2_vif_datasets()

    assert all({'pvals', 'pvals_adj'} <= set(df.columns) for df in per_study)
    from kosmic.meta_analysis.direction import annotate_direction
    out = annotate_direction(
        pd.DataFrame({'names': ['PW_A', 'PW_B']}), per_study).set_index('names')
    assert bool(out.loc['PW_A', 'direction_conflict'])
    assert not out.loc['PW_B', 'direction_conflict']
