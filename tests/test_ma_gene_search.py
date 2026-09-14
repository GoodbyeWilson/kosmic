"""The meta-analysis gene search must work on the results table too.

'Find gene' lives under the volcano on the Plots tab. It now follows
you to the Results tab and sits under the table there -- so it has to
select a row, otherwise it reads as a search box that does nothing on
the tab most of the reading happens in.
"""
from __future__ import annotations

import os

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from kosmic.gui.meta_analysis.pages.gene_ma.results import ResultsTab  # noqa: E402
from kosmic.gui.meta_analysis.pages.gene_ma_page import GeneMAPage  # noqa: E402


@pytest.fixture(scope="module")
def app():
    inst = QApplication.instance() or QApplication([])
    yield inst


@pytest.fixture
def results_df():
    return pd.DataFrame({
        "names": ["ACTA2", "COL1A1", "COL1A2", "POSTN"],
        "logfoldchanges": [0.5, 1.2, 0.9, 3.1],
        "se": [0.2, 0.2, 0.2, 0.3],
        "pvals_pooled": [0.2, 0.001, 0.01, 1e-8],
        "fdr": [0.3, 0.004, 0.02, 1e-6],
        "n_studies": [2, 2, 2, 2],
    })


def _selected(tab):
    rows = tab.table.selectionModel().selectedRows()
    if not rows:
        return None
    item = tab.table.item(rows[0].row(), 0)
    return item.text() if item else None


# -- ResultsTab.select_matching ---------------------------------------

def test_select_matching_exact(app, results_df):
    tab = ResultsTab()
    tab.populate(results_df, ["reml"])
    assert tab.select_matching("POSTN") is True
    assert _selected(tab) == "POSTN"


def test_select_matching_is_case_insensitive(app, results_df):
    tab = ResultsTab()
    tab.populate(results_df, ["reml"])
    assert tab.select_matching("postn") is True
    assert _selected(tab) == "POSTN"


def test_select_matching_prefers_prefix_over_substring(app, results_df):
    """'COL1' must land on COL1A1, not on whichever row is highest."""
    tab = ResultsTab()
    tab.populate(results_df, ["reml"])
    assert tab.select_matching("COL1") is True
    assert _selected(tab) == "COL1A1"


def test_select_matching_exact_beats_prefix(app, results_df):
    tab = ResultsTab()
    tab.populate(results_df, ["reml"])
    assert tab.select_matching("COL1A2") is True
    assert _selected(tab) == "COL1A2"


def test_select_matching_no_hit_leaves_selection(app, results_df):
    tab = ResultsTab()
    tab.populate(results_df, ["reml"])
    tab.select_matching("ACTA2")
    assert tab.select_matching("ZZZZ") is False
    assert _selected(tab) == "ACTA2"


def test_select_matching_blank_is_a_no_op(app, results_df):
    tab = ResultsTab()
    tab.populate(results_df, ["reml"])
    assert tab.select_matching("") is False
    assert tab.select_matching("   ") is False


def test_select_matching_on_empty_table(app):
    tab = ResultsTab()
    assert tab.select_matching("POSTN") is False


def test_select_matching_does_not_reemit(app, results_df):
    """Typing must not fire gene_selected -- that would jump tabs."""
    tab = ResultsTab()
    tab.populate(results_df, ["reml"])
    seen = []
    tab.gene_selected.connect(seen.append)
    tab.select_matching("POSTN")
    assert seen == []


# -- the strip's placement on the page --------------------------------

def _results_index(page):
    for i in range(page.tabs.count()):
        if page.tabs.widget(i) is page._results_tab:
            return i
    raise AssertionError("no results tab")


def test_strip_starts_under_the_plots(app):
    page = GeneMAPage()
    assert page._plot_controls.indexOf(page._search_strip) == 0
    assert page._search_strip.parent() is page.tabs.widget(0)


def test_strip_follows_you_to_the_results_tab(app):
    page = GeneMAPage()
    page.tabs.setCurrentIndex(_results_index(page))
    assert page._search_strip.parent() is page._results_tab
    assert page._plot_controls.indexOf(page._search_strip) == -1


def test_strip_comes_back_from_the_results_tab(app):
    """One widget moving means it must survive the round trip."""
    page = GeneMAPage()
    page.tabs.setCurrentIndex(_results_index(page))
    page.tabs.setCurrentIndex(0)
    assert page._plot_controls.indexOf(page._search_strip) == 0
    assert page._gene_search.parent() is page._search_strip


def test_strip_parks_with_the_plots_on_other_tabs(app):
    """Venn has nothing for the search to act on."""
    page = GeneMAPage()
    venn = next(i for i in range(page.tabs.count())
                if page.tabs.tabText(i) == "Venn")
    page.tabs.setCurrentIndex(venn)
    assert page._plot_controls.indexOf(page._search_strip) == 0


def test_search_box_belongs_to_the_strip(app):
    """So moving the strip carries the search box with it."""
    page = GeneMAPage()
    assert page._gene_search.parent() is page._search_strip


# -- the selected-gene readout is a legend on the volcano --------------

def test_readout_floats_on_the_volcano(app):
    """Not in the strip: it describes the plot, so it sits on the plot."""
    page = GeneMAPage()
    assert page._gene_info_label.parent() is page._volcano_plot
    assert page._gene_info_label.objectName() == "volcano_gene_legend"


def test_readout_hidden_until_a_gene_is_picked(app):
    page = GeneMAPage()
    assert page._gene_info_label.text() == ""
    assert page._gene_info_label.isHidden()


def test_readout_shows_and_hides_with_content(app):
    page = GeneMAPage()
    page._set_gene_legend("GABRE<br>log2FC: 1.797")
    assert not page._gene_info_label.isHidden()
    page._set_gene_legend("")
    assert page._gene_info_label.isHidden()


def test_readout_clears_the_axis_and_title(app):
    """(margin, margin) put it across the plot title -- anchor to the vb."""
    page = GeneMAPage()
    page.resize(1200, 800)
    page.show()
    page._set_gene_legend("GABRE<br>log2FC: 1.797 [0.963, 2.631]")
    app.processEvents()
    page._place_gene_legend()
    vb = page._volcano_plot.getPlotItem().vb.geometry()
    pos = page._gene_info_label.pos()
    assert pos.x() >= int(vb.left())
    assert pos.y() >= int(vb.top())
    page.hide()


def test_placement_survives_being_set_off_tab(app):
    """isVisible() is False off-tab; placement must not be skipped."""
    page = GeneMAPage()
    page.resize(1200, 800)
    page.show()
    ridx = _results_index(page)
    page.tabs.setCurrentIndex(ridx)
    page._set_gene_legend("GABRE<br>log2FC: 1.797")
    app.processEvents()
    assert page._gene_info_label.pos().x() > 0
    page.hide()


def test_plot_only_controls_stay_with_the_plots(app):
    page = GeneMAPage()
    plots = page.tabs.widget(0)
    assert page._yaxis_combo.parent() is plots
    assert page._olink_cvd_check.parent() is plots
    assert page._olink_explore_check.parent() is plots


def test_typing_selects_a_results_row(app, results_df):
    page = GeneMAPage()
    page._results_tab.populate(results_df, ["reml"])
    page._gene_search.setText("COL1")
    assert _selected(page._results_tab) == "COL1A1"
