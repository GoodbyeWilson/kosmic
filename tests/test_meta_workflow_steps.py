"""Discovery gained the Enrichment and Validation steps its mode card
always advertised, and enrichment is now the DE page rather than a
cut-down copy of it.

The mode chooser has always described Discovery as four sub-steps
(pooling, results, enrichment, validation) while STEPS_EXPLORATORY
implemented one. These tests pin the flow to what is advertised, and
pin the shared-code reuse that stopped the two enrichment UIs drifting.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from kosmic.gui.de_analysis.pages.enrichment_page import EnrichmentPage  # noqa: E402
from kosmic.gui.meta_analysis.pages.enrichment_step import (  # noqa: E402
    MetaEnrichmentPage,
)
from kosmic.gui.meta_analysis.pages.validation_page import ValidationPage  # noqa: E402
from kosmic.gui.meta_analysis.workspace import MetaAnalysisWorkspace  # noqa: E402


@pytest.fixture(scope="module")
def app():
    yield QApplication.instance() or QApplication([])


@pytest.fixture
def meta_df():
    rng = np.random.default_rng(0)
    n = 200
    return pd.DataFrame({
        "names": [f"G{i}" for i in range(n)],
        "logfoldchanges": rng.normal(0, 1.2, n),
        "se": rng.uniform(0.1, 0.5, n),
        "z_stat": rng.normal(0, 2, n),
        "pvals_pooled": rng.uniform(0, 1, n),
        "fdr": rng.uniform(0, 1, n),
        "n_studies": 2,
    })


# -- the advertised flow ----------------------------------------------

def test_discovery_has_the_advertised_steps(app):
    steps = [t for t, _ in MetaAnalysisWorkspace.steps_for_mode("exploratory")]
    assert steps == ["Select Studies", "Choose Mode",
                     "Discovery Meta-Analysis", "Enrichment", "Validation",
                     "Methods"]


def test_step_map_matches_step_list(app):
    steps = MetaAnalysisWorkspace.steps_for_mode("exploratory")
    assert len(steps) == len(MetaAnalysisWorkspace._EXPLORATORY_MAP)


def test_page_name_matches_the_nav_label(app):
    """The nav said Discovery Meta-Analysis; the log said Consensus."""
    ws = MetaAnalysisWorkspace()
    assert ws._stack_page_name(ws.PAGE_GENE_MA) == "Discovery Meta-Analysis"
    assert ws._stack_page_name(ws.PAGE_ENRICHMENT) == "Enrichment"
    assert ws._stack_page_name(ws.PAGE_VALIDATION) == "Validation"


def test_existing_page_indices_are_unchanged(app):
    """New pages were appended; saved projects reference the old ones."""
    assert MetaAnalysisWorkspace.PAGE_SELECT == 0
    assert MetaAnalysisWorkspace.PAGE_GENE_MA == 4
    assert MetaAnalysisWorkspace.PAGE_METHODS_RES == 6


# -- enrichment is the DE page, not a copy ----------------------------

@pytest.fixture
def enrich(app):
    """Hold the wrapper: a parentless QWidget with no Python reference is
    collected, and Qt then deletes the children under it."""
    page = MetaEnrichmentPage()
    yield page


def test_meta_enrichment_hosts_the_de_page(enrich):
    assert isinstance(enrich, EnrichmentPage)


def test_meta_enrichment_offers_both_methods(enrich):
    """The old meta card was ORA-only; GSEA was the substantive gap."""
    ep = enrich
    methods = [ep._method_combo.itemData(i)
               for i in range(ep._method_combo.count())]
    assert set(methods) == {"ora", "gsea"}


def test_meta_enrichment_has_the_full_library_list(enrich):
    """Seven hardcoded Enrichr libraries became the DE page's eleven."""
    assert enrich._library_combo.count() == 11


def test_meta_ranks_gsea_by_pooled_z(enrich):
    assert enrich.rank_column == "z_stat"


def test_de_still_ranks_gsea_by_log2fc(app):
    """The seam must not change the DE page's own behaviour."""
    assert EnrichmentPage.rank_column == "logfoldchanges"


def test_pooled_columns_are_translated(app, meta_df):
    """Pooling calls it fdr; the DE page reads pvals_adj."""
    p = MetaEnrichmentPage()
    p.set_meta_results(meta_df)
    de = p.ws.de_results
    assert "pvals_adj" in de.columns
    assert (de["pvals_adj"] == meta_df["fdr"]).all()
    assert p.ws.gene_de_genome_wide is not None


def test_meta_enrichment_clears(app, meta_df):
    p = MetaEnrichmentPage()
    p.set_meta_results(meta_df)
    p.set_meta_results(None)
    assert p.ws.de_results is None


# -- direction filter (the meta-only feature, now shared) -------------

def test_direction_filter_splits_up_and_down(app, meta_df):
    p = MetaEnrichmentPage()
    p.set_meta_results(meta_df)
    ep = p
    counts = {}
    for i in range(ep._direction_combo.count()):
        ep._direction_combo.setCurrentIndex(i)
        study, _, _ = ep._get_sig_genes()
        lfc = meta_df.set_index("names").loc[study, "logfoldchanges"]
        counts[ep._direction_combo.currentData()] = (
            int((lfc > 0).sum()), int((lfc < 0).sum()))
    assert counts["up"][1] == 0
    assert counts["down"][0] == 0
    assert counts["both"] == (counts["up"][0], counts["down"][1])


def test_direction_filter_leaves_background_alone(app, meta_df):
    """The ORA denominator must stay the whole tested universe."""
    p = MetaEnrichmentPage()
    p.set_meta_results(meta_df)
    ep = p
    backgrounds = []
    for i in range(ep._direction_combo.count()):
        ep._direction_combo.setCurrentIndex(i)
        _, bg, _ = ep._get_sig_genes()
        backgrounds.append(len(bg))
    assert len(set(backgrounds)) == 1


# -- validation step ---------------------------------------------------

def test_validation_has_all_four_routes(app):
    p = ValidationPage()
    labels = [p.tabs.tabText(i) for i in range(p.tabs.count())]
    assert labels == ["Cross-Dataset Reproducibility", "LOO Validation",
                      "Olink Panels", "GWAS Overlap"]


def test_repro_needs_two_studies_loo_needs_a_run(app, meta_df):
    p = ValidationPage()
    assert not p._repro_btn.isEnabled()
    assert not p._loo_btn.isEnabled()

    p.set_studies([{"name": "S1", "df": meta_df},
                   {"name": "S2", "df": meta_df}])
    assert p._repro_btn.isEnabled()
    assert not p._loo_btn.isEnabled()
    assert "pooling run" in p._status_label.text()


def test_loo_guard_counts_pooled_studies_not_loaded_ones(app, meta_df):
    """LOO folds over the studies the run used, not everything loaded."""
    p = ValidationPage()
    p.set_studies([{"name": "S1", "df": meta_df},
                   {"name": "S2", "df": meta_df}])
    p.set_run_context({"meta_df": meta_df, "method_keys": ["dl"],
                       "params": {"calibration": "analytical"},
                       "de_dfs": [meta_df] * 3,
                       "labels": ["S1", "S2", "S3"]})
    assert p._loo_btn.isEnabled()


def test_loo_result_goes_back_to_the_owning_page(app):
    """GeneMAPage keeps the per-mode cache; Validation hands it back."""
    ws = MetaAnalysisWorkspace()
    ws._on_loo_complete({"avg_replication_pct": 50.0})
    assert ws._gene_ma_page._loo_result is not None


def test_pooling_run_feeds_both_downstream_steps(app, meta_df):
    ws = MetaAnalysisWorkspace()
    g = ws._gene_ma_page
    g._meta_df = meta_df
    g._method_keys = ["dl"]
    g._last_consensus_params = {"calibration": "analytical",
                                "method_keys": ["dl"]}
    g._last_consensus_de_dfs = [meta_df] * 3
    g._last_consensus_labels = ["S1", "S2", "S3"]
    ws._on_gene_ma_analysis_complete()
    assert ws._enrichment_page.ws.de_results is not None
    assert ws._validation_page._loo_btn.isEnabled()


def test_moved_tabs_are_gone_from_the_pooling_page(app):
    ws = MetaAnalysisWorkspace()
    labels = [ws._gene_ma_page.tabs.tabText(i)
              for i in range(ws._gene_ma_page.tabs.count())]
    for gone in ("Enrichment", "Cross-Dataset Reproducibility",
                 "LOO Validation"):
        assert gone not in labels


# -- enrichment sidebar (shared by DE and meta) -----------------------

def test_run_and_export_sit_outside_the_accordion(enrich):
    """Run was inside the Settings group, so collapsing it hid Run."""
    from PyQt6.QtWidgets import QPushButton
    from kosmic.gui.shared.widgets import StageSummaryCard

    in_cards = set()
    for card in enrich.findChildren(StageSummaryCard):
        in_cards.update(id(b) for b in card.findChildren(QPushButton))
    outside = [b.text() for b in enrich.sidebar.findChildren(QPushButton)
               if id(b) not in in_cards]
    assert "Run Enrichment" in outside
    assert "Export CSV..." in outside


def test_method_library_and_algorithm_share_one_card(enrich):
    """Picking GSEA rewrites the library list, so they are one decision."""
    from kosmic.gui.shared.widgets import StageSummaryCard

    holder = None
    for card in enrich.findChildren(StageSummaryCard):
        if enrich._method_combo in card.findChildren(type(enrich._method_combo)):
            holder = card
    assert holder is not None
    widgets = holder.findChildren(type(enrich._method_combo))
    assert enrich._library_combo in widgets
    assert enrich._algo_combo in widgets


def test_method_still_constrains_the_library_list(enrich):
    enrich._method_combo.setCurrentIndex(0)          # ORA
    ora_n = enrich._library_combo.count()
    enrich._method_combo.setCurrentIndex(1)          # GSEA
    gsea_n = enrich._library_combo.count()
    assert ora_n == 11 and gsea_n == 8
    assert not enrich._fdr_spin.isEnabled()          # GSEA is threshold-free
    enrich._method_combo.setCurrentIndex(0)
    assert enrich._fdr_spin.isEnabled()


def _summary(card):
    box = card._summary_box
    return [box.itemAt(i).widget().text() for i in range(box.count())
            if box.itemAt(i).widget() is not None]


def test_every_card_summarises_itself(enrich):
    """A collapsed card with no summary is worse than no card."""
    from kosmic.gui.shared.widgets import StageSummaryCard

    for card in enrich.findChildren(StageSummaryCard):
        lines = _summary(card)
        assert lines and any(ln.strip() for ln in lines)


def test_algorithm_shows_without_the_page_being_visible(enrich):
    """isVisible() is False off-screen, which hid it from the summary."""
    enrich._method_combo.setCurrentIndex(0)
    for i in range(enrich._library_combo.count()):
        if enrich._library_combo.itemData(i) == "GO_BP":
            enrich._library_combo.setCurrentIndex(i)
            break
    assert "elim" in " ".join(_summary(enrich._terms_card))


def test_meta_context_carries_the_project_dir(enrich, tmp_path):
    """GSEA persists results via ws.project_dir; it was never set."""
    enrich.set_project_folder(tmp_path)
    assert enrich.ws.project_dir == tmp_path


# -- nothing may exceed the card's width ------------------------------

def test_hint_labels_wrap_by_default(app):
    """A non-wrapping QLabel reports its whole string as its minimum."""
    from kosmic.gui.shared.widgets import HintLabel

    label = HintLabel("A fairly long sentence of inline help text that "
                      "would otherwise demand a very wide container.")
    assert label.wordWrap()


def test_no_control_demands_a_wildly_oversized_card(enrich):
    """Catch the pathology, not font-metric noise.

    A HintLabel with word wrap off reported its entire sentence as its
    minimum width -- 1,560px inside a 274px card, 5.7x over -- and
    dragged every sibling control out of the container.

    The bound is deliberately loose (2x) because the headless Qt font is
    roughly twice the width of the real UI font, so exact widths here
    mean nothing. Anything approaching 2x is a widget that cannot wrap,
    which is the bug; a few percent over is measurement noise.
    """
    from PyQt6.QtWidgets import QWidget

    enrich.resize(1400, 900)
    enrich.show()
    try:
        for key in ("genes", "method", "terms"):
            enrich._accordion.toggle(key)
            enrich.grab()  # force a real layout pass
            card = getattr(enrich, f"_{key}_card")
            budget = card._settings_scroll.viewport().width()
            for w in card._settings_scroll.widget().findChildren(QWidget):
                need = w.minimumSizeHint().width()
                assert need < budget * 2, (
                    f"{key}: {w.__class__.__name__} wants {need}px "
                    f"in a {budget}px card")
    finally:
        enrich.hide()
