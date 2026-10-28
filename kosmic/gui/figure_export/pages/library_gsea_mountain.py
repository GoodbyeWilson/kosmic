# Library GSEA running-enrichment ("mountain") plot page.
#
# Draws the running-enrichment curve for one gene set from the last library
# GSEA run on the DE Enrichment tab or the Meta-Analysis Enrichment step (see
# _enrichment.py). The curves are held in memory by that run
# ('gsea_library_details'), so this figure is generated in the same session
# as the GSEA -- it is not reloaded from disk.

from __future__ import annotations

from kosmic.gui.figure_export.pages._enrichment import EnrichmentFigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls
from kosmic.gui.shared.theme import NoScrollComboBox


class _MountainControls(FigureControls):
    """Figure size plus a gene-set selector."""

    def __init__(self):
        super().__init__()
        self.term_combo = NoScrollComboBox()
        self.term_combo.currentIndexChanged.connect(self.changed)
        self.add_row("Gene set:", self.term_combo)


class LibraryGseaMountainPage(EnrichmentFigurePage):
    TITLE = "Library GSEA (running enrichment)"
    CHECKS_DE_STALENESS = True
    KIND = 'gsea'
    UNAVAILABLE_MESSAGE = (
        "Run preranked GSEA against a library on the DE or Meta-Analysis "
        "Enrichment step first -- this figure draws the running-enrichment "
        "curve from that run."
    )

    def _build_controls(self) -> FigureControls:
        ctrl = self._add_source_row(_MountainControls())
        # Each source has its own gene sets.
        ctrl.source.currentIndexChanged.connect(self._sync_terms)
        return ctrl

    # -- data access ----------------------------------------------------

    def _has_run(self, src, is_meta: bool) -> bool:
        return bool(getattr(src, 'gsea_library_details', None))

    def _details(self):
        src = self._source()
        return getattr(src, 'gsea_library_details', None) or None

    def _ordered_terms(self):
        """Gene sets with a stored curve, most down-regulated (lowest NES) first
        so the default selection is the top depleted set (e.g. OXPHOS)."""
        det = self._details()
        if not det:
            return []
        return [t for t, _ in sorted(
            det.items(), key=lambda kv: kv[1].get('nes', 0.0))]

    # -- lifecycle ------------------------------------------------------

    def on_activated(self):
        self._sync_sources()
        self._sync_terms()
        self.invalidate()  # curves live in memory and may have changed
        super().on_activated()

    def _sync_terms(self):
        terms = self._ordered_terms()
        combo = self._controls.term_combo
        from kosmic.reference.pathways.libraries import term_label
        current = combo.currentData()
        combo.blockSignals(True)
        combo.clear()
        for t in terms:                      # shown without the prefix; data is the name
            combo.addItem(term_label(t), t)
        if current in terms:
            combo.setCurrentIndex(terms.index(current))
        combo.blockSignals(False)

    # -- render ---------------------------------------------------------

    def _make_render_func(self):
        src = self._source()
        det = self._details()
        terms = self._ordered_terms()
        if not det or not terms:
            return None
        term = self._controls.term_combo.currentData() or terms[0]
        if term not in det:
            term = terms[0]
        d = det[term]

        nes = d.get('nes')
        fdr = None
        summary = getattr(src, 'gsea_library_results', None)
        if summary is not None and len(summary):
            hit = summary[summary['names'].astype(str) == term]
            if len(hit):
                import pandas as pd
                fdr = pd.to_numeric(hit['pvals_adj'], errors='coerce').iloc[0]
        dl = getattr(src, 'disease_label', None)
        cl = getattr(src, 'control_label', None)
        res_curve = d['RES']
        hits = d['hits']
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()
        def _render(rc=res_curve, h=hits, t=term, nes=nes, fdr=fdr, dl=dl, cl=cl,
                    f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.gsea_mountain import (
                create_gsea_mountain_plot,
            )
            return create_gsea_mountain_plot(
                rc, h, t, nes=nes, fdr=fdr, disease_label=dl, control_label=cl,
                figsize=sz, font_sizes=f)
        return _render

    def _default_export_basename(self) -> str:
        prefix = 'meta_' if self._is_meta() else ''
        return f"{prefix}library_gsea_running_enrichment"
