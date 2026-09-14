# Library GSEA running-enrichment ("mountain") plot page.
#
# Draws the running-enrichment curve for one gene set from the last library
# GSEA run on the DE Enrichment tab. The curves are held in memory by that
# run ('de_ws.gsea_library_details'), so this figure is generated in the same
# session as the GSEA -- it is not reloaded from disk.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls
from kosmic.gui.shared.theme import NoScrollComboBox


class _MountainControls(FigureControls):
    """Figure size plus a gene-set selector."""

    def __init__(self):
        super().__init__()
        self.term_combo = NoScrollComboBox()
        self.term_combo.currentIndexChanged.connect(self.changed)
        self.add_row("Gene set:", self.term_combo)


class LibraryGseaMountainPage(FigurePage):
    TITLE = "Library GSEA (running enrichment)"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run preranked GSEA against a library on the DE Enrichment tab first -- "
        "this figure draws the running-enrichment curve from that run."
    )

    def _build_controls(self) -> FigureControls:
        return _MountainControls()

    # -- data access ----------------------------------------------------

    def _details(self):
        de = self.workspace.de_ws
        if de is None:
            return None
        return getattr(de, 'gsea_library_details', None) or None

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
        self._sync_terms()
        self.invalidate()  # curves live in memory and may have changed
        super().on_activated()

    def _sync_terms(self):
        terms = self._ordered_terms()
        combo = self._controls.term_combo
        current = combo.currentText()
        combo.blockSignals(True)
        combo.clear()
        combo.addItems(terms)
        if current in terms:
            combo.setCurrentText(current)
        combo.blockSignals(False)

    def dependencies_met(self) -> bool:
        return bool(self._ordered_terms())

    # -- render ---------------------------------------------------------

    def _make_render_func(self):
        det = self._details()
        terms = self._ordered_terms()
        if not det or not terms:
            return None
        term = self._controls.term_combo.currentText() or terms[0]
        if term not in det:
            term = terms[0]
        d = det[term]

        de = self.workspace.de_ws
        nes = d.get('nes')
        fdr = None
        summary = getattr(de, 'gsea_library_results', None)
        if summary is not None and len(summary):
            hit = summary[summary['names'].astype(str) == term]
            if len(hit):
                import pandas as pd
                fdr = pd.to_numeric(hit['pvals_adj'], errors='coerce').iloc[0]
        dl = getattr(de, 'disease_label', None)
        cl = getattr(de, 'control_label', None)
        res_curve = d['RES']
        hits = d['hits']
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()
        from kosmic.gui.shared.theme import get_color
        theme_colors = {'plot_disease': get_color('plot_disease'),
                        'fg': get_color('plot_fg')}

        def _render(rc=res_curve, h=hits, t=term, nes=nes, fdr=fdr, dl=dl, cl=cl,
                    f=font_sizes, sz=figsize, tc=theme_colors):
            from kosmic.visualisation.de.gsea_mountain import (
                create_gsea_mountain_plot,
            )
            return create_gsea_mountain_plot(
                rc, h, t, nes=nes, fdr=fdr, disease_label=dl, control_label=cl,
                figsize=sz, font_sizes=f, theme_colors=tc)
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(de.project_dir)

    def _default_export_basename(self) -> str:
        return "library_gsea_running_enrichment"
