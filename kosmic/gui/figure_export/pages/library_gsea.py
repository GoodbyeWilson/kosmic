# Library preranked GSEA (fgsea) enrichment page.
#
# A matplotlib mirror of the Enrichment tab's library NES bar chart
# (Hallmark, WikiPathways, GO, ...), from a study's DE or the Meta-Analysis
# workspace's Enrichment step (see _enrichment.py). It never recomputes: it
# renders the summary of the last library run, held in memory or, for a
# study, saved in its folder.
# The set-count rule matches the in-app plot (all sets for a small library,
# else top per direction) so the figure mirrors what the Enrichment tab shows.

from __future__ import annotations

from kosmic.gui.figure_export.pages._enrichment import EnrichmentFigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class LibraryGseaPage(EnrichmentFigurePage):
    TITLE = "Library Enrichment (GSEA)"
    CHECKS_DE_STALENESS = True
    KIND = 'gsea'
    UNAVAILABLE_MESSAGE = (
        "Run preranked GSEA against a library (Hallmark, WikiPathways, ...) on "
        "the DE or Meta-Analysis Enrichment step first -- this figure mirrors "
        "those results."
    )

    def _build_controls(self) -> FigureControls:
        return FigureControls(default_w=10, default_h=8)

    def _make_render_func(self):
        src = self._source()
        res = getattr(src, 'gsea_library_results', None) if src is not None else None
        if res is None or not len(res):
            return None
        d = getattr(src, 'disease_label', None)
        c = getattr(src, 'control_label', None)
        lib = getattr(src, 'gsea_library_name', None)
        if lib:
            from kosmic.reference.pathways.libraries import library_label
            lib = library_label(lib)
        contrast = self._contrast(src)
        title = f'{lib}: {contrast}' if lib and contrast else None
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()

        def _render(res=res, d=d, c=c, t=title, f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.fgsea_plot import create_fgsea_nes_plot
            return create_fgsea_nes_plot(
                res, disease_label=d, control_label=c, title=t,
                figsize=sz, font_sizes=f)
        return _render

    def _default_export_basename(self) -> str:
        prefix = 'meta_' if self._is_meta() else ''
        return f"{prefix}library_enrichment_gsea"
