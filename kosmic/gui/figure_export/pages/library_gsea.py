# Library preranked GSEA (fgsea) enrichment page.
#
# A matplotlib mirror of the Enrichment tab's library NES bar chart
# (Hallmark, WikiPathways, GO, ...), from the DE workspace or the
# Meta-Analysis workspace's Enrichment step (see _enrichment.py). It never
# recomputes: it renders the in-memory summary of the last library run and,
# for the DE workspace, falls back to the results saved to disk by that run.
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
        return self._add_source_row(FigureControls())

    @staticmethod
    def _saved_results(de):
        """The DE workspace's library GSEA summary saved by its last run, or
        None."""
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_pathway_scoring_dir
        from kosmic.de.fgsea import GSEA_LIBRARY_RESULTS_FILE, load_fgsea_results
        return load_fgsea_results(
            de_pathway_scoring_dir(de.project_dir) / GSEA_LIBRARY_RESULTS_FILE)

    def _results(self, src, is_meta: bool):
        """The summary to render: the in-memory result, else (DE only) the
        one saved to disk. None if neither exists."""
        res = getattr(src, 'gsea_library_results', None)
        if res is not None and len(res) > 0:
            return res
        return None if is_meta else self._saved_results(src)

    def _has_run(self, src, is_meta: bool) -> bool:
        return src is not None and self._results(src, is_meta) is not None

    def _make_render_func(self):
        src = self._source()
        if src is None:
            return None
        res = self._results(src, self._is_meta())
        d = getattr(src, 'disease_label', None)
        c = getattr(src, 'control_label', None)
        lib = getattr(src, 'gsea_library_name', None)
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
