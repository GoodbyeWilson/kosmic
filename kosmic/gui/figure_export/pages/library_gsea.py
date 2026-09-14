# Library preranked GSEA (fgsea) enrichment page.
#
# A matplotlib mirror of the DE Enrichment tab's library NES bar chart
# (Hallmark, WikiPathways, GO, ...). It never recomputes: it renders the
# in-memory summary from the last library run ('de_ws.gsea_library_results')
# and falls back to the results saved to disk by that run. The set-count rule
# matches the in-app plot (all sets for a small library, else top per
# direction) so the figure mirrors what the DE window shows.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class LibraryGseaPage(FigurePage):
    TITLE = "Library Enrichment (GSEA)"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run preranked GSEA against a library (Hallmark, WikiPathways, ...) on "
        "the DE Enrichment tab first -- this figure mirrors those results."
    )

    def _build_controls(self) -> FigureControls:
        return FigureControls()

    def _results_path(self) -> Optional[Path]:
        """On-disk library GSEA summary saved by the last run, or None."""
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_pathway_scoring_dir
        from kosmic.de.fgsea import GSEA_LIBRARY_RESULTS_FILE
        return de_pathway_scoring_dir(de.project_dir) / GSEA_LIBRARY_RESULTS_FILE

    def _results(self):
        """Library GSEA summary to render: the in-memory result from this
        session's run, else the one saved to disk. None if neither exists."""
        de = self.workspace.de_ws
        if de is None:
            return None
        res = getattr(de, 'gsea_library_results', None)
        if res is not None and len(res) > 0:
            return res
        from kosmic.de.fgsea import load_fgsea_results
        path = self._results_path()
        return load_fgsea_results(path) if path else None

    def dependencies_met(self) -> bool:
        return self._results() is not None

    def _make_render_func(self):
        res = self._results()
        if res is None:
            return None
        de = self.workspace.de_ws
        d = getattr(de, 'disease_label', None)
        c = getattr(de, 'control_label', None)
        lib = getattr(de, 'gsea_library_name', None)
        title = None
        if lib and d and c:
            title = f'{lib}: {d} vs {c}'
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()
        from kosmic.gui.shared.theme import get_color
        theme_colors = {'plot_disease': get_color('plot_disease'),
                        'plot_ns': get_color('fg_secondary')}

        def _render(res=res, d=d, c=c, t=title, f=font_sizes, sz=figsize,
                    tc=theme_colors):
            from kosmic.visualisation.de.fgsea_plot import create_fgsea_nes_plot
            return create_fgsea_nes_plot(
                res, disease_label=d, control_label=c, title=t,
                figsize=sz, font_sizes=f, theme_colors=tc)
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(de.project_dir)

    def _default_export_basename(self) -> str:
        return "library_enrichment_gsea"
