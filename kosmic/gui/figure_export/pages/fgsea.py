# Preranked GSEA (fgsea) enrichment page.
#
# A matplotlib mirror of the DE window's NES bar chart. It never recomputes:
# it renders the in-memory fgsea summary from the last run
# ('de_ws.fgsea_results') and falls back to the results saved to disk by that
# run, so the figure always reflects the latest actual GSEA (against the
# genome-wide gene ranking) rather than a recomputed approximation.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class FgseaPage(FigurePage):
    TITLE = "Pathway Enrichment (GSEA)"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run preranked GSEA (fgsea) on the Pathway DE page first -- this figure "
        "mirrors those enrichment results."
    )

    def _build_controls(self) -> FigureControls:
        return FigureControls()

    def _results_path(self) -> Optional[Path]:
        """On-disk fgsea summary saved by the last run, or None."""
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_pathway_scoring_dir
        from kosmic.de.fgsea import FGSEA_RESULTS_FILE
        return de_pathway_scoring_dir(de.project_dir) / FGSEA_RESULTS_FILE

    def _results(self):
        """fgsea summary to render: the in-memory result from this session's run,
        else the one saved to disk. None if neither exists."""
        de = self.workspace.de_ws
        if de is None:
            return None
        res = getattr(de, 'fgsea_results', None)
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
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()
        from kosmic.gui.shared.theme import get_color
        theme_colors = {'plot_disease': get_color('plot_disease'),
                        'plot_ns': get_color('fg_secondary')}

        def _render(res=res, d=d, c=c, f=font_sizes, sz=figsize, tc=theme_colors):
            from kosmic.visualisation.de.fgsea_plot import create_fgsea_nes_plot
            return create_fgsea_nes_plot(
                res, disease_label=d, control_label=c,
                figsize=sz, font_sizes=f, theme_colors=tc)
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(de.project_dir)

    def _default_export_basename(self) -> str:
        return "pathway_enrichment_gsea"
