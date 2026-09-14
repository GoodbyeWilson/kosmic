# Enrichment bar plot page (top -log10(P) for enriched terms).
#
# Reads from 'workspace.de_ws.enrichment_results' (populated after running
# discovery-mode enrichment in the DE workspace). Plot styling reads theme
# colours.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class EnrichmentBarPage(FigurePage):
    TITLE = "Enrichment Bar"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run pathway enrichment in the DE workflow to view this figure."
    )

    def _build_controls(self) -> FigureControls:
        return FigureControls(default_w=10, default_h=8)

    def dependencies_met(self) -> bool:
        de = self.workspace.de_ws
        if de is None:
            return False
        enr = getattr(de, 'enrichment_results', None)
        return enr is not None and not enr.empty

    def _make_render_func(self):
        if not self.dependencies_met():
            return None
        de = self.workspace.de_ws
        df = de.enrichment_results.copy()
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()

        from kosmic.gui.shared.theme import get_color
        colors = dict(
            bg_color=get_color('plot_bg'),
            fg_color=get_color('plot_fg'),
            bar_color=get_color('plot_control'),
            sig_line_color=get_color('error'),
        )

        def _render(df=df, f=font_sizes, sz=figsize, colors=colors):
            from kosmic.visualisation.de.enrichment import create_enrichment_bar_plot
            return create_enrichment_bar_plot(
                df, figsize=sz, font_sizes=f, **colors,
            )
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(de.project_dir)

    def _default_export_basename(self) -> str:
        return "enrichment_bar"
