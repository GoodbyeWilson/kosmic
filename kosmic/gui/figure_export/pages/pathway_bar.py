# Pathway coverage bar plot page (significant up/down per pathway).
#
# Reads from 'workspace.de_ws' (DE results + pathway_coverage).

from __future__ import annotations

from pathlib import Path
from typing import Optional

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class PathwayBarPage(FigurePage):
    TITLE = "Pathway Coverage Bar"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run DE with selected gene sets to view this figure."
    )

    def _build_controls(self) -> FigureControls:
        return FigureControls(default_w=14, default_h=8)

    def dependencies_met(self) -> bool:
        de = self.workspace.de_ws
        if de is None:
            return False
        if getattr(de, 'de_results', None) is None:
            return False
        return bool(getattr(de, 'pathway_coverage', None))

    def _make_render_func(self):
        if not self.dependencies_met():
            return None
        de = self.workspace.de_ws
        dr = de.de_results.copy()
        pc = de.pathway_coverage
        d = de.disease_label
        c = de.control_label
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()

        def _render(dr=dr, pc=pc, d=d, c=c, f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.pathway_plots import create_pathway_bar_plot
            return create_pathway_bar_plot(dr, pc, d, c, font_sizes=f, figsize=sz)
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(de.project_dir)

    def _default_export_basename(self) -> str:
        return "pathway_coverage_bar"
