# Pathway activity bubble plot page.
#
# Reads from 'workspace.de_ws.pathway_de_results' (requires Effect_Size
# column, i.e. pathway-level scoring must have run -- not just gene-level DE).

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import QLineEdit

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class _BubbleControls(FigureControls):
    def __init__(self, default_w: int = 12, default_h: int = 10):
        super().__init__(default_w, default_h)

        self.title_edit = QLineEdit()
        self.title_edit.setPlaceholderText("Auto from gene set name")
        self.title_edit.textChanged.connect(self.changed)
        self.add_row("Title:", self.title_edit)


class PathwayBubblePage(FigurePage):
    TITLE = "Pathway Bubble"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run pathway-level DE scoring to view this figure."
    )

    def _build_controls(self) -> _BubbleControls:
        return _BubbleControls()

    def dependencies_met(self) -> bool:
        de = self.workspace.de_ws
        if de is None:
            return False
        pdr = getattr(de, 'pathway_de_results', None)
        if pdr is None:
            return False
        # bubble_plot.create_bubble_plot accepts either column name.
        return ('Effect_Size' in pdr.columns
                or 'logfoldchanges' in pdr.columns)

    def _make_render_func(self):
        if not self.dependencies_met():
            return None
        de = self.workspace.de_ws
        sd = de.pathway_de_results.copy()
        d = de.disease_label
        c = de.control_label
        title_override = self._controls.title_edit.text().strip()
        gl = getattr(de, 'geneset_label', "") or ""
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()

        def _render(sd=sd, d=d, c=c, ct=title_override, gl=gl,
                    f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.bubble_plot import create_bubble_plot
            if ct:
                title = ct
            elif 'Metabolic' in gl:
                title = f'Metabolic Pathway Activity: {d} vs {c}'
            else:
                title = f'Gene Set Activity: {d} vs {c}'
            y_label = 'Metabolic Pathways' if 'Metabolic' in gl else 'Gene Sets'
            return create_bubble_plot(
                sd, d, c, plot_title=title, y_label=y_label,
                font_sizes=f, figsize=sz or (12, 10),
            )
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(de.project_dir)

    def _default_export_basename(self) -> str:
        return "pathway_bubble"
