# Pathway activity bubble plot page.
#
# Reads the study's pathway DE results (FigurePage._de; requires an
# Effect_Size column, i.e. pathway-level scoring must have run -- not just
# gene-level DE).

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import QLineEdit

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class _BubbleControls(FigureControls):
    def __init__(self, default_w: int = 8, default_h: int = 6):
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

    STALENESS_ATTR = 'pathway_de_results'

    def dependencies_met(self) -> bool:
        de = self._de()
        return de is not None and de.available('pathway_de_results')

    def _make_render_func(self):
        de = self._de()
        pdr = de.pathway_de_results if de is not None else None
        # bubble_plot.create_bubble_plot accepts either column name.
        if pdr is None or not {'Effect_Size', 'logfoldchanges'} & set(pdr.columns):
            return None
        sd = pdr.copy()
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
                font_sizes=f, figsize=sz or (8, 6),
            )
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        return self._de_export_dir()

    def _default_export_basename(self) -> str:
        return "pathway_bubble"
