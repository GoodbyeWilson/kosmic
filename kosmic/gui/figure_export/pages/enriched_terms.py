# Enriched terms page: the top terms of the last enrichment run, as a dot
# plot (fold enrichment, gene count, significance) or a bar plot.
#
# Reads from 'workspace.de_ws.enrichment_results' (populated after running
# discovery-mode enrichment in the DE workspace).

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import QComboBox, QSpinBox

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls

VIEW_DOT, VIEW_BAR = "Dot", "Bar"


class _EnrichedTermsControls(FigureControls):
    def __init__(self):
        super().__init__()
        self.view = QComboBox()
        self.view.addItems([VIEW_DOT, VIEW_BAR])
        self.view.currentIndexChanged.connect(self.changed)
        self.add_row("View:", self.view)

        self.top_n = QSpinBox()
        self.top_n.setRange(3, 50)
        self.top_n.setValue(20)
        self.top_n.valueChanged.connect(self.changed)
        self.add_row("Terms:", self.top_n)


class EnrichedTermsPage(FigurePage):
    TITLE = "Enriched Terms"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run pathway enrichment in the DE workflow to view this figure."
    )

    def _build_controls(self) -> _EnrichedTermsControls:
        return _EnrichedTermsControls()

    def dependencies_met(self) -> bool:
        de = self.workspace.de_ws
        if de is None:
            return False
        enr = getattr(de, 'enrichment_results', None)
        return enr is not None and not enr.empty

    def _make_render_func(self):
        if not self.dependencies_met():
            return None
        df = self.workspace.de_ws.enrichment_results.copy()
        label = df.attrs.get('significance_label', 'FDR')
        dot = self._controls.view.currentText() == VIEW_DOT
        top_n = self._controls.top_n.value()
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()

        def _render(df=df, label=label, dot=dot, top_n=top_n, f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.enrichment import (
                create_enrichment_bar_plot, create_enrichment_dot_plot,
            )
            make = create_enrichment_dot_plot if dot else create_enrichment_bar_plot
            return make(df, top_n=top_n, figsize=sz, font_sizes=f,
                        significance_label=label)
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(de.project_dir)

    def _default_export_basename(self) -> str:
        return f"enriched_terms_{self._controls.view.currentText().lower()}"
