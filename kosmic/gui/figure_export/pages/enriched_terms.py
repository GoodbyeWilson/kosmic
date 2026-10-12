# Enriched terms page: the top terms of the last enrichment run, as a dot
# plot (fold enrichment, gene count, significance) or a bar plot.
#
# Reads 'enrichment_results' from the DE workspace's Enrichment page or the
# Meta-Analysis workspace's Enrichment step (see _enrichment.py).

from __future__ import annotations

from PyQt6.QtWidgets import QComboBox, QSpinBox

from kosmic.gui.figure_export.pages._enrichment import EnrichmentFigurePage
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


class EnrichedTermsPage(EnrichmentFigurePage):
    TITLE = "Enriched Terms"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run enrichment (GO or an Enrichr library) in the DE or Meta-Analysis "
        "workspace to view this figure."
    )

    def _build_controls(self) -> _EnrichedTermsControls:
        return self._add_source_row(_EnrichedTermsControls())

    def _make_render_func(self):
        src = self._source()
        if src is None:
            return None
        df = src.enrichment_results.copy()
        label = df.attrs.get('significance_label', 'FDR')
        contrast = self._contrast(src)
        title = f'Top enriched terms: {contrast}' if contrast else 'Top enriched terms'
        dot = self._controls.view.currentText() == VIEW_DOT
        top_n = self._controls.top_n.value()
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()

        def _render(df=df, label=label, title=title, dot=dot, top_n=top_n,
                    f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.enrichment import (
                create_enrichment_bar_plot, create_enrichment_dot_plot,
            )
            make = create_enrichment_dot_plot if dot else create_enrichment_bar_plot
            return make(df, top_n=top_n, figsize=sz, font_sizes=f, title=title,
                        significance_label=label)
        return _render

    def _default_export_basename(self) -> str:
        prefix = 'meta_' if self._is_meta() else ''
        return f"{prefix}enriched_terms_{self._controls.view.currentText().lower()}"
