# Enriched term links page: the top enriched terms and the DE genes behind
# them, as a chord diagram or a gene x term grid.
#
# Reads 'workspace.de_ws.enrichment_results' (the last enrichment run in the
# DE workspace: GO or an Enrichr library) and 'de_ws.de_results' for the
# genes' log2 fold changes.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import QCheckBox, QComboBox, QSpinBox

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls

VIEW_CHORD, VIEW_GRID = "Chord", "Grid"
ORDER_LFC, ORDER_TERM = "By fold change", "By term"


class _TermLinksControls(FigureControls):
    def __init__(self):
        super().__init__()

        self.view = QComboBox()
        self.view.addItems([VIEW_CHORD, VIEW_GRID])
        self.view.currentIndexChanged.connect(self.changed)
        self.add_row("View:", self.view)

        self.order = QComboBox()
        self.order.addItems([ORDER_LFC, ORDER_TERM])
        self.order.currentIndexChanged.connect(self.changed)
        self.add_row("Order genes:", self.order)

        self.n_terms = QSpinBox()
        self.n_terms.setRange(2, 12)
        self.n_terms.setValue(8)
        self.n_terms.valueChanged.connect(self.changed)
        self.add_row("Terms:", self.n_terms)

        self.per_term = QSpinBox()
        self.per_term.setRange(1, 20)
        self.per_term.setValue(8)
        self.per_term.setToolTip(
            "Genes shown per term: those with the largest |log2 fold change|.")
        self.per_term.valueChanged.connect(self.changed)
        self.add_row("Genes per term:", self.per_term)

        self.exclude_mt = QCheckBox("Exclude mitochondrial genes (MT-)")
        self.exclude_mt.setToolTip(
            "Leave MT- genes out of the figure. In single-nucleus data their\n"
            "fold changes are largely technical (ambient or carry-over RNA).\n"
            "Only the genes drawn change; the enrichment result does not.")
        self.exclude_mt.toggled.connect(self.changed)
        self.add_row("", self.exclude_mt)


class TermLinksPage(FigurePage):
    TITLE = "Enriched Term Links"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run enrichment (GO or an Enrichr library) in the DE workflow to view "
        "this figure."
    )

    def _build_controls(self) -> _TermLinksControls:
        return _TermLinksControls()

    def dependencies_met(self) -> bool:
        de = self.workspace.de_ws
        if de is None:
            return False
        enr = getattr(de, 'enrichment_results', None)
        dr = getattr(de, 'de_results', None)
        return (enr is not None and not enr.empty and 'Genes' in enr.columns
                and dr is not None and not dr.empty)

    def _make_render_func(self):
        if not self.dependencies_met():
            return None
        de = self.workspace.de_ws
        enr = de.enrichment_results.copy()
        dr = de.de_results[['names', 'logfoldchanges']].copy()
        label = enr.attrs.get('significance_label', 'FDR')
        d, c = de.disease_label, de.control_label
        ctrl = self._controls
        chord = ctrl.view.currentText() == VIEW_CHORD
        order = 'lfc' if ctrl.order.currentText() == ORDER_LFC else 'term'
        n_terms, per_term = ctrl.n_terms.value(), ctrl.per_term.value()
        exclude = ('mitochondrial',) if ctrl.exclude_mt.isChecked() else ()
        figsize = ctrl.figsize.get_figsize()
        font_sizes = self._font_sizes()

        def _render(enr=enr, dr=dr, label=label, d=d, c=c, chord=chord, order=order,
                    n_terms=n_terms, per_term=per_term, exclude=exclude,
                    sz=figsize, f=font_sizes):
            from kosmic.visualisation.de.term_links import (
                create_term_chord, create_term_grid, select_term_genes,
            )
            members, lfc, n_passing = select_term_genes(
                enr, dr, n_terms, per_term, exclude_groups=exclude)
            title = f'Top enriched terms: {d} vs {c}'
            if not n_passing:
                title += f'\n(no term at {label} < 0.05)'
            make = create_term_chord if chord else create_term_grid
            return make(members, lfc, order=order, title=title, figsize=sz,
                        font_sizes=f)
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(de.project_dir)

    def _default_export_basename(self) -> str:
        view = self._controls.view.currentText().lower()
        return f"enriched_term_links_{view}"
