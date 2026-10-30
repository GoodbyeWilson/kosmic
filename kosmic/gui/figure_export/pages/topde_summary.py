# DE genes by pathway: up / down DE gene counts per gene set (optionally
# significant genes only), with each gene set's coverage when known.
#
# Reads the study's DE results (FigurePage._de): the DE results,
# pathway_gene_sets and pathway_coverage.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import QCheckBox

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls
from kosmic import DEFAULT_FDR


class _TopDEControls(FigureControls):
    def __init__(self, default_w: int = 12, default_h: int = 8):
        super().__init__(default_w, default_h)

        self.sig_only = QCheckBox("Significant genes only")
        self.sig_only.setChecked(True)
        self.sig_only.toggled.connect(self.changed)
        self.add_row("", self.sig_only)
        self.add_lfc_gate()


class TopDESummaryPage(FigurePage):
    TITLE = "DE Genes by Pathway"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run DE with selected gene sets to view this figure."
    )

    def _build_controls(self) -> _TopDEControls:
        return _TopDEControls()

    def dependencies_met(self) -> bool:
        de = self._de()
        return (de is not None and de.available('de_results')
                and de.available('pathway_gene_sets'))

    def _make_render_func(self):
        de = self._de()
        if de is None or de.de_results is None or not de.pathway_gene_sets:
            return None
        dr = de.de_results.copy()
        gs = dict(de.pathway_gene_sets)
        cov = getattr(de, 'pathway_coverage', None) or None
        d = de.disease_label
        c = de.control_label
        sig_only = self._controls.sig_only.isChecked()
        gate = self._controls.lfc_gate.value()
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()

        def _render(dr=dr, gs=gs, cov=cov, d=d, c=c, so=sig_only, gate=gate,
                    f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.top_de_plots import create_pathway_summary_chart
            if so:
                subset = dr[(dr['pvals_adj'] < DEFAULT_FDR)
                            & (dr['logfoldchanges'].abs() > gate)]
                label = 'Significant'
            else:
                subset = dr
                label = 'All Genes'
            return create_pathway_summary_chart(
                subset, gs, d, c, label=label, font_sizes=f, figsize=sz,
                pathway_coverage=cov,
            )
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        return self._de_export_dir()

    def _default_export_basename(self) -> str:
        suffix = "sig" if self._controls.sig_only.isChecked() else "all"
        return f"top_de_summary_{suffix}"
