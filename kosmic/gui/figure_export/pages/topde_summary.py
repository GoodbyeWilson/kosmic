# Top DE summary chart page (gene counts per pathway, optionally sig-only).
#
# Reads from 'workspace.de_ws' (DE results + pathway_gene_sets).

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import QCheckBox

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls
from kosmic import DEFAULT_FDR, DEFAULT_LFC_THRESHOLD


class _TopDEControls(FigureControls):
    def __init__(self, default_w: int = 12, default_h: int = 8):
        super().__init__(default_w, default_h)

        self.sig_only = QCheckBox("Significant genes only")
        self.sig_only.setChecked(True)
        self.sig_only.toggled.connect(self.changed)
        self.add_row("", self.sig_only)


class TopDESummaryPage(FigurePage):
    TITLE = "Top DE Summary"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run DE with selected gene sets to view this figure."
    )

    def _build_controls(self) -> _TopDEControls:
        return _TopDEControls()

    def dependencies_met(self) -> bool:
        de = self.workspace.de_ws
        if de is None:
            return False
        if getattr(de, 'de_results', None) is None:
            return False
        return bool(getattr(de, 'pathway_gene_sets', None))

    def _make_render_func(self):
        if not self.dependencies_met():
            return None
        de = self.workspace.de_ws
        dr = de.de_results.copy()
        gs = dict(de.pathway_gene_sets)
        d = de.disease_label
        c = de.control_label
        sig_only = self._controls.sig_only.isChecked()
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()

        def _render(dr=dr, gs=gs, d=d, c=c, so=sig_only,
                    f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.top_de_plots import create_pathway_summary_chart
            if so:
                subset = dr[(dr['pvals_adj'] < DEFAULT_FDR)
                            & (dr['logfoldchanges'].abs() > DEFAULT_LFC_THRESHOLD)]
                label = 'Significant'
            else:
                subset = dr
                label = 'All Genes'
            return create_pathway_summary_chart(
                subset, gs, d, c, label=label, font_sizes=f, figsize=sz,
            )
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(de.project_dir)

    def _default_export_basename(self) -> str:
        suffix = "sig" if self._controls.sig_only.isChecked() else "all"
        return f"top_de_summary_{suffix}"
