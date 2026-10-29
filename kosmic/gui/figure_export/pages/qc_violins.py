# QC violin plot page (genes/cell, total counts, MT%).
#
# Reads the obs of the study's scRNA dataset (FigurePage._adata).

from __future__ import annotations

from pathlib import Path
from typing import Optional

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class QCViolinsPage(FigurePage):
    TITLE = "QC Violins"
    UNAVAILABLE_MESSAGE = (
        "Compute QC metrics in the scRNA workflow to view this figure."
    )

    def _build_controls(self) -> FigureControls:
        return FigureControls(default_w=12, default_h=6)

    def dependencies_met(self) -> bool:
        adata = self._adata()
        if adata is None:
            return False
        return 'n_genes_by_counts' in adata.obs.columns

    def _make_render_func(self):
        adata = self._adata()
        if adata is None or 'n_genes_by_counts' not in adata.obs.columns:
            return None
        obs = adata.obs
        figsize = self._controls.figsize.get_figsize()
        font_sizes = self._font_sizes()

        def _render(obs=obs, sz=figsize, f=font_sizes):
            from kosmic.visualisation.scrna.qc_plots import create_qc_violin_plots
            return create_qc_violin_plots(obs, figsize=sz, dark_mode=False, font_sizes=f)
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        return self._scrna_export_dir()

    def _default_export_basename(self) -> str:
        return "qc_violins"
