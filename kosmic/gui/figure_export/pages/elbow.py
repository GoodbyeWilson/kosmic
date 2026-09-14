# PCA elbow plot page.
#
# Reads variance ratios from 'workspace.scrna_ws.current_adata.uns['pca']'.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class PCAElbowPage(FigurePage):
    TITLE = "PCA Elbow"
    UNAVAILABLE_MESSAGE = (
        "Run PCA in the scRNA workflow to view this figure."
    )

    def _build_controls(self) -> FigureControls:
        return FigureControls(default_w=8, default_h=6)

    def _adata(self):
        scrna = self.workspace.scrna_ws
        if scrna is None:
            return None
        return getattr(scrna, 'current_adata', None)

    def dependencies_met(self) -> bool:
        adata = self._adata()
        if adata is None:
            return False
        pca = adata.uns.get('pca') if hasattr(adata, 'uns') else None
        return bool(pca and 'variance_ratio' in pca)

    def _make_render_func(self):
        adata = self._adata()
        if not self.dependencies_met():
            return None

        import numpy as np
        variance_ratio = np.asarray(
            adata.uns['pca']['variance_ratio'], dtype=float)

        def _render(variance_ratio=variance_ratio):
            from kosmic.visualisation.scrna.pca_plots import (
                find_elbow, create_elbow_plot,
            )
            suggested = find_elbow(variance_ratio)
            return create_elbow_plot(variance_ratio, suggested_pcs=suggested)
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        scrna = self.workspace.scrna_ws
        if scrna is None or not getattr(scrna, 'current_project_dir', None):
            return None
        from kosmic.paths import scrna_figures_dir
        return scrna_figures_dir(scrna.current_project_dir)

    def _default_export_basename(self) -> str:
        return "pca_elbow"
