# Highly Variable Genes page (scanpy dispersion plot).
#
# Reads from 'workspace.scrna_ws.current_adata' after HVG selection has run.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class VariableGenesPage(FigurePage):
    TITLE = "Variable Genes"
    UNAVAILABLE_MESSAGE = (
        "Compute highly variable genes in the scRNA workflow to view this figure."
    )

    def _build_controls(self) -> FigureControls:
        return FigureControls(default_w=10, default_h=8)

    def _adata(self):
        scrna = self.workspace.scrna_ws
        if scrna is None:
            return None
        return getattr(scrna, 'current_adata', None)

    def dependencies_met(self) -> bool:
        adata = self._adata()
        if adata is None:
            return False
        return ('highly_variable' in adata.var.columns
                and 'hvg' in adata.uns)

    def _make_render_func(self):
        adata = self._adata()
        if not self.dependencies_met():
            return None
        figsize = self._controls.figsize.get_figsize()

        def _render(adata=adata, sz=figsize):
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import scanpy as sc
            sc.pl.highly_variable_genes(adata, show=False)
            fig = plt.gcf()
            fig.set_size_inches(*sz)
            return fig
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        scrna = self.workspace.scrna_ws
        if scrna is None or not getattr(scrna, 'current_project_dir', None):
            return None
        from kosmic.paths import scrna_figures_dir
        return scrna_figures_dir(scrna.current_project_dir)

    def _default_export_basename(self) -> str:
        return "variable_genes"
