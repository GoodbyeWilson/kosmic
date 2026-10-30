# Highly Variable Genes page: mean against normalised variance per gene.
#
# Reads the study's scRNA dataset (FigurePage._adata) after HVG selection
# has run.

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
        return FigureControls(default_w=7, default_h=5)

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
        var = adata.var
        flavor = adata.uns['hvg'].get('flavor', 'seurat')
        font_sizes = self._font_sizes()

        def _render(var=var, flavor=flavor, sz=figsize, f=font_sizes):
            from kosmic.visualisation.scrna.hvg_plot import create_hvg_plot
            return create_hvg_plot(var, flavor, figsize=sz, font_sizes=f)
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        return self._scrna_export_dir()

    def _default_export_basename(self) -> str:
        return "variable_genes"
