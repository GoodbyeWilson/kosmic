# Dataset summary multi-panel figure page.
#
# Reads from 'workspace.de_ws' (DE results + pathway_gene_sets +
# pathway_coverage + h5ad_path + geneset_label).

from __future__ import annotations

from pathlib import Path
from typing import Optional

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class DatasetSummaryPage(FigurePage):
    TITLE = "Dataset Summary"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run DE with selected gene sets to view this figure."
    )

    def _build_controls(self) -> FigureControls:
        return FigureControls(default_w=16, default_h=10)

    def dependencies_met(self) -> bool:
        de = self.workspace.de_ws
        if de is None:
            return False
        return (getattr(de, 'de_results', None) is not None
                and bool(getattr(de, 'pathway_gene_sets', None))
                and bool(getattr(de, 'pathway_coverage', None)))

    def _make_render_func(self):
        if not self.dependencies_met():
            return None
        de = self.workspace.de_ws
        dr = de.de_results.copy()
        gs = dict(de.pathway_gene_sets)
        pc = de.pathway_coverage
        h5ad = getattr(de, 'h5ad_path', None)
        dn = Path(h5ad).stem if h5ad else "Dataset"
        d = de.disease_label
        c = de.control_label
        gl = getattr(de, 'geneset_label', "")
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()

        def _render(dr=dr, gs=gs, pc=pc, dn=dn, d=d, c=c, gl=gl,
                    f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.dataset_summary import (
                compute_dataset_summary, create_summary_figure,
            )
            sdf = compute_dataset_summary(dr, gs, pc)
            return create_summary_figure(
                sdf, dn, d, c, gl, font_sizes=f, figsize=sz,
            )
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(de.project_dir)

    def _default_export_basename(self) -> str:
        de = self.workspace.de_ws
        h5ad = getattr(de, 'h5ad_path', None) if de else None
        stem = Path(h5ad).stem if h5ad else "dataset"
        return f"dataset_summary_{stem}"
