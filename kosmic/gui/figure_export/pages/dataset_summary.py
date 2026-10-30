# Dataset summary multi-panel figure page.
#
# Reads the study's DE results (FigurePage._de): the DE results,
# pathway_gene_sets, pathway_coverage, h5ad_path and geneset_label.

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
        ctrl = FigureControls(default_w=16, default_h=10)
        ctrl.add_lfc_gate()
        return ctrl

    def dependencies_met(self) -> bool:
        de = self._de()
        return (de is not None and de.available('de_results')
                and de.available('pathway_gene_sets')
                and de.available('pathway_coverage'))

    def _make_render_func(self):
        de = self._de()
        if (de is None or de.de_results is None or not de.pathway_gene_sets
                or not de.pathway_coverage):
            return None
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
        gate = self._controls.lfc_gate.value()

        def _render(dr=dr, gs=gs, pc=pc, dn=dn, d=d, c=c, gl=gl, gate=gate,
                    f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.dataset_summary import (
                compute_dataset_summary, create_summary_figure,
            )
            sdf = compute_dataset_summary(dr, gs, pc, lfc_threshold=gate)
            return create_summary_figure(
                sdf, dn, d, c, gl, font_sizes=f, figsize=sz,
            )
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        return self._de_export_dir()

    def _default_export_basename(self) -> str:
        de = self._de()
        h5ad = getattr(de, 'h5ad_path', None) if de else None
        stem = Path(h5ad).stem if h5ad else "dataset"
        return f"dataset_summary_{stem}"
