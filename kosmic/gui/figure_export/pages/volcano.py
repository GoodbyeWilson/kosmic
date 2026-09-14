# Volcano plot page: gene-level DE results.
#
# Reads from 'workspace.de_ws' (DEWorkspace). Disabled until DE has run.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import QComboBox, QDoubleSpinBox, QSpinBox

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class _VolcanoControls(FigureControls):
    def __init__(self, default_w: int = 8, default_h: int = 8):
        super().__init__(default_w, default_h)

        self.max_labels = QSpinBox()
        self.max_labels.setRange(0, 50)
        self.max_labels.setValue(10)
        self.max_labels.valueChanged.connect(self.changed)
        self.add_row("Max labels:", self.max_labels)

        self.point_size = QSpinBox()
        self.point_size.setRange(1, 30)
        self.point_size.setValue(9)
        self.point_size.valueChanged.connect(self.changed)
        self.add_row("Point size:", self.point_size)

        self.color_scheme = QComboBox()
        self.color_scheme.addItems(
            ["Red / Blue", "Orange / Teal", "Magenta / Cyan"])
        self.color_scheme.currentIndexChanged.connect(self.changed)
        self.add_row("Colors:", self.color_scheme)

        self.pval_type = QComboBox()
        self.pval_type.addItems(["Adjusted (FDR)", "Nominal"])
        self.pval_type.currentIndexChanged.connect(self.changed)
        self.add_row("Y-axis p-value:", self.pval_type)

        # 0 = automatic range.
        self.x_limit = QDoubleSpinBox()
        self.x_limit.setRange(0.0, 50.0)
        self.x_limit.setSingleStep(0.5)
        self.x_limit.setDecimals(1)
        self.x_limit.setValue(0.0)
        self.x_limit.setSpecialValueText("Auto")
        self.x_limit.valueChanged.connect(self.changed)
        self.add_row("X limit (+/- log2FC):", self.x_limit)

        self.y_limit = QDoubleSpinBox()
        self.y_limit.setRange(0.0, 1000.0)
        self.y_limit.setSingleStep(1.0)
        self.y_limit.setDecimals(1)
        self.y_limit.setValue(0.0)
        self.y_limit.setSpecialValueText("Auto")
        self.y_limit.valueChanged.connect(self.changed)
        self.add_row("Y limit (-log10 p):", self.y_limit)


class VolcanoPage(FigurePage):
    TITLE = "Volcano Plot"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = "Run a Differential Expression analysis to view this figure."

    def _build_controls(self) -> _VolcanoControls:
        return _VolcanoControls()

    def dependencies_met(self) -> bool:
        de = self.workspace.de_ws
        return de is not None and getattr(de, 'de_results', None) is not None

    def _make_render_func(self):
        if not self.dependencies_met():
            return None

        de = self.workspace.de_ws
        ctrl = self._controls
        dr = de.de_results.copy()
        d = de.disease_label
        c = de.control_label
        dn = Path(de.h5ad_path).stem if getattr(de, 'h5ad_path', None) else "Dataset"
        ml = ctrl.max_labels.value()
        ps = ctrl.point_size.value() ** 2
        cs = ctrl.color_scheme.currentIndex()
        pcol = 'pvals' if ctrl.pval_type.currentIndex() == 1 else 'pvals_adj'
        xlim = ctrl.x_limit.value() or None
        ylim = ctrl.y_limit.value() or None
        figsize = ctrl.figsize.get_figsize() or (8, 8)
        font_sizes = self._font_sizes()

        gene_highlights = set()
        cov = getattr(de, 'pathway_coverage', None) or {}
        for info in cov.values():
            gene_highlights.update(info.get('genes', []))

        def _render(dr=dr, d=d, c=c, dn=dn, ml=ml, ps=ps, cs=cs, pcol=pcol,
                    xlim=xlim, ylim=ylim, figsize=figsize, font_sizes=font_sizes,
                    gene_highlights=gene_highlights):
            from kosmic.visualisation.de.volcano import create_volcano_plot
            return create_volcano_plot(
                dr, d, c, dataset_name=dn, pval_col=pcol,
                genes_of_interest=gene_highlights or None,
                max_labels=ml, point_size=ps, color_scheme=cs,
                x_limit=xlim, y_limit=ylim,
                font_sizes=font_sizes, figsize=figsize,
            )
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(de.project_dir)
