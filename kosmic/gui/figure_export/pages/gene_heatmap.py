# Gene heatmap page (per-pathway sample x gene heatmap).
#
# Reads from 'workspace.de_ws' (current_adata + pathway_gene_sets +
# sample_col + condition_col + min_cells + disease/control_label).

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QComboBox, QHBoxLayout, QPushButton, QSizePolicy, QSpinBox,
)

from kosmic.de.de_analysis import tested_gene_set
from kosmic.gui.de_analysis.plot_settings import HEATMAP_CMAPS
from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls
from kosmic.gui.shared.widgets import SecondaryLabel


class _HeatmapControls(FigureControls):
    pathway_changed = pyqtSignal()

    def __init__(self, default_w: int = 10, default_h: int = 8):
        super().__init__(default_w, default_h)

        # Pathway selector with prev/next nav and a counter label.
        pw_row = QHBoxLayout()
        pw_row.setSpacing(4)

        self._prev_btn = QPushButton("<")
        self._prev_btn.setFixedWidth(28)
        self._prev_btn.clicked.connect(self._prev)
        pw_row.addWidget(self._prev_btn)

        self.pathway_combo = QComboBox()
        self.pathway_combo.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.pathway_combo.currentIndexChanged.connect(self._on_pathway_changed)
        pw_row.addWidget(self.pathway_combo)

        self._next_btn = QPushButton(">")
        self._next_btn.setFixedWidth(28)
        self._next_btn.clicked.connect(self._next)
        pw_row.addWidget(self._next_btn)

        from PyQt6.QtWidgets import QWidget
        pw_container = QWidget()
        pw_container.setLayout(pw_row)
        self.add_row("Pathway:", pw_container)

        self._pw_label = SecondaryLabel("")
        self._pw_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.add_row("", self._pw_label)

        self.cmap = QComboBox()
        self.cmap.addItems(HEATMAP_CMAPS)
        self.cmap.currentIndexChanged.connect(self.changed)
        self.add_row("Colormap:", self.cmap)

        self.vmin = QSpinBox()
        self.vmin.setRange(0, 49)
        self.vmin.setValue(5)
        self.vmin.setSuffix("%")
        self.vmin.valueChanged.connect(self.changed)
        self.add_row("Min percentile:", self.vmin)

        self.vmax = QSpinBox()
        self.vmax.setRange(51, 100)
        self.vmax.setValue(95)
        self.vmax.setSuffix("%")
        self.vmax.valueChanged.connect(self.changed)
        self.add_row("Max percentile:", self.vmax)

    def _on_pathway_changed(self):
        self._update_nav_state()
        self.changed.emit()

    def _update_nav_state(self):
        idx = self.pathway_combo.currentIndex()
        n = self.pathway_combo.count()
        self._prev_btn.setEnabled(idx > 0)
        self._next_btn.setEnabled(0 <= idx < n - 1)
        self._pw_label.setText(f"{idx + 1} / {n}" if n > 0 else "")

    def _prev(self):
        i = self.pathway_combo.currentIndex()
        if i > 0:
            self.pathway_combo.setCurrentIndex(i - 1)

    def _next(self):
        i = self.pathway_combo.currentIndex()
        if i < self.pathway_combo.count() - 1:
            self.pathway_combo.setCurrentIndex(i + 1)

    def populate_pathways(self, pathways):
        prev = self.pathway_combo.currentText()
        self.pathway_combo.blockSignals(True)
        self.pathway_combo.clear()
        for pw in pathways:
            self.pathway_combo.addItem(pw)
        self.pathway_combo.blockSignals(False)
        if prev:
            idx = self.pathway_combo.findText(prev)
            if idx >= 0:
                self.pathway_combo.setCurrentIndex(idx)
        self._update_nav_state()


class GeneHeatmapPage(FigurePage):
    TITLE = "Gene Heatmap"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Load an adata with sample/condition assignments and selected gene sets "
        "to view this figure."
    )

    def __init__(self, workspace):
        # Lazy cache, keyed on data_version + heatmap params.
        self._sample_df = None
        self._available_pathways: dict = {}
        self._prep_version: int = -1
        self._prep_signature: tuple = ()
        super().__init__(workspace)

    def _build_controls(self) -> _HeatmapControls:
        return _HeatmapControls()

    def dependencies_met(self) -> bool:
        de = self.workspace.de_ws
        if de is None:
            return False
        return (getattr(de, 'current_adata', None) is not None
                and bool(getattr(de, 'pathway_gene_sets', None))
                and bool(getattr(de, 'sample_col', None))
                and bool(getattr(de, 'condition_col', None)))

    def _prep_inputs(self):
        de = self.workspace.de_ws
        return (
            de.current_adata,
            de.pathway_gene_sets,
            de.sample_col,
            de.condition_col,
            getattr(de, 'min_cells', 3),
        )

    def _ensure_prepared(self):
        """Recompute sample_df + available_pathways if upstream changed."""
        if not self.dependencies_met():
            self._sample_df = None
            self._available_pathways = {}
            return
        adata, gs, sample_col, cond_col, min_cells = self._prep_inputs()
        de = self.workspace.de_ws
        allowed = tested_gene_set(getattr(de, 'de_results', None))
        signature = (id(adata), sample_col, cond_col, min_cells, len(gs),
                     id(getattr(de, 'de_results', None)),
                     len(allowed) if allowed is not None else -1)
        version = self.workspace.data_version()
        if (self._sample_df is not None
                and self._prep_version == version
                and self._prep_signature == signature):
            return
        try:
            from kosmic.visualisation.de.gene_heatmap import prepare_heatmap_data
            # Restrict to the DE-tested gene set (shared with the in-app heatmap)
            # so genes filtered out of DE, e.g. one-donor genes, never show.
            sd, avail = prepare_heatmap_data(
                adata, gs, sample_col, cond_col, min_cells=min_cells,
                allowed_genes=allowed)
            self._sample_df = sd
            self._available_pathways = avail
            self._prep_version = version
            self._prep_signature = signature
        except Exception as e:
            self.log_message.emit(f"Heatmap prep failed: {e}")
            self._sample_df = None
            self._available_pathways = {}

    def on_activated(self):
        self._ensure_prepared()
        if self._available_pathways:
            self._controls.populate_pathways(list(self._available_pathways.keys()))
        super().on_activated()

    def invalidate(self):
        super().invalidate()
        self._sample_df = None
        self._available_pathways = {}
        self._prep_version = -1
        self._prep_signature = ()

    def _make_render_func(self):
        if not self.dependencies_met():
            return None
        if self._sample_df is None or not self._available_pathways:
            return None
        ctrl = self._controls
        pw_name = ctrl.pathway_combo.currentText()
        if not pw_name:
            return None
        genes = self._available_pathways.get(pw_name)
        if not genes:
            return None

        de = self.workspace.de_ws
        sd = self._sample_df
        d = de.disease_label
        c = de.control_label
        cmap = ctrl.cmap.currentText()
        vmin = ctrl.vmin.value()
        vmax = ctrl.vmax.value()
        font_sizes = self._font_sizes()
        figsize = ctrl.figsize.get_figsize()

        def _render(sd=sd, genes=genes, pn=pw_name, c=c, d=d, cmap=cmap,
                    vmin=vmin, vmax=vmax, f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.gene_heatmap import create_gene_heatmap
            return create_gene_heatmap(
                sd, genes, pn, c, d, cmap=cmap,
                vmin_pct=vmin, vmax_pct=vmax,
                font_sizes=f, figsize=sz,
            )
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(de.project_dir)

    def _default_export_basename(self) -> str:
        pw = self._controls.pathway_combo.currentText() if self._controls else ""
        safe = pw.replace('/', '_').replace(' ', '_') if pw else "pathway"
        return f"heatmap_{safe}"
