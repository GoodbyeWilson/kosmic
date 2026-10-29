# Gene heatmap page (per-pathway sample x gene heatmap).
#
# Reads the study's DE run (FigurePage._de: pathway_gene_sets, sample_col,
# condition_col, min_cells, counts_source, disease/control_label) and the
# counts of the gene sets' genes (StudySource.de_counts). Preparing the
# per-sample table reads those counts, from the study's file when the DE
# workspace does not hold it, so it runs on the render thread; the pathway
# selector is filled once the first preparation is drawn.

from __future__ import annotations

import threading
from pathlib import Path
from typing import Optional

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QComboBox, QHBoxLayout, QPushButton, QSizePolicy, QSpinBox,
)

from kosmic.de.de_analysis import tested_gene_set
from kosmic.gui.de_analysis.plot_settings import HEATMAP_CMAPS
from kosmic.visualisation.de.gene_heatmap import HOUSE_CMAP
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
        self.cmap.addItems([HOUSE_CMAP] + HEATMAP_CMAPS)
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
        # (signature, sample_df, available_pathways) of the last preparation.
        # Written on the render thread, read on the GUI thread.
        self._prepared: Optional[tuple] = None
        self._prep_lock = threading.Lock()
        super().__init__(workspace)

    def _build_controls(self) -> _HeatmapControls:
        return _HeatmapControls()

    def dependencies_met(self) -> bool:
        de = self._de()
        return (de is not None and de.available('pathway_gene_sets')
                and de.available('sample_col') and de.available('condition_col')
                and (de.current_adata is not None or de.h5ad_path is not None))

    def _inputs(self) -> Optional[dict]:
        """What the preparation needs, taken on the GUI thread; None when the
        figure is not available."""
        if not self.dependencies_met():
            return None
        de = self._de()
        gs = de.pathway_gene_sets
        if not gs:
            return None
        # Restrict to the DE-tested gene set (shared with the in-app heatmap)
        # so genes filtered out of DE, e.g. one-donor genes, never show.
        allowed = tested_gene_set(de.de_results)
        counts_source = de.counts_source or 'raw'
        inputs = dict(gene_sets=dict(gs), sample_col=de.sample_col,
                      condition_col=de.condition_col, min_cells=de.min_cells,
                      allowed=allowed, counts_source=counts_source)
        inputs['signature'] = (self.workspace.data_version(), inputs['sample_col'],
                               inputs['condition_col'], inputs['min_cells'],
                               tuple(sorted(gs)), counts_source,
                               len(allowed) if allowed is not None else -1)
        return inputs

    def _prepare(self, inputs: dict):
        """(sample_df, available_pathways) for 'inputs', computed once per
        signature. Reads the counts of the gene sets' genes: the DE
        workspace's dataset when it holds the study with the run's counts,
        else those genes from the study's file. Call from a render thread."""
        with self._prep_lock:
            if self._prepared is not None and self._prepared[0] == inputs['signature']:
                return self._prepared[1], self._prepared[2]
        from kosmic.scrna.inspect.detection import detect_species, format_gene_for_species
        from kosmic.visualisation.de.gene_heatmap import prepare_heatmap_data
        src = self._study_source()
        names = src.var_names() or []
        species = detect_species(names) if names else None
        present = set(names)
        genes = list(dict.fromkeys(
            g for gene_list in inputs['gene_sets'].values()
            for g in (format_gene_for_species(x, species or 'human') for x in gene_list)
            if g in present))
        adata = src.de_counts(genes, inputs['counts_source'])
        if adata is None:
            return None, {}
        sd, avail = prepare_heatmap_data(
            adata, inputs['gene_sets'], inputs['sample_col'], inputs['condition_col'],
            min_cells=inputs['min_cells'], allowed_genes=inputs['allowed'],
            species=species)
        with self._prep_lock:
            self._prepared = (inputs['signature'], sd, avail)
        return sd, avail

    def _sync_pathways(self):
        """List the pathways of the last preparation in the selector, without
        drawing again (the render already drew the one selected)."""
        with self._prep_lock:
            avail = self._prepared[2] if self._prepared is not None else {}
        combo = self._controls.pathway_combo
        names = list(avail)
        if names == [combo.itemText(i) for i in range(combo.count())]:
            return
        self._controls.blockSignals(True)
        self._controls.populate_pathways(names)
        self._controls.blockSignals(False)

    def _on_render_finished(self, png_bytes: bytes):
        super()._on_render_finished(png_bytes)
        self._sync_pathways()

    def invalidate(self):
        super().invalidate()
        with self._prep_lock:
            self._prepared = None

    def _make_render_func(self):
        inputs = self._inputs()
        if inputs is None:
            return None
        de = self._de()
        ctrl = self._controls
        pw_name = ctrl.pathway_combo.currentText()
        d = de.disease_label
        c = de.control_label
        cmap = ctrl.cmap.currentText()
        vmin = ctrl.vmin.value()
        vmax = ctrl.vmax.value()
        font_sizes = self._font_sizes()
        figsize = ctrl.figsize.get_figsize()

        def _render(inputs=inputs, pn=pw_name, c=c, d=d, cmap=cmap,
                    vmin=vmin, vmax=vmax, f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.gene_heatmap import create_gene_heatmap
            sd, avail = self._prepare(inputs)
            # The pathway selected, else the first, before the selector is
            # filled (it is filled from this preparation once drawn).
            name = pn if pn in avail else next(iter(avail), None)
            if sd is None or name is None:
                return None
            return create_gene_heatmap(
                sd, avail[name], name, c, d, cmap=cmap,
                vmin_pct=vmin, vmax_pct=vmax,
                font_sizes=f, figsize=sz,
            )
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        return self._de_export_dir()

    def _default_export_basename(self) -> str:
        pw = self._controls.pathway_combo.currentText() if self._controls else ""
        safe = pw.replace('/', '_').replace(' ', '_') if pw else "pathway"
        return f"heatmap_{safe}"
