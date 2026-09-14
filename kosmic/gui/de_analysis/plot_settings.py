# Plot Settings Dialog
# Global plot configuration with per-plot-type tabs.
# Auto-detects which plot is active and applies changes in-place.
# And it persists settings via QSettings.

from PyQt6.QtCore import QSettings
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QTabWidget,
    QWidget, QLabel, QComboBox, QPushButton,
    QFormLayout, QSpinBox,
)

from kosmic.gui.shared.theme import get_color, NoScrollDoubleSpinBox
from kosmic import DEFAULT_BASE_FONT_SIZE, DEFAULT_EXPORT_DPI
from kosmic.gui.shared.widgets import SecondaryLabel


# ── Colormap choices ──────────────────────────────────────────────────

HEATMAP_CMAPS = [
    'viridis', 'RdBu_r', 'coolwarm', 'plasma', 'inferno',
    'magma', 'cividis', 'PiYG', 'BrBG', 'seismic',
]

_PAGE_GENE_DE = 4
_PAGE_PATHWAY_DE = 5

TAB_GENERAL = 0
TAB_HEATMAP = 1
TAB_VOLCANO = 2


class PlotSettingsDialog(QDialog):
    """Non-modal plot settings dialog with Apply button."""

    def __init__(self, de_workspace, parent=None):
        super().__init__(parent)
        self.ws = de_workspace
        self.settings = QSettings("KOSMIC", "KOSMIC")
        self.setWindowTitle("Plot Settings")
        self.setMinimumWidth(420)
        self.setMinimumHeight(350)

        self._setup_ui()
        self._load_settings()
        self._auto_select_tab()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------

    def _setup_ui(self):
        layout = QVBoxLayout(self)

        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_general_tab(), "General")
        self.tabs.addTab(self._build_heatmap_tab(), "Heatmap")
        self.tabs.addTab(self._build_volcano_tab(), "Volcano")
        layout.addWidget(self.tabs)

        # Button row
        btn_row = QHBoxLayout()
        btn_row.addStretch()

        self.apply_btn = QPushButton("Apply")
        self.apply_btn.clicked.connect(self._apply)
        btn_row.addWidget(self.apply_btn)

        self.close_btn = QPushButton("Close")
        self.close_btn.clicked.connect(self.close)
        btn_row.addWidget(self.close_btn)

        layout.addLayout(btn_row)

    # ── General ───────────────────────────────────────────────────────

    def _build_general_tab(self):
        w = QWidget()
        form = QFormLayout(w)
        form.setContentsMargins(12, 12, 12, 12)

        self.export_dpi_spin = QSpinBox()
        self.export_dpi_spin.setRange(72, 600)
        self.export_dpi_spin.setValue(150)
        self.export_dpi_spin.setSuffix(" dpi")
        self.export_dpi_spin.setToolTip(
            "Resolution for all exported PNG images.\n"
            "150 = screen quality, 300 = publication quality, 600 = high-res print")
        form.addRow("Export DPI:", self.export_dpi_spin)

        self.base_font_spin = QSpinBox()
        self.base_font_spin.setRange(6, 24)
        self.base_font_spin.setValue(11)
        self.base_font_spin.setSuffix(" pt")
        self.base_font_spin.setToolTip(
            "Base font size for all plot text.\n"
            "All other sizes scale from this:\n"
            "  Title = base × 1.4\n"
            "  Axis labels = base × 1.1\n"
            "  Ticks / legend = base × 0.9\n"
            "  Annotations = base × 0.8")
        form.addRow("Base font size:", self.base_font_spin)

        form.addRow(QLabel(""))  # spacer

        # Show the computed sizes as a live preview
        self.font_preview = SecondaryLabel()
        self.font_preview.setWordWrap(True)
        self._update_font_preview()
        self.base_font_spin.valueChanged.connect(self._update_font_preview)
        form.addRow(self.font_preview)

        return w

    def _update_font_preview(self):
        from kosmic.visualisation import default_font_sizes
        fs = default_font_sizes(self.base_font_spin.value())
        self.font_preview.setText(
            f"Typography scale (base = {fs['base']} pt):\n"
            f"  Title:          {fs['title']:.0f} pt\n"
            f"  Axis labels:    {fs['axis_label']:.0f} pt\n"
            f"  Ticks / legend: {fs['tick']:.0f} pt\n"
            f"  Colorbar:       {fs['colorbar_label']:.0f} pt\n"
            f"  Annotations:    {fs['annotation']:.0f} pt\n\n"
            f"Applies to all static plots (heatmaps, volcano, bar, dotplots)."
        )

    # ── Heatmap ───────────────────────────────────────────────────────

    def _build_heatmap_tab(self):
        w = QWidget()
        form = QFormLayout(w)
        form.setContentsMargins(12, 12, 12, 12)

        self.hm_cmap = QComboBox()
        self.hm_cmap.addItems(HEATMAP_CMAPS)
        form.addRow("Colormap:", self.hm_cmap)

        self.hm_vmin = NoScrollDoubleSpinBox()
        self.hm_vmin.setRange(0, 50)
        self.hm_vmin.setValue(5)
        self.hm_vmin.setSuffix("%")
        form.addRow("Clip lower:", self.hm_vmin)

        self.hm_vmax = NoScrollDoubleSpinBox()
        self.hm_vmax.setRange(50, 100)
        self.hm_vmax.setValue(95)
        self.hm_vmax.setSuffix("%")
        form.addRow("Clip upper:", self.hm_vmax)

        return w

    # ── Volcano ───────────────────────────────────────────────────────

    def _build_volcano_tab(self):
        w = QWidget()
        form = QFormLayout(w)
        form.setContentsMargins(12, 12, 12, 12)

        self.vol_sig_color = QComboBox()
        self.vol_sig_color.addItems([
            'red / blue', 'orange / teal', 'magenta / cyan',
        ])
        form.addRow("Sig. colors:", self.vol_sig_color)

        self.vol_point_size = QSpinBox()
        self.vol_point_size.setRange(3, 20)
        self.vol_point_size.setValue(9)
        form.addRow("Point size:", self.vol_point_size)

        self.vol_max_labels = QSpinBox()
        self.vol_max_labels.setRange(0, 50)
        self.vol_max_labels.setValue(10)
        form.addRow("Max gene labels:", self.vol_max_labels)

        return w

    # ------------------------------------------------------------------
    # Auto-select tab based on active plot
    # ------------------------------------------------------------------

    def _auto_select_tab(self):
        """Select the settings tab matching the currently visible plot."""
        if self.ws is None:
            return

        page_idx = self.ws.stack.currentIndex()

        if page_idx == _PAGE_GENE_DE:
            # Check which results sub-tab is active
            gene_de = self.ws.gene_de_page
            if hasattr(gene_de, '_results_tabs'):
                sub = gene_de._results_tabs.currentIndex()
                tab_name = gene_de._results_tabs.tabText(sub)
                if 'Volcano' in tab_name:
                    self.tabs.setCurrentIndex(TAB_VOLCANO)
                elif 'Heatmap' in tab_name:
                    self.tabs.setCurrentIndex(TAB_HEATMAP)
                else:
                    self.tabs.setCurrentIndex(TAB_GENERAL)
            else:
                self.tabs.setCurrentIndex(TAB_GENERAL)

        elif page_idx == _PAGE_PATHWAY_DE:
            self.tabs.setCurrentIndex(TAB_GENERAL)
        else:
            self.tabs.setCurrentIndex(TAB_GENERAL)

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def _load_settings(self):
        """Load saved plot preferences from QSettings."""
        s = self.settings
        # General
        self.export_dpi_spin.setValue(
            int(s.value("plot/export_dpi", DEFAULT_EXPORT_DPI)))
        self.base_font_spin.setValue(
            int(s.value("plot/base_font_size", DEFAULT_BASE_FONT_SIZE)))

        # Heatmap
        cmap = s.value("plot/heatmap_cmap", "viridis")
        idx = self.hm_cmap.findText(cmap)
        if idx >= 0:
            self.hm_cmap.setCurrentIndex(idx)
        self.hm_vmin.setValue(float(s.value("plot/heatmap_vmin", 5)))
        self.hm_vmax.setValue(float(s.value("plot/heatmap_vmax", 95)))

        # Volcano
        self.vol_point_size.setValue(
            int(s.value("plot/volcano_point_size", 9)))
        self.vol_max_labels.setValue(
            int(s.value("plot/volcano_max_labels", 10)))
        color_idx = int(s.value("plot/volcano_color_scheme", 0))
        self.vol_sig_color.setCurrentIndex(
            min(color_idx, self.vol_sig_color.count() - 1))

    def _save_settings(self):
        """Persist current values to QSettings."""
        s = self.settings
        s.setValue("plot/export_dpi", self.export_dpi_spin.value())
        s.setValue("plot/base_font_size", self.base_font_spin.value())

        s.setValue("plot/heatmap_cmap", self.hm_cmap.currentText())
        s.setValue("plot/heatmap_vmin", self.hm_vmin.value())
        s.setValue("plot/heatmap_vmax", self.hm_vmax.value())

        s.setValue("plot/volcano_point_size", self.vol_point_size.value())
        s.setValue("plot/volcano_max_labels", self.vol_max_labels.value())
        s.setValue("plot/volcano_color_scheme",
                   self.vol_sig_color.currentIndex())

    # ------------------------------------------------------------------
    # Apply
    # ------------------------------------------------------------------

    def _apply(self):
        """Save settings and regenerate the currently visible plot."""
        self._save_settings()
        self.settings.sync()  # flush to disk before plot code reads back

        current_tab = self.tabs.currentIndex()

        tab_names = {TAB_GENERAL: 'General', TAB_HEATMAP: 'Heatmap',
                     TAB_VOLCANO: 'Volcano'}
        self.ws.log_message.emit(
            f"Applying {tab_names.get(current_tab, '?')} settings...")

        try:
            if current_tab == TAB_HEATMAP:
                self._apply_heatmap()
            elif current_tab == TAB_VOLCANO:
                self._apply_volcano()
            else:
                # General tab (base font, DPI) — refresh all active plots
                self._apply_general()
        except Exception as e:
            import traceback
            self.ws.log_message.emit(f"Plot settings apply error: {e}")
            self.ws.log_message.emit(traceback.format_exc())

    def _apply_general(self):
        """Refresh all active plots after base font / DPI change."""
        self.ws.log_message.emit(
            f"General settings saved (base={self.base_font_spin.value()}pt, "
            f"DPI={self.export_dpi_spin.value()}).  Refreshing plots...")
        try:
            self._apply_volcano()
        except Exception as e:
            self.ws.log_message.emit(f"Volcano refresh skipped: {e}")
        try:
            self._apply_heatmap()
        except Exception as e:
            self.ws.log_message.emit(f"Heatmap refresh skipped: {e}")

    def _apply_heatmap(self):
        """Regenerate the gene heatmap with current settings."""
        page = self.ws.gene_de_page
        pw = page.heatmap_pathway_combo.currentText()
        self.ws.log_message.emit(
            f"Heatmap apply: pathway='{pw}', cmap={self.hm_cmap.currentText()}, "
            f"clip={self.hm_vmin.value()}-{self.hm_vmax.value()}%, "
            f"sample_df={'yes' if self.ws.sample_df is not None else 'no'}")
        if pw:
            page._generate_heatmap(pw)
        else:
            self.ws.log_message.emit("No pathway selected for heatmap")

    def _apply_volcano(self):
        """Regenerate the volcano plot with current settings."""
        from kosmic.gui.shared.theme import style_pg_plot, get_font_sizes
        page = self.ws.gene_de_page
        if self.ws.de_results is not None:
            # Re-set volcano data to pick up new point size / colors / fonts
            pval_threshold = page.pval_filter.value()
            fc_threshold = page.fc_filter.value()
            pathway_genes = set()
            if self.ws.pathway_coverage:
                for info in self.ws.pathway_coverage.values():
                    pathway_genes.update(info['genes'])
            page.interactive_volcano.set_data(
                self.ws.de_results,
                pval_threshold=pval_threshold,
                fc_threshold=fc_threshold,
                pathway_genes=pathway_genes if pathway_genes else None,
            )
            # Re-apply theme AFTER set_data (clear() resets axes)
            style_pg_plot(page.interactive_volcano,
                          title='Volcano Plot',
                          left_label='-Log10 Adjusted P-value (FDR)')
            fs = get_font_sizes()
            page.interactive_volcano.setLabel(
                'bottom', 'Log2 Fold Change',
                color=get_color('fg_secondary'),
                size=f'{fs["axis_label"]:.0f}pt')
            self.ws.log_message.emit(
                f"Volcano updated: base={fs['base']}pt, "
                f"point_size={self.vol_point_size.value()}, "
                f"color_scheme={self.vol_sig_color.currentText()}")
        else:
            self.ws.log_message.emit("No DE results to update volcano plot")

