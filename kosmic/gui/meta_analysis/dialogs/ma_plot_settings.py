# Meta-Analysis Plot Settings
# Non-modal dialog for meta-analysis plot options (Y-axis, significance
# toggles, Olink overlays). Opened from View > Plot Settings when the
# meta-analysis workspace is active.
#
# Updates the GeneMAPage plots live via direct widget references.

from __future__ import annotations

from PyQt6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QFormLayout, QGroupBox,
    QHBoxLayout, QPushButton, QVBoxLayout,
)


class MAPlotSettingsDialog(QDialog):
    """Non-modal dialog wrapping the meta-analysis plot controls."""

    def __init__(self, gene_ma_page, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Meta-Analysis Plot Settings")
        self.setMinimumWidth(360)
        self._cp = gene_ma_page

        layout = QVBoxLayout(self)

        # Y-axis
        yaxis_group = QGroupBox("Volcano Plot")
        yg = QFormLayout(yaxis_group)
        self._yaxis_combo = QComboBox()
        self._yaxis_combo.addItems(["-log10(FDR)", "Z-statistic"])
        # Sync with current state
        self._yaxis_combo.setCurrentText(
            gene_ma_page._yaxis_combo.currentText())
        self._yaxis_combo.currentIndexChanged.connect(self._apply_yaxis)
        yg.addRow("Y-axis:", self._yaxis_combo)
        layout.addWidget(yaxis_group)

        # Olink overlays
        olink_group = QGroupBox("Olink Panel Overlays")
        og = QVBoxLayout(olink_group)
        self._olink_cvd = QCheckBox("CVD III (amber rings)")
        self._olink_explore = QCheckBox("Explore 3072 (cyan rings)")
        self._olink_cvd.setChecked(
            gene_ma_page._olink_cvd_check.isChecked())
        self._olink_explore.setChecked(
            gene_ma_page._olink_explore_check.isChecked())
        self._olink_cvd.toggled.connect(self._apply_olink)
        self._olink_explore.toggled.connect(self._apply_olink)
        og.addWidget(self._olink_cvd)
        og.addWidget(self._olink_explore)
        layout.addWidget(olink_group)

        # Close button
        btn_row = QHBoxLayout()
        btn_row.addStretch()
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

    def _apply_yaxis(self):
        self._cp._yaxis_combo.setCurrentText(
            self._yaxis_combo.currentText())

    def _apply_olink(self):
        self._cp._olink_cvd_check.setChecked(self._olink_cvd.isChecked())
        self._cp._olink_explore_check.setChecked(self._olink_explore.isChecked())
