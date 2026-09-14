# Add Data dialog: Project-hosted intake (ADR-001 Phase E). Hosts the
# GEO import and assemble-from-parts tools for the active study; local
# file imports stay simple file dialogs. On completion the host
# refreshes the project scan and the study's Dataset screen.
from __future__ import annotations

from pathlib import Path

from PyQt6.QtWidgets import QDialog, QHBoxLayout, QTabWidget, QVBoxLayout

from kosmic.gui.intake.assemble_page import AssemblePage
from kosmic.gui.intake.geo_page import GeoImportPage
from kosmic.gui.shared.widgets import CaptionLabel, SecondaryButton


class AddDataDialog(QDialog):
    """Intake dialog for one study: Download from GEO | Assemble from parts.

    'study_dir' is the target study folder; 'log_cb' receives log
    lines; 'on_data_changed' fires whenever files land in the study
    (raw downloads or an assembled h5ad) so the host can refresh.
    """

    def __init__(self, study_dir: Path, log_cb, on_data_changed,
                 parent=None, start: str = 'geo'):
        super().__init__(parent)
        self._study_dir = Path(study_dir)
        self._changed = False

        def _mark_changed(*_a):
            self._changed = True
            on_data_changed()

        self.setWindowTitle(f"Add Data -- {self._study_dir.name}")
        self.resize(980, 680)

        layout = QVBoxLayout(self)
        layout.setSpacing(8)
        layout.addWidget(CaptionLabel(
            f"Data added here is registered to the study "
            f"'{self._study_dir.name}'. Conversion of downloaded raw "
            f"files happens on the study's Dataset screen."))

        tabs = QTabWidget()
        self._geo_page = GeoImportPage(
            study_dir_provider=lambda: self._study_dir,
            log_cb=log_cb,
            on_raw_data_changed=_mark_changed,
        )
        tabs.addTab(self._geo_page, "Download from GEO")
        self._assemble_page = AssemblePage(
            study_dir_provider=lambda: self._study_dir,
            log_cb=log_cb,
            on_h5ad_ready=_mark_changed,
        )
        tabs.addTab(self._assemble_page, "Assemble from Parts")
        tabs.setCurrentIndex(1 if start == 'assemble' else 0)
        layout.addWidget(tabs, 1)

        close_row = QHBoxLayout()
        close_row.addStretch()
        close_btn = SecondaryButton("Close")
        close_btn.clicked.connect(self.accept)
        close_row.addWidget(close_btn)
        layout.addLayout(close_row)

    @property
    def data_changed(self) -> bool:
        """True when any download/assembly landed files in the study."""
        return self._changed
