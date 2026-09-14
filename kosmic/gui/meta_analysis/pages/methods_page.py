# Discovery step 6: the auto-generated record of what actually ran.
#
# This was a "View Methods" button in the pooling page's sidebar that
# opened a modal dialog. Two problems with that: what it prints spans
# every study and every step -- including Enrichment and Validation,
# which are their own steps -- so it did not belong to the pooling page;
# and a modal meant you could not read it beside the results it
# describes. The DE workspace already ends on a Methods step, so this
# matches it.
#
# Nothing here computes anything. Every line comes from the provenance
# sidecars written as each stage ran, so what you read is what happened,
# not a description of what the settings currently say.
from __future__ import annotations

from typing import Optional

from PyQt6.QtGui import QFont, QGuiApplication
from PyQt6.QtWidgets import (
    QHBoxLayout, QPushButton, QTextEdit, QWidget,
)

from kosmic.gui.shared.widgets import SecondaryLabel, SimplePage


class MetaMethodsPage(SimplePage):
    """Combined cross-study methods, rendered from provenance."""

    help_id = "meta/methods"

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._project_folder = None
        self._study_names = None
        self._setup_ui()

    def _setup_ui(self):
        lay = self.body_layout

        self._status = SecondaryLabel(
            "Set a project directory and run an analysis to see the "
            "recorded methods.")
        self._status.setWordWrap(True)
        lay.addWidget(self._status)

        self._view = QTextEdit()
        self._view.setReadOnly(True)
        # No wrapping: the rendered stages are aligned key/value lines
        # and rewrapping them makes the parameter columns unreadable.
        self._view.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self._view.setFont(QFont("Consolas", 10))
        lay.addWidget(self._view, 1)

        row = QHBoxLayout()
        self._copy_btn = QPushButton("Copy to clipboard")
        self._copy_btn.clicked.connect(self._copy_text)
        self._copy_btn.setEnabled(False)
        row.addWidget(self._copy_btn)

        self._refresh_btn = QPushButton("Refresh")
        self._refresh_btn.setToolTip(
            "Re-read the provenance sidecars. Useful after running a "
            "step in another workspace.")
        self._refresh_btn.clicked.connect(self.refresh)
        row.addWidget(self._refresh_btn)
        row.addStretch()
        lay.addLayout(row)

    # -- inputs ---------------------------------------------------------

    def set_project(self, project_folder, study_names=None):
        self._project_folder = project_folder
        self._study_names = list(study_names) if study_names else None
        self.refresh()

    def on_activated(self):
        """Re-read whenever the step is opened, so it is never stale."""
        self.refresh()

    def refresh(self):
        if not self._project_folder:
            self._view.clear()
            self._copy_btn.setEnabled(False)
            self._status.setText(
                "Set a project directory to see the recorded methods.")
            return

        from kosmic.meta_analysis.io import build_combined_methods
        text = build_combined_methods(self._project_folder, self._study_names)
        if not text:
            self._view.clear()
            self._copy_btn.setEnabled(False)
            self._status.setText(
                "No recorded methods yet. Each stage writes its own entry "
                "as it runs, so this fills in as you work through the "
                "studies and the analysis.")
            return

        self._view.setPlainText(text)
        self._copy_btn.setEnabled(True)
        n_studies = text.count("===== Study:")
        self._status.setText(
            f"Recorded methods across {n_studies} "
            f"stud{'y' if n_studies == 1 else 'ies'} plus the "
            f"meta-analysis. Copy-pasteable.")

    def _copy_text(self):
        text = self._view.toPlainText()
        if not text:
            return
        QGuiApplication.clipboard().setText(text)
        self._status.setText("Copied to clipboard.")
