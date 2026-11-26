# Methods Comparison workspace (ADR-017), embedded in AppWindow's stack.
#
# Operates on the project, like Meta-Analysis: it compares merged
# pseudobulk on the shared atlas with REML, SumRank, Stouffer and vote
# counting for one counts source and cell type. Figures reads the
# comparison it holds through 'figure_inputs()', or else the newest saved
# comparison under 'methods_comparison/'.

from pathlib import Path
from typing import Optional

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QStackedWidget, QVBoxLayout, QWidget


class ComparisonWorkspace(QWidget):
    """The Methods Comparison workspace: one page, no workflow steps."""

    status_message = pyqtSignal(str)
    log_message = pyqtSignal(str)
    project_directory_changed = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        from kosmic.gui.comparison.pages.compare_page import ComparePage

        self.progress_bar = None
        self.project_folder: Optional[Path] = None
        self.bundle = None   # the ComparisonBundle shown, once one is run or loaded

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.stack = QStackedWidget()
        self.compare_page = ComparePage(self)
        self.stack.addWidget(self.compare_page)
        lay.addWidget(self.stack)

    def set_project_directory_external(self, folder) -> None:
        """Set the project from AppWindow.project_changed."""
        if not folder or not Path(folder).is_dir():
            return
        self.project_folder = Path(folder)
        self.bundle = None
        self.project_directory_changed.emit(str(folder))
        self.compare_page.set_project(self.project_folder)

    def on_activated(self) -> None:
        """Called by AppWindow each time the workspace is shown."""

    def figure_inputs(self) -> dict:
        """What Figures reads from this workspace: the comparison shown and
        the project it belongs to (ADR-017 §6)."""
        return {'bundle': self.bundle, 'project_folder': self.project_folder}
