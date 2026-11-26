# The Methods Comparison workspace's one page (ADR-017): choose a counts
# source and cell type, the methods and the call rule, check the merged
# design, run the comparison and read the summary. The controls are added
# in the next build step; this is the page's frame.

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QVBoxLayout

from kosmic.gui.shared.widgets import HintLabel, SidebarPage


class ComparePage(SidebarPage):
    """Merged pseudobulk against the meta-analysis methods."""

    help_id = "comparison/index"

    def __init__(self, workspace):
        super().__init__()
        self.workspace = workspace
        self.sidebar_layout.addStretch()
        lay = QVBoxLayout(self.content)
        lay.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.placeholder = HintLabel(
            "Open a project to compare merged pseudobulk with meta-analysis.")
        lay.addWidget(self.placeholder)

    def set_project(self, project):
        self.placeholder.setText(
            "Open a project to compare merged pseudobulk with meta-analysis."
            if project is None else
            f"Project: {project.name}. The comparison controls are not built yet.")
