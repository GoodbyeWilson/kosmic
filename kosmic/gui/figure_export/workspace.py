# FigureExportWorkspace -- sibling to scRNA / DE / Meta. Holds references
# to the scRNA + DE workspaces; figure pages read state from those
# workspaces directly. Callers bump 'data_version' when underlying state
# changes so cached pixmaps invalidate.
from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QFrame, QHBoxLayout, QListWidget, QListWidgetItem,
    QStackedWidget, QWidget,
)

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.shared import borderless
from kosmic.gui.figure_export.pages.dataset_summary import DatasetSummaryPage
from kosmic.gui.figure_export.pages.elbow import PCAElbowPage
from kosmic.gui.figure_export.pages.enrichment_bar import EnrichmentBarPage
from kosmic.gui.figure_export.pages.fgsea import FgseaPage
from kosmic.gui.figure_export.pages.gene_heatmap import GeneHeatmapPage
from kosmic.gui.figure_export.pages.library_gsea import LibraryGseaPage
from kosmic.gui.figure_export.pages.library_gsea_mountain import LibraryGseaMountainPage
from kosmic.gui.figure_export.pages.patient_dotplot import PatientDotplotPage
from kosmic.gui.figure_export.pages.pathway_bar import PathwayBarPage
from kosmic.gui.figure_export.pages.pathway_bubble import PathwayBubblePage
from kosmic.gui.figure_export.pages.qc_violins import QCViolinsPage
from kosmic.gui.figure_export.pages.topde_summary import TopDESummaryPage
from kosmic.gui.figure_export.pages.umap import UMAPPage
from kosmic.gui.figure_export.pages.variable_genes import VariableGenesPage
from kosmic.gui.figure_export.pages.volcano import VolcanoPage


# (group_label, page_class)
_PAGE_REGISTRY: list[tuple[str, type[FigurePage]]] = [
    # Embeddings
    ("Embeddings", UMAPPage),
    # Single-cell QC
    ("Single-cell QC", QCViolinsPage),
    ("Single-cell QC", VariableGenesPage),
    # Clustering
    ("Clustering", PCAElbowPage),
    # Differential Expression
    ("Differential Expression", VolcanoPage),
    ("Differential Expression", GeneHeatmapPage),
    ("Differential Expression", PathwayBarPage),
    ("Differential Expression", TopDESummaryPage),
    ("Differential Expression", DatasetSummaryPage),
    # Pathway Activity
    ("Pathway Activity", PathwayBubblePage),
    ("Pathway Activity", PatientDotplotPage),
    # Enrichment
    ("Enrichment", EnrichmentBarPage),
    ("Enrichment", FgseaPage),
    ("Enrichment", LibraryGseaPage),
    ("Enrichment", LibraryGseaMountainPage),
]


class FigureExportWorkspace(QWidget):
    """Top-level Figure Export workspace."""

    status_message = pyqtSignal(str)
    log_message = pyqtSignal(str)

    def __init__(self, scrna_workspace=None, de_workspace=None):
        super().__init__()
        self.scrna_ws = scrna_workspace
        self.de_ws = de_workspace
        self._data_version = 0
        self._pages: list[FigurePage] = []

        self._setup_ui()

    # -- Setters --------------------------------------------------------

    def set_scrna_workspace(self, ws):
        self.scrna_ws = ws
        self.bump_data_version()

    def set_de_workspace(self, ws):
        self.de_ws = ws
        self.bump_data_version()

    def bump_data_version(self):
        """Call after upstream workspace state changes; invalidates page caches."""
        self._data_version += 1
        for p in self._pages:
            p.invalidate()

    def data_version(self) -> int:
        return self._data_version

    # -- UI -------------------------------------------------------------

    def _setup_ui(self):
        layout = borderless(QHBoxLayout, self)

        # Left page-list pane
        self._list = QListWidget()
        self._list.setObjectName("figure_export_nav")
        self._list.setFixedWidth(220)
        self._list.setFrameShape(QFrame.Shape.NoFrame)
        self._list.itemSelectionChanged.connect(self._on_selection_changed)

        self._stack = QStackedWidget()

        current_group: Optional[str] = None
        for group_label, page_cls in _PAGE_REGISTRY:
            if group_label != current_group:
                hdr = QListWidgetItem(group_label.upper())
                hdr.setFlags(Qt.ItemFlag.NoItemFlags)
                font = hdr.font()
                font.setBold(True)
                font.setPointSize(max(8, font.pointSize() - 1))
                hdr.setFont(font)
                self._list.addItem(hdr)
                current_group = group_label

            page = page_cls(workspace=self)
            page.log_message.connect(self.log_message.emit)
            self._stack.addWidget(page)
            self._pages.append(page)

            stack_index = self._stack.count() - 1
            row = QListWidgetItem(f"  {page.TITLE}")
            row.setData(Qt.ItemDataRole.UserRole, stack_index)
            self._list.addItem(row)

        layout.addWidget(self._list)
        layout.addWidget(self._stack, stretch=1)

        # Select the first selectable row
        for i in range(self._list.count()):
            item = self._list.item(i)
            if item.flags() & Qt.ItemFlag.ItemIsSelectable:
                self._list.setCurrentRow(i)
                break

    def _on_selection_changed(self):
        item = self._list.currentItem()
        if item is None:
            return
        idx = item.data(Qt.ItemDataRole.UserRole)
        if idx is None:
            return
        self._stack.setCurrentIndex(int(idx))
        page = self._stack.currentWidget()
        if isinstance(page, FigurePage):
            page.on_activated()

    # -- Theming --------------------------------------------------------

    def refresh_theme(self):
        for p in self._pages:
            p.refresh_theme()

    # -- Activation hook -- called when this workspace becomes visible --

    def on_activated(self):
        """Called by main.py when the user switches to this workspace."""
        page = self._stack.currentWidget()
        if isinstance(page, FigurePage):
            page.on_activated()
