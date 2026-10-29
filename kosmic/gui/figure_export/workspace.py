# FigureExportWorkspace -- sibling to scRNA / DE / Meta. Holds references
# to the scRNA + DE workspaces; figure pages read state from those
# workspaces directly. Callers bump 'data_version' when underlying state
# changes so cached pixmaps invalidate.
#
# The figure list shows one section at a time -- the figures of the scRNA,
# DE or Meta-Analysis workspace, chosen with a switch above it -- and greys
# the figures whose data is not available, naming the step to run in the
# tooltip. main.py opens the section of the workspace the user came from.
from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QBrush, QColor
from PyQt6.QtWidgets import (
    QButtonGroup, QCheckBox, QFrame, QHBoxLayout, QListWidget,
    QListWidgetItem, QPushButton, QStackedWidget, QVBoxLayout, QWidget,
)

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.shared import borderless
from kosmic.gui.shared.theme import get_color
from kosmic.gui.figure_export.pages.cell_types import (
    CompositionPage, MarkerDotPlotPage, MarkerScoreUMAPPage,
)
from kosmic.gui.figure_export.pages.dataset_summary import DatasetSummaryPage
from kosmic.gui.figure_export.pages.elbow import PCAElbowPage
from kosmic.gui.figure_export.pages.enriched_terms import EnrichedTermsPage
from kosmic.gui.figure_export.pages.fgsea import FgseaPage
from kosmic.gui.figure_export.pages.gene_heatmap import GeneHeatmapPage
from kosmic.gui.figure_export.pages.leading_edge import LeadingEdgeHeatmapPage
from kosmic.gui.figure_export.pages.library_gsea import LibraryGseaPage
from kosmic.gui.figure_export.pages.library_gsea_mountain import LibraryGseaMountainPage
from kosmic.gui.figure_export.pages.meta_figures import (
    LeaveOneOutPage, MetaVolcanoPage, PathwayForestPage, PathwayMetaVolcanoPage,
    StudyForestPage, VoteLevelsPage,
)
from kosmic.gui.figure_export.pages.patient_dotplot import PatientDotplotPage
from kosmic.gui.figure_export.pages.pathway_bubble import PathwayBubblePage
from kosmic.gui.figure_export.pages.qc_violins import QCViolinsPage
from kosmic.gui.figure_export.pages.term_links import TermLinksPage
from kosmic.gui.figure_export.pages.topde_summary import TopDESummaryPage
from kosmic.gui.figure_export.pages.umap import UMAPPage
from kosmic.gui.figure_export.pages.variable_genes import VariableGenesPage
from kosmic.gui.figure_export.pages.volcano import VolcanoPage


# (group_label, page_class)
_PAGE_REGISTRY: list[tuple[str, type[FigurePage]]] = [
    # Embeddings
    ("Embeddings", UMAPPage),
    # Cell types
    ("Cell Types", MarkerDotPlotPage),
    ("Cell Types", MarkerScoreUMAPPage),
    ("Cell Types", CompositionPage),
    # Single-cell QC
    ("Single-cell QC", QCViolinsPage),
    ("Single-cell QC", VariableGenesPage),
    # Clustering
    ("Clustering", PCAElbowPage),
    # Differential Expression
    ("Differential Expression", VolcanoPage),
    ("Differential Expression", GeneHeatmapPage),
    ("Differential Expression", TopDESummaryPage),
    ("Differential Expression", DatasetSummaryPage),
    # Pathway Activity
    ("Pathway Activity", PathwayBubblePage),
    ("Pathway Activity", PatientDotplotPage),
    # Enrichment
    ("Enrichment", EnrichedTermsPage),
    ("Enrichment", TermLinksPage),
    ("Enrichment", FgseaPage),
    ("Enrichment", LibraryGseaPage),
    ("Enrichment", LibraryGseaMountainPage),
    ("Enrichment", LeadingEdgeHeatmapPage),
    # Meta-Analysis (ADR-010)
    ("Meta-Analysis", MetaVolcanoPage),
    ("Meta-Analysis", PathwayMetaVolcanoPage),
    ("Meta-Analysis", StudyForestPage),
    ("Meta-Analysis", PathwayForestPage),
    ("Meta-Analysis", LeaveOneOutPage),
    ("Meta-Analysis", VoteLevelsPage),
]

# Sections of the figure list: (key, switch label, groups). Each group, and
# so each figure, is in one section. Enrichment sits with DE, where most
# enrichment runs are made; its figures also draw a Meta-Analysis run.
_SECTIONS: list[tuple[str, str, tuple[str, ...]]] = [
    ("scrna", "scRNA", ("Embeddings", "Cell Types", "Single-cell QC", "Clustering")),
    ("de", "DE", ("Differential Expression", "Pathway Activity", "Enrichment")),
    ("meta", "Meta", ("Meta-Analysis",)),
]
_SECTION_OF_GROUP = {g: key for key, _label, groups in _SECTIONS for g in groups}


class FigureExportWorkspace(QWidget):
    """Top-level Figure Export workspace."""

    status_message = pyqtSignal(str)
    log_message = pyqtSignal(str)

    def __init__(self, scrna_workspace=None, de_workspace=None, meta_workspace=None):
        super().__init__()
        self.scrna_ws = scrna_workspace
        self.de_ws = de_workspace
        # Read through its figure_inputs() (ADR-010).
        self.meta_ws = meta_workspace
        self._data_version = 0
        self._pages: list[FigurePage] = []
        # (section, group, page) for every figure, in registry order.
        self._entries: list[tuple[str, str, FigurePage]] = []
        self._section = _SECTIONS[0][0]

        self._setup_ui()

    # -- Setters --------------------------------------------------------

    def set_scrna_workspace(self, ws):
        self.scrna_ws = ws
        self.bump_data_version()

    def set_de_workspace(self, ws):
        self.de_ws = ws
        self.bump_data_version()

    def set_meta_workspace(self, ws):
        self.meta_ws = ws
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

        # Left pane: section switch, availability filter, figure list.
        nav = QWidget()
        nav.setFixedWidth(220)
        nav_lay = QVBoxLayout(nav)
        nav_lay.setContentsMargins(8, 8, 8, 0)
        nav_lay.setSpacing(6)

        switch = QHBoxLayout()
        switch.setSpacing(6)
        self._section_group = QButtonGroup(self)
        self._section_group.setExclusive(True)
        self._section_buttons: dict[str, QPushButton] = {}
        for key, label, _groups in _SECTIONS:
            btn = QPushButton(label)
            btn.setProperty("role", "filter_chip")
            btn.setCheckable(True)
            btn.setChecked(key == self._section)
            btn.clicked.connect(lambda _c=False, k=key: self.show_section(k))
            self._section_group.addButton(btn)
            self._section_buttons[key] = btn
            switch.addWidget(btn)
        switch.addStretch(1)
        nav_lay.addLayout(switch)

        self._available_only = QCheckBox("Show available only")
        self._available_only.setToolTip(
            "Hide the figures whose data is not loaded or not yet computed.")
        self._available_only.toggled.connect(self._rebuild_list)
        nav_lay.addWidget(self._available_only)

        self._list = QListWidget()
        self._list.setObjectName("figure_export_nav")
        self._list.setFrameShape(QFrame.Shape.NoFrame)
        self._list.itemSelectionChanged.connect(self._on_selection_changed)
        nav_lay.addWidget(self._list, 1)

        self._stack = QStackedWidget()
        for group_label, page_cls in _PAGE_REGISTRY:
            page = page_cls(workspace=self)
            page.log_message.connect(self.log_message.emit)
            self._stack.addWidget(page)
            self._pages.append(page)
            self._entries.append((_SECTION_OF_GROUP[group_label], group_label, page))

        layout.addWidget(nav)
        layout.addWidget(self._stack, stretch=1)
        self._rebuild_list()

    # -- Sections and availability ---------------------------------------

    def show_section(self, key: Optional[str]) -> None:
        """Show the figures of one section ('scrna', 'de' or 'meta'); None
        or an unknown key keeps the current one."""
        if key not in self._section_buttons:
            return
        self._section_buttons[key].setChecked(True)
        if key != self._section:
            self._section = key
            self._rebuild_list()

    def current_section(self) -> str:
        return self._section

    @staticmethod
    def _available(page: FigurePage) -> bool:
        try:
            return bool(page.dependencies_met())
        except Exception:  # noqa: BLE001 -- a failing check means unavailable
            return False

    def _rebuild_list(self) -> None:
        """Fill the list with the current section's figures: unavailable ones
        greyed with the step to run as tooltip, or left out when 'Show
        available only' is ticked. Keeps the open figure selected when it
        is still listed, else opens the first available one."""
        current = self._stack.currentWidget()
        hide_unavailable = self._available_only.isChecked()
        grey = QBrush(QColor(get_color('fg_tertiary')))

        self._list.blockSignals(True)
        self._list.clear()
        rows: list[tuple[QListWidgetItem, FigurePage, bool]] = []
        group = None
        for section, group_label, page in self._entries:
            if section != self._section:
                continue
            available = self._available(page)
            if hide_unavailable and not available:
                continue
            if group_label != group:
                hdr = QListWidgetItem(group_label.upper())
                hdr.setFlags(Qt.ItemFlag.NoItemFlags)
                font = hdr.font()
                font.setBold(True)
                font.setPointSize(max(8, font.pointSize() - 1))
                hdr.setFont(font)
                self._list.addItem(hdr)
                group = group_label
            row = QListWidgetItem(f"  {page.TITLE}")
            row.setData(Qt.ItemDataRole.UserRole, self._stack.indexOf(page))
            if not available:
                row.setForeground(grey)
                row.setToolTip(page.UNAVAILABLE_MESSAGE)
            self._list.addItem(row)
            rows.append((row, page, available))

        pick = (next((r for r, p, _a in rows if p is current), None)
                or next((r for r, _p, a in rows if a), None)
                or (rows[0][0] if rows else None))
        if pick is not None:
            self._list.setCurrentItem(pick)
        self._list.blockSignals(False)
        if pick is not None and self._stack.widget(
                int(pick.data(Qt.ItemDataRole.UserRole))) is not current:
            self._on_selection_changed()

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
        self._rebuild_list()           # the grey of unavailable figures

    # -- Activation hook -- called when this workspace becomes visible --

    def on_activated(self):
        """Called by main.py when the user switches to this workspace.
        Availability may have changed since the list was last filled."""
        self._rebuild_list()
        page = self._stack.currentWidget()
        if isinstance(page, FigurePage):
            page.on_activated()
