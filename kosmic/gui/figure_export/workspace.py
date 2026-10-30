# FigureExportWorkspace -- sibling to scRNA / DE / Meta. Holds references
# to the scRNA, DE and Meta-Analysis workspaces; the scRNA and DE figure
# pages read them through the study source (below), the Meta figures
# through the Meta workspace's figure_inputs(). Callers bump 'data_version'
# when underlying state changes so cached pixmaps invalidate.
#
# The figure list shows one section at a time -- the figures of the scRNA,
# DE or Meta-Analysis workspace, chosen with a switch above it -- and greys
# the figures whose data is not available, naming the step to run in the
# tooltip. main.py opens the section of the workspace the user came from.
#
# The scRNA and DE figures show one study, named above the switch, through
# 'self.source' (ADR-014): a workspace's data when it holds the study, else
# the study's folder. main.py sets it to the active study, or to the study
# of the Project page's 'Open in Figures'. When
# the study has more than one DE run (one per cell type and method, say),
# the DE section shows a Run selector above the list.
#
# The Meta figures draw 'meta_inputs()': the Meta-Analysis workspace's run
# when it holds one, else a run saved in the project's meta_analysis/
# folder (kosmic.meta_analysis.saved), chosen with the same Run selector in
# the Meta section.
from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QBrush, QColor
from PyQt6.QtWidgets import (
    QButtonGroup, QCheckBox, QComboBox, QFrame, QHBoxLayout, QLabel, QListWidget,
    QListWidgetItem, QPushButton, QStackedWidget, QVBoxLayout, QWidget,
)

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.study_source import StudySource
from kosmic.gui.shared import borderless
from kosmic.gui.shared.theme import get_color
from kosmic.gui.shared.widgets import SecondaryLabel
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
    # The same enrichment figures, drawing the pooled run
    ("Pooled Enrichment", EnrichedTermsPage),
    ("Pooled Enrichment", TermLinksPage),
    ("Pooled Enrichment", LibraryGseaPage),
    ("Pooled Enrichment", LibraryGseaMountainPage),
    ("Pooled Enrichment", LeadingEdgeHeatmapPage),
]

# Sections of the figure list: (key, switch label, groups). Each group is
# in one section. The enrichment figures are listed in both DE (a study's
# run) and Meta (the pooled run); the section decides which they draw.
_SECTIONS: list[tuple[str, str, tuple[str, ...]]] = [
    ("scrna", "scRNA", ("Embeddings", "Cell Types", "Single-cell QC", "Clustering")),
    ("de", "DE", ("Differential Expression", "Pathway Activity", "Enrichment")),
    ("meta", "Meta", ("Meta-Analysis", "Pooled Enrichment")),
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
        # The study the scRNA and DE figures show, and where their data
        # comes from: a workspace that holds it, else its folder (ADR-014).
        self.source = StudySource(lambda: self.scrna_ws, lambda: self.de_ws, parent=self)
        self.source.changed.connect(self._on_source_changed)
        # The project whose saved meta-analysis runs the Meta figures can
        # draw, the saved run chosen (None: the Meta workspace's, else the
        # first saved), the runs found when the list was last filled, and
        # (run, files signature, inputs) of the saved run last read.
        self._project: Optional[Path] = None
        self._meta_run = None
        self._meta_runs: Optional[list] = None
        self._meta_saved: Optional[tuple] = None

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

    def set_study(self, study_dir) -> None:
        """Show the scRNA and DE figures of 'study_dir' (a study folder)."""
        self.source.set_study(study_dir)

    def set_project(self, folder) -> None:
        """The project whose saved meta-analysis runs the Meta figures can
        draw."""
        new = Path(folder) if folder else None
        if new == self._project:
            return
        self._project = new
        self._meta_run, self._meta_runs, self._meta_saved = None, None, None
        self.bump_data_version()

    # -- The Meta figures' run ----------------------------------------------

    def _held_meta_inputs(self) -> Optional[dict]:
        """The Meta-Analysis workspace's figure_inputs() when it holds a
        pooling run, else None."""
        inputs = self.meta_ws.figure_inputs() if self.meta_ws is not None else None
        if inputs and (inputs.get('gene_results') is not None
                       or inputs.get('pathway_results') is not None):
            return inputs
        return None

    def _project_folder(self) -> Optional[Path]:
        if self._project is not None:
            return self._project
        study = self.source.study()
        return study.parent if study is not None else None

    def meta_runs(self) -> list:
        """The project's saved meta-analysis runs (kosmic.meta_analysis.saved)."""
        if self._meta_runs is None:
            from kosmic.meta_analysis.saved import list_meta_runs
            project = self._project_folder()
            self._meta_runs = list_meta_runs(project) if project is not None else []
        return self._meta_runs

    def meta_inputs(self):
        """What the Meta figures draw: the Meta-Analysis workspace's run
        unless a saved run was chosen, else the first saved run, read from
        the project folder (a mapping with the same keys)."""
        held = self._held_meta_inputs()
        if self._meta_run is None and held is not None:
            return held
        runs = self.meta_runs()
        run = self._meta_run if self._meta_run is not None else (runs[0] if runs else None)
        if run is None:
            return self.meta_ws.figure_inputs() if self.meta_ws is not None else None
        try:
            sig = tuple(sorted((p.name, p.stat().st_mtime_ns) for p in run.folder.iterdir()))
        except OSError:
            sig = None
        if self._meta_saved is None or self._meta_saved[:2] != (run, sig):
            from kosmic.meta_analysis.saved import figure_inputs
            self._meta_saved = (run, sig, figure_inputs(run))
        return self._meta_saved[2]

    def _set_meta_run(self, run) -> None:
        if run == self._meta_run:
            return
        self._meta_run = run
        self._on_source_changed()

    def _on_source_changed(self):
        # New data for the same figures: previews cached against the old
        # objects must not be reused (a new object can reuse an old id).
        self.bump_data_version()
        if hasattr(self, '_list'):
            self._show_study()
            self._rebuild_list()
            page = self._stack.currentWidget()
            if self.isVisible() and isinstance(page, FigurePage):
                page.on_activated()

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

        self._study_label = SecondaryLabel("")
        self._study_label.setWordWrap(True)
        nav_lay.addWidget(self._study_label)
        self._show_study()

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

        # The run the DE or Meta figures draw (two runs or more).
        self._run_row = QWidget()
        run_lay = borderless(QHBoxLayout, self._run_row)
        run_lay.setSpacing(6)
        run_lay.addWidget(QLabel("Run:"))
        self._run_combo = QComboBox()
        self._run_combo.currentIndexChanged.connect(self._on_run_chosen)
        run_lay.addWidget(self._run_combo, 1)
        self._run_row.setVisible(False)
        nav_lay.addWidget(self._run_row)

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
            page.study_section = _SECTION_OF_GROUP[group_label]
            self._stack.addWidget(page)
            self._pages.append(page)
            self._entries.append((page.study_section, group_label, page))

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

    def _show_study(self) -> None:
        """Name the study the scRNA and DE figures show."""
        study = self.source.study()
        self._study_label.setText(f"Study: {study.name}" if study else "No study selected")
        self._study_label.setToolTip(str(study) if study else "")

    def _sync_runs(self) -> None:
        """Fill the Run selector for the section shown. DE: the DE
        workspace's run (when it holds one for this study), then each run
        saved in the study's folder. Meta: the Meta-Analysis workspace's
        run (when it holds one), then each run saved in the project. Shown
        when there is a choice."""
        choices: list[tuple[str, object]] = []
        if self._section == 'meta':
            self._meta_runs = None                       # look again
            if self._held_meta_inputs() is not None:
                choices.append(("Run in the Meta-Analysis workspace", None))
            choices += [(run.label, run) for run in self.meta_runs()]
            chosen = self._meta_run
            tip = ("The meta-analysis run the Meta figures draw: the one in the\n"
                   "Meta-Analysis workspace, or one saved in the project.")
        else:
            live = self.source._live_de()
            if live is not None and getattr(live, 'de_results', None) is not None:
                choices.append(("Run in the DE workspace", None))
            choices += [(run.label, run) for run in self.source.de_runs()]
            chosen = self.source.de_run()
            tip = ("The DE run the DE figures draw: the one in the DE workspace, or\n"
                   "one saved in the study's folder (one per cell type and method).")
        self._run_combo.setToolTip(tip)
        self._run_combo.blockSignals(True)
        self._run_combo.clear()
        for label, run in choices:
            self._run_combo.addItem(label, run)
        index = next((i for i, (_l, run) in enumerate(choices) if run == chosen), 0)
        self._run_combo.setCurrentIndex(index)
        self._run_combo.blockSignals(False)
        self._run_row.setVisible(self._section in ('de', 'meta') and len(choices) > 1)

    def _on_run_chosen(self, index: int) -> None:
        if index < 0:
            return
        run = self._run_combo.itemData(index)
        if self._section == 'meta':
            self._set_meta_run(run)
        else:
            self.source.set_de_run(run)

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
        self._sync_runs()
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
                row.setToolTip(page.unavailable_reason())
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
