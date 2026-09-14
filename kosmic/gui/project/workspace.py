# Project Workspace.
#
# The entry point for the app: pick or create a project folder, add
# studies, get data into them, optionally build a shared atlas, and
# check everything is ready for analysis. Five workflow steps (ADR-002)
# drive a step-specific sidebar over one content pane: the studies
# table + details panel, or the recent-projects view when no project
# is open. Step completion comes from 'kosmic.project_workflow', which
# reads what is on disk.
#
# Pure UI; AppWindow is the only consumer of its signals and owns the
# explorer pane that lists the steps.
from __future__ import annotations

from datetime import datetime
from pathlib import Path

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor, QCursor, QFont
from PyQt6.QtWidgets import (
    QAbstractItemView, QFileDialog, QFrame, QHBoxLayout, QHeaderView,
    QInputDialog, QLabel, QLineEdit, QListWidget, QListWidgetItem, QMenu,
    QMessageBox, QPushButton, QStackedWidget, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget,
)

from kosmic.gui.shared.icon_provider import make_icon
from kosmic.gui.shared.run_worker import run_worker
from kosmic.gui.shared.widgets.base_worker import BaseWorker
from kosmic.gui.shared.theme import get_color
from kosmic.gui.shared.widgets import (
    CaptionLabel, HeaderLabel, HintLabel, InfoPanel, PrimaryButton,
    SecondaryButton, SecondaryLabel, SectionHeader, SidebarPage,
)
from kosmic.manifest import study_overview
from kosmic.paths import MASTER_ACCESSION, list_studies
from kosmic.project_workflow import (
    STEP_ATLAS, STEP_PROJECT, STEP_REVIEW, STEPS,
    ProjectProgress, project_progress, study_kind,
)


_STATUS_COLUMNS = (
    ("raw_data",   "Raw"),
    ("processed",  "Processed"),
    ("de_results", "DE"),
    ("pseudobulk", "Pseudobulk"),
)

# Deterministic per-row identity dots for the study table (cycled).
_DOT_COLORS = ('#9b59b6', '#2ecc71', '#3498db', '#f39c12',
               '#1abc9c', '#e84393', '#e74c3c', '#16a085')


def _dot_icon(color_hex: str):
    """Small filled-circle icon used as a study identity marker."""
    from PyQt6.QtGui import QIcon, QPainter, QPixmap
    pm = QPixmap(12, 12)
    pm.fill(Qt.GlobalColor.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setBrush(QColor(color_hex))
    p.setPen(Qt.PenStyle.NoPen)
    p.drawEllipse(1, 1, 10, 10)
    p.end()
    return QIcon(pm)


class _GeneOverlapWorker(BaseWorker):
    """Measure how many gene names the project's studies have in common.

    Reads each study's ``var`` index only -- no expression data -- so this
    stays cheap even where the matrices run to tens of gigabytes. Runs off
    the GUI thread anyway because it touches every study's file.
    """

    def __init__(self, project_dir, parent=None):
        super().__init__(parent)
        self._project_dir = str(project_dir)

    def _run(self):
        from kosmic.scrna.inspect.gene_names import (
            load_hgnc_lookup, load_locus_groups,
        )
        from kosmic.scrna.load.gene_overlap import project_gene_overlap

        self.progress.emit("Reading gene names...")
        lookup, _ = load_hgnc_lookup()
        groups = load_locus_groups()
        protein_coding = {sym for sym, group in groups.items()
                          if group == "protein-coding gene"}
        return project_gene_overlap(
            self._project_dir, lookup=lookup, protein_coding=protein_coding)


class _ProjectPage(SidebarPage):
    """Sidebar (per-step panels) + content (welcome | studies)."""

    SIDEBAR_WIDTH = 300
    SIDEBAR_MIN_WIDTH = 260


class ProjectWorkspace(QWidget):
    """KOSMIC Project hub: five workflow steps over the studies table."""

    WORKFLOW_STEPS = STEPS
    # F1 opens the help page for the step you are on (AppWindow reads
    # help_id off the widget; a property serves it per step).
    _STEP_HELP = ("project/project", "project/studies", "project/data",
                  "project/atlas", "project/review")

    @property
    def help_id(self) -> str:
        return self._STEP_HELP[self._current_step]

    open_requested = pyqtSignal(str)        # absolute path to open
    create_requested = pyqtSignal(str)      # absolute path to scaffold
    recent_project_clicked = pyqtSignal(str)  # absolute path from recents list
    study_selected = pyqtSignal(str)        # accession picked from study list
    studies_add_requested = pyqtSignal(list)  # accession names to scaffold
    download_requested = pyqtSignal()       # jump to scRNA Load Data tab
    refresh_requested = pyqtSignal()
    study_action_requested = pyqtSignal(str, str)  # (accession, action). action in:
    # 'set_active', 'download_raw', 'process', 'reprocess', 'open_folder',
    # 'remove', 'propagate_labels'.
    combine_requested = pyqtSignal()        # open the Combine modal
    open_in_requested = pyqtSignal(str, str)  # (accession, mode) mode in
    # 'scrna' | 'de' | 'meta'; AppWindow activates the study + workspace.
    import_requested = pyqtSignal(str, str)  # (accession, route) for the
    # study selected in the table; route in 'h5ad' | 'rds' | 'csv' | 'geo'
    # | 'assemble'. AppWindow makes that study active, copies the file into
    # raw_data/ and converts it into processed_data/ from here, so a study
    # is analysis-ready without a detour through another workspace.
    step_changed = pyqtSignal(int)      # the selected workflow step
    progress_changed = pyqtSignal()     # step completion / status recomputed

    _CONTENT_WELCOME = 0
    _CONTENT_STUDIES = 1

    def __init__(self, last_directory: str = "", parent=None):
        super().__init__(parent)
        self._last_directory = last_directory or str(Path.home())
        self._project_dir: Path | None = None
        self._active_accession: str | None = None
        self._details_accession: str | None = None
        self._row_kinds: dict[str, tuple[str, str | None]] = {}
        self._progress = ProjectProgress()
        self._current_step = STEP_PROJECT
        self._recent_paths: list[str] = []
        self._build_ui()

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        self._page = _ProjectPage(self)
        outer.addWidget(self._page)

        # Sidebar: one panel per workflow step, swapped by on_sidebar_step.
        self._side_stack = QStackedWidget()
        for builder in (self._build_project_panel,
                        self._build_studies_panel,
                        self._build_data_panel,
                        self._build_atlas_panel,
                        self._build_review_panel):
            self._side_stack.addWidget(builder())
        self._page.sidebar_layout.addWidget(self._side_stack)
        self._page.sidebar_layout.addStretch()

        # Content: welcome (no project) | studies table + details.
        self._content_stack = QStackedWidget()
        self._content_stack.addWidget(self._build_welcome_pane())
        self._content_stack.addWidget(self._build_studies_pane())
        self._page.content_layout.addWidget(self._content_stack)

        self._content_stack.setCurrentIndex(self._CONTENT_WELCOME)
        self._side_stack.setCurrentIndex(STEP_PROJECT)
        self.refresh_theme()

    @staticmethod
    def _panel(title: str, hint: str) -> tuple[QWidget, QVBoxLayout]:
        """A sidebar panel with the step's header + one-line hint."""
        panel = QWidget()
        lay = QVBoxLayout(panel)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(10)
        lay.addWidget(SectionHeader(title.upper()))
        lay.addWidget(HintLabel(hint))
        # The stack is as tall as its tallest panel; without top
        # alignment a shorter panel spreads its widgets down it.
        lay.setAlignment(Qt.AlignmentFlag.AlignTop)
        return panel, lay

    # -- step 1: Project ------------------------------------------------

    def _build_project_panel(self) -> QWidget:
        panel, lay = self._panel(
            "Project", "A project is a folder with one subfolder per "
            "study. KOSMIC reopens the last one you used.")

        self._project_name_label = HeaderLabel("No project open")
        lay.addWidget(self._project_name_label)
        self._project_path_label = CaptionLabel("")
        self._project_path_label.setWordWrap(True)
        lay.addWidget(self._project_path_label)

        self._open_btn = PrimaryButton("  Open Project...")
        self._open_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._open_btn.clicked.connect(self._on_open_clicked)
        lay.addWidget(self._open_btn)
        self._create_btn = SecondaryButton("Create New Project...")
        self._create_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._create_btn.clicked.connect(self._on_create_clicked)
        lay.addWidget(self._create_btn)

        self._side_recents_header = SecondaryLabel("Recent projects")
        self._side_recents_header.hide()
        lay.addWidget(self._side_recents_header)
        self._side_recents = QWidget()
        self._side_recents_layout = QVBoxLayout(self._side_recents)
        self._side_recents_layout.setContentsMargins(0, 0, 0, 0)
        self._side_recents_layout.setSpacing(2)
        lay.addWidget(self._side_recents)
        return panel

    # -- step 2: Studies ------------------------------------------------

    def _build_studies_panel(self) -> QWidget:
        panel, lay = self._panel(
            "Studies", "Add each dataset you want to compare as a study, "
            "named by its accession. This creates the folders; data "
            "comes next.")

        add_btn = PrimaryButton("+ Add Study")
        add_btn.setToolTip("Scaffold study folders from accession IDs.")
        add_btn.clicked.connect(self._on_add_study_clicked)
        lay.addWidget(add_btn)

        self._remove_btn = SecondaryButton("Remove selected study...")
        self._remove_btn.setToolTip(
            "Remove the study selected in the table from the project.")
        self._remove_btn.clicked.connect(self._on_remove_clicked)
        lay.addWidget(self._remove_btn)

        self._studies_caption = CaptionLabel("")
        self._studies_caption.setWordWrap(True)
        lay.addWidget(self._studies_caption)
        return panel

    # -- step 3: Data ---------------------------------------------------

    def _build_data_panel(self) -> QWidget:
        panel, lay = self._panel(
            "Data", "Select a study in the table, then import or download "
            "its data. KOSMIC keeps the file in raw_data/ and converts it "
            "into the working dataset in processed_data/.")

        self._data_target_label = SecondaryLabel("")
        self._data_target_label.setWordWrap(True)
        lay.addWidget(self._data_target_label)

        import_btn = PrimaryButton("Add Data...")
        import_btn.setToolTip(
            "Add data to the selected study: import a local file, download "
            "from GEO, or assemble from parts.")
        import_menu = QMenu(import_btn)
        import_menu.addAction(
            "Import h5ad...", lambda: self._emit_import('h5ad'))
        import_menu.addAction(
            "Import Seurat RDS / Robj...",
            lambda: self._emit_import('rds'))
        import_menu.addAction(
            "Import CSV / TSV...", lambda: self._emit_import('csv'))
        import_menu.addSeparator()
        import_menu.addAction(
            "Download from GEO...", lambda: self._emit_import('geo'))
        import_menu.addAction(
            "Assemble from parts...",
            lambda: self._emit_import('assemble'))
        import_btn.setMenu(import_menu)
        self._import_btn = import_btn
        lay.addWidget(import_btn)

        self._data_missing_label = CaptionLabel("")
        self._data_missing_label.setWordWrap(True)
        lay.addWidget(self._data_missing_label)
        return panel

    # -- step 4: Shared Atlas -------------------------------------------

    def _build_atlas_panel(self) -> QWidget:
        panel, lay = self._panel(
            "Shared Atlas", "Optional. Combine the studies into one "
            "atlas, cluster and annotate it once in scRNA, then push "
            "those cell types back so they mean the same thing in every "
            "study.")

        combine_btn = SecondaryButton("Create shared atlas...")
        combine_btn.setToolTip(
            "Merge the project's studies into one master dataset for "
            "shared clustering and annotation.")
        combine_btn.clicked.connect(self.combine_requested.emit)
        lay.addWidget(combine_btn)

        # Combining is an inner join on gene names, so the atlas carries
        # only what every study has. That figure is not predictable from
        # the per-study counts -- a study on an older annotation build
        # can share barely two thirds of its genes with the rest -- so
        # show it beside the button that acts on it.
        self._overlap_label = CaptionLabel("")
        self._overlap_label.setWordWrap(True)
        self._overlap_label.hide()
        lay.addWidget(self._overlap_label)

        self._propagate_btn = SecondaryButton("Propagate labels to studies...")
        self._propagate_btn.setToolTip(
            "Copy the atlas's cell-type annotation onto each study as "
            "cell_type_atlas.")
        self._propagate_btn.clicked.connect(
            lambda: self.study_action_requested.emit(
                MASTER_ACCESSION, "propagate_labels"))
        lay.addWidget(self._propagate_btn)

        self._atlas_caption = CaptionLabel("")
        self._atlas_caption.setWordWrap(True)
        lay.addWidget(self._atlas_caption)
        return panel

    # -- step 5: Review -------------------------------------------------

    def _build_review_panel(self) -> QWidget:
        panel, lay = self._panel(
            "Review", "Meta-analysis pools studies that have DE results "
            "and pseudobulk counts. Filter the table to see who is "
            "behind.")

        self._filter_list = QListWidget()
        self._filter_list.setProperty("role", "borderless")
        self._filter_list.setFixedHeight(150)
        for key, label in (('all', "All studies"),
                           ('raw_data', "Have raw data"),
                           ('processed', "Have processed data"),
                           ('de_results', "Have DE results"),
                           ('pseudobulk', "Have pseudobulk counts")):
            item = QListWidgetItem(label)
            item.setData(Qt.ItemDataRole.UserRole, key)
            self._filter_list.addItem(item)
        self._filter_list.setCurrentRow(0)
        self._filter_list.currentRowChanged.connect(
            lambda _r: self._apply_study_search())
        lay.addWidget(self._filter_list)

        self._meta_btn = PrimaryButton("Go to Meta-Analysis")
        self._meta_btn.clicked.connect(self._on_go_to_meta)
        lay.addWidget(self._meta_btn)

        self._review_caption = CaptionLabel("")
        self._review_caption.setWordWrap(True)
        lay.addWidget(self._review_caption)
        return panel

    # -- content: welcome -----------------------------------------------

    def _build_welcome_pane(self) -> QWidget:
        page = QWidget()
        outer = QVBoxLayout(page)
        outer.setContentsMargins(40, 60, 40, 40)
        outer.setSpacing(14)
        outer.addStretch(2)

        title = QLabel("Open a project to begin")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title_font = QFont("Segoe UI", 16)
        title_font.setBold(True)
        title.setFont(title_font)
        outer.addWidget(title)

        subtitle = SecondaryLabel(
            "A KOSMIC project is a folder containing one subfolder per "
            "study. Open one, create one, or pick a recent project."
        )
        subtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        subtitle.setWordWrap(True)
        outer.addWidget(subtitle)
        outer.addSpacing(8)

        btn_row = QHBoxLayout()
        btn_row.addStretch()
        open_btn = PrimaryButton("  Open Project...")
        open_btn.setMinimumHeight(36)
        open_btn.setFixedWidth(200)
        open_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        open_btn.clicked.connect(self._on_open_clicked)
        btn_row.addWidget(open_btn)
        self._welcome_open_btn = open_btn
        create_btn = SecondaryButton("Create New Project...")
        create_btn.setMinimumHeight(36)
        create_btn.setFixedWidth(220)
        create_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        create_btn.clicked.connect(self._on_create_clicked)
        btn_row.addWidget(create_btn)
        btn_row.addStretch()
        outer.addLayout(btn_row)

        outer.addSpacing(16)

        self._recents_header = SecondaryLabel("Recent projects")
        self._recents_header.setAlignment(Qt.AlignmentFlag.AlignCenter)
        rh_font = self._recents_header.font()
        rh_font.setBold(True)
        self._recents_header.setFont(rh_font)
        self._recents_header.hide()
        outer.addWidget(self._recents_header)

        self._recents_container = QWidget()
        self._recents_layout = QVBoxLayout(self._recents_container)
        self._recents_layout.setContentsMargins(0, 4, 0, 0)
        self._recents_layout.setSpacing(2)
        self._recents_layout.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        outer.addWidget(self._recents_container)

        outer.addStretch(3)
        return page

    # -- content: studies table + details ------------------------------

    def _build_studies_pane(self) -> QWidget:
        right = QWidget()
        layout = QVBoxLayout(right)
        layout.setContentsMargins(16, 20, 24, 16)
        layout.setSpacing(12)

        head_row = QHBoxLayout()
        head_row.setSpacing(10)
        title_col = QVBoxLayout()
        title_col.setSpacing(2)
        as_title = QLabel("Studies")
        as_title.setProperty("role", "page_header")
        title_col.addWidget(as_title)
        title_col.addWidget(CaptionLabel(
            "The datasets in this project. Each study is processed on its "
            "own, then the results are pooled in Meta-Analysis."))
        head_row.addLayout(title_col)
        head_row.addStretch()
        self._refresh_btn = SecondaryButton("Refresh")
        self._refresh_btn.clicked.connect(self.refresh_requested.emit)
        head_row.addWidget(self._refresh_btn)
        self._search_edit = QLineEdit()
        self._search_edit.setPlaceholderText("Search studies...")
        self._search_edit.setFixedWidth(220)
        self._search_edit.textChanged.connect(self._apply_study_search)
        head_row.addWidget(self._search_edit)
        layout.addLayout(head_row)

        # Columns: Study | Active | Cells | Genes | asset flags |
        # Conditions | Modified. Everything you can do to a row lives in
        # the details panel below or the step sidebar; the table itself
        # is for choosing. The blank-headed column carries the "Active"
        # chip, kept out of the Study column because its text drives
        # search and accession lookup.
        self._EXTRA_COLUMNS = ("", "Cells", "Genes")
        self._TAIL_COLUMNS = ("Conditions", "Modified")
        n_cols = 1 + len(self._EXTRA_COLUMNS) + len(_STATUS_COLUMNS) \
            + len(self._TAIL_COLUMNS)
        self._study_table = QTableWidget(0, n_cols)
        self._study_table.setHorizontalHeaderLabels(
            ["Study"] + list(self._EXTRA_COLUMNS)
            + [label for _, label in _STATUS_COLUMNS]
            + list(self._TAIL_COLUMNS)
        )
        self._study_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
        )
        self._study_table.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection
        )
        self._study_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self._study_table.verticalHeader().setVisible(False)
        self._study_table.verticalHeader().setDefaultSectionSize(32)
        self._study_table.itemSelectionChanged.connect(
            self._on_study_row_selected)
        header = self._study_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for i in range(1, n_cols):
            header.setSectionResizeMode(i, QHeaderView.ResizeMode.ResizeToContents)
        layout.addWidget(self._study_table, 1)

        # Inline study-details panel (Study Manager view, ADR-001)
        layout.addWidget(self._build_study_details_panel())

        self._empty_studies_label = SecondaryLabel(
            "No studies yet. Use '+ Add Study' in the Studies step to "
            "create a folder for each dataset."
        )
        self._empty_studies_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self._empty_studies_label)
        self._empty_studies_label.hide()
        return right

    def _build_study_details_panel(self) -> QWidget:
        panel = QFrame()
        panel.setObjectName("info_panel")
        self._details_panel = panel
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(16, 12, 16, 14)
        layout.setSpacing(10)

        head = QHBoxLayout()
        head.setSpacing(10)
        self._details_dot = QLabel()
        head.addWidget(self._details_dot)
        self._details_name = QLabel("")
        self._details_name.setObjectName("hero_dataset_name")
        head.addWidget(self._details_name)
        # What kind of folder this is: study, atlas or subset. They
        # look identical in the table and on disk, and a subset row
        # sitting beside its parent is the first thing that confuses
        # a new user.
        self._details_kind = CaptionLabel("")
        head.addWidget(self._details_kind)
        self._details_active_badge = QLabel("Active")
        self._details_active_badge.setProperty("role", "active_badge")
        self._details_active_badge.hide()
        head.addWidget(self._details_active_badge)
        head.addStretch()
        folder_btn = SecondaryButton("Open folder")
        folder_btn.setToolTip("Show this study's folder in the file manager.")
        folder_btn.clicked.connect(
            lambda: self.study_action_requested.emit(
                self._details_accession or "", "open_folder"))
        head.addWidget(folder_btn)
        for label, mode in (("Open in scRNA", "scrna"),
                            ("Open in DE", "de"),
                            ("Open in Meta", "meta")):
            btn = SecondaryButton(label)
            btn.clicked.connect(
                lambda _c=False, m=mode:
                self.open_in_requested.emit(self._details_accession or "", m))
            head.addWidget(btn)
        layout.addLayout(head)

        panels = QHBoxLayout()
        panels.setSpacing(10)

        self._ov_panel = InfoPanel("STUDY OVERVIEW")
        self._ov_panel.add_row('id', "Study ID")
        self._ov_panel.add_row('dataset', "Dataset")
        self._ov_panel.add_row('modified', "Last modified")
        self._ov_panel.add_row('status', "Status")
        panels.addWidget(self._ov_panel, 1)

        self._assets_panel = InfoPanel("DATA ASSETS")
        self._assets_panel.add_row('raw', "Raw data")
        self._assets_panel.add_row('processed', "Processed data")
        self._assets_panel.add_row('de', "DE results")
        self._assets_panel.add_row('pseudobulk', "Pseudobulk")
        panels.addWidget(self._assets_panel, 1)

        self._mapping_panel = InfoPanel("SAMPLE & CONDITION MAPPING")
        self._mapping_panel.add_row('sample', "Sample column")
        self._mapping_panel.add_row('condition', "Condition column")
        self._mapping_panel.add_row('control', "Control")
        self._mapping_panel.add_row('disease', "Disease")
        self._mapping_panel.add_row('excluded', "Exclusions")
        panels.addWidget(self._mapping_panel, 1)

        self._next_panel = InfoPanel("NEXT STEP")
        self._next_panel.add_row('verdict', "Status")
        self._next_panel.add_row('missing', "Missing")
        self._next_panel.add_row('action', "Do next")
        panels.addWidget(self._next_panel, 1)

        layout.addLayout(panels)
        panel.hide()
        return panel

    # ------------------------------------------------------------------
    # Workflow steps (AppWindow's explorer pane talks to these)
    # ------------------------------------------------------------------

    @property
    def selected_accession(self) -> str | None:
        """The study highlighted in the table (shown in the details panel)."""
        return self._details_accession

    @property
    def current_step(self) -> int:
        return self._current_step

    def on_sidebar_step(self, index: int) -> None:
        """Explorer step clicked: swap the sidebar to that step's panel."""
        if self._project_dir is None:
            index = STEP_PROJECT   # nothing else applies without a project
        index = max(STEP_PROJECT, min(STEP_REVIEW, index))
        self._current_step = index
        self._side_stack.setCurrentIndex(index)
        self.step_changed.emit(index)

    def iter_completed_steps(self):
        """Yield ``(step_index, done)`` for every workflow step."""
        return enumerate(self._progress.complete)

    def compute_sidebar_status(self) -> tuple[str, str]:
        """``(text, state)`` for the explorer STATUS panel."""
        return self._progress.status

    # ------------------------------------------------------------------
    # Public API: state switching + refresh
    # ------------------------------------------------------------------

    def show_empty(self) -> None:
        self._project_dir = None
        self._project_name_label.setText("No project open")
        self._project_path_label.setText("")
        self._content_stack.setCurrentIndex(self._CONTENT_WELCOME)
        self._progress = project_progress(None)
        self._refresh_step_panels()
        self.on_sidebar_step(STEP_PROJECT)
        self.progress_changed.emit()

    def show_project(self, project_dir: Path, active_accession: str | None) -> None:
        """Render ``project_dir`` with optional active study.

        On a change of project the selected step jumps to the first
        thing left to do; a refresh of the same project keeps the step
        the user is on.
        """
        new_project = self._project_dir != Path(project_dir)
        self._project_dir = Path(project_dir)
        self._render_header(project_dir)
        self._render_study_table(project_dir, active_accession)
        self._content_stack.setCurrentIndex(self._CONTENT_STUDIES)
        if new_project:
            self.on_sidebar_step(self._first_incomplete_step())
        self.progress_changed.emit()
        self._start_gene_overlap(project_dir)

    def _first_incomplete_step(self) -> int:
        # The atlas step is optional, so it never holds the cursor.
        for i, done in enumerate(self._progress.complete):
            if not done and i != STEP_ATLAS:
                return i
        return STEP_REVIEW

    def update_active_study(self, active_accession: str | None) -> None:
        """Refresh the active-row markers (bold font + details badge).

        Also moves the row selection onto the newly active study, so
        "Set as active" visibly does something and the details panel below
        follows it. Selection stays a separate idea -- clicking another row
        still inspects that study without changing which one is active --
        but the two coincide until you deliberately look elsewhere.

        Skips the filesystem scan that ``show_project`` runs, so this is
        cheap to call on every study-click.
        """
        if self._project_dir is None:
            return

        self._active_accession = active_accession
        active_font = QFont()
        active_font.setBold(True)
        plain_font = QFont()
        plain_font.setBold(False)

        self._study_table.blockSignals(True)
        for row in range(self._study_table.rowCount()):
            item = self._study_table.item(row, 0)
            if item is None:
                continue
            accession = item.data(Qt.ItemDataRole.UserRole) or item.text()
            item.setFont(active_font if accession == active_accession
                         else plain_font)
            self._study_table.setItem(
                row, 1, self._active_chip(accession == active_accession))
        self._study_table.blockSignals(False)

        # Move the highlight onto the active row: without this "Set as
        # active" only toggles a bold font, which reads as nothing having
        # happened. _select_study_row triggers the details panel to follow.
        if active_accession:
            self._select_study_row(active_accession)

        if self._details_accession:
            self._details_active_badge.setVisible(
                self._details_accession == active_accession)
            self._ov_panel.set_value(
                'status',
                "Active" if self._details_accession == active_accession
                else "Inactive",
                ok=True if self._details_accession == active_accession
                else None)
        self._refresh_data_panel()

    def refresh_theme(self):
        icon = make_icon("folder", "#ffffff", 16)
        self._open_btn.setIcon(icon)
        self._welcome_open_btn.setIcon(icon)

    def set_last_directory(self, directory: str) -> None:
        """Remember a directory to default the next file dialog to."""
        self._last_directory = directory or str(Path.home())

    def set_recent_projects(self, paths) -> None:
        """Render the recent-projects lists (welcome pane + Project step)."""
        valid = [str(p) for p in paths if p and Path(p).is_dir()]
        self._recent_paths = valid
        for header, layout in ((self._recents_header, self._recents_layout),
                               (self._side_recents_header,
                                self._side_recents_layout)):
            while layout.count():
                item = layout.takeAt(0)
                w = item.widget()
                if w is not None:
                    w.deleteLater()
            header.setVisible(bool(valid))
            for path in valid:
                btn = QPushButton(Path(path).name)
                btn.setToolTip(path)
                btn.setFlat(True)
                btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
                btn.clicked.connect(
                    lambda _c=False, p=path:
                    self.recent_project_clicked.emit(p))
                layout.addWidget(btn)

    # ------------------------------------------------------------------
    # Per-step sidebar state
    # ------------------------------------------------------------------

    def _refresh_step_panels(self) -> None:
        p = self._progress
        # Studies
        kinds = [k for k, _ in self._row_kinds.values()]
        parts = [f"{p.n_studies} {'study' if p.n_studies == 1 else 'studies'}"]
        n_sub = kinds.count("subset")
        if n_sub:
            parts.append(f"{n_sub} {'subset' if n_sub == 1 else 'subsets'}")
        if "atlas" in kinds:
            parts.append("shared atlas")
        self._studies_caption.setText(
            " · ".join(parts) if p.project_open else "")
        self._remove_btn.setEnabled(bool(self._row_kinds))
        # Data
        self._refresh_data_panel()
        # Atlas
        self._propagate_btn.setEnabled(p.has_atlas)
        if not p.project_open:
            self._atlas_caption.setText("")
        elif p.has_atlas:
            self._atlas_caption.setText(
                "Atlas built. Cluster and annotate _master in scRNA, then "
                "propagate its labels.")
        elif p.n_studies < 2:
            self._atlas_caption.setText(
                "Needs at least two studies with processed data.")
        else:
            self._atlas_caption.setText("No shared atlas yet.")
        # Review
        self._meta_btn.setEnabled(p.complete[STEP_REVIEW])
        if p.project_open:
            self._review_caption.setText(
                f"{p.n_ready} of {p.n_rows} ready for meta-analysis.")
        else:
            self._review_caption.setText("")

    def _refresh_data_panel(self) -> None:
        p = self._progress
        selected = self._details_accession
        if not p.project_open:
            self._data_target_label.setText("")
            self._data_missing_label.setText("")
            self._import_btn.setEnabled(False)
            return
        if selected:
            self._data_target_label.setText(f"Selected study: {selected}")
            self._import_btn.setEnabled(True)
        else:
            self._data_target_label.setText(
                "Select a study in the table to add data to it.")
            self._import_btn.setEnabled(False)
        if p.missing_data:
            names = ", ".join(p.missing_data)
            self._data_missing_label.setText(f"Without data: {names}")
        elif p.n_studies:
            self._data_missing_label.setText("Every study has data.")
        else:
            self._data_missing_label.setText("")

    # ------------------------------------------------------------------
    # Details panel
    # ------------------------------------------------------------------

    def _on_study_row_selected(self):
        rows = self._study_table.selectionModel().selectedRows()
        if not rows:
            return
        item = self._study_table.item(rows[0].row(), 0)
        if item is None:
            return
        accession = item.data(Qt.ItemDataRole.UserRole) or item.text()
        self._render_study_details(accession)
        self._refresh_data_panel()

    def _render_study_details(self, accession: str) -> None:
        project_dir = self._project_dir
        if project_dir is None or not accession:
            self._details_panel.hide()
            return
        study_dir = Path(project_dir) / accession
        overview = study_overview(study_dir)
        status = overview["status"]
        dataset = overview["dataset"] or {}
        semantics = overview["semantics"] or {}

        self._details_accession = accession
        accessions = list_studies(project_dir)
        idx = accessions.index(accession) if accession in accessions else 0
        self._details_dot.setPixmap(
            _dot_icon(_DOT_COLORS[idx % len(_DOT_COLORS)]).pixmap(12, 12))
        self._details_name.setText(accession)
        kind, parent = self._row_kinds.get(accession) or study_kind(study_dir)
        self._details_kind.setText(
            "Shared atlas" if kind == "atlas"
            else f"Subset of {parent}" if kind == "subset" else "")
        self._details_active_badge.setVisible(
            accession == self._active_accession)

        self._ov_panel.set_value('id', accession)
        if dataset:
            self._ov_panel.set_value(
                'dataset',
                f"{dataset.get('n_cells', 0):,} cells × "
                f"{dataset.get('n_genes', 0):,} genes",
                tooltip=dataset.get('file', ''))
        else:
            self._ov_panel.set_value('dataset', "not recorded")
        try:
            mtime = datetime.fromtimestamp(study_dir.stat().st_mtime)
            self._ov_panel.set_value(
                'modified', mtime.strftime("%d %b %Y %H:%M"))
        except OSError:
            self._ov_panel.set_value('modified', "?")
        is_active = accession == self._active_accession
        self._ov_panel.set_value(
            'status', "Active" if is_active else "Inactive",
            ok=True if is_active else None)

        asset_paths = {
            'raw': "raw_data/", 'processed': "processed_data/",
            'de': "results/de_analysis/", 'pseudobulk': "results/de_analysis/",
        }
        for key, flag in (('raw', status["raw_data"]),
                          ('processed', status["processed"]),
                          ('de', status["de_results"]),
                          ('pseudobulk', status["pseudobulk"])):
            self._assets_panel.set_value(
                key, "Available" if flag else "None",
                ok=True if flag else None,
                tooltip=str(study_dir / asset_paths[key]))

        role_map = semantics.get('role_map', {})
        controls = [v for v, r in role_map.items() if r == 'control']
        diseases = [v for v, r in role_map.items() if r == 'disease']
        excludes = [v for v, r in role_map.items() if r == 'exclude']

        def _sem(key, value):
            self._mapping_panel.set_value(
                key, str(value) if value else "not set",
                ok=True if value else None)

        _sem('sample', semantics.get('sample_column'))
        _sem('condition', semantics.get('condition_column'))
        _sem('control', ", ".join(controls))
        _sem('disease', ", ".join(diseases))
        self._mapping_panel.set_value(
            'excluded', ", ".join(excludes) if excludes else "None")

        self._render_next_step(kind, status, semantics)
        self._details_panel.show()

    def _render_next_step(self, kind: str, status: dict, semantics: dict) -> None:
        """The per-study version of the workflow: what this row needs next."""
        missing = [label for (key, label), flag in zip(
            _STATUS_COLUMNS,
            (status["raw_data"], status["processed"],
             status["de_results"], status["pseudobulk"])) if not flag]
        self._next_panel.set_value(
            'missing', ", ".join(missing) if missing else "None",
            ok=None if missing else True)

        if kind == "atlas":
            verdict, action, ok = ("Shared atlas", (
                "Cluster and annotate in scRNA, then propagate labels"
                if status["processed"] else "Not built"), None)
        elif not status["raw_data"] and not status["processed"]:
            verdict, action, ok = "No data", "Add Data (Data step)", False
        elif not status["processed"]:
            verdict, action, ok = ("Raw only",
                                   "Convert: Open in scRNA → Load Data", False)
        elif not semantics:
            verdict, action, ok = ("Not configured",
                                   "Open in scRNA → Inspect: set sample, "
                                   "condition and roles", False)
        elif not (status["de_results"] and status["pseudobulk"]):
            verdict, action, ok = ("Processed",
                                   "Open in DE and run differential "
                                   "expression", None)
        else:
            verdict, action, ok = ("Ready for meta-analysis",
                                   "Include in Meta-Analysis", True)
        self._next_panel.set_value('verdict', verdict, ok=ok)
        self._next_panel.set_value('action', action)

    def _apply_study_search(self):
        """Row visibility = search match AND readiness-filter match."""
        needle = self._search_edit.text().strip().lower()
        filter_item = self._filter_list.currentItem()
        filter_key = (filter_item.data(Qt.ItemDataRole.UserRole)
                      if filter_item else 'all')
        for row in range(self._study_table.rowCount()):
            item = self._study_table.item(row, 0)
            if item is None:
                continue
            name = item.text().lower()
            status = item.data(Qt.ItemDataRole.UserRole + 1) or {}
            visible = ((not needle or needle in name)
                       and (filter_key == 'all'
                            or bool(status.get(filter_key))))
            self._study_table.setRowHidden(row, not visible)

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _start_gene_overlap(self, project_dir: Path) -> None:
        """Measure the shared gene set in the background."""
        self._overlap_label.setText("Counting shared genes...")
        self._overlap_label.show()
        self._overlap_worker = run_worker(
            _GeneOverlapWorker(project_dir, parent=self),
            on_finished=self._on_gene_overlap,
            on_failed=lambda _msg: self._overlap_label.hide(),
        )

    def _on_gene_overlap(self, result: dict) -> None:
        """Show the shared gene count, and what harmonising would add.

        Two studies is the minimum for an intersection to mean anything;
        below that the panel would just restate the study's own gene count,
        so it stays hidden.
        """
        counted = [s for s in result["studies"].values()
                   if s["n_genes"] is not None]
        if len(counted) < 2:
            self._overlap_label.hide()
            return

        lines = [f"{result['shared']:,} genes shared across "
                 f"{len(counted)} studies"]
        pc = result.get("shared_protein_coding")
        if pc is not None:
            lines.append(f"{pc:,} protein-coding")

        harmonised = result.get("harmonised") or {}
        gain = harmonised.get("shared", 0) - result["shared"]
        if gain > 0:
            lines.append(f"+{gain:,} more if gene names are harmonised")

        unreadable = [name for name, s in result["studies"].items()
                      if s["error"]]
        if unreadable:
            lines.append(f"({len(unreadable)} not readable: "
                         f"{', '.join(sorted(unreadable))})")

        self._overlap_label.setText("\n".join(lines))
        self._overlap_label.setToolTip("\n".join(
            f"{name}: {s['n_genes']:,} genes"
            + (f", {s['n_protein_coding']:,} protein-coding"
               if s["n_protein_coding"] is not None else "")
            for name, s in sorted(result["studies"].items())
            if s["n_genes"] is not None))
        self._overlap_label.show()

    def _render_header(self, project_dir: Path) -> None:
        self._project_name_label.setText(project_dir.name)
        self._project_name_label.setToolTip(str(project_dir))
        try:
            mtime = datetime.fromtimestamp(project_dir.stat().st_mtime)
            mtime_str = mtime.strftime("%d %b %Y %H:%M")
        except OSError:
            mtime_str = "unknown"
        self._project_path_label.setText(
            f"{project_dir}\nLast modified {mtime_str}")

    def _render_study_table(
        self, project_dir: Path, active_accession: str | None
    ) -> None:
        accessions = list_studies(project_dir)
        self._active_accession = active_accession

        self._study_table.blockSignals(True)
        self._study_table.setRowCount(0)

        active_font = QFont()
        active_font.setBold(True)

        rows: dict[str, dict] = {}
        self._row_kinds = {}
        for i, accession in enumerate(accessions):
            row = self._study_table.rowCount()
            self._study_table.insertRow(row)
            is_active = accession == active_accession

            study_item = QTableWidgetItem(accession)
            study_item.setData(Qt.ItemDataRole.UserRole, accession)
            study_item.setIcon(_dot_icon(_DOT_COLORS[i % len(_DOT_COLORS)]))
            if is_active:
                study_item.setFont(active_font)
            self._study_table.setItem(row, 0, study_item)

            study_dir = project_dir / accession
            overview = study_overview(study_dir)
            status = overview["status"]
            study_item.setData(Qt.ItemDataRole.UserRole + 1, dict(status))
            kind = study_kind(study_dir)
            self._row_kinds[accession] = kind
            rows[accession] = {"kind": kind[0], "status": status}
            if kind[0] == "atlas":
                study_item.setToolTip("Shared atlas built by Combine studies")
            elif kind[0] == "subset":
                study_item.setToolTip(f"Cell-type subset of {kind[1]}")
            dataset = overview["dataset"] or {}
            semantics = overview["semantics"] or {}

            self._study_table.setItem(row, 1, self._active_chip(is_active))

            cells = (f"{dataset.get('n_cells', 0):,}"
                     if dataset.get('n_cells') else "")
            genes = (f"{dataset.get('n_genes', 0):,}"
                     if dataset.get('n_genes') else "")
            for col_idx, text in ((2, cells), (3, genes)):
                cell = QTableWidgetItem(text)
                cell.setTextAlignment(Qt.AlignmentFlag.AlignRight
                                      | Qt.AlignmentFlag.AlignVCenter)
                self._study_table.setItem(row, col_idx, cell)

            offset = 1 + len(self._EXTRA_COLUMNS)
            for col_idx, (key, _label) in enumerate(_STATUS_COLUMNS,
                                                    start=offset):
                cell = QTableWidgetItem("✓" if status[key] else "")
                cell.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self._study_table.setItem(row, col_idx, cell)

            role_map = semantics.get('role_map', {})
            roles = sorted({r for r in role_map.values() if r != 'exclude'})
            conditions = (f"{len(roles)} ({', '.join(roles)})"
                          if roles else "")
            try:
                mtime = datetime.fromtimestamp(study_dir.stat().st_mtime)
                modified = mtime.strftime("%d %b %Y %H:%M")
            except OSError:
                modified = ""
            tail_start = offset + len(_STATUS_COLUMNS)
            self._study_table.setItem(
                row, tail_start, QTableWidgetItem(conditions))
            self._study_table.setItem(
                row, tail_start + 1, QTableWidgetItem(modified))


        self._study_table.blockSignals(False)
        self._empty_studies_label.setVisible(not accessions)
        self._apply_study_search()

        # Readiness filter counts (Review step)
        counts = {'all': len(accessions), 'raw_data': 0, 'processed': 0,
                  'de_results': 0, 'pseudobulk': 0}
        for r in rows.values():
            for key in ('raw_data', 'processed', 'de_results', 'pseudobulk'):
                if r["status"].get(key):
                    counts[key] += 1
        base_labels = {'all': "All studies", 'raw_data': "Have raw data",
                       'processed': "Have processed data",
                       'de_results': "Have DE results",
                       'pseudobulk': "Have pseudobulk counts"}
        for idx in range(self._filter_list.count()):
            fitem = self._filter_list.item(idx)
            key = fitem.data(Qt.ItemDataRole.UserRole)
            fitem.setText(f"{base_labels[key]}  ({counts[key]})")

        self._progress = project_progress(rows)
        self._refresh_step_panels()

        # Select the active (or first) row so the details panel shows
        target = active_accession if active_accession in accessions else (
            accessions[0] if accessions else None)
        if target:
            self._study_table.selectRow(accessions.index(target))
        else:
            self._details_panel.hide()

    @staticmethod
    def _parse_accessions(raw: str) -> list[str]:
        """Parse pasted text into a clean accession list.

        Splits on newlines and commas; trims whitespace; rejects entries
        with path separators or that match a reserved name.
        """
        candidates: list[str] = []
        for chunk in raw.replace(",", "\n").splitlines():
            name = chunk.strip()
            if not name or name.startswith("."):
                continue
            if "/" in name or "\\" in name:
                continue
            if name == "meta_analysis":
                continue
            if name not in candidates:
                candidates.append(name)
        return candidates

    def _on_add_study_clicked(self):
        text, ok = QInputDialog.getMultiLineText(
            self, "Add Studies",
            "Accession IDs, one per line (e.g. GSE183852):")
        if not ok:
            return
        accessions = self._parse_accessions(text)
        if accessions:
            self.studies_add_requested.emit(accessions)

    def _emit_import(self, route: str) -> None:
        if self._details_accession:
            self.import_requested.emit(self._details_accession, route)

    def _on_remove_clicked(self):
        if self._details_accession:
            self.study_action_requested.emit(self._details_accession, "remove")

    def _on_go_to_meta(self):
        # Meta works across the project, but the AppWindow slot wants a
        # study to name; any row will do -- prefer the active one.
        accession = self._active_accession or next(iter(self._row_kinds), "")
        if accession:
            self.open_in_requested.emit(accession, "meta")

    @staticmethod
    def _active_chip(is_active: bool) -> QTableWidgetItem:
        """The 'Active' chip shown beside the study name."""
        item = QTableWidgetItem("Active" if is_active else "")
        item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        if is_active:
            font = QFont()
            font.setBold(True)
            item.setFont(font)
            item.setForeground(QColor(get_color('accent_primary')))
        return item

    def _select_study_row(self, accession: str) -> None:
        """Select a study's table row (shows it in the details panel)."""
        for row in range(self._study_table.rowCount()):
            item = self._study_table.item(row, 0)
            if item is None:
                continue
            if (item.data(Qt.ItemDataRole.UserRole) or item.text()) == accession:
                self._study_table.selectRow(row)
                return

    def _on_open_clicked(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Open KOSMIC Project", self._last_directory
        )
        if folder:
            self.open_requested.emit(folder)

    def _on_create_clicked(self):
        parent = QFileDialog.getExistingDirectory(
            self, "Choose Parent Folder for New Project", self._last_directory
        )
        if not parent:
            return
        name, ok = QInputDialog.getText(
            self, "New Project", "Project name:"
        )
        name = name.strip() if ok else ""
        if not name:
            return
        path = Path(parent) / name
        if path.exists() and any(path.iterdir()):
            reply = QMessageBox.question(
                self,
                "Folder Not Empty",
                f"'{path.name}' already exists and contains files. "
                "Use it as a KOSMIC project anyway?",
            )
            if reply != QMessageBox.StandardButton.Yes:
                return
        self.create_requested.emit(str(path))
