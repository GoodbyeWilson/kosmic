# Data Inspection & Preparation Tab
# Inspect h5ad contents, explore metadata, map standard columns, assign
# sample conditions, merge external metadata.

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QLineEdit,
    QPushButton, QTableWidgetItem, QHeaderView, QComboBox,
    QFrame, QFileDialog, QDialog,
    QListWidget, QListWidgetItem, QScrollArea,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor
from pathlib import Path

import numpy as np
import pandas as pd

from kosmic.gui.shared.theme import NoScrollComboBox, get_color
from kosmic.gui.shared.icon_provider import make_icon
from kosmic.gui.shared.widgets import (
    BaseWorker, CopyableTableWidget, SettingsGroup, StatusLabel,
    SidebarTabbedPage,
    SecondaryLabel, SecondaryButton, CaptionLabel, IconTile, InfoPanel,
    PrimaryButton, SectionHeader, StatBlock,
)
from kosmic.gui.shared import dialogs, run_worker


class _H5adSaveWorker(BaseWorker):
    """Background-thread 'write_h5ad'. 'finished_ok' payload is the path written."""

    def __init__(self, adata, output_path: str):
        super().__init__()
        self.adata = adata
        self.output_path = str(output_path)

    def _run(self):
        self.progress.emit(f"Saving to {Path(self.output_path).name}...")
        self.adata.write_h5ad(self.output_path)
        return self.output_path


class ConvertTab(SidebarTabbedPage):
    """Tab for inspecting and preparing h5ad data.

    A sidebar page rather than a bare tab strip: DATA SUMMARY and CURRENT
    MAPPING describe the dataset as a whole, so they stay in view while you
    move between Setup, Samples and Tools. Previously the summary lived in
    the Setup tab's right column and a partial copy of the mapping lived in
    the Samples tab's, which meant neither was visible from Tools and the
    mapping was stated twice.
    """

    help_id = "scrna/inspect"

    log_message = pyqtSignal(str)  # routed to OutputPanel

    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.project_dir = None
        self.adata = None
        self.summary = {}
        self.h5ad_path = None

        self.setup_ui()

    def _build_sidebar(self):
        """Dataset-level context, visible from every inner tab.

        Setup configures the mapping, Samples verifies it and Tools acts on
        it -- all three are about the same dataset, so the summary and the
        resolved mapping belong beside the tab strip rather than inside one
        of the tabs.
        """
        side = self.sidebar_layout

        side.addWidget(SectionHeader("DATA SUMMARY"))
        self._summary_panel = InfoPanel("")
        self._summary_panel.add_row('cells', "Cells")
        self._summary_panel.add_row('genes', "Genes")
        self._summary_panel.add_row('meta', "Metadata columns")
        self._summary_panel.add_row(
            'raw', "Raw counts",
            "Original counts snapshot used by DE and pathway scoring.")
        self._summary_panel.add_row(
            'sparsity', "Sparsity",
            "Fraction of zero entries. High sparsity is normal for scRNA.")
        side.addWidget(self._summary_panel)

        side.addSpacing(8)

        side.addWidget(SectionHeader("CURRENT MAPPING"))
        self._mapping_view = InfoPanel("")
        self._mapping_view.add_row('sample', "Sample")
        self._mapping_view.add_row('condition', "Condition")
        self._mapping_view.add_row('cell_type', "Cell type")
        self._mapping_view.add_row('control', "Control")
        self._mapping_view.add_row('disease', "Disease")
        self._mapping_view.add_row(
            'excluded', "Excluded",
            "Condition values set to Exclude. Their cells stay in the "
            "dataset but are left out of the disease-vs-control contrast.")
        side.addWidget(self._mapping_view)

        self._edit_link = QPushButton("Edit in Setup  ›")
        self._edit_link.setProperty("role", "link_accent")
        self._edit_link.setCursor(Qt.CursorShape.PointingHandCursor)
        self._edit_link.clicked.connect(lambda: self._tabs.setCurrentIndex(0))
        side.addWidget(self._edit_link,
                       alignment=Qt.AlignmentFlag.AlignLeft)

        side.addStretch()

    def setup_ui(self):
        self._build_sidebar()
        # SidebarTabbedPage already installs the QTabWidget on the content
        # side; '_tabs' is the name the rest of this file uses for it.
        self._tabs = self.tabs

        # ── Setup tab ─────────────────────────────────────────────
        setup_scroll = QScrollArea()
        setup_scroll.setWidgetResizable(True)
        setup_scroll.setFrameShape(QFrame.Shape.NoFrame)
        setup_page = QWidget()
        page_layout = QVBoxLayout(setup_page)
        page_layout.setContentsMargins(16, 12, 16, 12)
        page_layout.setSpacing(12)

        st_header = QHBoxLayout()
        st_col = QVBoxLayout()
        st_col.setSpacing(2)
        st_title = QLabel("Setup")
        st_title.setProperty("role", "header")
        st_col.addWidget(st_title)
        st_col.addWidget(CaptionLabel(
            "Confirm how KOSMIC should interpret this dataset before "
            "analysis."))
        st_header.addLayout(st_col)
        st_header.addStretch()
        self.save_setup_btn = PrimaryButton("Save && continue")
        self.save_setup_btn.setToolTip(
            "Write the mapping and roles to the h5ad, then continue to "
            "the Samples overview.")
        self.save_setup_btn.clicked.connect(self._save_setup)
        self.save_setup_btn.setEnabled(False)
        self.save_setup_btn.setMinimumHeight(32)
        st_header.addWidget(self.save_setup_btn,
                            alignment=Qt.AlignmentFlag.AlignTop)
        page_layout.addLayout(st_header)

        # The summary that used to sit in this tab's right column is now
        # in the sidebar, where every tab can see it.
        body = QVBoxLayout()
        body.setSpacing(12)
        self._setup_overview_section(body)
        self._setup_sample_setup_section(body)
        body.addStretch()
        page_layout.addLayout(body)

        setup_scroll.setWidget(setup_page)
        self._tabs.addTab(setup_scroll, "Setup")

        # ── Samples tab ───────────────────────────────────────────
        # Samples overview: stat tiles + grouped per-sample QC table +
        # column-mapping side panel (one row per sample; full obs is
        # per-cell -- millions of rows).
        samples_page = QWidget()
        samples_layout = QVBoxLayout(samples_page)
        samples_layout.setContentsMargins(16, 12, 16, 12)
        samples_layout.setSpacing(10)

        sm_header = QHBoxLayout()
        sm_title_col = QVBoxLayout()
        sm_title_col.setSpacing(2)
        sm_title = QLabel("Samples Overview")
        sm_title.setProperty("role", "header")
        sm_title_col.addWidget(sm_title)
        sm_title_col.addWidget(CaptionLabel(
            "Review detected samples and cell counts. Verify column "
            "mappings and edit roles if needed."))
        sm_header.addLayout(sm_title_col)
        sm_header.addStretch()
        edit_mapping_btn = SecondaryButton("Edit Column Mapping")
        edit_mapping_btn.clicked.connect(
            lambda: self._tabs.setCurrentIndex(0))
        sm_header.addWidget(edit_mapping_btn,
                            alignment=Qt.AlignmentFlag.AlignTop)
        samples_layout.addLayout(sm_header)

        tiles_row = QHBoxLayout()
        tiles_row.setSpacing(10)
        self._sample_tiles = {}
        for key, icon, caption in (("cells", "cells", "Cells"),
                                   ("genes", "dna", "Genes"),
                                   ("samples", "layers", "Samples"),
                                   ("conditions", "split", "Conditions")):
            tile = QFrame()
            tile.setObjectName("info_panel")
            tl = QHBoxLayout(tile)
            tl.setContentsMargins(14, 10, 14, 10)
            tl.setSpacing(12)
            tl.addWidget(IconTile(icon, size=40, icon_size=18))
            block = StatBlock(caption)
            tl.addWidget(block)
            tl.addStretch()
            tiles_row.addWidget(tile, 1)
            self._sample_tiles[key] = block
        samples_layout.addLayout(tiles_row)

        controls_row = QHBoxLayout()
        controls_row.setSpacing(8)
        controls_row.addWidget(QLabel("Group by"))
        self._samples_group_combo = NoScrollComboBox()
        self._samples_group_combo.addItems(["Condition", "None"])
        self._samples_group_combo.currentTextChanged.connect(
            lambda _t: self._render_samples_table())
        controls_row.addWidget(self._samples_group_combo)
        controls_row.addSpacing(8)
        controls_row.addWidget(QLabel("Show"))
        self._samples_show_combo = NoScrollComboBox()
        self._samples_show_combo.addItems(
            ["All samples", "Disease only", "Control only"])
        self._samples_show_combo.currentTextChanged.connect(
            lambda _t: self._render_samples_table())
        controls_row.addWidget(self._samples_show_combo)
        controls_row.addSpacing(8)
        self._samples_search = QLineEdit()
        self._samples_search.setPlaceholderText("Search samples...")
        self._samples_search.setFixedWidth(220)
        self._samples_search.textChanged.connect(
            lambda _t: self._render_samples_table())
        controls_row.addWidget(self._samples_search)
        controls_row.addStretch()
        export_btn = SecondaryButton("Export Table")
        export_btn.setToolTip("Save the per-sample table as CSV.")
        export_btn.clicked.connect(self._export_samples_table)
        controls_row.addWidget(export_btn)
        samples_layout.addLayout(controls_row)

        # The mapping report this tab used to carry is in the sidebar now,
        # so the table gets the full width.
        self.samples_table = CopyableTableWidget()
        self.samples_table.setAlternatingRowColors(False)
        self.samples_table.horizontalHeader().setStretchLastSection(True)
        samples_layout.addWidget(self.samples_table, 1)
        self._tabs.addTab(samples_page, "Samples")

        # ── Tools tab ─────────────────────────────────────────────
        # Pre-analysis cell filter. Sole occupant, so its content sits
        # directly on the page rather than inside a collapsible group.
        tools_scroll = QScrollArea()
        tools_scroll.setWidgetResizable(True)
        tools_scroll.setFrameShape(QFrame.Shape.NoFrame)
        tools_page = QWidget()
        tools_layout = QVBoxLayout(tools_page)
        tools_layout.setContentsMargins(16, 12, 16, 12)
        tools_layout.setSpacing(10)
        t_header = QVBoxLayout()
        t_header.setSpacing(2)
        t_title = QLabel("Tools")
        t_title.setProperty("role", "header")
        t_header.addWidget(t_title)
        t_header.addWidget(CaptionLabel(
            "Remove cells you do not want, before anything is computed on "
            "them."))
        tools_layout.addLayout(t_header)
        self._build_subset_section(tools_layout)
        tools_layout.addStretch()
        tools_scroll.setWidget(tools_page)
        self._tabs.addTab(tools_scroll, "Tools")

        # Status sits in the sidebar with the rest of the dataset context.
        self.status_label = QLabel("Select a project folder to begin")
        self.status_label.setProperty("role", "padded_status")
        self.status_label.setWordWrap(True)
        self.sidebar_layout.addWidget(self.status_label)

    @staticmethod
    def _card_title_row(number: str, title: str):
        """Numbered step badge + small-caps card title."""
        row = QHBoxLayout()
        row.setSpacing(10)
        badge = QLabel(number)
        badge.setProperty("role", "step_badge")
        badge.setFixedSize(22, 22)
        badge.setAlignment(Qt.AlignmentFlag.AlignCenter)
        row.addWidget(badge)
        lab = QLabel(title)
        lab.setProperty("role", "panel_title")
        row.addWidget(lab)
        row.addStretch()
        return row

    def _setup_overview_section(self, parent_layout):
        """Dataset mapping card: the three column designations."""
        map_card = QFrame()
        map_card.setObjectName("info_panel")
        ml = QVBoxLayout(map_card)
        ml.setContentsMargins(16, 14, 16, 16)
        ml.setSpacing(10)
        ml.addLayout(self._card_title_row("1", "DATASET MAPPING"))
        map_cap = CaptionLabel(
            "Tell KOSMIC which columns identify the sample, the condition, "
            "and the cell type.")
        map_cap.setWordWrap(True)
        ml.addWidget(map_cap)

        grid = QGridLayout()
        grid.setHorizontalSpacing(16)
        grid.setVerticalSpacing(10)
        grid.setColumnStretch(2, 1)
        self._mapping_combos = {}
        self._mapping_infos = {}
        for row, (designation, label) in enumerate((
                ('Sample', "Sample"),
                ('Condition', "Condition"),
                ('Cell type', "Cell type"))):
            grid.addWidget(QLabel(label), row, 0)
            combo = NoScrollComboBox()
            combo.setFixedWidth(370)
            combo.currentTextChanged.connect(
                lambda col, d=designation:
                self._on_mapping_combo_changed(d, col))
            grid.addWidget(combo, row, 1)
            info = CaptionLabel("")
            grid.addWidget(info, row, 2)
            self._mapping_combos[designation] = combo
            self._mapping_infos[designation] = info
        ml.addLayout(grid)

        # The three combos above are shortcuts into the full column list;
        # that list is one row per obs column (29 of them here) against the
        # Samples tab's one row per sample, so it is a different question
        # and belongs in its own window rather than expanded inline.
        browse_row = QHBoxLayout()
        self._browse_columns_btn = SecondaryButton("Browse all columns...")
        self._browse_columns_btn.setToolTip(
            "Open the full list of metadata columns: type, number of "
            "distinct values, example values, and which one is designated "
            "as sample / condition / cell type.")
        self._browse_columns_btn.clicked.connect(self._open_columns_dialog)
        browse_row.addWidget(self._browse_columns_btn)
        browse_row.addStretch()
        ml.addLayout(browse_row)

        parent_layout.addWidget(map_card)

    def _build_subset_section(self, subset_layout):
        """Build the cell filter directly onto the Tools page.

        It used to sit in a collapsible group box, which made sense when
        Tools held two things. As the only occupant, the group added a
        header repeating the tab's own title and a collapse control for
        nothing.

        Distinct from the Subset workflow step: this is pre-analysis
        triage (e.g. keep only the author's endothelial cells) that
        preserves existing embeddings and annotations; the Subset step
        later derives a new analysis dataset from computed clusters.
        """

        subset_info = SecondaryLabel(
            "Everything is kept by default. Untick what you do not want, "
            "and it is removed before anything is computed on it -- no QC, "
            "no normalisation, no clustering, and no influence on the "
            "embedding.\n\n"
            "Not the Subset step: that comes after clustering and writes a "
            "new dataset. This comes first and edits the working one.")
        subset_info.setWordWrap(True)
        subset_layout.addWidget(subset_info)

        # This deletes cells from the working h5ad and rewrites it. The
        # original is only recoverable by re-importing from raw_data/, so
        # say so before the controls rather than in the confirm dialog
        # after the user has already decided.
        subset_warning = StatusLabel(
            "Destructive: unticked cells are deleted from the working "
            "dataset and the file is rewritten. Getting them back means "
            "re-importing the study from raw_data/.",
            state='warning')
        subset_warning.setWordWrap(True)
        subset_layout.addWidget(subset_warning)

        subset_col_layout = QHBoxLayout()
        subset_col_layout.addWidget(QLabel("Filter by column:"))
        self.subset_column_combo = NoScrollComboBox()
        self.subset_column_combo.addItem("-- Select column --")
        self.subset_column_combo.currentIndexChanged.connect(self._on_subset_column_changed)
        subset_col_layout.addWidget(self.subset_column_combo, 1)
        subset_layout.addLayout(subset_col_layout)

        # Checkboxes rather than multi-select highlighting: highlighting
        # does not say which side of the line a row is on, and the
        # consequence here is deletion. A ticked row is kept.
        subset_layout.addWidget(SecondaryLabel(
            "Untick anything you want deleted. Ticked values are kept."))
        self.subset_values_list = QListWidget()
        self.subset_values_list.setSelectionMode(
            QListWidget.SelectionMode.NoSelection)
        self.subset_values_list.itemChanged.connect(
            self._refresh_subset_preview)
        subset_layout.addWidget(self.subset_values_list)

        self.subset_preview_label = SecondaryLabel("")
        self.subset_preview_label.setWordWrap(True)
        subset_layout.addWidget(self.subset_preview_label)

        subset_btn_layout = QHBoxLayout()
        self.select_all_subset_btn = QPushButton("Select All")
        self.select_all_subset_btn.clicked.connect(self._select_all_subset_values)
        self.select_all_subset_btn.setEnabled(False)
        subset_btn_layout.addWidget(self.select_all_subset_btn)

        self.deselect_all_subset_btn = QPushButton("Deselect All")
        self.deselect_all_subset_btn.clicked.connect(self._deselect_all_subset_values)
        self.deselect_all_subset_btn.setEnabled(False)
        subset_btn_layout.addWidget(self.deselect_all_subset_btn)

        self.apply_subset_btn = QPushButton("Delete unticked cells")
        self.apply_subset_btn.clicked.connect(self._apply_subset)
        self.apply_subset_btn.setEnabled(False)
        self.apply_subset_btn.setToolTip(
            "Delete every cell whose value you have unticked, and rewrite "
            "the working dataset.")
        subset_btn_layout.addWidget(self.apply_subset_btn)
        subset_btn_layout.addStretch()
        subset_layout.addLayout(subset_btn_layout)

        self.subset_status_label = StatusLabel("", state='success')
        subset_layout.addWidget(self.subset_status_label)

    def _setup_sample_setup_section(self, parent_layout):
        """Condition roles card + completion banner + column expander.

        Role assignment writes adata.obs['_role'] for downstream DE / CC perm."""
        role_card = QFrame()
        role_card.setObjectName("info_panel")
        role_layout = QVBoxLayout(role_card)
        role_layout.setContentsMargins(16, 14, 16, 16)
        role_layout.setSpacing(10)
        role_layout.addLayout(self._card_title_row("2", "CONDITION ROLES"))
        role_cap = CaptionLabel(
            "Assign each condition value a role for downstream analysis.")
        role_cap.setWordWrap(True)
        role_layout.addWidget(role_cap)

        self.role_table = CopyableTableWidget()
        self.role_table.setColumnCount(4)
        self.role_table.setHorizontalHeaderLabels(
            ["Condition value", "# Samples", "# Cells", "Role"])
        for col in (0, 1, 2):
            self.role_table.horizontalHeader().setSectionResizeMode(
                col, QHeaderView.ResizeMode.Stretch)
        # Fixed width: ResizeToContents ignores cell widgets and would
        # clip the role dropdowns.
        self.role_table.horizontalHeader().setSectionResizeMode(
            3, QHeaderView.ResizeMode.Fixed)
        self.role_table.setColumnWidth(3, 150)
        self.role_table.verticalHeader().setVisible(False)
        self.role_table.verticalHeader().setDefaultSectionSize(40)
        # Column 0 is renamable (see _refresh_role_table); the other
        # columns are made non-editable per item.
        self.role_table.setEditTriggers(
            CopyableTableWidget.EditTrigger.DoubleClicked
            | CopyableTableWidget.EditTrigger.EditKeyPressed)
        self.role_table.setSelectionMode(
            CopyableTableWidget.SelectionMode.SingleSelection)
        self.role_table.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.role_table.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        # Renaming a condition value has to reach the banner (so it reads
        # as unsaved) and the per-sample dropdown (so it offers the new
        # name). Guarded against firing while the table is being filled.
        self.role_table.itemChanged.connect(self._on_role_value_edited)
        role_layout.addWidget(self.role_table)
        self._role_rows = []
        self._populating_roles = False

        self.role_status_label = StatusLabel("", state='success')
        self.role_status_label.setWordWrap(True)
        role_layout.addWidget(self.role_status_label)

        # The table above assigns a role to each *value* of the condition
        # column. A dataset that arrived without one has no values to
        # assign, so the labels have to be built per sample first -- the
        # same job, one level down, and rare enough to belong behind a
        # button rather than as a second card competing for attention.
        assign_row = QHBoxLayout()
        self._sample_assign_btn = SecondaryButton("Set condition per sample...")
        self._sample_assign_btn.setToolTip(
            "Label each sample individually. Only needed when the dataset "
            "has no condition column; the labels you enter become the "
            "values this card assigns roles to.")
        self._sample_assign_btn.clicked.connect(self._open_sample_assign_dialog)
        assign_row.addWidget(self._sample_assign_btn)
        assign_row.addStretch()
        role_layout.addLayout(assign_row)

        parent_layout.addWidget(role_card)

        # Completion banner: green when the mapping is analysis-ready.
        self.mapping_banner = QFrame()
        self.mapping_banner.setProperty("role", "banner_info")
        bl = QHBoxLayout(self.mapping_banner)
        bl.setContentsMargins(14, 10, 14, 10)
        bl.setSpacing(10)
        self._banner_icon = QLabel()
        self._banner_icon.setFixedWidth(18)
        bl.addWidget(self._banner_icon)
        self._banner_text = QLabel("")
        bl.addWidget(self._banner_text)
        bl.addStretch()
        parent_layout.addWidget(self.mapping_banner)

        # Full column table, demoted behind a card-style expander bar.
        # The column browser lives in a dialog built once and reused. The
        # table itself must outlive any single showing: the Designate-as
        # combos in its last column are the mapping mechanism, and
        # '_on_mapping_combo_changed' writes through to them.
        self.meta_table = CopyableTableWidget()
        self.meta_table.setColumnCount(5)
        self.meta_table.setHorizontalHeaderLabels(
            ["Column", "Type", "Unique", "Sample Values", "Designate as"])
        hdr = self.meta_table.horizontalHeader()
        hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        hdr.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        hdr.setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)
        self.meta_table.setAlternatingRowColors(False)
        self.meta_table.verticalHeader().setVisible(False)

        self._columns_dialog = QDialog(self)
        self._columns_dialog.setWindowTitle("Metadata columns")
        self._columns_dialog.resize(900, 560)
        dlg_layout = QVBoxLayout(self._columns_dialog)
        dlg_layout.setContentsMargins(16, 16, 16, 16)
        dlg_layout.setSpacing(10)
        self._columns_caption = CaptionLabel(
            "One row per metadata column. Use 'Designate as' to tell KOSMIC "
            "which column identifies the sample, the condition and the cell "
            "type; the picks appear back on the Setup card.")
        self._columns_caption.setWordWrap(True)
        dlg_layout.addWidget(self._columns_caption)
        dlg_layout.addWidget(self.meta_table, 1)
        close_row = QHBoxLayout()
        close_row.addStretch()
        close_btn = SecondaryButton("Close")
        close_btn.clicked.connect(self._columns_dialog.accept)
        close_row.addWidget(close_btn)
        dlg_layout.addLayout(close_row)

        sample_assign_group = SettingsGroup("Sample Condition Assignment", collapsible=True, expanded=True)
        sample_assign_layout = QVBoxLayout()
        sample_assign_layout.setContentsMargins(0, 0, 0, 0)
        sample_assign_layout.setSpacing(3)
        sample_assign_group.add_layout(sample_assign_layout)

        sample_assign_info = SecondaryLabel("Assign condition labels to individual samples:")
        sample_assign_layout.addWidget(sample_assign_info)

        condition_labels_layout = QHBoxLayout()
        condition_labels_layout.addWidget(QLabel("Condition labels:"))

        self.condition_label_1 = QLineEdit("Control")
        self.condition_label_1.setPlaceholderText("e.g., Control, Sham, Healthy")
        self.condition_label_1.textChanged.connect(self._update_condition_options)
        condition_labels_layout.addWidget(self.condition_label_1)

        self.condition_label_2 = QLineEdit("Disease")
        self.condition_label_2.setPlaceholderText("e.g., Disease, MCAO, HF")
        self.condition_label_2.textChanged.connect(self._update_condition_options)
        condition_labels_layout.addWidget(self.condition_label_2)

        condition_labels_layout.addStretch()
        self._condition_labels_row = QWidget()
        self._condition_labels_row.setLayout(condition_labels_layout)
        sample_assign_layout.addWidget(self._condition_labels_row)

        # Optional "show metadata" column: surfaces dominant value per
        # sample so the user can see subtype/tissue/sex while assigning.
        meta_show_layout = QHBoxLayout()
        meta_show_layout.addWidget(QLabel("Show metadata column:"))
        self.meta_show_combo = NoScrollComboBox()
        self.meta_show_combo.addItem("(none)")
        self.meta_show_combo.setToolTip(
            "Pick an obs column to display its dominant value "
            "per sample alongside the condition assignment. "
            "Useful for sub-cohort / etiology / tissue source.")
        self.meta_show_combo.currentIndexChanged.connect(
            self._populate_sample_assignment_table)
        meta_show_layout.addWidget(self.meta_show_combo, 1)
        sample_assign_layout.addLayout(meta_show_layout)

        self.sample_assign_table = CopyableTableWidget()
        self.sample_assign_table.setColumnCount(2)
        self.sample_assign_table.setHorizontalHeaderLabels(["Sample", "Condition"])
        self.sample_assign_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.sample_assign_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        sample_assign_layout.addWidget(self.sample_assign_table, 1)

        self.apply_sample_conditions_btn = QPushButton("Apply Sample Conditions")
        self.apply_sample_conditions_btn.clicked.connect(self._apply_sample_conditions)
        self.apply_sample_conditions_btn.setEnabled(False)
        self.apply_sample_conditions_btn.setFixedHeight(28)
        sample_assign_layout.addWidget(self.apply_sample_conditions_btn, 0)

        # Built once and reused: the table's per-sample combos are read
        # back by _apply_sample_conditions, so the widgets have to outlive
        # any single showing.
        self._sample_assign_group = sample_assign_group
        self._sample_assign_dialog = QDialog(self)
        self._sample_assign_dialog.setWindowTitle("Condition per sample")
        self._sample_assign_dialog.resize(720, 560)
        sa_layout = QVBoxLayout(self._sample_assign_dialog)
        sa_layout.setContentsMargins(16, 16, 16, 16)
        sa_layout.setSpacing(10)
        self._sample_assign_hint = CaptionLabel("")
        self._sample_assign_hint.setWordWrap(True)
        sa_layout.addWidget(self._sample_assign_hint)
        sa_layout.addWidget(sample_assign_group, 1)
        sa_close_row = QHBoxLayout()
        sa_close_row.addStretch()
        sa_close = SecondaryButton("Close")
        sa_close.clicked.connect(self._sample_assign_dialog.accept)
        sa_close_row.addWidget(sa_close)
        sa_layout.addLayout(sa_close_row)
        self._update_sample_assign_visibility()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def set_project_directory(self, directory: str):
        self.project_dir = Path(directory)

    def reset_state(self):
        """Clear all cached state (called on raw data reload)."""
        self.adata = None
        self.summary = {}
        self.h5ad_path = None
        self._last_adata_version = -1
        for key in ('cells', 'genes', 'meta', 'raw', 'sparsity'):
            self._summary_panel.set_value(key, "")
        self.meta_table.setRowCount(0)
        self.samples_table.setRowCount(0)
        self.status_label.setText("Select a project folder to begin")

    def set_data(self, adata, file_path=None):
        """Receive adata broadcast from Load Data tab via ScRNAWorkspace.set_adata()."""
        from kosmic.scrna.load.converters import build_adata_summary

        self.adata = adata
        self.h5ad_path = Path(file_path) if file_path else None
        self.summary = build_adata_summary(adata, file_path or "")
        self._update_display()
        self._populate_role_rows()
        self._populate_samples_table()
        self._update_sample_assign_visibility()
        self.status_label.setText(
            f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes"
        )

    def on_tab_activated(self):
        """Called when tab becomes visible — pull latest adata from workspace."""
        ws = self.main_window
        if ws.current_adata is None:
            self.status_label.setText("No data loaded — go to Load Data first")
            return
        version = getattr(ws, '_adata_version', 0)
        data_changed = version != getattr(self, '_last_adata_version', -1)
        if data_changed:
            self._last_adata_version = version
            self.set_data(ws.current_adata, ws.current_h5ad_path)
        if self.adata is not None:
            has_sample = 'sample' in self.adata.obs.columns
            has_condition = ('condition' in self.adata.obs.columns
                             or '_role' in self.adata.obs.columns)
            if has_sample and 'condition' in self.adata.obs.columns:
                ws.mark_step_complete(2)
            # Land on the view that matters: the Samples overview once
            # the mapping is done, Setup while it still needs doing.
            if data_changed:
                self._tabs.setCurrentIndex(
                    1 if (has_sample and has_condition) else 0)

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def _update_display(self):
        """Update all display elements with loaded data."""
        if not self.summary:
            return

        s = self.summary
        self._summary_panel.set_value('cells', f"{s['n_cells']:,}")
        self._summary_panel.set_value(
            'genes', f"{s.get('n_vars', s.get('n_genes', 0)):,}")
        self._summary_panel.set_value('meta', str(len(s['obs_columns'])))
        self._summary_panel.set_value(
            'raw', "Stored" if s['has_raw'] else "Not stored",
            ok=True if s['has_raw'] else None)
        self._summary_panel.set_value(
            'sparsity',
            f"{s['sparsity']:.1%}" if s.get('sparsity') else "dense")

        # Metadata table
        self.meta_table.setRowCount(0)
        self.subset_column_combo.clear()
        self.subset_values_list.clear()
        self.meta_show_combo.blockSignals(True)
        self.meta_show_combo.clear()
        self.subset_column_combo.addItem("-- Select column --")
        self.meta_show_combo.addItem("(none)")

        for col in self._sorted_obs_columns(s['obs_summary']):
            col_info = s['obs_summary'][col]
            row = self.meta_table.rowCount()
            self.meta_table.insertRow(row)

            self.meta_table.setItem(row, 0, QTableWidgetItem(col))
            self.meta_table.setItem(row, 1, QTableWidgetItem(col_info['dtype']))
            self.meta_table.setItem(row, 2, QTableWidgetItem(str(col_info['n_unique'])))
            self.meta_table.setItem(row, 3, QTableWidgetItem(", ".join(col_info['sample_values'])))

            # Per-row "Designate as" combo: column browsing and mapping share one widget.
            designate_combo = QComboBox()
            designate_combo.addItems(['—', 'Condition', 'Sample', 'Cell type'])
            designate_combo.currentTextChanged.connect(
                lambda role, c=col: self._on_designate_changed(c, role))
            self.meta_table.setCellWidget(row, 4, designate_combo)

            self.subset_column_combo.addItem(col)
            self.meta_show_combo.addItem(col)
        self.meta_show_combo.blockSignals(False)

        self._auto_detect_designations()
        self.save_setup_btn.setEnabled(True)
        self._browse_columns_btn.setText(
            f"Browse all {len(self.summary['obs_columns'])} columns...")
        self._refresh_mapping_editor()
        self._populate_sample_assignment_table()

    def _auto_detect_designations(self):
        """
        Set the per-row Designate-as combos by name heuristics.
        Standardised columns ('condition' / 'sample' / 'cell_type')
        win when present; otherwise substring-match against common
        biological-column patterns.
        """
        if not self.summary or self.adata is None:
            return

        obs_cols = list(self.summary['obs_summary'].keys())
        cols_lower = {c.lower(): c for c in obs_cols}
        designations = {}

        for std_name, role in (('condition', 'Condition'),
                               ('sample', 'Sample'),
                               ('cell_type', 'Cell type')):
            if std_name in obs_cols:
                designations[role] = std_name

        if 'Condition' not in designations:
            for pat in ('condition', 'group', 'disease', 'status', 'treatment'):
                for cl, c in cols_lower.items():
                    if pat in cl:
                        designations['Condition'] = c
                        break
                if 'Condition' in designations:
                    break
        if 'Sample' not in designations:
            for pat in ('donor_id', 'donor', 'patient', 'subject',
                        'sample_id', 'sample', 'orig.ident', 'orig_ident',
                        'biosample', 'individual', 'specimen'):
                for cl, c in cols_lower.items():
                    if pat in cl:
                        designations['Sample'] = c
                        break
                if 'Sample' in designations:
                    break
        if 'Cell type' not in designations:
            for pat in ('cell_type', 'celltype', 'cell.type',
                        'names', 'cluster', 'annotation'):
                for cl, c in cols_lower.items():
                    if pat in cl or cl == pat:
                        designations['Cell type'] = c
                        break
                if 'Cell type' in designations:
                    break

        for r in range(self.meta_table.rowCount()):
            item = self.meta_table.item(r, 0)
            combo = self.meta_table.cellWidget(r, 4)
            if item is None or combo is None:
                continue
            target = '—'
            for role, designated_col in designations.items():
                if item.text() == designated_col:
                    target = role
                    break
            combo.blockSignals(True)
            combo.setCurrentText(target)
            combo.blockSignals(False)

    @staticmethod
    def _sorted_obs_columns(obs_summary):
        """
        Sort obs columns: standardised first, then user-set bio/clinical
        patterns, then 'other', then scanpy-style computed columns last.
        Pattern lists live in 'kosmic.scrna.inspect.batch' (centralised config);
        the sort itself stays here -- it's a small UX wrapper.
        """
        from kosmic.scrna.inspect.batch import (
            STANDARDISED_OBS_NAMES, USER_OBS_PATTERNS, COMPUTED_OBS_PATTERNS,
        )
        standardised_set = set(STANDARDISED_OBS_NAMES)

        def priority(col):
            cl = str(col).lower()
            if col in standardised_set:
                return (0, cl)
            if any(p in cl for p in COMPUTED_OBS_PATTERNS):
                return (3, cl)
            if any(p in cl for p in USER_OBS_PATTERNS):
                return (1, cl)
            return (2, cl)

        return sorted(obs_summary.keys(), key=priority)

    def _on_designate_changed(self, col_name, role):
        """
        Enforce uniqueness of Condition/Sample/Cell type and live-refresh
        the Role Assignment table when Condition changes.
        """
        if role != '—':
            for r in range(self.meta_table.rowCount()):
                item = self.meta_table.item(r, 0)
                combo = self.meta_table.cellWidget(r, 4)
                if item is None or combo is None:
                    continue
                if item.text() == col_name:
                    continue
                if combo.currentText() == role:
                    combo.blockSignals(True)
                    combo.setCurrentText('—')
                    combo.blockSignals(False)
        self._populate_role_rows()
        self._refresh_mapping_editor()

    def _designated_col(self, designation: str):
        """Return the obs column currently designated as *designation*
        in the metadata table, or None."""
        for r in range(self.meta_table.rowCount()):
            item = self.meta_table.item(r, 0)
            combo = self.meta_table.cellWidget(r, 4)
            if item is None or combo is None:
                continue
            if combo.currentText() == designation:
                return item.text()
        return None

    def _refresh_mapping_editor(self):
        """Sync the Dataset Mapping combos + captions from the table
        designations (the table stays the single mapping mechanism)."""
        if self.adata is None or not hasattr(self, '_mapping_combos'):
            return
        obs = self.adata.obs
        cols = [str(c) for c in obs.columns if obs[c].nunique() <= 200]
        self._mapping_sync = True
        try:
            for designation, combo in self._mapping_combos.items():
                col = (self._designated_condition_col()
                       if designation == 'Condition'
                       else self._designated_col(designation))
                combo.blockSignals(True)
                combo.clear()
                combo.addItem("")
                combo.addItems(cols)
                if col:
                    idx = combo.findText(str(col))
                    if idx >= 0:
                        combo.setCurrentIndex(idx)
                combo.blockSignals(False)
                info = self._mapping_infos[designation]
                if col and col in obs.columns:
                    n = obs[col].nunique()
                    if designation == 'Sample':
                        info.setText(f"{n} unique samples")
                    elif designation == 'Cell type':
                        info.setText(f"{n} cell types")
                    else:
                        values = [str(v) for v in obs[col].unique()[:3]]
                        info.setText(
                            f"{n} values · " + " · ".join(values)
                            + (" ..." if n > 3 else ""))
                else:
                    info.setText("not set")
        finally:
            self._mapping_sync = False

    def _designated_condition_col(self):
        """
        Return the obs column name currently designated as
        'Condition' in the metadata table, or fall back to the
        standardised 'condition' if it already exists.
        """
        for r in range(self.meta_table.rowCount()):
            item = self.meta_table.item(r, 0)
            combo = self.meta_table.cellWidget(r, 4)
            if item is None or combo is None:
                continue
            if combo.currentText() == 'Condition':
                return item.text()
        if (self.adata is not None
                and 'condition' in self.adata.obs.columns):
            return 'condition'
        return None

    # ------------------------------------------------------------------
    # Subset
    # ------------------------------------------------------------------

    def _on_subset_column_changed(self, index):
        """Handle subset column selection - populate values list."""
        self.subset_values_list.clear()
        self.subset_status_label.setText("")

        if index <= 0 or self.adata is None:
            self.select_all_subset_btn.setEnabled(False)
            self.deselect_all_subset_btn.setEnabled(False)
            self.apply_subset_btn.setEnabled(False)
            return

        col_name = self.subset_column_combo.currentText()
        if col_name not in self.adata.obs.columns:
            return

        value_counts = self.adata.obs[col_name].value_counts()

        self.subset_values_list.blockSignals(True)
        for value, count in value_counts.items():
            item = QListWidgetItem(f"{value}  ({count:,} cells)")
            item.setData(Qt.ItemDataRole.UserRole, value)
            item.setData(Qt.ItemDataRole.UserRole + 1, int(count))
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            # Everything starts kept, so an accidental Apply is a no-op
            # rather than a deletion.
            item.setCheckState(Qt.CheckState.Checked)
            self.subset_values_list.addItem(item)
        self.subset_values_list.blockSignals(False)

        self.select_all_subset_btn.setEnabled(True)
        self.deselect_all_subset_btn.setEnabled(True)
        self._refresh_subset_preview()

    def _subset_checked_items(self):
        """List items currently ticked (i.e. to keep)."""
        return [self.subset_values_list.item(i)
                for i in range(self.subset_values_list.count())
                if self.subset_values_list.item(i).checkState()
                == Qt.CheckState.Checked]

    def _refresh_subset_preview(self):
        """Say how many cells stay and how many go, before anything runs."""
        total_items = self.subset_values_list.count()
        if not total_items or self.adata is None:
            self.subset_preview_label.setText("")
            self.apply_subset_btn.setEnabled(False)
            return

        keep = self._subset_checked_items()
        n_keep = sum(item.data(Qt.ItemDataRole.UserRole + 1) or 0
                     for item in keep)
        total = int(self.adata.n_obs)
        n_drop = total - n_keep

        if not keep:
            self.subset_preview_label.setText(
                "Everything unticked — that would delete every cell.")
            self.apply_subset_btn.setEnabled(False)
            return
        if n_drop <= 0:
            self.subset_preview_label.setText(
                f"Nothing unticked — all {total:,} cells would be kept.")
            self.apply_subset_btn.setEnabled(False)
            return

        pct = 100.0 * n_drop / total if total else 0.0
        self.subset_preview_label.setText(
            f"Keeping {n_keep:,} cells ({len(keep)} of {total_items} values). "
            f"Deleting {n_drop:,} cells ({pct:.0f}%).")
        self.apply_subset_btn.setEnabled(True)

    def _set_all_subset_checks(self, state):
        self.subset_values_list.blockSignals(True)
        for i in range(self.subset_values_list.count()):
            self.subset_values_list.item(i).setCheckState(state)
        self.subset_values_list.blockSignals(False)
        self._refresh_subset_preview()

    def _select_all_subset_values(self):
        self._set_all_subset_checks(Qt.CheckState.Checked)

    def _deselect_all_subset_values(self):
        self._set_all_subset_checks(Qt.CheckState.Unchecked)

    def _apply_subset(self):
        """Apply subset filter to keep only selected values."""
        if self.adata is None:
            return

        col_name = self.subset_column_combo.currentText()
        if col_name == "-- Select column --" or col_name not in self.adata.obs.columns:
            return

        selected_values = [item.data(Qt.ItemDataRole.UserRole)
                           for item in self._subset_checked_items()]

        if not selected_values:
            dialogs.warning(
                self, "Everything unticked",
                "Leave at least one value ticked. Unticking them all would "
                "delete every cell.")
            return

        if len(selected_values) == self.subset_values_list.count():
            dialogs.info(self, "Nothing to remove",
                         "Nothing is unticked, so no cells would be "
                         "deleted.")
            return

        original_cells = self.adata.n_obs
        mask = self.adata.obs[col_name].isin(selected_values)
        new_cells = mask.sum()

        removed = original_cells - new_cells
        pct = 100.0 * removed / max(original_cells, 1)
        kept_preview = ", ".join(str(v) for v in selected_values[:5])
        if len(selected_values) > 5:
            kept_preview += f" (+{len(selected_values) - 5} more)"
        if not dialogs.confirm(
                self, "Delete cells?",
                f"This deletes {removed:,} of {original_cells:,} cells "
                f"({pct:.0f}%), leaving {new_cells:,}." + "\n\n"
                + f"Keeping: {kept_preview}" + "\n\n"
                + "The working dataset is rewritten. Annotations and "
                "embeddings are kept for the remaining cells; the "
                "deleted cells are only recoverable by re-importing "
                "from raw_data/." + "\n\nContinue?"):
            return

        self.adata = self.adata[mask].copy()

        if self.h5ad_path:
            self.adata.write_h5ad(self.h5ad_path)
            self.log_message.emit(f"Subsetted data saved: {new_cells:,} cells")

        self.summary['n_cells'] = self.adata.n_obs
        self._summary_panel.set_value(
            'cells', f"{self.summary['n_cells']:,}")
        self.subset_status_label.setText(f"✓ Subsetted to {new_cells:,} cells")
        self.subset_status_label.set_state('success')

        self._on_subset_column_changed(self.subset_column_combo.currentIndex())
        self._populate_sample_assignment_table()

        if self.main_window is not None and self.h5ad_path:
            self.main_window.set_adata(self.adata, str(self.h5ad_path))
            # This deletes cells from the working h5ad permanently, and
            # every other data-changing step records itself -- without this
            # the study's methods would describe an analysis of cells that
            # are no longer in the file.
            if hasattr(self.main_window, 'record_provenance'):
                self.main_window.record_provenance('filter_dataset', {
                    'column': col_name,
                    'kept_values': [str(v) for v in selected_values],
                    'n_cells_before': int(original_cells),
                    'n_cells_after': int(new_cells),
                })

        dialogs.info(
            self, "Subset Applied",
            f"Data subsetted successfully!\n\n"
            f"Original: {original_cells:,} cells\n"
            f"After subset: {new_cells:,} cells\n"
            f"Removed: {original_cells - new_cells:,} cells"
        )

    # ------------------------------------------------------------------
    # Column mapping
    # ------------------------------------------------------------------

    def _save_setup(self):
        """Save column designations and role assignment to adata in one h5ad write."""
        from kosmic.scrna.inspect.roles import resolve_roles

        if self.adata is None:
            return

        # ── Step 1: resolve column designations ──────────────────────
        role_to_obscol = {
            'Condition': 'condition',
            'Sample': 'sample',
            'Cell type': 'cell_type',
        }
        designations = {}
        for r in range(self.meta_table.rowCount()):
            item = self.meta_table.item(r, 0)
            combo = self.meta_table.cellWidget(r, 4)
            if item is None or combo is None:
                continue
            role = combo.currentText()
            if role in role_to_obscol:
                designations[role] = item.text()

        changes = []
        for role, src_col in designations.items():
            std_col = role_to_obscol[role]
            if src_col not in self.adata.obs.columns:
                continue
            if std_col not in self.adata.obs.columns or src_col != std_col:
                self.adata.obs[std_col] = self.adata.obs[src_col]
                changes.append(f"{std_col} <- {src_col}")

        # ── Step 2: resolve role assignment ──────────────────────────
        cond_col_name = self._designated_condition_col()
        role_changes = []
        if cond_col_name is not None and self._role_rows:
            # A renamed condition value has to be rewritten in the obs
            # column before roles are resolved, and the role map keyed by
            # the new name -- otherwise resolve_roles matches nothing.
            renames = self._collect_condition_renames()
            if renames:
                col = self.adata.obs[cond_col_name].astype(str)
                self.adata.obs[cond_col_name] = col.replace(renames)
                role_changes.append(
                    f"{cond_col_name} renamed: "
                    + ", ".join(f"{old_v} -> {new_v}"
                                for old_v, new_v in sorted(renames.items())))

            role_map = {}
            for value, combo in self._role_rows:
                role = combo.currentText().lower()
                if role not in ('disease', 'control', 'exclude'):
                    continue
                role_map[renames.get(value, value)] = role

            n_dis = sum(1 for r in role_map.values() if r == 'disease')
            n_ctl = sum(1 for r in role_map.values() if r == 'control')
            if n_dis >= 1 and n_ctl >= 1:
                cond_values = self.adata.obs[cond_col_name].values
                is_disease, two_group = resolve_roles(cond_values, role_map)
                roles = np.where(
                    is_disease, 'disease',
                    np.where(two_group, 'control', 'exclude'))
                self.adata.obs['_role'] = pd.Categorical(
                    roles, categories=['control', 'disease', 'exclude'])
                self.adata.uns['role_map'] = dict(role_map)
                # Recorded so the DE Setup page can sync exactly, not by name pattern.
                self.adata.uns['role_condition_col'] = str(cond_col_name)
                role_changes.append(
                    f"_role ({n_dis} disease, {n_ctl} control)")

        # ── Step 3: write once ───────────────────────────────────────
        if not changes and not role_changes:
            self.status_label.setText("Nothing to save — no changes detected.")
            return

        # h5ad write is multi-second I/O -- run on a worker.
        self._pending_save_summary = changes + role_changes
        self.save_setup_btn.setEnabled(False)
        self.status_label.setText("Saving...")
        self._save_worker = _H5adSaveWorker(self.adata, self.h5ad_path)
        run_worker(
            self._save_worker,
            on_finished=self._on_save_setup_finished,
            on_failed=self._on_save_setup_failed,
            on_progress=lambda m: self.status_label.setText(m),
        )

    def _on_save_setup_finished(self, _path: str):
        self.save_setup_btn.setEnabled(True)
        summary = getattr(self, '_pending_save_summary', [])
        self.status_label.setText(
            "Saved: " + ", ".join(summary) if summary else "Saved.")
        self.main_window.set_adata(self.adata, str(self.h5ad_path))
        if hasattr(self.main_window, 'record_provenance'):
            self.main_window.record_provenance('setup', {
                'condition_col': self.adata.uns.get('role_condition_col'),
                'role_map': dict(self.adata.uns.get('role_map', {})),
                'changes': list(getattr(self, '_pending_save_summary', [])),
            })
        self._write_study_manifest()
        self.main_window.mark_step_complete(2)
        self._populate_role_rows()
        self._update_sample_assign_visibility()
        # Saving is what creates obs['_role'], and sample_overview prefers
        # it over the raw condition column -- so the Samples table built at
        # load time predates the roles it is meant to verify. Rebuild it
        # before landing there, or the control/disease colouring only
        # appears after a restart. This also refreshes the sidebar's
        # CURRENT MAPPING via _refresh_mapping_panel.
        self._populate_samples_table()
        # Contextual continue: setup saved, land on the verify view.
        self._tabs.setCurrentIndex(1)

    def _write_study_manifest(self):
        """Refresh the study manifest cache (ADR-001; this save flow is
        the manifest's only writer). Best-effort: never blocks a save."""
        if self.adata is None or not self.h5ad_path:
            return
        try:
            from kosmic.manifest import write_manifest
            h5ad = Path(self.h5ad_path)
            study_dir = (h5ad.parent.parent
                         if h5ad.parent.name == 'processed_data'
                         else h5ad.parent)
            obs_cols = self.adata.obs.columns
            write_manifest(
                study_dir,
                dataset={
                    'file': h5ad.name,
                    'n_cells': int(self.adata.n_obs),
                    'n_genes': int(self.adata.n_vars),
                },
                semantics={
                    'sample_column': 'sample' if 'sample' in obs_cols else None,
                    'condition_column': self.adata.uns.get('role_condition_col')
                                        or ('condition' if 'condition' in obs_cols else None),
                    'cell_type_column': 'cell_type' if 'cell_type' in obs_cols else None,
                    'role_map': dict(self.adata.uns.get('role_map', {})),
                },
            )
        except Exception as e:  # cache refresh must never block a save
            self.status_label.setText(f"Saved (manifest refresh failed: {e})")

    def _on_save_setup_failed(self, message: str):
        self.save_setup_btn.setEnabled(True)
        self.status_label.setText("Save failed.")
        dialogs.warning(
            self, "Save failed",
            f"Changes applied in memory but h5ad save failed:\n{message}")

    def _populate_samples_table(self):
        """Recompute the Samples overview (tiles + per-sample table) from
        'batch.sample_overview'; rendering honours the filter controls."""
        if not hasattr(self, 'samples_table'):
            return
        from kosmic.scrna.inspect.batch import sample_overview

        overview = sample_overview(self.adata)
        self._samples_df = overview['table']
        self._sample_tiles['cells'].set_value(f"{overview['n_cells']:,}")
        self._sample_tiles['genes'].set_value(f"{overview['n_genes']:,}")
        self._sample_tiles['samples'].set_value(str(overview['n_samples']))
        self._sample_tiles['conditions'].set_value(
            str(overview['n_conditions']))
        self._render_samples_table()
        self._refresh_mapping_panel()

    _SAMPLE_METRICS = (
        ('n_cells', "Cells", "{:,.0f}"),
        ('pct_doublets', "Doublets (%)", "{:.1f}"),
        ('pct_mito', "Mitochondrial (%)", "{:.1f}"),
        ('median_genes', "Median genes/cell", "{:,.0f}"),
    )

    def _samples_view_df(self):
        """The per-sample table filtered by the Show + search controls."""
        import pandas as pd
        df = getattr(self, '_samples_df', None)
        if df is None or df.empty:
            return pd.DataFrame()
        show = self._samples_show_combo.currentText()
        if 'condition' in df.columns and show != "All samples":
            wanted = 'disease' if show == "Disease only" else 'control'
            df = df[df['condition'].str.lower() == wanted]
        needle = self._samples_search.text().strip().lower()
        if needle:
            df = df[df['sample'].str.lower().str.contains(needle)]
        return df

    def _render_samples_table(self):
        """Render the filtered per-sample table, optionally grouped by
        condition with subtotal header rows and a Total footer."""
        table = self.samples_table
        table.setRowCount(0)
        df = self._samples_view_df()
        if df.empty:
            table.setColumnCount(1)
            table.setHorizontalHeaderLabels(
                ["(no sample column)" if self.adata is not None else ""])
            return

        has_condition = 'condition' in df.columns
        metrics = [(key, label, fmt) for key, label, fmt in
                   self._SAMPLE_METRICS if key in df.columns]
        headers = (["Sample", "Cells"]
                   + (["Condition"] if has_condition else [])
                   + [label for key, label, _f in metrics if key != 'n_cells'])
        table.setColumnCount(len(headers))
        table.setHorizontalHeaderLabels(headers)

        group_by_condition = (has_condition and
                              self._samples_group_combo.currentText() == "Condition")
        total_cells = int(self._samples_df['n_cells'].sum())

        def add_data_row(rec):
            row = table.rowCount()
            table.insertRow(row)
            table.setItem(row, 0, QTableWidgetItem(str(rec['sample'])))
            table.setItem(row, 1, QTableWidgetItem(f"{rec['n_cells']:,.0f}"))
            col = 2
            if has_condition:
                value = str(rec['condition'])
                low = value.lower()
                item = QTableWidgetItem(
                    value.capitalize() if low in ('disease', 'control')
                    else value)
                if low in ('disease', 'control'):
                    item.setForeground(QColor(get_color(f"plot_{low}")))
                    font = item.font()
                    font.setBold(True)
                    item.setFont(font)
                table.setItem(row, col, item)
                col += 1
            for key, _label, fmt in metrics:
                if key == 'n_cells':
                    continue
                value = rec.get(key)
                text = fmt.format(value) if pd.notna(value) else ""
                table.setItem(row, col, QTableWidgetItem(text))
                col += 1

        def add_group_header(name, sub):
            row = table.rowCount()
            table.insertRow(row)
            n_cells = int(sub['n_cells'].sum())
            pct = n_cells / total_cells * 100 if total_cells else 0
            low = str(name).lower()
            colour = (get_color('plot_disease') if low == 'disease'
                      else get_color('plot_control') if low == 'control'
                      else get_color('fg_secondary'))
            head = QTableWidgetItem(
                f"{str(name).capitalize()}  ({len(sub)} samples)")
            font = head.font()
            font.setBold(True)
            head.setFont(font)
            head.setForeground(QColor(colour))
            table.setItem(row, 0, head)
            sub_item = QTableWidgetItem(f"{n_cells:,} cells ({pct:.1f}%)")
            sub_item.setForeground(QColor(colour))
            table.setItem(row, table.columnCount() - 1, sub_item)

        if group_by_condition:
            order = sorted(
                df['condition'].unique(),
                key=lambda v: {'disease': 0, 'control': 1}.get(
                    str(v).lower(), 2))
            for name in order:
                sub = df[df['condition'] == name]
                add_group_header(name, sub)
                for _idx, rec in sub.iterrows():
                    add_data_row(rec)
        else:
            for _idx, rec in df.iterrows():
                add_data_row(rec)

        # Total footer
        row = table.rowCount()
        table.insertRow(row)
        total_item = QTableWidgetItem("Total")
        font = total_item.font()
        font.setBold(True)
        total_item.setFont(font)
        table.setItem(row, 0, total_item)
        shown = int(df['n_cells'].sum())
        pct = shown / total_cells * 100 if total_cells else 0
        shown_item = QTableWidgetItem(f"{shown:,} cells ({pct:.0f}%)")
        shown_item.setFont(font)
        table.setItem(row, 1, shown_item)
        table.resizeColumnsToContents()

    def _export_samples_table(self):
        """Save the filtered per-sample table as CSV."""
        df = self._samples_view_df()
        if df.empty:
            dialogs.info(self, "Nothing to Export",
                         "No sample rows to export yet.")
            return
        start = str(self.h5ad_path.parent) if self.h5ad_path else ""
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Samples Table", start + "/samples_overview.csv",
            "CSV Files (*.csv)")
        if path:
            df.to_csv(path, index=False)
            self.status_label.setText(f"Exported samples table: {path}")

    # ------------------------------------------------------------------
    # Column Mapping side panel (mirrors the Setup designations)
    # ------------------------------------------------------------------

    def _refresh_mapping_panel(self):
        """Fill the Samples tab's read-only mapping report from what
        Setup has established (obs columns + uns role_map)."""
        if self.adata is None or not hasattr(self, '_mapping_view'):
            return
        obs_cols = self.adata.obs.columns

        def _fill(key, value):
            self._mapping_view.set_value(
                key, str(value) if value else "not set",
                ok=True if value else None)

        _fill('sample', 'sample' if 'sample' in obs_cols else None)
        _fill('condition', self._designated_condition_col())
        _fill('cell_type', 'cell_type' if 'cell_type' in obs_cols else None)
        role_map = dict(self.adata.uns.get('role_map', {}))
        controls = [v for v, r in role_map.items() if r == 'control']
        diseases = [v for v, r in role_map.items() if r == 'disease']
        _fill('control', ", ".join(controls))
        _fill('disease', ", ".join(diseases))
        # Excluded arms are easy to set and then forget about, and nothing
        # downstream mentions them again -- their samples simply never
        # appear in a contrast. Naming them here keeps the omission visible.
        excluded = [v for v, r in role_map.items() if r == 'exclude']
        self._mapping_view.set_value(
            'excluded', ", ".join(excluded) if excluded else "none",
            ok=None if excluded else True)

    def _on_mapping_combo_changed(self, designation: str, col: str):
        """Write a Dataset Mapping pick through to the column table's
        Designate-as combo (the single mapping mechanism); Save commits
        to the h5ad."""
        if getattr(self, '_mapping_sync', False) or not col:
            return
        for r in range(self.meta_table.rowCount()):
            item = self.meta_table.item(r, 0)
            combo = self.meta_table.cellWidget(r, 4)
            if item is None or combo is None:
                continue
            if item.text() == col:
                combo.setCurrentText(designation)
                break

    def _update_sample_assign_visibility(self):
        """Say whether per-sample labelling is needed or merely available.

        When 'obs['condition']' already holds values the roles card has
        everything it needs, and labelling per sample would overwrite it --
        so the button stays reachable (one bad donor is a legitimate reason
        to relabel) but the hint says it is not required.
        """
        btn = getattr(self, '_sample_assign_btn', None)
        if btn is None:
            return

        # Resolve the *designated* condition column, not the literal name
        # 'condition': GSE292067's is 'group', and hardcoding the name told
        # the user their dataset had no condition column when it plainly did.
        needed = True
        col = self._sample_assign_condition_col()
        if col is not None:
            n_unique = (self.adata.obs[col].astype(str)
                        .replace('', pd.NA).dropna().nunique())
            needed = n_unique == 0
        self._sample_assign_col_name = col

        hint = getattr(self, '_sample_assign_hint', None)
        if needed:
            btn.setText("Set condition per sample...")
            if hint is not None:
                hint.setText(
                    "No condition column is designated yet. Name the two "
                    "groups, assign each sample to one, then Apply -- the "
                    "labels become the values you assign roles to back on "
                    "the Setup card.")
        else:
            btn.setText("Edit condition per sample...")
            if hint is not None:
                n = len(self._condition_options())
                hint.setText(
                    f"Condition column '{col}' already has {n} value(s), so "
                    "this is optional. The dropdown offers all of them; "
                    "changing one rewrites that column for every cell of "
                    "the sample. Use Condition Roles on the Setup card to "
                    "mark a whole value as excluded.")

    # ------------------------------------------------------------------
    # Role assignment
    # ------------------------------------------------------------------

    def _make_role_combo(self, current: str, enabled: bool = True):
        """Standard role dropdown; coloured text marks the chosen role."""
        combo = NoScrollComboBox()
        combo.addItems(
            ['Exclude', 'Control', 'Disease'] if enabled else ['Exclude'])
        combo.setCurrentText(current)
        combo.setEnabled(enabled)
        combo.setToolTip(
            "Excluded samples are dropped from DE and meta-analysis.")
        combo.currentTextChanged.connect(
            lambda _t: self._refresh_setup_completion())
        return combo

    def _set_role_cell(self, row: int, combo):
        self.role_table.setCellWidget(row, 3, combo)

    def _populate_role_rows(self):
        """
        Fill the roles card from the obs column currently designated as
        Condition. Pre-populates from 'adata.uns['role_map']' when present;
        otherwise guesses each value via DISEASE/CONTROL_VALUE_PATTERNS.
        """
        self._populating_roles = True
        self.role_table.setRowCount(0)
        self._role_rows = []
        self._size_role_table()
        if self.adata is None:
            self.role_status_label.setText("")
            self._populating_roles = False
            return

        cond_col_name = self._designated_condition_col()
        if cond_col_name is None or cond_col_name not in self.adata.obs.columns:
            self.role_status_label.setText("")
            self._populating_roles = False
            self._refresh_setup_completion()
            return
        col = self.adata.obs[cond_col_name]
        # Lookup is case-insensitive; canonical case is preserved.
        existing_map = {}
        uns_map = self.adata.uns.get('role_map') if hasattr(self.adata, 'uns') else None
        if isinstance(uns_map, dict):
            existing_map = {str(k).lower(): str(v) for k, v in uns_map.items()
                            if v in ('disease', 'control', 'exclude')}

        col_str = col.astype(str)
        unique_vals = sorted(col.dropna().astype(str).unique().tolist())
        if not unique_vals:
            self.role_status_label.setText(
                "Condition column has no values to assign.")
            self._populating_roles = False
            self._refresh_setup_completion()
            return

        # Cap rows; aggregate the rest into a single "other" row.
        WARN_THRESHOLD = 20
        too_many = len(unique_vals) > WARN_THRESHOLD
        if too_many:
            counts = col_str.value_counts()
            top_vals = counts.head(WARN_THRESHOLD).index.tolist()
            other_n = int(counts.iloc[WARN_THRESHOLD:].sum())
            unique_vals = top_vals
        else:
            other_n = 0

        from kosmic.scrna.inspect.roles import (
            CONTROL_VALUE_PATTERNS, DISEASE_VALUE_PATTERNS,
        )

        def _guess_role(value: str) -> str:
            v = value.lower().strip()
            for p in CONTROL_VALUE_PATTERNS:
                if v == p or p in v:
                    return 'control'
            for p in DISEASE_VALUE_PATTERNS:
                if v == p or p in v:
                    return 'disease'
            return 'exclude'

        sample_col = self._designated_col('Sample') or (
            'sample' if 'sample' in self.adata.obs.columns else None)

        unrecognised = []
        rows = []
        for val in unique_vals:
            mask = col_str == val
            n_cells = int(mask.sum())
            if sample_col:
                n_samples = int(self.adata.obs.loc[mask, sample_col].nunique())
                samples_text = f"{n_samples} samples"
            else:
                samples_text = "—"
            existing_role = existing_map.get(val.lower())
            if existing_role is None:
                existing_role = _guess_role(val)
                if existing_role == 'exclude':
                    unrecognised.append(val)
            rows.append((val, samples_text, f"{n_cells:,} cells",
                         existing_role))
        # Control first, then Disease, then Exclude -- reads as the
        # comparison it defines.
        order = {'control': 0, 'disease': 1, 'exclude': 2}
        rows.sort(key=lambda r: (order.get(r[3], 3), r[0].lower()))

        self.role_table.setRowCount(len(rows) + (1 if other_n else 0))
        for i, (val, samples_text, cells_text, role) in enumerate(rows):
            # Editable: renaming a condition value here is the only way to
            # change it. The original is kept in UserRole so _save_setup can
            # build the old -> new map; nothing is written until Save.
            value_item = QTableWidgetItem(val)
            value_item.setData(Qt.ItemDataRole.UserRole, val)
            value_item.setFlags(value_item.flags() | Qt.ItemFlag.ItemIsEditable)
            value_item.setToolTip(
                "Double-click to rename this condition value. The new name "
                "replaces the old one throughout the condition column when "
                "you press Save & continue.")
            self.role_table.setItem(i, 0, value_item)
            for col, text in ((1, samples_text), (2, cells_text)):
                cell = QTableWidgetItem(text)
                cell.setFlags(cell.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.role_table.setItem(i, col, cell)
            combo = self._make_role_combo(role.capitalize())
            self._set_role_cell(i, combo)
            self._role_rows.append((val, combo))
        if other_n:
            i = len(rows)
            self.role_table.setItem(
                i, 0, QTableWidgetItem(f"(other — {other_n} more values)"))
            self.role_table.setItem(i, 1, QTableWidgetItem(""))
            self.role_table.setItem(i, 2, QTableWidgetItem(""))
            self._set_role_cell(
                i, self._make_role_combo('Exclude', enabled=False))
        self._size_role_table()
        self._populating_roles = False

        if too_many:
            self.role_status_label.setText(
                f"Condition column has {len(col_str.unique())} unique "
                f"values — sure this is the right column? Top "
                f"{WARN_THRESHOLD} shown, the rest default to Exclude.")
            self.role_status_label.set_state('warning')
        elif unrecognised:
            # The auto-guess dropped these to Exclude; make it loud, not
            # silent -- excluded samples vanish from DE and meta-analysis.
            preview = ", ".join(unrecognised[:4])
            if len(unrecognised) > 4:
                preview += f" (+{len(unrecognised) - 4} more)"
            self.role_status_label.setText(
                f"{len(unrecognised)} value(s) not recognised and defaulted "
                f"to Exclude: {preview}. Excluded samples are dropped from "
                f"DE — review before saving.")
            self.role_status_label.set_state('warning')
        else:
            self.role_status_label.setText("")
            self.role_status_label.set_state('info')
        self._refresh_setup_completion()

    def _on_role_value_edited(self, item):
        """A condition value was renamed in the roles table."""
        if getattr(self, '_populating_roles', False) or item.column() != 0:
            return
        self._refresh_setup_completion()
        # The per-sample dialog reads _condition_options() when it opens,
        # so it will pick the new name up; refresh it now if it is showing.
        dialog = getattr(self, '_sample_assign_dialog', None)
        if dialog is not None and dialog.isVisible():
            self._populate_sample_assignment_table()

    def _size_role_table(self):
        """Content-driven height: the table shows all rows, no scroll."""
        n = self.role_table.rowCount()
        header_h = self.role_table.horizontalHeader().sizeHint().height()
        row_h = self.role_table.verticalHeader().defaultSectionSize()
        self.role_table.setFixedHeight(header_h + max(n, 1) * row_h + 4)

    def _staged_renames(self) -> dict:
        """Renames typed into the Condition value column: old -> new.

        Pure read, no dialogs -- safe to call from anything that needs to
        know the pending names (the per-sample dropdown, the banner).
        """
        renames = {}
        for i in range(self.role_table.rowCount()):
            item = self.role_table.item(i, 0)
            if item is None:
                continue
            original = item.data(Qt.ItemDataRole.UserRole)
            current = item.text().strip()
            if original and current and current != str(original):
                renames[str(original)] = current
        return renames

    def _collect_condition_renames(self) -> dict:
        """Renames to apply at save, having refused any that would merge.

        Skips blanks, unchanged names, and any rename onto a value another
        row already holds -- merging two conditions is not a rename, and
        doing it silently would fold two arms together.
        """
        renames = {}
        names = {}
        for i in range(self.role_table.rowCount()):
            item = self.role_table.item(i, 0)
            if item is not None:
                names[i] = item.text().strip()

        for i in range(self.role_table.rowCount()):
            item = self.role_table.item(i, 0)
            if item is None:
                continue
            original = item.data(Qt.ItemDataRole.UserRole)
            current = names[i]
            if not original or not current or current == original:
                continue
            if current in {n for j, n in names.items() if j != i}:
                dialogs.warning(
                    self, "Name already in use",
                    f"'{current}' is already a condition value. Renaming "
                    f"'{original}' to it would merge the two groups, which "
                    "this does not do. Pick a different name.")
                item.setText(str(original))
                continue
            renames[str(original)] = current
        return renames

    def _open_sample_assign_dialog(self):
        """Show the per-sample condition labelling dialog."""
        self._update_sample_assign_visibility()
        self._populate_sample_assignment_table()
        self._sample_assign_dialog.show()
        self._sample_assign_dialog.raise_()
        self._sample_assign_dialog.activateWindow()

    def _open_columns_dialog(self):
        """Show the metadata-column browser."""
        self._columns_dialog.show()
        self._columns_dialog.raise_()
        self._columns_dialog.activateWindow()

    def _set_banner(self, state: str, text: str):
        """Swap the completion banner between success and info styling."""
        self.mapping_banner.setProperty("role", f"banner_{state}")
        for w in (self.mapping_banner, self._banner_icon, self._banner_text):
            style = w.style()
            if style is not None:
                style.unpolish(w)
                style.polish(w)
        icon = 'circle-check' if state == 'success' else 'circle-dot'
        colour = get_color('success' if state == 'success' else 'fg_secondary')
        self._banner_icon.setPixmap(make_icon(icon, colour, 16).pixmap(16, 16))
        self._banner_text.setText(text)

    def _pending_role_map(self) -> dict:
        """The role map the form currently describes, keyed by final name."""
        renames = self._staged_renames()
        out = {}
        for value, combo in self._role_rows:
            role = combo.currentText().lower()
            if role in ('disease', 'control', 'exclude'):
                out[renames.get(value, value)] = role
        return out

    def _setup_is_dirty(self) -> bool:
        """True when the form differs from what is written in the h5ad.

        Computed rather than tracked with a flag, so it stays right however
        the form was changed -- a role dropdown, a renamed value, or a
        different condition column.
        """
        if self.adata is None:
            return False
        if self._staged_renames():
            return True
        saved = self.adata.uns.get('role_map')
        saved = {str(k): str(v) for k, v in saved.items()} if isinstance(saved, dict) else {}
        if self._pending_role_map() != saved:
            return True
        saved_col = self.adata.uns.get('role_condition_col')
        current_col = self._designated_condition_col()
        return bool(saved_col) and str(saved_col) != str(current_col)

    def _refresh_setup_completion(self):
        """Drive the banner under the roles card.

        It used to report only that the form was filled in -- and since the
        roles are pre-filled with guesses that was true the moment the page
        opened, so it read "complete" before anything had been saved. It
        now distinguishes incomplete, unsaved, and saved.
        """
        if self.adata is None:
            self.mapping_banner.setVisible(False)
            return
        self.mapping_banner.setVisible(True)

        has_sample = (self._designated_col('Sample') is not None
                      or 'sample' in self.adata.obs.columns)
        role_map = self._pending_role_map()
        roles = list(role_map.values())

        if not self._role_rows:
            self._set_banner(
                'info', "Select a condition column to assign roles")
        elif not has_sample:
            self._set_banner(
                'info', "Designate a sample column under Dataset Mapping")
        elif 'control' not in roles or 'disease' not in roles:
            self._set_banner(
                'info', "Assign at least one Control and one Disease value")
        elif self._setup_is_dirty():
            n_ctl = roles.count('control')
            n_dis = roles.count('disease')
            n_exc = roles.count('exclude')
            excl = f", {n_exc} excluded" if n_exc else ""
            self._set_banner(
                'info',
                f"Unsaved: {n_ctl} control, {n_dis} disease{excl}. "
                "Press Save & continue to write this to the dataset.")
        else:
            n_exc = roles.count('exclude')
            excl = f", {n_exc} excluded" if n_exc else ""
            self._set_banner(
                'success',
                f"Saved: {roles.count('control')} control, "
                f"{roles.count('disease')} disease{excl}.")

        self.save_setup_btn.setEnabled(self._setup_is_dirty())

    # ------------------------------------------------------------------
    # Sample condition assignment
    # ------------------------------------------------------------------

    def _infer_condition_from_sample_name(self, sample_name: str) -> str:
        """Infer condition (disease/control) from sample name patterns."""
        sample_lower = str(sample_name).lower()

        disease_patterns = ['chf', 'hf', 'disease', 'patient', 'case', 'tumor',
                           'cancer', 'treated', 'sick', 'affected', 'mutant']
        for pattern in disease_patterns:
            if pattern in sample_lower:
                return 'disease'

        control_patterns = ['normal', 'control', 'ctrl', 'healthy', 'wt', 'wildtype',
                           'wild_type', 'sham', 'untreated', 'donor']
        for pattern in control_patterns:
            if pattern in sample_lower:
                return 'control'

        return ''

    def _condition_options(self) -> list:
        """Values the per-sample dropdown offers.

        Every distinct value already in the condition column, not the first
        two. GSE292067's 'group' has three (Donor, NICM, DoxCM); offering
        only two meant a DoxCM sample's combo could not hold its own value,
        so it silently displayed the first option instead and Apply would
        have rewritten those cells -- collapsing a three-arm design to two.

        Falls back to the two free-text labels only when there is no
        condition column to read, which is the case this dialog exists for.
        """
        values = []
        col = self._sample_assign_condition_col()
        if col is not None:
            # Pending renames are staged in the roles table, not yet in obs.
            # Offering the old names here would let you assign a sample to a
            # value that is about to stop existing.
            renames = self._staged_renames()
            values = sorted(
                {renames.get(str(v), str(v))
                 for v in self.adata.obs[col].unique()
                 if str(v) not in ('', 'nan', 'None')})
        if values:
            return values
        label1 = self.condition_label_1.text().strip() or "Control"
        label2 = self.condition_label_2.text().strip() or "Disease"
        return [label1, label2]

    def _sample_assign_condition_col(self):
        """The condition column this dialog reads, or None."""
        if self.adata is None:
            return None
        designated = self._designated_condition_col()
        if designated and designated in self.adata.obs.columns:
            return designated
        for candidate in ('condition', 'Condition', 'group', 'Group',
                          'disease', 'Disease', 'status', 'Status'):
            if candidate in self.adata.obs.columns:
                return candidate
        return None

    def _populate_sample_assignment_table(self):
        """Populate the sample assignment table with unique sample values."""
        self.sample_assign_table.setRowCount(0)
        self.apply_sample_conditions_btn.setEnabled(False)

        if self.adata is None:
            return

        sample_col = None
        sample_candidates = [
            'sample', 'Sample', 'orig.ident', 'orig_ident',
            'biosample_id', 'donor_id', 'patient_id', 'subject_id',
            'donor', 'patient', 'subject', 'sample_id', 'SampleID',
            'individual', 'specimen', 'batch'
        ]
        for col in sample_candidates:
            if col in self.adata.obs.columns:
                sample_col = col
                break

        if sample_col is None:
            return

        unique_samples = self.adata.obs[sample_col].unique()
        if len(unique_samples) <= 1:
            return

        # Auto-detect condition labels from the existing condition column.
        existing_conditions = {}
        condition_col = None
        for candidate in ['condition', 'Condition', 'group', 'Group',
                          'disease', 'Disease', 'status', 'Status']:
            if candidate in self.adata.obs.columns:
                condition_col = candidate
                break

        if condition_col is not None:
            # Deliberately not copied into the two label boxes: a column
            # with three or more values does not fit in two, and writing
            # the first two there is what lost the third.
            # Map each sample to its existing condition.
            for sample in unique_samples:
                mask = self.adata.obs[sample_col] == sample
                conditions = self.adata.obs.loc[mask, condition_col].unique()
                if len(conditions) == 1:
                    existing_conditions[sample] = str(conditions[0])

        options = self._condition_options()
        # The free-text labels only build a column that does not exist yet;
        # with one already present they would just mislead.
        has_column = self._sample_assign_condition_col() is not None
        self._condition_labels_row.setVisible(not has_column)

        # Optional third column: dominant value of a user-picked obs column per sample.
        meta_col = self.meta_show_combo.currentText()
        show_meta = (meta_col not in ('(none)', '')
                     and meta_col in self.adata.obs.columns)

        if show_meta:
            self.sample_assign_table.setColumnCount(3)
            self.sample_assign_table.setHorizontalHeaderLabels(
                ["Sample", f"{meta_col}", "Condition"])
            self.sample_assign_table.horizontalHeader().setSectionResizeMode(
                0, QHeaderView.ResizeMode.Stretch)
            self.sample_assign_table.horizontalHeader().setSectionResizeMode(
                1, QHeaderView.ResizeMode.Stretch)
            self.sample_assign_table.horizontalHeader().setSectionResizeMode(
                2, QHeaderView.ResizeMode.Stretch)
        else:
            self.sample_assign_table.setColumnCount(2)
            self.sample_assign_table.setHorizontalHeaderLabels(
                ["Sample", "Condition"])
            self.sample_assign_table.horizontalHeader().setSectionResizeMode(
                0, QHeaderView.ResizeMode.Stretch)
            self.sample_assign_table.horizontalHeader().setSectionResizeMode(
                1, QHeaderView.ResizeMode.Stretch)

        # Compute dominant meta value per sample once (cheap)
        meta_per_sample = {}
        if show_meta:
            for sample in unique_samples:
                mask = self.adata.obs[sample_col] == sample
                vals = self.adata.obs.loc[mask, meta_col]
                vals = vals[~vals.astype(str).isin(['', 'nan'])]
                if len(vals) == 0:
                    meta_per_sample[sample] = ''
                    continue
                vc = vals.value_counts()
                dom = str(vc.index[0])
                pct = 100.0 * vc.iloc[0] / len(vals)
                # Show "value (98%)" if not unanimous, else just the value.
                if pct >= 99.5:
                    meta_per_sample[sample] = dom
                else:
                    meta_per_sample[sample] = f"{dom} ({pct:.0f}%)"

        for sample in sorted(unique_samples):
            row = self.sample_assign_table.rowCount()
            self.sample_assign_table.insertRow(row)

            sample_item = QTableWidgetItem(str(sample))
            sample_item.setFlags(sample_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.sample_assign_table.setItem(row, 0, sample_item)

            if show_meta:
                meta_item = QTableWidgetItem(meta_per_sample.get(sample, ''))
                meta_item.setFlags(
                    meta_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.sample_assign_table.setItem(row, 1, meta_item)
                cond_col_idx = 2
            else:
                cond_col_idx = 1

            condition_combo = NoScrollComboBox()
            condition_combo.addItems(options)

            if sample in existing_conditions and existing_conditions[sample] not in ['', 'nan']:
                condition_combo.setCurrentText(existing_conditions[sample])
            else:
                inferred = self._infer_condition_from_sample_name(sample)
                condition_combo.setCurrentText(inferred)

            self.sample_assign_table.setCellWidget(
                row, cond_col_idx, condition_combo)

        self.apply_sample_conditions_btn.setEnabled(True)

    def _update_condition_options(self):
        """Update the condition dropdown options when custom labels change."""
        if self.sample_assign_table.rowCount() == 0:
            return

        options = self._condition_options()
        cond_col_idx = self.sample_assign_table.columnCount() - 1
        for row in range(self.sample_assign_table.rowCount()):
            combo = self.sample_assign_table.cellWidget(row, cond_col_idx)
            if combo:
                current_text = combo.currentText()
                combo.clear()
                combo.addItems(options)
                if current_text in options:
                    combo.setCurrentText(current_text)

    def _apply_sample_conditions(self):
        """Apply the sample condition assignments to create/update the condition column."""
        if self.adata is None:
            return

        sample_col = None
        sample_candidates = [
            'sample', 'Sample', 'orig.ident', 'orig_ident',
            'biosample_id', 'donor_id', 'patient_id', 'subject_id',
            'donor', 'patient', 'subject', 'sample_id', 'SampleID',
            'individual', 'specimen', 'batch'
        ]
        for col in sample_candidates:
            if col in self.adata.obs.columns:
                sample_col = col
                break

        if sample_col is None:
            dialogs.warning(self, "No Sample Column",
                "No sample column is designated yet.\n\n"
                "Pick one under Dataset Mapping at the top of the Setup "
                "tab, then try again.")
            return

        cond_col_idx = self.sample_assign_table.columnCount() - 1

        sample_to_condition = {}
        missing_conditions = []
        for row in range(self.sample_assign_table.rowCount()):
            sample_name = self.sample_assign_table.item(row, 0).text()
            condition_combo = self.sample_assign_table.cellWidget(
                row, cond_col_idx)
            condition = condition_combo.currentText().strip()

            if not condition:
                missing_conditions.append(sample_name)
            else:
                sample_to_condition[sample_name] = condition

        if missing_conditions:
            dialogs.warning(self, "Missing Conditions",
                "Please specify conditions for:\n" + "\n".join(missing_conditions))
            return

        preview_lines = []
        for sample, condition in sorted(sample_to_condition.items()):
            mask = self.adata.obs[sample_col].astype(str) == str(sample)
            count = mask.sum()
            preview_lines.append(f"  {sample} ({count:,} cells) → {condition}")

        empty_sample_mask = (
            self.adata.obs[sample_col].isna() |
            (self.adata.obs[sample_col].astype(str) == '') |
            (self.adata.obs[sample_col].astype(str) == 'nan')
        )
        empty_count = empty_sample_mask.sum()

        msg = "This will apply the following condition mappings:\n\n"
        msg += "\n".join(preview_lines)
        if empty_count > 0:
            msg += f"\n\n⚠️ {empty_count:,} cells have no sample identifier and will get empty condition."

        if not dialogs.confirm(self, "Confirm Condition Assignment", msg + "\n\nProceed?"):
            return

        sample_values = self.adata.obs[sample_col].astype(str)
        self.adata.obs['condition'] = sample_values.map(sample_to_condition)
        self.adata.obs['condition'] = self.adata.obs['condition'].fillna('')

        self.adata.write_h5ad(self.h5ad_path)
        self._populate_role_rows()
        self._update_sample_assign_visibility()

        condition_counts = self.adata.obs['condition'].value_counts()
        self.summary['obs_summary']['condition'] = {
            'dtype': 'object',
            'n_unique': len(condition_counts),
            'sample_values': list(condition_counts.index[:5]),
            'value_counts': {str(k): int(v) for k, v in condition_counts.head(10).items()}
        }

        self._update_display()
        self.main_window.set_adata(self.adata, str(self.h5ad_path))
        if hasattr(self.main_window, 'record_provenance'):
            self.main_window.record_provenance('sample_conditions', {
                'condition_col': 'condition',
                'values': {str(k): int(v) for k, v
                           in self.adata.obs['condition'].value_counts().items()},
            })

        result_counts = self.adata.obs['condition'].value_counts()
        result_summary = "\n".join([f"  {cond}: {count:,} cells" for cond, count in result_counts.items()])
        dialogs.info(self, "Conditions Applied",
            f"Condition assignment complete!\n\n{result_summary}\n\nSaved to file.")
        self.status_label.setText("Sample conditions applied and saved")

    # ------------------------------------------------------------------
    # External metadata
    # ------------------------------------------------------------------

    def _load_h5ad_file(self, path: str):
        """Load an h5ad file and update display."""
        try:
            import anndata
            adata = anndata.read_h5ad(path)
            self.set_data(adata, path)
        except Exception as e:
            self.status_label.setText(f"Failed to reload: {e}")

    def refresh_theme(self):
        """No-op: subset_status_label colour is set dynamically per-state."""
        pass
