# DE Analysis — Load Data Page
# Load h5ad file and configure columns for DE analysis.
# Styled to match the scRNA Load Data page.

from pathlib import Path

from PyQt6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QGroupBox,
    QFileDialog, QFrame, QMenu,
    QDialog, QCheckBox,
)
from PyQt6.QtCore import Qt

from kosmic.gui.shared.theme import NoScrollComboBox
from kosmic.gui.shared.widgets import (
    SimplePage, CaptionLabel, DatasetOverview, InfoPanel, PrimaryButton,
    SecondaryButton,
)
from kosmic.paths import list_studies, processed_data_dir, study_status
from kosmic.gui.shared.widgets import BaseWorker
from kosmic.gui.shared import dialogs, run_worker


class _CountsLoadWorker(BaseWorker):
    """Load a study for DE: obs, var and the counts as X, nothing else.

    Emits 'finished_ok' with '(adata, path, message)' like the intake
    loader, so the Setup page's handler is shared.
    """

    def __init__(self, h5ad_path: str):
        super().__init__()
        self.h5ad_path = h5ad_path

    def _run(self):
        from kosmic.scrna.load.h5ad_meta import read_counts_adata
        adata = read_counts_adata(self.h5ad_path)
        msg = (f"{adata.n_obs:,} cells x {adata.n_vars:,} genes "
               f"(counts from {adata.uns.get('counts_loaded_from', '?')})")
        return (adata, self.h5ad_path, msg)


class SetupPage(SimplePage):
    """Load Data page: h5ad loading and column configuration."""

    help_id = "de/setup"
    SCROLLABLE = True
    CONTENT_MARGINS = (40, 30, 40, 30)
    CONTENT_SPACING = 0

    def __init__(self, workspace):
        super().__init__()
        self.ws = workspace
        self._load_worker = None
        self._setup_ui()

    def _setup_ui(self):
        # The dataset screen is the shared scaffold; this page only
        # supplies its stats, actions and panels.
        ov = DatasetOverview(
            "Dataset",
            "Review the study's dataset and condition roles before running DE.",
            stats=(("cells", "Cells"), ("genes", "Genes"),
                   ("samples", "Samples")))
        self._overview = ov
        self._hero_stats = ov.stats
        self._study_card = ov.hero

        # ADR-001: DE consumes the project's study; direct file loading
        # is the escape hatch, not the workflow.
        self._change_btn = ov.add_action(SecondaryButton("Change study..."))
        self._change_btn.setToolTip("Switch to another study in this project.")
        self._change_menu = QMenu(self._change_btn)
        self._change_menu.aboutToShow.connect(self._populate_change_menu)
        self._change_btn.setMenu(self._change_menu)
        self._import_btn = ov.add_action(
            SecondaryButton("Import external dataset..."))
        self._import_btn.setToolTip(
            "Load an h5ad from outside the project (one-off analyses).")
        self._import_btn.clicked.connect(self._load_h5ad)

        self._dataset_panel = ov.add_panel(InfoPanel("DATASET"))
        self._dataset_panel.add_row('source', "Source")
        self._dataset_panel.add_row('file', "File")
        self._dataset_panel.add_row(
            'sample_col', "Sample column",
            "Column with patient/sample IDs -- each unique value becomes "
            "one pseudobulk replicate.")

        self._roles_panel = ov.add_panel(InfoPanel("CONDITION ROLES"))
        self._roles_panel.add_row('condition_col', "Condition column")
        self._roles_panel.add_row('control', "Control group")
        self._roles_panel.add_row('disease', "Disease group")
        self._panels_inner = ov.panels_widget
        self._panels_inner.hide()

        layout = ov.extra_layout

        # Shown only when roles are already set. It replaces the column
        # grid rather than sitting beside it: the panels above already
        # state every value, so a second, greyed-out copy of the same
        # four facts said nothing and implied they were editable here.
        self._edit_roles_btn = SecondaryButton("Edit sample and conditions in scRNA →")
        self._edit_roles_btn.setToolTip(
            "Sample column and condition roles are set in "
            "scRNA → Inspect → Setup. Click to jump there.")
        self._edit_roles_btn.clicked.connect(self.ws.open_scrna_requested.emit)
        self._edit_roles_btn.hide()
        layout.addWidget(self._edit_roles_btn,
                         alignment=Qt.AlignmentFlag.AlignHCenter)

        # === No registered dataset / no study state ===
        self._no_dataset_frame = QFrame()
        self._no_dataset_frame.setProperty("role", "drop_zone")
        nd_layout = QVBoxLayout(self._no_dataset_frame)
        nd_layout.setContentsMargins(20, 18, 20, 18)
        nd_layout.setSpacing(6)
        self._no_dataset_label = QLabel(
            "No processed dataset is registered for this study.")
        self._no_dataset_label.setProperty("role", "file_row_name")
        nd_layout.addWidget(self._no_dataset_label)
        self._no_dataset_caption = CaptionLabel(
            "Run the scRNA pipeline on this study first, or import an "
            "external h5ad.")
        self._no_dataset_caption.setWordWrap(True)
        nd_layout.addWidget(self._no_dataset_caption)
        nd_btn_row = QHBoxLayout()
        nd_import_btn = SecondaryButton("Import external dataset...")
        nd_import_btn.clicked.connect(self._load_h5ad)
        nd_btn_row.addWidget(nd_import_btn)
        nd_btn_row.addStretch()
        nd_layout.addLayout(nd_btn_row)
        self._no_dataset_frame.hide()
        layout.addWidget(self._no_dataset_frame)

        # === scRNA pipeline not ready for DE yet ===
        self._not_ready_frame = QFrame()
        self._not_ready_frame.setProperty("role", "drop_zone")
        nr_layout = QVBoxLayout(self._not_ready_frame)
        nr_layout.setContentsMargins(20, 18, 20, 18)
        nr_layout.setSpacing(6)
        self._not_ready_label = QLabel("scRNA pipeline isn't ready for DE yet.")
        self._not_ready_label.setProperty("role", "file_row_name")
        nr_layout.addWidget(self._not_ready_label)
        self._not_ready_caption = CaptionLabel("")
        self._not_ready_caption.setWordWrap(True)
        nr_layout.addWidget(self._not_ready_caption)
        nr_btn_row = QHBoxLayout()
        nr_goto_btn = SecondaryButton("Go to scRNA →")
        nr_goto_btn.clicked.connect(self.ws.open_scrna_requested.emit)
        nr_btn_row.addWidget(nr_goto_btn)
        nr_btn_row.addStretch()
        nr_layout.addLayout(nr_btn_row)
        self._not_ready_frame.hide()
        layout.addWidget(self._not_ready_frame)

        # === Column Configuration (centered card) ===
        # Only meaningful once a dataset is loaded and turns out to
        # have no roles recorded; hidden until then.
        config_frame = QFrame()
        self._config_frame = config_frame
        config_frame.hide()
        config_frame.setMaximumWidth(500)
        config_frame.setProperty("role", "panel")
        config_layout = QVBoxLayout(config_frame)
        config_layout.setContentsMargins(16, 16, 16, 16)
        config_layout.setSpacing(10)

        title_row = QHBoxLayout()
        config_title = QLabel("Column Configuration")
        config_title.setProperty("role", "header")
        title_row.addWidget(config_title)
        title_row.addStretch()
        config_layout.addLayout(title_row)

        config_sub = CaptionLabel("Assign which columns identify samples and conditions.")
        config_sub.setWordWrap(True)
        config_layout.addWidget(config_sub)

        grid = QGridLayout()
        grid.setColumnStretch(0, 0)
        grid.setColumnStretch(1, 1)
        grid.setVerticalSpacing(8)
        grid.setHorizontalSpacing(12)

        lbl = QLabel("Sample column:")
        grid.addWidget(lbl, 0, 0)
        self.sample_col_combo = NoScrollComboBox()
        self.sample_col_combo.setMinimumWidth(180)
        self.sample_col_combo.setToolTip("Column with patient/sample IDs.\nEach unique value becomes one pseudobulk replicate.")
        self.sample_col_combo.currentTextChanged.connect(self._on_sample_col_changed)
        grid.addWidget(self.sample_col_combo, 0, 1)

        lbl = QLabel("Condition column:")
        grid.addWidget(lbl, 1, 0)
        self.condition_col_combo = NoScrollComboBox()
        self.condition_col_combo.setToolTip("Column with experimental conditions.")
        self.condition_col_combo.currentTextChanged.connect(self._on_condition_col_changed)
        grid.addWidget(self.condition_col_combo, 1, 1)

        lbl = QLabel("Control group:")
        grid.addWidget(lbl, 2, 0)
        self.control_combo = NoScrollComboBox()
        self.control_combo.setToolTip("Which condition value represents the control/baseline group.")
        self.control_combo.currentTextChanged.connect(self._on_control_changed)
        grid.addWidget(self.control_combo, 2, 1)

        lbl = QLabel("Disease group:")
        grid.addWidget(lbl, 3, 0)
        self.disease_combo = NoScrollComboBox()
        self.disease_combo.setToolTip("Which condition value represents the disease/treatment group.")
        self.disease_combo.currentTextChanged.connect(self._on_disease_changed)
        grid.addWidget(self.disease_combo, 3, 1)

        config_layout.addLayout(grid)

        self.combine_conditions_btn = QPushButton("Combine Conditions...")
        self.combine_conditions_btn.setObjectName("ghost_button")
        self.combine_conditions_btn.setToolTip("Merge multiple conditions into groups (e.g., DCM + HCM -> HF)")
        self.combine_conditions_btn.clicked.connect(self._show_combine_conditions_dialog)
        self.combine_conditions_btn.setEnabled(False)
        config_layout.addWidget(self.combine_conditions_btn)

        layout.addWidget(config_frame, alignment=Qt.AlignmentFlag.AlignHCenter)

        # Same placement as scRNA's Dataset page: the column sits a
        # little above centre, with the next step bottom-right.
        body = self.body_layout
        body.addStretch(2)
        center_row = QHBoxLayout()
        center_row.addStretch(1)
        center_row.addWidget(ov)
        center_row.addStretch(1)
        body.addLayout(center_row)
        body.addStretch(3)

        nav_row = QHBoxLayout()
        nav_row.addStretch()
        self._continue_btn = PrimaryButton("Continue to Choose Analysis  →")
        self._continue_btn.setMinimumHeight(32)
        self._continue_btn.setToolTip("Go to the next step: pick scoring or discovery.")
        self._continue_btn.clicked.connect(lambda: self.ws.switch_tab(1))
        self._continue_btn.hide()
        nav_row.addWidget(self._continue_btn)
        body.addLayout(nav_row)

    # ------------------------------------------------------------------
    # Activation
    # ------------------------------------------------------------------

    def on_activated(self):
        """Called when page becomes visible."""
        self._refresh_study_card()
        if self.ws.current_adata is not None:
            self._refresh_from_adata()
        else:
            self._auto_load_registered_dataset()

    def _registered_h5ad(self):
        """The study's working h5ad: the manifest's file if it still
        exists, else the newest in processed_data/. None if nothing."""
        from kosmic.manifest import study_overview
        proj = self.ws.project_dir
        if not proj:
            return None
        processed = processed_data_dir(Path(proj))
        dataset = study_overview(Path(proj))["dataset"] or {}
        if dataset.get('file'):
            named = processed / dataset['file']
            if named.is_file():
                return named
        h5ads = sorted(processed.glob('*.h5ad'), key=lambda p: p.stat().st_mtime)
        return h5ads[-1] if h5ads else None

    def _auto_load_registered_dataset(self):
        """Load the study's h5ad in the background, the way scRNA's
        Load Data step does, so arriving from the Project page shows
        the dataset rather than an empty form waiting for a click."""
        if self._load_worker is not None and self._load_worker.isRunning():
            return
        path = self._registered_h5ad()
        if path is None:
            return
        self._start_load(str(path))

    def show_loaded_state(self, label: str, detail: str) -> None:
        """Refresh column combos and show the loaded status on the study card."""
        self._not_ready_frame.hide()
        self._refresh_from_adata()
        self._refresh_study_card()
        self._set_study_status(label, 'success', tooltip=detail)
        self._overview.set_ready(True)
        self._continue_btn.show()

    def hide_loaded_state(self) -> None:
        """Clear the loaded status (called when adata is cleared)."""
        self._not_ready_frame.hide()
        self._overview.set_ready(False)
        self._continue_btn.hide()
        self._refresh_study_card()

    def show_not_ready_state(self, reason: str) -> None:
        """scRNA has data but 'de_readiness' failed -- show why, with a
        way to jump back to scRNA to finish the pipeline."""
        self._study_card.hide()
        self._panels_inner.hide()
        self._no_dataset_frame.hide()
        self._overview.set_ready(False)
        self._continue_btn.hide()
        self._not_ready_caption.setText(reason)
        self._not_ready_frame.show()

    # ------------------------------------------------------------------
    # Study card (ADR-001)
    # ------------------------------------------------------------------

    def _set_study_status(self, text: str, state: str, tooltip: str = "") -> None:
        self._overview.set_status(text, state, tooltip)

    def _refresh_study_card(self):
        """Reflect the active study's registered dataset (manifest +
        disk flags) on the study card, or show the empty state."""
        from kosmic.manifest import study_overview

        proj = self.ws.project_dir
        adata = self.ws.current_adata

        if adata is not None:
            self._not_ready_frame.hide()

        if not proj:
            if adata is not None:
                # External dataset loaded with no study context
                name = (Path(self.ws.h5ad_path).stem
                        if getattr(self.ws, 'h5ad_path', None) else "External dataset")
                self._overview.set_name(name)
                self._overview.set_subtitle("External dataset")
                self._set_study_status(
                    f"Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes",
                    'success')
                self._change_btn.setEnabled(False)
                self._study_card.show()
                self._panels_inner.show()
                self._no_dataset_frame.hide()
            else:
                self._no_dataset_label.setText("No study selected.")
                self._no_dataset_caption.setText(
                    "Open a project and pick a study in the Project "
                    "workspace, or import an external h5ad.")
                self._study_card.hide()
                self._panels_inner.hide()
                self._no_dataset_frame.show()
            return

        study_dir = Path(proj)
        self._overview.set_name(study_dir.name)
        self._change_btn.setEnabled(True)
        overview = study_overview(study_dir)
        dataset = overview["dataset"] or {}

        if adata is not None:
            file_name = (Path(self.ws.h5ad_path).name
                         if getattr(self.ws, 'h5ad_path', None)
                         else dataset.get('file', ''))
            self._overview.set_subtitle(
                f"Processed dataset · {file_name}" if file_name
                else "Processed dataset")
            self._set_study_status(
                f"✓ Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes",
                'success')
            self._study_card.show()
            self._panels_inner.show()
            self._no_dataset_frame.hide()
        elif overview["status"]["processed"]:
            if dataset:
                self._overview.set_subtitle(
                    f"Processed dataset · {dataset.get('file', '')}")
                self._set_study_status(
                    f"{dataset.get('n_cells', 0):,} cells × "
                    f"{dataset.get('n_genes', 0):,} genes · not loaded yet",
                    'info')
            else:
                self._overview.set_subtitle("Processed dataset present")
                self._set_study_status("Not loaded yet", 'info')
            self._study_card.show()
            self._panels_inner.hide()
            self._no_dataset_frame.hide()
            self._config_frame.hide()
        else:
            self._no_dataset_label.setText(
                "No processed dataset is registered for this study.")
            self._no_dataset_caption.setText(
                "Run the scRNA pipeline on this study first, or import an "
                "external h5ad.")
            self._study_card.hide()
            self._panels_inner.hide()
            self._config_frame.hide()
            self._no_dataset_frame.show()

    def _populate_change_menu(self):
        """Fill the Change-study menu with the project's studies (and any
        other h5ads within the current study)."""
        menu = self._change_menu
        menu.clear()
        proj = self.ws.project_dir
        root = Path(proj).parent if proj else None
        current = Path(proj).name if proj else None

        studies = list_studies(root) if root else []
        if not studies:
            act = menu.addAction("No studies in project")
            act.setEnabled(False)
        for acc in studies:
            st = study_status(root / acc)
            if acc == current:
                act = menu.addAction(f"{acc}  (current)")
                act.setEnabled(False)
            elif st["processed"]:
                act = menu.addAction(acc)
                act.triggered.connect(
                    lambda _c=False, a=acc:
                    self.ws.study_change_requested.emit(a))
            else:
                act = menu.addAction(f"{acc}  (no processed data)")
                act.setEnabled(False)

        # Other h5ads within the current study
        if proj:
            pd_dir = processed_data_dir(Path(proj))
            cur_file = str(self.ws.h5ad_path) if getattr(self.ws, 'h5ad_path', None) else ""
            others = ([f for f in sorted(pd_dir.glob("*.h5ad"))
                       if str(f) != cur_file]
                      if pd_dir.exists() else [])
            if others:
                menu.addSeparator()
                header = menu.addAction("Other datasets in this study")
                header.setEnabled(False)
                for f in others:
                    act = menu.addAction(f"  {f.name}")
                    act.triggered.connect(
                        lambda _c=False, p=str(f): self._do_load_h5ad(p))

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_h5ad(self):
        start_dir = ""
        if self.ws.project_dir:
            start_dir = str(self.ws.project_dir)
        elif self.ws.settings.value("last_directory", ""):
            start_dir = self.ws.settings.value("last_directory", "")

        path, _ = QFileDialog.getOpenFileName(
            self, "Load h5ad File", start_dir,
            "AnnData files (*.h5ad);;All files (*)",
        )
        if path:
            self._do_load_h5ad(path)

    def _do_load_h5ad(self, path):
        """User picked a file: load it (in the background)."""
        self._start_load(str(path))

    def _start_load(self, path: str) -> None:
        if self._load_worker is not None and self._load_worker.isRunning():
            return
        self._study_card.show()
        self._no_dataset_frame.hide()
        self._not_ready_frame.hide()
        self._config_frame.hide()
        self._set_study_status(f"Loading {Path(path).name}...", 'info')
        self.ws.status_message.emit(f"DE: loading {Path(path).name}...")
        # Counts, obs and var only: differential expression never reads
        # the normalised matrix, the embeddings or the graph, and loading
        # them doubled the memory of every study (the atlas did not fit).
        self._load_worker = _CountsLoadWorker(path)
        run_worker(
            self._load_worker,
            on_finished=self._on_loaded,
            on_failed=self._on_load_failed,
        )

    def _on_loaded(self, payload) -> None:
        adata, path, _message = payload
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
        # Same coercion as DEWorkspace.set_adata: numeric condition
        # columns must be strings or downstream plot routines silently
        # produce empty masks against the string labels in the combos.
        self.ws._coerce_condition_columns_to_string(adata)
        self.ws.current_adata = adata
        self.ws.h5ad_path = path
        # The file on disk is the project's source of truth (ADR-001):
        # a dataset loaded from it must not be silently replaced by
        # whatever scRNA happens to hold in memory.
        self.ws._manual_load = True
        self.ws.dataset_loaded.emit(str(path))

        # Reset stale state from previous dataset
        self.ws.de_results = None
        self.ws.significant_genes = None
        self.ws.pathway_coverage = None
        self.ws.sample_df = None
        self.ws.pathway_de_results = None
        self.ws.parent_stats_df = None
        self.ws.enrichment_results = None

        self.show_loaded_state(
            label=f"Data loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes",
            detail=str(path),
        )
        self.ws.log_message.emit(
            f"Loaded h5ad: {path} ({adata.n_obs:,} cells x {adata.n_vars:,} genes)")
        self.ws.status_message.emit(f"DE: Loaded {path} ({adata.n_obs:,} cells)")
        self.ws.mark_step_complete(0)  # Step 0 = Load Data

    def _on_load_failed(self, message: str) -> None:
        self._study_card.show()
        self._no_dataset_frame.hide()
        self._not_ready_frame.hide()
        self._set_study_status(f"Error: {message}", 'error')
        self.ws.log_message.emit(f"Error loading h5ad: {message}")

    @staticmethod
    def _sample_column_candidates(adata):
        """obs columns that could plausibly identify a donor.

        The old list was every column with under 100 values, which on a
        processed study is 34 entries including 'cell_type', 'study',
        'scrublet' and half a dozen clustering resolutions. Choosing one
        of those silently redefines what a pseudobulk replicate is --
        pick 'cell_type' and your "donors" become cell types.

        The decisive test is roles: 'create_pseudobulk' takes one role
        per sample group, so a column whose values straddle disease and
        control cannot be a sample column. Anything that fails that is
        not a judgement call, it is wrong.
        """
        if adata is None:
            return set()
        obs = adata.obs
        roles = (obs['_role'].astype(str) if '_role' in obs.columns else None)

        keep = set()
        for col in obs.columns:
            if col.startswith('_'):
                continue
            series = obs[col]
            if not (str(series.dtype) in ('category', 'object')
                    or hasattr(series, 'cat')):
                continue
            n_unique = series.nunique()
            if n_unique < 2 or n_unique > 100:
                continue
            if roles is not None:
                # Every value of a real sample column sits in one arm.
                if roles.groupby(series.astype(str), observed=True).nunique().max() > 1:
                    continue
            keep.add(col)
        return keep

    def _refresh_from_adata(self):
        """Populate all combos from current adata."""
        adata = self.ws.current_adata
        if adata is None:
            return

        # Populate column combos
        self.sample_col_combo.blockSignals(True)
        self.condition_col_combo.blockSignals(True)
        self.sample_col_combo.clear()
        self.condition_col_combo.clear()

        sample_candidates = self._sample_column_candidates(adata)
        for col in adata.obs.columns:
            n_unique = adata.obs[col].nunique()
            if n_unique <= 100:
                self.condition_col_combo.addItem(col)
                if col in sample_candidates:
                    self.sample_col_combo.addItem(col)

        # Restore previous selection if workspace already has one,
        # otherwise auto-select from candidates
        restored_sample = False
        if self.ws.sample_col:
            idx = self.sample_col_combo.findText(self.ws.sample_col)
            if idx >= 0:
                self.sample_col_combo.setCurrentIndex(idx)
                restored_sample = True

        if not restored_sample:
            sample_candidates = [
                'donor_id', 'sample', 'sample_id', 'Sample', 'patient',
                'Patient', 'donor', 'Donor', 'subject', 'Subject',
                'patient_id', 'PatientID', 'subject_id', 'individual',
                'biosample_id', 'orig.ident',
            ]
            for candidate in sample_candidates:
                idx = self.sample_col_combo.findText(candidate)
                if idx >= 0:
                    self.sample_col_combo.setCurrentIndex(idx)
                    break

        restored_condition = False
        # Priority order:
        #   1. '_role' -- the canonical role-aware column written by
        #      Inspect (values are already 'control' / 'disease' /
        #      'exclude'). Picking this gives plots meaningful labels
        #      instead of whatever the upstream condition column held
        #      (e.g. Int64 0/1 in the our Koenig file).
        #   2. 'role_condition_col' -- the Inspect-tab source column
        #      the role_map was built against. Useful if the user
        #      wants to see the original encoding.
        #   3. The workspace's last-used condition_col.
        #   4. A list of common name candidates.
        adata = self.ws.current_adata
        if adata is not None and '_role' in adata.obs.columns:
            idx = self.condition_col_combo.findText('_role')
            if idx >= 0:
                self.condition_col_combo.setCurrentIndex(idx)
                restored_condition = True

        if not restored_condition:
            inspect_cond_col = (
                adata.uns.get('role_condition_col')
                if adata is not None else None
            )
            if inspect_cond_col:
                idx = self.condition_col_combo.findText(str(inspect_cond_col))
                if idx >= 0:
                    self.condition_col_combo.setCurrentIndex(idx)
                    restored_condition = True

        if not restored_condition and self.ws.condition_col:
            idx = self.condition_col_combo.findText(self.ws.condition_col)
            if idx >= 0:
                self.condition_col_combo.setCurrentIndex(idx)
                restored_condition = True

        if not restored_condition:
            condition_candidates = [
                'condition', 'Condition', 'disease', 'Disease',
                'disease_status', 'diagnosis', 'Diagnosis',
                'group', 'Group', 'status', 'Status',
                'treatment', 'Treatment', 'cell_type', 'phenotype',
            ]
            for candidate in condition_candidates:
                idx = self.condition_col_combo.findText(candidate)
                if idx >= 0:
                    self.condition_col_combo.setCurrentIndex(idx)
                    break

        self.sample_col_combo.blockSignals(False)
        self.condition_col_combo.blockSignals(False)

        # Sync to workspace
        self.ws.sample_col = self.sample_col_combo.currentText()
        self.ws.condition_col = self.condition_col_combo.currentText()

        # Trigger condition column change to populate control/disease
        self._on_condition_col_changed(self.condition_col_combo.currentText())

        self._apply_role_lock_state(adata)
        self._update_info_panels()

    # ------------------------------------------------------------------
    # Condition-role lock (read-only when scRNA already assigned roles)
    # ------------------------------------------------------------------

    def _apply_role_lock_state(self, adata) -> None:
        """Lock condition/control/disease to read-only when 'adata'
        already carries an authoritative role assignment from scRNA's
        Inspect tab ('role_map' or a resolved '_role' column).

        Editing the role assignment is scRNA's job (Inspect -> Setup);
        letting DE silently redefine it here is what used to make it
        unclear whether a DE-side edit should propagate back to scRNA.
        Locking removes the question -- fix it at the source instead.
        An external dataset with no role source stays fully editable,
        since there's nowhere else to set it.

        The lock is scoped to the canonical column, not the adata as a
        whole: "Combine Conditions" still works to build a DE-local
        derived grouping (it's not scRNA's role assignment, so there's
        nothing to defer to), and once that derived column is selected,
        Control/Disease unlock for it specifically.
        """
        canonical_col = None
        role_map = adata.uns.get('role_map') if adata is not None else None
        if adata is not None:
            if '_role' in adata.obs.columns:
                canonical_col = '_role'
            elif isinstance(role_map, dict) and role_map:
                canonical_col = adata.uns.get('role_condition_col')

        has_role_source = canonical_col is not None
        self._role_locked = has_role_source
        current_col = self.condition_col_combo.currentText()
        roles_locked_here = has_role_source and current_col == canonical_col

        # With roles set there is nothing here to choose -- the sample
        # column and the condition roles were both decided in Inspect --
        # so the whole grid goes rather than being shown greyed out.
        # The DATASET / CONDITION ROLES panels above carry the values.
        if hasattr(self, '_config_frame'):
            self._config_frame.setVisible(not has_role_source)

        self.sample_col_combo.setEnabled(not has_role_source)
        self.condition_col_combo.setEnabled(not has_role_source)
        self.control_combo.setEnabled(not roles_locked_here)
        self.disease_combo.setEnabled(not roles_locked_here)

        lock_tip = ("Set in scRNA → Inspect → Setup. Click \"Edit in "
                    "scRNA\" to change, or use \"Combine Conditions\" to "
                    "build a grouping for this analysis." if has_role_source else "")
        self.condition_col_combo.setToolTip(
            lock_tip or "Column with experimental conditions.")
        self.sample_col_combo.setToolTip(
            lock_tip or "Column with patient/sample IDs.\n"
            "Each unique value becomes one pseudobulk replicate.")
        role_tip = lock_tip if roles_locked_here else ""
        self.control_combo.setToolTip(
            role_tip or "Which condition value represents the control/baseline group.")
        self.disease_combo.setToolTip(
            role_tip or "Which condition value represents the disease/treatment group.")
        self._edit_roles_btn.setVisible(has_role_source)

    def _update_info_panels(self) -> None:
        """Fill the DATASET / CONDITION ROLES panels from current state."""
        adata = self.ws.current_adata
        if adata is None:
            return

        proj = self.ws.project_dir
        h5ad_path = getattr(self.ws, 'h5ad_path', None)
        self._dataset_panel.set_value(
            'source', Path(proj).name if proj else "External dataset")
        file_name = Path(h5ad_path).name if h5ad_path else "?"
        self._dataset_panel.set_value(
            'file', file_name, tooltip=str(h5ad_path) if h5ad_path else None)

        sample_col = self.ws.sample_col
        if sample_col and sample_col in adata.obs.columns:
            n_samples = adata.obs[sample_col].nunique()
            self._dataset_panel.set_value(
                'sample_col', f"{sample_col} ({n_samples:,})", ok=True)
        else:
            self._dataset_panel.set_value('sample_col', "Not set", ok=False)

        condition_col = self.ws.condition_col
        self._roles_panel.set_value(
            'condition_col', condition_col or "Not set", ok=bool(condition_col))
        if condition_col and condition_col in adata.obs.columns:
            counts = adata.obs[condition_col].astype(str).value_counts()
            control, disease = self.ws.control_label, self.ws.disease_label
            if control:
                n = int(counts.get(str(control), 0))
                self._roles_panel.set_value(
                    'control', f"{control} ({n:,} cells)", ok=n > 0)
            else:
                self._roles_panel.set_value('control', "Not set", ok=False)
            if disease:
                n = int(counts.get(str(disease), 0))
                self._roles_panel.set_value(
                    'disease', f"{disease} ({n:,} cells)", ok=n > 0)
            else:
                self._roles_panel.set_value('disease', "Not set", ok=False)
        else:
            self._roles_panel.set_value('control', "—")
            self._roles_panel.set_value('disease', "—")

        self._hero_stats['cells'].set_value(f"{adata.n_obs:,}")
        self._hero_stats['genes'].set_value(f"{adata.n_vars:,}")
        if sample_col and sample_col in adata.obs.columns:
            self._hero_stats['samples'].set_value(f"{adata.obs[sample_col].nunique():,}")
        else:
            self._hero_stats['samples'].set_value("—")

    # ------------------------------------------------------------------
    # Column config sync
    # ------------------------------------------------------------------

    def _on_sample_col_changed(self, text):
        self.ws.sample_col = text
        self._update_info_panels()

    def _on_condition_col_changed(self, col_name):
        self.ws.condition_col = col_name
        if self.ws.current_adata is None or not col_name:
            return

        self.control_combo.clear()
        self.disease_combo.clear()

        unique_values = self.ws.current_adata.obs[col_name].unique().tolist()
        for val in unique_values:
            self.control_combo.addItem(str(val))
            self.disease_combo.addItem(str(val))

        self.combine_conditions_btn.setEnabled(len(unique_values) > 2)

        # Special case: '_role' already holds role names directly,
        # so just pair "control" / "disease" with the matching combo
        if col_name == '_role':
            for val in unique_values:
                if str(val).lower() == 'control':
                    idx = self.control_combo.findText(str(val))
                    if idx >= 0:
                        self.control_combo.setCurrentIndex(idx)
                elif str(val).lower() == 'disease':
                    idx = self.disease_combo.findText(str(val))
                    if idx >= 0:
                        self.disease_combo.setCurrentIndex(idx)
            self._apply_role_lock_state(self.ws.current_adata)
            return

        # Prefer adata.uns['role_map'] when the Inspect tab has set it.
        # Falls back to string-pattern matching only when the dataset
        role_map = self.ws.current_adata.uns.get('role_map') if self.ws.current_adata is not None else None
        if isinstance(role_map, dict) and role_map:
            from kosmic.scrna.inspect.roles import roles_to_labels
            disease_label, control_label = roles_to_labels(role_map)
            if control_label is not None:
                idx = self.control_combo.findText(str(control_label))
                if idx >= 0:
                    self.control_combo.setCurrentIndex(idx)
            if disease_label is not None:
                idx = self.disease_combo.findText(str(disease_label))
                if idx >= 0:
                    self.disease_combo.setCurrentIndex(idx)
            self._apply_role_lock_state(self.ws.current_adata)
            return

        # Fallback: string-pattern auto-detect when no role_map.
        control_patterns = ['control', 'ctrl', 'sham', 'healthy', 'donor', 'normal', 'wt', 'nf']
        disease_patterns = ['disease', 'tumor', 'dcm', 'hf', 'mcao', 'treated', 'patient', 'nicm', 'cm', 'hcm']

        for i, val in enumerate(unique_values):
            val_lower = str(val).lower()
            for pattern in control_patterns:
                if pattern in val_lower:
                    self.control_combo.setCurrentIndex(i)
                    break
            for pattern in disease_patterns:
                if pattern in val_lower:
                    self.disease_combo.setCurrentIndex(i)
                    break

        self._apply_role_lock_state(self.ws.current_adata)

    def _on_control_changed(self, text):
        self.ws.control_label = text
        self._update_info_panels()

    def _on_disease_changed(self, text):
        self.ws.disease_label = text
        self._update_info_panels()

    # ------------------------------------------------------------------
    # Combine conditions dialog
    # ------------------------------------------------------------------

    def _show_combine_conditions_dialog(self):
        if self.ws.current_adata is None:
            return

        condition_col = self.condition_col_combo.currentText()
        if not condition_col:
            dialogs.warning(self, "Error", "Select a condition column first")
            return

        unique_conditions = sorted(self.ws.current_adata.obs[condition_col].unique().tolist())

        dialog = QDialog(self)
        dialog.setWindowTitle("Combine Conditions")
        dialog.setMinimumWidth(400)
        dlayout = QVBoxLayout(dialog)

        instructions = QLabel(
            "Select conditions to combine into a single group.\n"
            "This creates a new column for your analysis."
        )
        instructions.setWordWrap(True)
        dlayout.addWidget(instructions)

        dlayout.addWidget(QLabel(f"Current conditions in '{condition_col}':"))
        checkbox_group = QGroupBox()
        checkbox_layout = QVBoxLayout(checkbox_group)
        checkboxes = {}
        for cond in unique_conditions:
            cb = QCheckBox(f"{cond} ({(self.ws.current_adata.obs[condition_col] == cond).sum():,} cells)")
            checkboxes[cond] = cb
            checkbox_layout.addWidget(cb)
        dlayout.addWidget(checkbox_group)

        name_layout = QHBoxLayout()
        name_layout.addWidget(QLabel("New group name:"))
        name_edit = QLineEdit()
        name_edit.setPlaceholderText("e.g., HF")
        name_layout.addWidget(name_edit)
        dlayout.addLayout(name_layout)

        col_layout = QHBoxLayout()
        col_layout.addWidget(QLabel("New column name:"))
        col_edit = QLineEdit()
        col_edit.setText(f"{condition_col}_grouped")
        col_layout.addWidget(col_edit)
        dlayout.addLayout(col_layout)

        btn_layout = QHBoxLayout()
        apply_btn = QPushButton("Apply")
        cancel_btn = QPushButton("Cancel")
        btn_layout.addWidget(apply_btn)
        btn_layout.addWidget(cancel_btn)
        dlayout.addLayout(btn_layout)

        def apply_grouping():
            selected = [cond for cond, cb in checkboxes.items() if cb.isChecked()]
            new_name = name_edit.text().strip()
            new_col = col_edit.text().strip()

            if len(selected) < 2:
                dialogs.warning(dialog, "Error", "Select at least 2 conditions to combine")
                return
            if not new_name:
                dialogs.warning(dialog, "Error", "Enter a name for the combined group")
                return
            if not new_col:
                dialogs.warning(dialog, "Error", "Enter a name for the new column")
                return

            mapping = {}
            for cond in unique_conditions:
                mapping[cond] = new_name if cond in selected else cond

            self.ws.current_adata.obs[new_col] = self.ws.current_adata.obs[condition_col].map(mapping)

            if self.condition_col_combo.findText(new_col) == -1:
                self.condition_col_combo.addItem(new_col)

            idx = self.condition_col_combo.findText(new_col)
            self.condition_col_combo.setCurrentIndex(idx)

            self.ws.log_message.emit(f"Created new column '{new_col}': {', '.join(selected)} -> {new_name}")

            # Save adata so the grouped column persists across sessions
            h5ad_path = getattr(self.ws, 'h5ad_path', None)
            if h5ad_path:
                try:
                    self.ws.current_adata.write_h5ad(h5ad_path)
                    self.ws.log_message.emit(f"Saved updated data to {h5ad_path}")
                except Exception as e:
                    self.ws.log_message.emit(f"Warning: could not save — {e}")

            dialog.accept()

        apply_btn.clicked.connect(apply_grouping)
        cancel_btn.clicked.connect(dialog.reject)

        dialog.exec()

    def refresh_theme(self):
        """No-op: every widget on this page is styled via theme roles."""
