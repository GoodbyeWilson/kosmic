# Meta-analysis workspace: 8-page workflow embedded in AppWindow's stack.

from pathlib import Path
from kosmic import DEFAULT_FDR, MIN_STUDIES
from kosmic.gui.shared import dialogs
from kosmic.gui.shared.widgets import CaptionLabel, PrimaryButton, SecondaryButton

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QCheckBox, QComboBox, QFileDialog,
    QFrame, QScrollArea, QStackedWidget,
)
from PyQt6.QtCore import QSettings, Qt, pyqtSignal

from kosmic.paths import meta_output_dir
from kosmic.meta_analysis.io import (
    load_de_results,
    require_se,
    discover_de_results,
    meta_selection_folder,
    split_cell_types,
)

# Choices in the Select Studies cell-type box that are not cell types.
_ALL_ENTRIES = '__all__'
_WHOLE_STUDY = '__whole__'
_EXTERNAL = '__external__'


class MetaAnalysisWorkspace(QWidget):
    """
    Embeddable workspace for meta-analysis visualization.

    Contains all UI, state, and analysis logic.  Can be embedded directly
    in any parent layout (e.g. AppWindow's QStackedWidget) without creating
    a hidden QMainWindow.
    """

    status_message = pyqtSignal(str)
    log_message = pyqtSignal(str)  # routed to OutputPanel
    step_completed = pyqtSignal(int)
    data_status = pyqtSignal(str, str)  # (text, color) for sidebar status panel
    mode_changed = pyqtSignal(str)  # consensus mode: "exploratory" or "hypothesis"
    tab_changed = pyqtSignal(int)   # sidebar index of the newly active sub-page
    project_directory_changed = pyqtSignal(str)  # syncs sidebar file tree

    def __init__(self, parent=None):
        super().__init__(parent)

        self.settings = QSettings("KOSMIC", "KOSMIC")

        # Sidebar progress bar (wired externally by main.py)
        self.progress_bar = None
        self._step_completed = {}

        # Dataset state (populated by the Select Studies page)
        self.datasets = []
        self.dataset_labels = []

        # Pathway-level state (loaded from pathway_de_results.csv)
        self.pathway_datasets = []
        self.pathway_dataset_labels = []

        # Project folder (set by AppWindow.project_changed)
        self._project_folder = ""

        # DE method picked in View > Meta-Analysis Settings; used by
        # _rescan_project to choose which {accession}_DE_*.csv files
        # to load. Initial value loaded from QSettings; main.py rewrites
        # it whenever the MA Settings dialog is closed.
        from kosmic.gui.meta_analysis.dialogs.ma_settings_dialog import load_de_method
        self.de_method = load_de_method(self.settings)

        self._setup_ui()
        self._load_settings()

    # Stack page constants -- index into self._master_stack.
    #: Discovered-entry count above which nothing is
    #: pre-selected -- see _rescan_project.
    _AUTOSELECT_LIMIT = 4

    PAGE_SELECT           = 0
    PAGE_CHOOSE_MODE      = 1
    PAGE_PATHWAY_EXPLORER = 2
    PAGE_PATHWAY_MA       = 3
    PAGE_GENE_MA          = 4
    PAGE_METHODS_CMP      = 5
    PAGE_METHODS_RES      = 6
    PAGE_ENRICHMENT       = 7
    PAGE_VALIDATION       = 8
    PAGE_METHODS          = 9

    # Sidebar step -> stack page maps. Each list[i] == stack page that
    # the ith sidebar workflow step maps to in this mode.
    _INITIAL_MAP     = (PAGE_SELECT, PAGE_CHOOSE_MODE)
    _EXPLORATORY_MAP = (PAGE_SELECT, PAGE_CHOOSE_MODE, PAGE_GENE_MA,
                        PAGE_ENRICHMENT, PAGE_VALIDATION, PAGE_METHODS)
    _HYPOTHESIS_MAP  = (PAGE_SELECT, PAGE_CHOOSE_MODE,
                        PAGE_PATHWAY_EXPLORER, PAGE_PATHWAY_MA,
                        PAGE_GENE_MA, PAGE_METHODS)
    _COMPARISON_MAP  = (PAGE_SELECT, PAGE_CHOOSE_MODE,
                        PAGE_METHODS_CMP, PAGE_METHODS_RES)

    # Sidebar step list per consensus mode.
    STEPS_INITIAL = (
        ("Select Studies", "Review and confirm loaded DE results."),
        ("Choose Mode", "Hypothesis-driven or Exploratory."),
    )
    STEPS_EXPLORATORY = STEPS_INITIAL + (
        ("Discovery Meta-Analysis", "Pool gene effects across all studies."),
        ("Enrichment", "Pathway enrichment of the pooled gene list."),
        ("Validation", "Stability of the result, and outside evidence for it."),
        ("Methods", "Auto-generated record of what actually ran; copy-pasteable."),
    )
    STEPS_HYPOTHESIS = STEPS_INITIAL + (
        ("Select Pathways", "Choose which pathways to investigate."),
        ("Pathway Meta-Analysis", "Pathway-level pooling across studies."),
        ("Focused Meta-Analysis", "Gene-level drill-down within significant pathways."),
        ("Methods", "Auto-generated record of what actually ran; copy-pasteable."),
    )
    STEPS_COMPARISON = STEPS_INITIAL + (
        ("Individual Methods", "Evaluate each method standalone: per-method volcanos, overlap, meta-discovery yield."),
        ("Consensus Evaluation", "Build and compare method combinations. Is pooling worth it?"),
    )

    @classmethod
    def steps_for_mode(cls, mode):
        """Return the sidebar step list for the given mode."""
        if mode == 'exploratory':
            return cls.STEPS_EXPLORATORY
        if mode == 'hypothesis':
            return cls.STEPS_HYPOTHESIS
        if mode == 'comparison':
            return cls.STEPS_COMPARISON
        return cls.STEPS_INITIAL

    def _setup_ui(self):
        """Setup the user interface. See PAGE_* constants for stack indices."""
        outer_layout = QVBoxLayout(self)
        outer_layout.setContentsMargins(0, 0, 0, 0)

        self._master_stack = QStackedWidget()
        outer_layout.addWidget(self._master_stack)

        # Page 0: Import Studies
        self._master_stack.addWidget(self._build_import_studies_page())

        # Page 2: Choose Mode.
        # The "exploratory" key for Discovery is preserved for backward
        # compatibility with saved projects.
        from kosmic.gui.shared.theme import get_color
        from kosmic.gui.shared.widgets import ModeChooserPage, ModeOption
        self._choose_mode_page = ModeChooserPage(
            title="Choose your analysis approach",
            subtitle=("Three routes through the same pooled studies. "
                      "You can switch modes at any time."),
            options=[
                ModeOption(
                    key="hypothesis", icon_name="flask-conical",
                    title="Hypothesis-driven",
                    subtitle=("Test specific pathways. Pool methods within "
                              "each pathway and drill down to genes. "
                              "Pathway-subset FDR for more power."),
                    steps=["Select pathways of interest",
                           "Pathway-level consensus MA",
                           "Drill into genes within significant pathways",
                           "Per-pathway gene volcano + forest plots"],
                    accent_color=get_color('accent_primary'),
                ),
                ModeOption(
                    key="exploratory", icon_name="microscope",
                    title="Discovery",
                    subtitle=("Find DE genes genome-wide. Single method or "
                              "consensus across methods. Volcano, "
                              "enrichment, external validation."),
                    steps=["Genome-wide pooling (single method or consensus)",
                           "Volcano + Venn + Results table",
                           "Pathway enrichment on result",
                           "Validation (proteomics / Olink / GWAS)"],
                    accent_color=get_color('error'),
                ),
                ModeOption(
                    key="comparison", icon_name="bar-chart-2",
                    title="Methods Comparison",
                    subtitle=("Benchmark methods. Which method to trust on "
                              "its own? Is the consensus worth doing? "
                              "Per-method LOO, AUC, meta-discovery yield, "
                              "cross-dataset agreement."),
                    steps=["Run all methods in parallel (shared DE)",
                           "Per-method LOO replication + held-out AUC",
                           "Per-method meta-discovery counts",
                           "Pairwise cross-dataset agreement"],
                    accent_color=get_color('accent_primary'),
                ),
            ],
        )
        self._choose_mode_page.mode_selected.connect(
            self._on_consensus_mode_selected)
        self._master_stack.addWidget(self._choose_mode_page)

        # Page 3: PathwayExplorer (hypothesis mode pathway selection,
        # shared widget; the DE workspace uses the same one).
        from kosmic.gui.shared.pathway_explorer import PathwayExplorer
        self._pathway_explorer = PathwayExplorer(
            parent=self, context='meta')
        self._pathway_explorer.pathways_changed.connect(
            self._on_pathway_explorer_changed)
        self._pathway_explorer.confirmed.connect(
            self._on_pathway_confirmed)
        self._master_stack.addWidget(self._pathway_explorer)

        # Page 4: Pathway Consensus MA (multi-method)
        from kosmic.gui.meta_analysis.pages.pathway_ma_page import (
            PathwayMAPage)
        self._pw_ma_page = PathwayMAPage()
        self._pw_ma_page.log_message.connect(self._log)
        self._pw_ma_page.analysis_complete.connect(
            self._on_pathway_ma_analysis_complete)
        self._master_stack.addWidget(self._pw_ma_page)

        # Page 5: Consensus MA
        from kosmic.gui.meta_analysis.pages.gene_ma_page import GeneMAPage
        self._gene_ma_page = GeneMAPage()
        self._gene_ma_page.log_message.connect(self._log)
        # Live progress -> sidebar status panel during long runs.
        # On completion analysis_complete fires _update_data_status(5),
        # which restores the data-driven "Done: N/M significant" text.
        self._gene_ma_page.sidebar_status.connect(self.data_status.emit)
        self._gene_ma_page.analysis_complete.connect(
            self._on_gene_ma_analysis_complete)
        self._master_stack.addWidget(self._gene_ma_page)

        # Forward PathwayExplorer's default gene set (loaded during
        # its __init__ before the signal was connected). Must happen
        # after _gene_ma_page is created.
        initial_pathways = self._pathway_explorer.get_flat_pathways()
        if initial_pathways:
            self._gene_ma_page.set_hypothesis_pathways(initial_pathways)

        # Page 6: Individual Methods (Methods Comparison mode, step 1)
        from kosmic.gui.meta_analysis.pages.methods_comparison_page import (
            MethodsComparisonPage)
        self._methods_comparison_page = MethodsComparisonPage()
        self._methods_comparison_page.log_message.connect(self._log)
        self._master_stack.addWidget(self._methods_comparison_page)

        # Page 7: Consensus Evaluation (Methods Comparison mode, step 2)
        from kosmic.gui.meta_analysis.pages.methods_comparison_results_page import (
            MethodsComparisonResultsPage)
        self._methods_comparison_results_page = MethodsComparisonResultsPage()
        self._methods_comparison_results_page.log_message.connect(self._log)
        self._master_stack.addWidget(self._methods_comparison_results_page)

        # Page 7: Enrichment (Discovery step 4). Hosts the DE
        # workspace's EnrichmentPage so both arms run the same code.
        from kosmic.gui.meta_analysis.pages.enrichment_step import (
            MetaEnrichmentPage)
        self._enrichment_page = MetaEnrichmentPage()
        self._enrichment_page.log_message.connect(self._log)
        self._master_stack.addWidget(self._enrichment_page)

        # Page 8: Validation (Discovery step 5).
        from kosmic.gui.meta_analysis.pages.validation_page import (
            ValidationPage)
        self._validation_page = ValidationPage()
        self._validation_page.log_message.connect(self._log)
        self._validation_page.loo_complete.connect(
            self._on_loo_complete)
        self._master_stack.addWidget(self._validation_page)

        # Page 9: Methods (last step, as in the DE workspace). Renders
        # the provenance every stage writes as it runs.
        from kosmic.gui.meta_analysis.pages.methods_page import (
            MetaMethodsPage)
        self._methods_page = MetaMethodsPage()
        self._master_stack.addWidget(self._methods_page)

        self._consensus_mode = None  # set by choose page

        # Initial status
        self.status_message.emit("Ready — set a project directory to begin")

    # ------------------------------------------------------------------
    # Gene-level pages
    # ------------------------------------------------------------------

    def _build_import_studies_page(self) -> QWidget:
        """Page 1: choose which project studies feed the meta-analysis.

        ADR-001: studies with DE results are discovered from the project
        and all selected by default; external CSVs are the escape hatch."""
        self._all_gene_datasets = []
        self._all_pathway_datasets = []
        self._excluded_studies = set()
        #: Above this many discovered entries, nothing is pre-selected.
        self._selection_touched = False

        page = QWidget()
        outer = QVBoxLayout(page)
        outer.setContentsMargins(0, 0, 0, 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)

        content = QWidget()
        lay = QVBoxLayout(content)
        lay.setContentsMargins(40, 30, 40, 30)
        lay.setSpacing(16)

        # No leading stretch: with one above and one below, the page
        # centred itself, so every change in content height moved the
        # whole block. Anchor to the top and let the trailing stretch
        # take the slack.
        lay.addSpacing(8)

        # Title
        title = QLabel("Select Studies")
        title.setProperty("role", "header")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lay.addWidget(title)

        self._select_summary = CaptionLabel(
            "Studies with differential-expression results are discovered "
            "from the project automatically."
        )
        self._select_summary.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._select_summary.setWordWrap(True)
        lay.addWidget(self._select_summary)

        lay.addSpacing(8)

        # === Eligible studies card ===
        data_frame = QFrame()
        data_frame.setObjectName("info_panel")
        data_frame.setMaximumWidth(560)
        data_frame.setMinimumWidth(520)
        data_layout = QVBoxLayout(data_frame)
        data_layout.setContentsMargins(16, 14, 16, 14)
        data_layout.setSpacing(8)

        head_row = QHBoxLayout()
        data_title = QLabel("STUDIES WITH DE RESULTS")
        data_title.setProperty("role", "panel_title")
        head_row.addWidget(data_title)
        head_row.addStretch()
        self._rescan_btn = SecondaryButton("Rescan")
        self._rescan_btn.setToolTip(
            "Look for DE results again, e.g. after running DE on a study.")
        self._rescan_btn.clicked.connect(self._rescan_project)
        head_row.addWidget(self._rescan_btn)
        self.add_btn = SecondaryButton("Import external result...")
        self.add_btn.setToolTip("Add a DE results CSV from outside this project.")
        self.add_btn.clicked.connect(self._add_datasets)
        head_row.addWidget(self.add_btn)
        data_layout.addLayout(head_row)

        # A per-cell-type run writes one accession per cell type, so this
        # list went from a handful of studies to a couple of dozen rows.
        # Ticking them by hand, or worse leaving them all ticked, is how
        # you end up pooling fibroblasts with endothelium. Choosing a cell
        # type selects that type in every study and nothing else.
        tools_row = QHBoxLayout()
        tools_row.setSpacing(6)
        tools_row.addWidget(QLabel("Cell type:"))
        self._cell_type_combo = QComboBox()
        self._cell_type_combo.setToolTip(
            "Choose a cell type to select it in every study that has it.\n"
            "Pooling is only meaningful within one cell type. Types found\n"
            f"in fewer than {MIN_STUDIES} studies cannot be pooled and are "
            "greyed out.")
        self._cell_type_combo.currentIndexChanged.connect(
            self._on_cell_type_chosen)
        tools_row.addWidget(self._cell_type_combo, 1)

        self._select_all_btn = SecondaryButton("Select all")
        self._select_all_btn.clicked.connect(lambda: self._set_all_shown(True))
        tools_row.addWidget(self._select_all_btn)

        self._select_none_btn = SecondaryButton("Clear")
        self._select_none_btn.clicked.connect(lambda: self._set_all_shown(False))
        tools_row.addWidget(self._select_none_btn)
        data_layout.addLayout(tools_row)

        # A fixed-height scroller, so filtering changes what is listed
        # and nothing else. Rows used to sit straight in the layout, so
        # hiding some shrank the panel and shifted everything on screen.
        checks_host = QWidget()
        self._study_checks_layout = QVBoxLayout(checks_host)
        self._study_checks_layout.setContentsMargins(0, 0, 0, 0)
        self._study_checks_layout.setSpacing(4)
        # AlignTop rather than a trailing spacer: _render_study_checklist
        # empties this layout item by item and would eat the spacer.
        self._study_checks_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        checks_scroll = QScrollArea()
        checks_scroll.setWidgetResizable(True)
        checks_scroll.setFrameShape(QFrame.Shape.NoFrame)
        checks_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        checks_scroll.setProperty("role", "thin_scroll")
        checks_scroll.setFixedHeight(280)
        checks_scroll.setWidget(checks_host)
        data_layout.addWidget(checks_scroll)

        self._no_eligible_label = CaptionLabel(
            "No DE results found in this project yet. Run Differential "
            "Expression on a study first, or import an external result.")
        self._no_eligible_label.setWordWrap(True)
        data_layout.addWidget(self._no_eligible_label)

        data_row = QHBoxLayout()
        data_row.addStretch()
        data_row.addWidget(data_frame)
        data_row.addStretch()
        lay.addLayout(data_row)

        lay.addSpacing(16)

        # Next button
        next_row = QHBoxLayout()
        next_row.addStretch()
        self._import_next_btn = PrimaryButton("Continue →")
        self._import_next_btn.setMinimumHeight(32)
        self._import_next_btn.setMinimumWidth(200)
        def _confirm_and_continue():
            if self.datasets:
                self.mark_step_complete(0)  # Select Studies
            self.switch_page(1)  # advance to Choose Mode
        self._import_next_btn.clicked.connect(_confirm_and_continue)
        next_row.addWidget(self._import_next_btn)
        next_row.addStretch()
        lay.addLayout(next_row)

        lay.addStretch(2)

        scroll.setWidget(content)
        outer.addWidget(scroll)
        return page

    def _log(self, message: str):
        """Add message to log."""
        self.log_message.emit(message)

    # ------------------------------------------------------------------
    # Project directory
    # ------------------------------------------------------------------

    def _rescan_project(self):
        """Rescan the project directory for DE results and auto-load them."""
        if not self._project_folder:
            return

        discovered = discover_de_results(Path(self._project_folder))
        self._scan_signature = self._results_signature(discovered)
        if not discovered:
            self._log(f"No DE results found in {self._project_folder}")
            self.status_message.emit("No DE results found")
            return

        # Auto-load all DE results (prefer the 3 per-method CSVs)
        de_files = [e for e in discovered if e['file_type'] == 'de_results']
        gene_files = [e for e in discovered if e['file_type'] == 'gene']

        # Prefer DE results; fall back to GeneStats files
        to_load = de_files if de_files else gene_files

        if not to_load:
            self._log(f"No loadable files found in {self._project_folder}")
            return

        # Group by accession — pick method matching user selection
        by_accession: dict[str, list] = {}
        for entry in to_load:
            acc = entry['accession']
            by_accession.setdefault(acc, []).append(entry)

        # Map the user-selected DE method label to a filename pattern.
        de_method_files = {
            "DESeq2": "deseq2",
            "Welch CPM": "welch_cpm",
            "Welch Raw": "welch_raw",
        }
        preferred = de_method_files.get(self.de_method, 'deseq2')

        selected = []
        for acc, entries in by_accession.items():
            if entries[0]['file_type'] == 'de_results':
                # Pick the user-selected method, fall back to whatever is available
                match = next((e for e in entries if e.get('de_method') == preferred), None)
                if match:
                    selected.append(match)
                else:
                    selected.append(entries[0])
                    self._log(f"  {acc}: {preferred} not found, using {entries[0].get('de_method', '?')}")
            else:
                # GeneStats: prefer Norm
                norm = next((e for e in entries if e.get('variant') == 'Norm'), None)
                selected.append(norm or entries[0])

        # Also collect pathway DE results keyed by analysis_folder
        pathway_de_files = [e for e in discovered if e['file_type'] == 'pathway_de']
        pw_de_by_folder: dict[str, list] = {}
        for entry in pathway_de_files:
            pw_de_by_folder.setdefault(entry['analysis_folder'], []).append(entry)

        # Rebuild the discovered pool; keep externally imported results
        externals = [d for d in getattr(self, '_all_gene_datasets', [])
                     if d.get('external')]
        self._all_gene_datasets = []
        self._all_pathway_datasets = []

        for entry in selected:
            filepath = str(entry["path"])
            try:
                df = load_de_results(filepath)
                require_se(df, filepath)
                # For Welch EB files, use moderated SE if available
                if entry.get('de_method') in ('welch_ttest_eb', 'welch_cpm_eb') and 'se_moderated' in df.columns:
                    valid_mod = df['se_moderated'].notna() & (df['se_moderated'] > 0)
                    if valid_mod.any():
                        df.loc[valid_mod, 'se'] = df.loc[valid_mod, 'se_moderated']
                name = entry["accession"]
                self._all_gene_datasets.append({
                    'path': filepath, 'name': name, 'df': df,
                    'study': entry['analysis_folder']})
                self._log(f"Loaded: {name} ({len(df)} genes)")
            except Exception as e:
                self._log(f"Failed to load {filepath}: {e}")
        by_folder: dict[str, list] = {}
        for entry in to_load:
            by_folder.setdefault(entry['analysis_folder'], []).append(
                entry['accession'])
        cell_types = split_cell_types(by_folder)
        for d in self._all_gene_datasets:
            d['cell_type'] = cell_types.get(d['name'])
        self._all_gene_datasets.extend(externals)

        # Load pathway DE datasets separately
        for folder, entries in pw_de_by_folder.items():
            # Prefer base pathway_de_results.csv over parent_
            base = [e for e in entries if 'parent_' not in e['path'].name]
            entry = base[0] if base else entries[0]
            filepath = str(entry["path"])
            try:
                pw_df = load_de_results(filepath, dataset_name=entry["accession"])
                require_se(pw_df, filepath)
                name = entry["accession"]
                self._all_pathway_datasets.append({'path': filepath, 'name': name, 'df': pw_df})
                self._log(f"Loaded pathway DE: {name} ({len(pw_df)} pathways)")
            except Exception as e:
                self._log(f"Failed to load pathway DE {filepath}: {e}")

        # Drop exclusions for studies that no longer exist
        names = {d['name'] for d in self._all_gene_datasets}
        self._excluded_studies &= names

        # Pre-ticking everything was fine when a project held three or
        # four studies. Once per-cell-type DE writes an accession per
        # cell type it is an actively wrong default: running as-loaded
        # would pool fibroblasts with endothelium. Above a handful of
        # entries, start from nothing and make the choice deliberate.
        if len(names) > self._AUTOSELECT_LIMIT and not self._selection_touched:
            self._excluded_studies = set(names)

        self._render_study_checklist()
        self._apply_study_selection()
        self._log(f"Scanned {self._project_folder}: "
                  f"{len(self._all_gene_datasets)} dataset(s) from "
                  f"{len(by_accession)} studies")

    @staticmethod
    def _results_signature(discovered) -> frozenset:
        """Which result files exist and when each last changed."""
        sig = set()
        for e in discovered:
            try:
                sig.add((str(e['path']), e['path'].stat().st_mtime_ns))
            except OSError:
                continue
        return frozenset(sig)

    def on_activated(self):
        """Rescan when DE results were added or rewritten since the last scan.

        DE run in this session writes its results while the Meta workspace
        is off screen; without this, Select Studies still says it found
        none. A rescan reloads every result table, so it runs only when the
        set of files or their modification times changed.
        """
        if not self._project_folder:
            return
        discovered = discover_de_results(Path(self._project_folder))
        if self._results_signature(discovered) != getattr(
                self, '_scan_signature', None):
            self._rescan_project()

    def set_project_directory_external(self, folder: str):
        """Set project directory from AppWindow.project_changed."""
        if not folder or not Path(folder).is_dir():
            return
        self._project_folder = folder
        self.project_directory_changed.emit(folder)
        self.settings.setValue("meta/project_dir", folder)
        Path(folder, "meta_analysis").mkdir(exist_ok=True)
        self._log(f"Project directory set: {folder}")
        self.status_message.emit(f"Project: {Path(folder).name}")
        self._rescan_project()

    def _get_output_dir(self) -> Path | None:
        """Return the meta-analysis output directory, or None if not set."""
        if not self._project_folder:
            return None
        out = meta_output_dir(self._project_folder)
        out.mkdir(exist_ok=True)
        return out

    # ------------------------------------------------------------------
    # Dataset management
    # ------------------------------------------------------------------

    def _add_datasets(self):
        """Import an external DE results CSV (escape hatch; project
        studies are discovered automatically)."""
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Import External DE Results",
            "",
            "CSV Files (*.csv);;All Files (*)"
        )

        for filepath in files:
            try:
                df = load_de_results(filepath)
                require_se(df, filepath)
                name = Path(filepath).stem
                existing = next(
                    (d for d in self._all_gene_datasets if d['path'] == filepath),
                    None)
                if existing is not None:
                    existing['df'] = df
                    self._log(f"Reloaded: {name} ({len(df)} genes)")
                else:
                    self._all_gene_datasets.append(
                        {'path': filepath, 'name': name, 'df': df,
                         'external': True})
                    self._log(f"Imported external result: {name} ({len(df)} genes)")
            except Exception as e:
                dialogs.warning(self, "Error", f"Failed to load {filepath}:\n{str(e)}")

        self._render_study_checklist()
        self._apply_study_selection()

    def _remove_external(self, name: str):
        """Remove an externally imported result from the pool."""
        self._all_gene_datasets = [
            d for d in self._all_gene_datasets
            if not (d.get('external') and d['name'] == name)]
        self._excluded_studies.discard(name)
        self._render_study_checklist()
        self._apply_study_selection()

    def _on_study_check_toggled(self, name: str, checked: bool):
        self._selection_touched = True
        if checked:
            self._excluded_studies.discard(name)
        else:
            self._excluded_studies.add(name)
        self._apply_study_selection()

    def _render_study_checklist(self):
        """Rebuild the eligible-studies checkbox rows from the pool."""
        while self._study_checks_layout.count():
            item = self._study_checks_layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

        for d in self._all_gene_datasets:
            row = QWidget()
            row_layout = QHBoxLayout(row)
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.setSpacing(8)
            label = f"{d['name']}  ({len(d['df']):,} genes)"
            if d.get('external'):
                label += "  · external"
            cb = QCheckBox(label)
            cb.setChecked(d['name'] not in self._excluded_studies)
            cb.toggled.connect(
                lambda checked, n=d['name']:
                self._on_study_check_toggled(n, checked))
            row.setProperty("study_name", d['name'])
            row.setProperty("group", self._group_of(d))
            row_layout.addWidget(cb)
            row_layout.addStretch()
            if d.get('external'):
                rm = QPushButton("Remove")
                rm.setProperty("role", "row_action")
                rm.clicked.connect(
                    lambda _c=False, n=d['name']: self._remove_external(n))
                row_layout.addWidget(rm)
            self._study_checks_layout.addWidget(row)

        self._no_eligible_label.setVisible(not self._all_gene_datasets)
        self._populate_cell_type_combo()
        self._apply_study_filter()

    @staticmethod
    def _group_of(d) -> str:
        """The cell-type box choice an entry belongs to."""
        if d.get('external'):
            return _EXTERNAL
        return d.get('cell_type') or _WHOLE_STUDY

    def _populate_cell_type_combo(self):
        """Offer each cell type once, with the number of studies that have it.

        Keeps the current choice across a rescan when it still exists.
        Choices found in fewer than MIN_STUDIES studies are shown but
        disabled, since the pooling step would refuse them.
        """
        combo = self._cell_type_combo
        previous = combo.currentData()
        counts: dict[str, int] = {}
        for d in self._all_gene_datasets:
            g = self._group_of(d)
            counts[g] = counts.get(g, 0) + 1

        combo.blockSignals(True)
        combo.clear()
        combo.addItem(f"All entries ({len(self._all_gene_datasets)})",
                      _ALL_ENTRIES)
        special = {_WHOLE_STUDY: "Whole study", _EXTERNAL: "External imports"}
        order = ([g for g in (_WHOLE_STUDY,) if g in counts]
                 + sorted(g for g in counts if g not in special)
                 + [g for g in (_EXTERNAL,) if g in counts])
        model = combo.model()
        for g in order:
            n = counts[g]
            unit = "study" if n == 1 else "studies"
            if g == _EXTERNAL:
                unit = "result" if n == 1 else "results"
            combo.addItem(f"{special.get(g, g)}  ({n} {unit})", g)
            if g != _EXTERNAL and n < MIN_STUDIES:
                model.item(combo.count() - 1).setEnabled(False)
        idx = combo.findData(previous)
        combo.setCurrentIndex(idx if idx >= 0 else 0)
        combo.blockSignals(False)

    def _on_cell_type_chosen(self, _index=None):
        """Show the chosen group and select exactly it."""
        group = self._cell_type_combo.currentData()
        if group not in (None, _ALL_ENTRIES):
            self._selection_touched = True
            self._excluded_studies = {
                d['name'] for d in self._all_gene_datasets
                if self._group_of(d) != group}
            for i in range(self._study_checks_layout.count()):
                w = self._study_checks_layout.itemAt(i).widget()
                if w is None:
                    continue
                for cb in w.findChildren(QCheckBox):
                    cb.blockSignals(True)
                    cb.setChecked(w.property("study_name")
                                  not in self._excluded_studies)
                    cb.blockSignals(False)
            self._apply_study_selection()
        self._apply_study_filter()

    def _visible_study_rows(self):
        """Rows currently passing the filter, as (widget, name) pairs."""
        rows = []
        for i in range(self._study_checks_layout.count()):
            item = self._study_checks_layout.itemAt(i)
            w = item.widget() if item else None
            if w is None or w.isHidden():
                continue
            rows.append((w, str(w.property("study_name") or "")))
        return rows

    def _apply_study_filter(self, _index=None):
        """Show only the chosen cell type's rows, and relabel the buttons.

        The buttons act on what is shown rather than on everything, so a
        Select all never quietly ticks rows of another cell type.
        """
        group = self._cell_type_combo.currentData()
        filtered = group not in (None, _ALL_ENTRIES)
        n_shown = 0
        for i in range(self._study_checks_layout.count()):
            item = self._study_checks_layout.itemAt(i)
            w = item.widget() if item else None
            if w is None:
                continue
            match = not filtered or w.property("group") == group
            w.setVisible(match)
            n_shown += int(match)

        self._select_all_btn.setText(
            "Select shown" if filtered else "Select all")
        self._select_none_btn.setText("Clear shown" if filtered else "Clear")
        for btn in (self._select_all_btn, self._select_none_btn):
            btn.setEnabled(n_shown > 0)

    def _set_all_shown(self, checked: bool):
        """Tick or untick every row currently visible."""
        self._selection_touched = True
        for w, name in self._visible_study_rows():
            for cb in w.findChildren(QCheckBox):
                cb.setChecked(checked)
            if checked:
                self._excluded_studies.discard(name)
            else:
                self._excluded_studies.add(name)
        self._apply_study_selection()

    def _apply_study_selection(self):
        """Rebuild the working dataset lists from the pool minus exclusions."""
        excluded = self._excluded_studies
        self.datasets = [d for d in self._all_gene_datasets
                         if d['name'] not in excluded]
        self.dataset_labels = [d['name'] for d in self.datasets]
        self.pathway_datasets = [d for d in self._all_pathway_datasets
                                 if d['name'] not in excluded]
        self.pathway_dataset_labels = [d['name'] for d in self.pathway_datasets]

        if hasattr(self, '_pw_ma_page') and self._pw_ma_page is not None:
            self._pw_ma_page.set_datasets(
                self.pathway_datasets,
                gene_datasets=self.datasets,
                project_folder=self._project_folder)

        n_all = len(self._all_gene_datasets)
        n_sel = len(self.datasets)
        if n_all == 0:
            self._select_summary.setText(
                "Studies with differential-expression results are "
                "discovered from the project automatically.")
        elif n_sel == 0:
            self._select_summary.setText(
                f"{n_all} entries found · none selected. Pooling is only "
                f"meaningful within one cell type, so choose one in the "
                f"Cell type box.")
        else:
            head = (f"{n_all} studies available for meta-analysis · all selected"
                    if n_sel == n_all
                    else f"{n_sel} of {n_all} available entries selected")
            self._select_summary.setText(
                f"{head} · results go to meta_analysis/{self.output_selection}/")
        self.status_message.emit(f"{n_sel} of {n_all} dataset(s) selected")














    # ------------------------------------------------------------------
    # Page switching (driven by sidebar workflow steps)
    # ------------------------------------------------------------------
    def current_stack_page(self) -> int:
        """Index of the currently visible page in the master stack."""
        return self._master_stack.currentIndex()

    @property
    def project_folder(self) -> str:
        """Currently selected parent folder containing study subfolders."""
        return self._project_folder

    @property
    def output_selection(self) -> str:
        """Folder under meta_analysis/ for the current selection (ADR-007).

        Named after the selected project results' cell type; external
        imports have none and do not count.
        """
        return meta_selection_folder(
            [d.get('cell_type') for d in self.datasets if not d.get('external')])

    @property
    def selected_study_folders(self) -> list:
        """Study folders of the selected project results, for the Methods step."""
        return sorted({d['study'] for d in self.datasets
                       if not d.get('external') and d.get('study')})

    @property
    def consensus_mode(self) -> str:
        """Active consensus mode: 'exploratory' / 'hypothesis' / 'comparison' / None."""
        return self._consensus_mode

    def iter_completed_steps(self):
        """Yield '(step_index, done)' for every workflow step."""
        return self._step_completed.items()

    @property
    def gene_ma_page(self):
        """The consensus meta-analysis page (used by the MA plot dialog)."""
        return self._gene_ma_page

    def mark_step_complete(self, step_index: int):
        """Mark a workflow step as complete (shows tick in sidebar)."""
        self._step_completed[step_index] = True
        self.step_completed.emit(step_index)
        # Refresh status text for current page
        self._update_data_status(self.current_stack_page())

    def _update_data_status(self, stack_page: int):
        """
        Emit contextual guidance for the sidebar status panel.

        Parameter is a STACK page index (not sidebar step index):
        0=Set Directory, 1=Import Studies, 2=Choose Mode,
        3=Pathway Explorer, 4=Pathway MA, 5=Gene MA (Consensus),
        6=Methods Comparison, 7=Methods Comparison Results.
        """
        n_ds = len(self.datasets)

        if stack_page == self.PAGE_SELECT:
            if n_ds > 0:
                self.data_status.emit(
                    f"{n_ds} studies. Click Continue.", '#4CAF50')
            else:
                self.data_status.emit(
                    "Set project directory first.", '#FF9800')

        elif stack_page == self.PAGE_CHOOSE_MODE:
            if self._consensus_mode:
                self.data_status.emit(
                    f"Mode: {self._consensus_mode}.", '#4CAF50')
            else:
                self.data_status.emit(
                    "Choose a mode.", '#2196F3')

        elif stack_page == self.PAGE_PATHWAY_EXPLORER:
            pw = getattr(self._pathway_explorer,
                         '_current_pathways', {})
            if pw:
                n_g = len({g for genes in pw.values() for g in genes})
                self.data_status.emit(
                    f"{len(pw)} pathways, {n_g} genes.", '#4CAF50')
            else:
                self.data_status.emit(
                    "Select a gene set.", '#2196F3')

        elif stack_page == self.PAGE_PATHWAY_MA:
            pw_df = getattr(self._pw_ma_page, '_meta_df', None)
            if pw_df is not None and len(pw_df) > 0:
                n_sig = int((pw_df['fdr'] < DEFAULT_FDR).sum())
                self.data_status.emit(
                    f"Done: {n_sig}/{len(pw_df)} significant.",
                    '#4CAF50')
            else:
                self.data_status.emit(
                    "Run pathway meta-analysis.", '#2196F3')

        elif stack_page == self.PAGE_GENE_MA:
            meta_df = getattr(self._gene_ma_page, '_meta_df', None)
            if meta_df is not None and len(meta_df) > 0:
                fdr = meta_df['fdr'].values.astype(float)
                n_sig = int((fdr < DEFAULT_FDR).sum())
                self.data_status.emit(
                    f"Done: {n_sig}/{len(meta_df)} significant.",
                    '#4CAF50')
            else:
                self.data_status.emit(
                    "Configure and run.", '#2196F3')

        elif stack_page == self.PAGE_METHODS_CMP:
            if self._methods_comparison_page.get_method_results() is not None:
                self.data_status.emit(
                    "Per-method discovery complete.", '#4CAF50')
            else:
                self.data_status.emit(
                    "Run per-method discovery.", '#2196F3')

        elif stack_page == self.PAGE_METHODS_RES:
            if self._methods_comparison_page.get_method_results() is None:
                self.data_status.emit(
                    "Run Individual Methods first.", '#FF9800')
            else:
                self.data_status.emit(
                    "Build and compare consensus configurations.",
                    '#2196F3')

    def get_page_map(self):
        """Return the active sidebar-step -> stack-page map."""
        if self._consensus_mode == 'exploratory':
            return self._EXPLORATORY_MAP
        elif self._consensus_mode == 'hypothesis':
            return self._HYPOTHESIS_MAP
        elif self._consensus_mode == 'comparison':
            return self._COMPARISON_MAP
        return self._INITIAL_MAP

    def _on_gene_ma_analysis_complete(self):
        """
        Refresh sidebar status + tick the matching workflow step.

        The sidebar step that maps to 'GeneMAPage' differs per mode;
        in Comparison mode it isn't in the flow at all.
        """
        self._update_data_status(self.PAGE_GENE_MA)
        page_map = self.get_page_map()
        if self.PAGE_GENE_MA in page_map:
            self.mark_step_complete(page_map.index(self.PAGE_GENE_MA))

        # Enrichment and Validation are downstream steps now, so they
        # are told the run finished rather than reading this page.
        context = self._gene_ma_page.get_run_context()
        self._enrichment_page.set_meta_results(context.get('meta_df'))
        self._validation_page.set_run_context(context)

    def _on_loo_complete(self, result):
        """GeneMAPage keeps the per-mode result cache, so hand it back."""
        self._gene_ma_page.cache_loo_result(result)

    def _on_pathway_ma_analysis_complete(self):
        """Refresh sidebar status + tick the Pathway MA step (hypothesis mode)."""
        self._update_data_status(self.PAGE_PATHWAY_MA)
        page_map = self.get_page_map()
        if self.PAGE_PATHWAY_MA in page_map:
            self.mark_step_complete(page_map.index(self.PAGE_PATHWAY_MA))

    def _stack_page_name(self, stack_page: int) -> str:
        names = {
            0: "Import Studies",
            1: "Choose Mode",
            2: "Select Pathways",
            3: "Pathway-Level Meta-Analysis",
            4: "Discovery Meta-Analysis",
            5: "Individual Methods",
            6: "Consensus Evaluation",
            7: "Enrichment",
            8: "Validation",
            9: "Methods",
        }
        return names.get(stack_page, "")

    def switch_page(self, sidebar_index: int):
        """
        Switch to the workflow step at the given sidebar index.

        The sidebar index is translated through the active mode's page
        map to a stack page index.
        """
        page_map = self.get_page_map()
        if not (0 <= sidebar_index < len(page_map)):
            return
        stack_page = page_map[sidebar_index]
        self._master_stack.setCurrentIndex(stack_page)

        self.status_message.emit(self._stack_page_name(stack_page))
        self._update_data_status(stack_page)
        self.tab_changed.emit(sidebar_index)

        # Page-specific activation
        if stack_page == self.PAGE_CHOOSE_MODE:
            self._choose_mode_page.on_activated()
        elif stack_page == self.PAGE_METHODS:
            # Re-read on every visit: stages are written as they run, and
            # a study may have been reprocessed in another workspace.
            self._methods_page.set_project(
                self._project_folder, self.selected_study_folders,
                selection=self.output_selection)
        elif stack_page == self.PAGE_ENRICHMENT:
            self._enrichment_page.progress_bar = self.progress_bar
            self._enrichment_page.set_project_folder(
                self._project_folder, selection=self._run_selection())
        elif stack_page == self.PAGE_VALIDATION:
            self._validation_page.progress_bar = self.progress_bar
            self._validation_page.set_project_folder(
                self._project_folder, selection=self._run_selection())
            self._validation_page.set_studies(self.datasets)
        elif stack_page == self.PAGE_PATHWAY_MA:
            pw_gene_sets = getattr(self._pathway_explorer,
                                   '_current_pathways', {})
            self._pw_ma_page.set_datasets(
                self.pathway_datasets,
                gene_datasets=self.datasets,
                project_folder=self._project_folder,
                pathway_gene_sets=pw_gene_sets,
                output_selection=self.output_selection)
        elif stack_page == self.PAGE_GENE_MA:
            self._gene_ma_page.progress_bar = self.progress_bar
            # In hypothesis mode, also pass the pathway MA results so the
            # consensus page can pull significant pathways + their gene sets.
            pw_results = None
            pw_gene_sets = None
            if self._consensus_mode == 'hypothesis':
                pw_results = self._pw_ma_page.meta_df
                pw_gene_sets = self._pw_ma_page.pathway_gene_sets
            self._gene_ma_page.set_datasets(
                self.datasets, self.dataset_labels, self._project_folder,
                pathway_ma_results=pw_results,
                pathway_gene_sets=pw_gene_sets,
                output_selection=self.output_selection)
            if self._consensus_mode:
                self._gene_ma_page.set_mode(self._consensus_mode)
        elif stack_page == self.PAGE_METHODS_CMP:
            self._methods_comparison_page.progress_bar = self.progress_bar
            self._methods_comparison_page.set_datasets(
                self.datasets, self.dataset_labels, self._project_folder)
            try:
                from kosmic.gui.meta_analysis.pseudobulk_discovery import (
                    find_pseudobulk_paths)
                pb_paths = find_pseudobulk_paths(
                    self.datasets, self._project_folder, log_cb=self._log)
                self._methods_comparison_page.set_pseudobulk_paths(pb_paths)
            except Exception as e:
                self._log(f"pseudobulk path discovery failed: {e}")
        elif stack_page == self.PAGE_METHODS_RES:
            # Pull per-method results AND LOO results forward so the
            # trade-off frontier can use both cheap-metric and LOO metrics.
            self._methods_comparison_results_page.progress_bar = self.progress_bar
            method_results = self._methods_comparison_page.get_method_results()
            loo_results = self._methods_comparison_page.get_loo_results()
            self._methods_comparison_results_page.set_method_results(
                method_results, project_folder=self._project_folder,
                loo_results=loo_results)

    def _run_selection(self):
        """Output folder of the last pooling run, else the current selection."""
        return self._gene_ma_page.run_output_selection or self.output_selection

    def _on_pathway_explorer_changed(self, pathways_dict):
        """Forward selected pathways to the consensus page."""
        self._gene_ma_page.set_hypothesis_pathways(pathways_dict)
        n_pw = len(pathways_dict)
        n_genes = len({g for genes in pathways_dict.values()
                       for g in genes})
        if n_pw > 0:
            self.data_status.emit(
                f"{n_pw} pathways, {n_genes} genes selected. "
                f"Click Confirm to proceed.", '#4CAF50')
        else:
            self.data_status.emit(
                "Select a pathway gene set.", '#2196F3')

    def _on_pathway_confirmed(self):
        """User confirmed pathway selection — tick step and advance."""
        self.mark_step_complete(2)  # Select Pathways (hypothesis mode)
        pw = self._pathway_explorer.get_flat_pathways()
        n_pw = len(pw)
        n_genes = len({g for genes in pw.values() for g in genes})
        self.data_status.emit(
            f"Pathways confirmed: {n_pw} pathways, {n_genes} genes.",
            '#4CAF50')
        self.switch_page(3)  # advance to Pathway MA

    def _on_consensus_mode_selected(self, mode: str):
        """Called when the user picks a mode card."""
        self._consensus_mode = mode
        self._gene_ma_page.set_mode(mode)
        self.mark_step_complete(1)  # mark "Choose Mode" complete

        # Tell the host to rebuild the sidebar with the new step list
        self.mode_changed.emit(mode)

        # Auto-advance to the next sidebar step (index 2) in the new map.
        self.switch_page(2)

    def refresh_theme(self):
        """Refresh widgets that don't re-style from the global stylesheet."""
        if hasattr(self, '_pw_ma_page'):
            self._pw_ma_page.refresh_theme()

    def _load_settings(self):
        """Load saved settings and auto-scan if a directory was saved."""
        saved_dir = self.settings.value("meta/project_dir", "")
        if saved_dir and Path(saved_dir).is_dir():
            self._project_folder = saved_dir
            self._rescan_project()

