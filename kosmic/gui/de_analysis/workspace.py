# DE Analysis Workspace
# Top-level workspace for differential expression analysis.
# Operates as a standalone mode in the NavigationRail, with its own
# sidebar workflow steps.
#
# Can receive adata from the scRNA pipeline or load standalone h5ad files.

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QStackedWidget, QFileDialog,
)
from PyQt6.QtCore import QSettings, pyqtSignal
from pathlib import Path

from kosmic import DE_MIN_CELLS, DE_MIN_COUNTS
from kosmic.gui.shared import borderless


class DEWorkspace(QWidget):
    """
    Embeddable workspace for DE analysis mode.

    Owns a QStackedWidget with workflow step pages.
    Communicates with AppWindow via signals.

    Supports two modes:
    - "scoring": select gene sets upfront, run Gene DE, then Pathway DE
    - "discovery": run genome-wide Gene DE, then enrichment analysis
    """

    status_message = pyqtSignal(str)
    log_message = pyqtSignal(str)
    step_completed = pyqtSignal(int)
    tab_changed = pyqtSignal(int)
    project_directory_changed = pyqtSignal(str)
    dataset_loaded = pyqtSignal(str)  # a study file was loaded here by hand (path)
    dataset_loading = pyqtSignal(str)  # about to load a study file by hand (path): release others first
    mode_changed = pyqtSignal(str)  # "scoring" or "discovery"
    study_change_requested = pyqtSignal(str)  # accession; AppWindow switches study
    run_changed = pyqtSignal()  # another gene DE run was chosen in the run bar (ADR-016)
    open_scrna_requested = pyqtSignal()  # AppWindow switches to scRNA -> Inspect

    # Sidebar step list per analysis mode.
    STEPS_INITIAL = (
        ("Load Data", "Load h5ad file, assign sample and condition columns."),
        ("Choose Analysis", "Select scoring or discovery mode."),
    )
    STEPS_SCORING = (
        ("Load Data", "Load h5ad file, assign sample and condition columns."),
        ("Choose Analysis", "Select scoring or discovery mode."),
        ("Select Gene Sets", "Choose pathways for differential expression analysis."),
        ("Gene Differential Expression", "Pseudobulk gene-level DE -- run first; the pathway analysis builds on its results."),
        ("Pathway Differential Expression", "Pathway-level activity scoring and testing."),
        ("Methods", "Auto-generated record of the analysis pipeline; copy-pasteable."),
    )
    STEPS_DISCOVERY = (
        ("Load Data", "Load h5ad file, assign sample and condition columns."),
        ("Choose Analysis", "Select scoring or discovery mode."),
        ("Gene Differential Expression", "Genome-wide pseudobulk differential expression."),
        ("Enrichment", "Pathway enrichment of significant genes."),
        ("Methods", "Auto-generated record of the analysis pipeline; copy-pasteable."),
    )

    # Sidebar step index -> stack page index. Stack pages after Set
    # Directory removal: 0=Setup, 1=Choose, 2=GeneSets, 3=GeneDE,
    # 4=PathDE, 5=Enrichment, 6=Methods. Scoring runs Gene DE (page 3)
    # before Pathway DE (page 4).
    _SCORING_MAP = (0, 1, 2, 3, 4, 6)
    _DISCOVERY_MAP = (0, 1, 3, 5, 6)
    _INITIAL_MAP = (0, 1)
    # Stack pages that show the run bar (ADR-016): Gene DE, Pathway DE,
    # Enrichment and Methods.
    _RUN_BAR_PAGES = (3, 4, 5, 6)

    @classmethod
    def steps_for_mode(cls, mode):
        """Return the sidebar step list for the given analysis_mode."""
        if mode == 'scoring':
            return cls.STEPS_SCORING
        if mode == 'discovery':
            return cls.STEPS_DISCOVERY
        return cls.STEPS_INITIAL

    def __init__(self, parent=None):
        super().__init__(parent)
        self.settings = QSettings("KOSMIC", "KOSMIC")
        self._step_completed = {}

        # ------------------------------------------------------------------
        # Shared state — all pages read/write these via self (workspace)
        # ------------------------------------------------------------------

        # Data
        self.current_adata = None
        self.h5ad_path = None
        self.project_dir = None
        # (resolved file path, count source) of the file the DE results were
        # made from, or None. Releasing the dataset keeps it, so reloading
        # the same file with the same counts keeps the results (ADR-014).
        self.results_dataset = None
        # Which counts a study is loaded with: None for the original
        # counts, or a layer ('decontX_counts'). Chosen on the Dataset
        # page; kept for the session, so each study loaded next uses it.
        self.counts_choice = None
        # scRNA cross-push bookkeeping: the scrna._adata_version last
        # synced into this workspace, and whether the user has since
        # loaded their own file (which then takes precedence -- an
        # auto cross-push must not clobber an explicit import).
        self._last_scrna_adata_version = None
        self._manual_load = False

        # Config
        self.sample_col = ''
        self.condition_col = ''
        self.control_label = ''
        self.disease_label = ''
        # From config.toml, not hardcoded: 100 silently overrode the
        # documented default of 10 and excluded whole cell types (7 of 13
        # testable on the DCM master, against 12 at the config value).
        self.min_cells = DE_MIN_CELLS
        self.min_counts = DE_MIN_COUNTS   # summed transcripts per donor; 0 = off
        # Gene universe for BH correction and the result table: None
        # (every gene the fit kept) or 'shared_atlas' (the genes the
        # project's shared atlas holds, so per-study results are tested
        # over the same list as the mega-analysis on the atlas).
        self.gene_universe = None
        self.de_method = 'deseq2'
        self.unit = 'sample'   # 'sample' (pseudobulk) or 'cell' (scanpy Wilcoxon)
        self.geneset_label = 'Metabolic -- Comprehensive (12 pathways)'

        # Empty until the Gene Sets page (Scoring mode) or the Enrichment
        # page (Discovery mode) populates it.
        self.pathway_gene_sets = {}
        self.current_hierarchy = None
        self.parent_gene_sets = None

        # The study's gene DE runs and the one the steps after gene DE use
        # (ADR-016): 'kosmic.de.batch.CellTypeRun's, a whole-study run (no
        # cell type) first. Chosen in the run bar.
        self.cell_type_runs = []
        self.run_index = -1
        self._run_cells = None   # ((dataset id, run index), the run's cells)

        # Results
        self.de_results = None
        self.significant_genes = None
        # (fdr, lfc_threshold, assumed) of the run shown: the rule its counts
        # and files use (ADR-015). None until a run is made or loaded.
        self.significance = None
        self.pathway_coverage = None
        self.sample_df = None
        self.pathway_de_results = None  # Meta-analysis compatible (same format as gene DE)
        self.parent_stats_df = None
        self.enrichment_results = None  # DataFrame from EnrichmentPage (discovery mode)
        self.fgsea_results = None  # preranked GSEA of the a-priori panel (Pathway DE page)
        self.gsea_library_results = None  # preranked GSEA of a library (Enrichment page)
        self.gsea_library_name = None  # display name of that library (for figure titles)
        self.gsea_library_details = None  # running-enrichment curves per set (mountain plot)

        # Analysis mode: None (not chosen), "scoring", or "discovery"
        self.analysis_mode = None

        # UI (wired externally by main.py)
        self.progress_bar = None

        self._setup_ui()

    @property
    def gse_accession(self) -> str:
        """Accession used to label exported files. Derived from the active study folder name."""
        return self.project_dir.name if self.project_dir else ''

    def _setup_ui(self):
        layout = borderless(QVBoxLayout, self)

        # The run bar above every page (ADR-016); shown by refresh_run_bar.
        from kosmic.gui.de_analysis.run_bar import RunBar
        self.run_bar = RunBar(self)
        self.run_bar.run_chosen.connect(self.select_run)
        self.run_bar.hide()
        layout.addWidget(self.run_bar)

        self.stack = QStackedWidget()

        # Import pages lazily to avoid circular imports
        from kosmic.gui.de_analysis.pages.setup_page import SetupPage
        from kosmic.gui.de_analysis.pages.gene_set_page import GeneSetPage
        from kosmic.gui.de_analysis.pages.gene_de_page import GeneDEPage
        from kosmic.gui.de_analysis.pages.pathway_de_page import PathwayDEPage
        from kosmic.gui.de_analysis.pages.enrichment_page import EnrichmentPage
        from kosmic.gui.shared.theme import get_color
        from kosmic.gui.shared.widgets import ModeChooserPage, ModeOption

        # Page 0: Load Data
        self.setup_page = SetupPage(self)
        self.stack.addWidget(self.setup_page)

        # Page 1: Choose Analysis Mode
        self.choose_page = ModeChooserPage(
            title="Choose your analysis approach",
            subtitle=("Both approaches run full-genome DE. "
                      "You can switch modes later."),
            options=[
                ModeOption(
                    key="scoring", icon_name="flask-conical",
                    title="Hypothesis-driven",
                    subtitle=("Test specific pathways or gene sets for "
                              "differential activity between conditions."),
                    steps=["Select gene sets",
                           "Gene DE (pseudobulk)",
                           "Pathway-level DE and scoring",
                           "Publication figures"],
                    accent_color=get_color('accent_primary'),
                ),
                ModeOption(
                    key="discovery", icon_name="microscope",
                    title="Discovery",
                    subtitle=("Find which biological processes are affected "
                              "without prior gene set selection."),
                    steps=["Gene DE (full genome)",
                           "Enrichment analysis",
                           "Publication figures"],
                    accent_color=get_color('error'),
                ),
            ],
        )
        self.choose_page.mode_selected.connect(self._on_mode_selected)
        self.stack.addWidget(self.choose_page)

        # Page 2: Select Gene Sets (scoring mode)
        self.gene_set_page = GeneSetPage(self)
        self.stack.addWidget(self.gene_set_page)

        # Page 3: Gene DE (both modes)
        self.gene_de_page = GeneDEPage(self)
        self.stack.addWidget(self.gene_de_page)

        # Page 4: Pathway DE (scoring mode)
        self.pathway_de_page = PathwayDEPage(self)
        self.stack.addWidget(self.pathway_de_page)

        # Page 5: Enrichment (discovery mode)
        self.enrichment_page = EnrichmentPage(self)
        self.stack.addWidget(self.enrichment_page)

        # Page 6: Methods (auto-generated record of the analysis chain)
        from kosmic.gui.de_analysis.pages.methods_page import MethodsPage
        self.methods_page = MethodsPage(self)
        self.stack.addWidget(self.methods_page)

        # Gene DE complete → mode-aware step marking + cascade
        self.gene_de_page.analysis_complete.connect(self._on_gene_de_complete)

        # Enrichment complete → mark step
        self.enrichment_page.enrichment_complete.connect(self._on_enrichment_complete)

        layout.addWidget(self.stack, 1)

    # ------------------------------------------------------------------
    # Navigation
    # ------------------------------------------------------------------

    def get_page_map(self):
        """Return the page map for the current analysis mode."""
        if self.analysis_mode == 'scoring':
            return self._SCORING_MAP
        elif self.analysis_mode == 'discovery':
            return self._DISCOVERY_MAP
        else:
            return self._INITIAL_MAP

    def switch_tab(self, index: int):
        """Switch to the given sidebar step index via the mode-specific page map."""
        page_map = self.get_page_map()
        if 0 <= index < len(page_map):
            page = page_map[index]
            if 0 <= page < self.stack.count():
                self.stack.setCurrentIndex(page)
                if (page in self._RUN_BAR_PAGES and not self.cell_type_runs
                        and self.current_adata is not None and self.de_results is None):
                    self.gene_de_page._try_load_previous_results()
                self.refresh_run_bar()
                self.tab_changed.emit(index)
                widget = self.stack.widget(page)
                if hasattr(widget, 'on_activated'):
                    widget.on_activated()

    def on_sidebar_step(self, index: int) -> None:
        self.switch_tab(index)

    # ------------------------------------------------------------------
    # The run the steps after gene DE use (ADR-016)
    # ------------------------------------------------------------------

    @property
    def shown_run(self):
        """The chosen 'CellTypeRun', or None."""
        if 0 <= self.run_index < len(self.cell_type_runs):
            return self.cell_type_runs[self.run_index]
        return None

    def set_runs(self, runs, select=None) -> None:
        """List 'runs' in the run bar and choose one: the run at index
        'select', else the first that ran. Emits 'run_changed' when a run
        is chosen."""
        self.cell_type_runs = list(runs)
        ok = [i for i, r in enumerate(self.cell_type_runs) if r.status == 'ok']
        if select is None or select not in ok:
            select = ok[0] if ok else -1
        self.run_index = -1
        self.run_bar.set_runs(self.cell_type_runs, select)
        if select >= 0:
            self.select_run(select)
        else:
            self.refresh_run_bar()

    def select_run(self, index: int) -> None:
        """Choose the run at 'index' for every page; a run that did not run
        cannot be chosen."""
        runs = self.cell_type_runs
        if not (0 <= index < len(runs)) or runs[index].status != 'ok':
            return
        if index == self.run_index:
            return
        self.run_index = index
        self._run_cells = None
        self.run_bar.set_runs(runs, index)
        # Results of the steps after gene DE belong to the run they were made
        # on; the pages reload the new run's from its folder.
        self._clear_downstream()
        # The pages apply the run (Gene DE sets de_results and its rule)
        # before the bar's line is written from them.
        self.run_changed.emit()
        self.refresh_run_bar()

    def _clear_downstream(self) -> None:
        self.pathway_de_results = None
        self.parent_stats_df = None
        self.pathway_score_data = None
        self.fgsea_results = None
        self.enrichment_results = None
        self.gsea_library_results = None
        self.gsea_library_name = None
        self.gsea_library_details = None

    def run_record(self) -> dict:
        """The chosen run's saved '.json' record ({} without one)."""
        from kosmic.de.saved import read_run_record
        run = self.shown_run
        if run is None or run.de_path is None:
            return {}
        return read_run_record(run.de_path) or {}

    def run_counts(self) -> str:
        """The counts folder of the chosen run, else of the loaded counts."""
        run = self.shown_run
        if run is not None and run.counts:
            return run.counts
        return self.gene_de_page.counts_folder()

    def run_folder(self):
        """Where the analyses on the chosen run save and find their outputs
        (ADR-016): 'pathway_scoring/<counts>/[<cell type>/]', or None
        without a study folder."""
        from kosmic.paths import de_downstream_dir
        if not self.project_dir:
            return None
        run = self.shown_run
        return de_downstream_dir(self.project_dir, self.run_counts(),
                                 run.cell_type if run is not None else None)

    def run_saved_folder(self):
        """The folder to read the chosen run's saved outputs from: its own,
        or for a whole-study run the flat folder of earlier versions
        (kosmic.de.saved.downstream_dir)."""
        from kosmic.de.saved import downstream_dir
        if not self.project_dir:
            return None
        run = self.shown_run
        return downstream_dir(self.project_dir, self.run_counts(),
                              run.cell_type if run is not None else None)

    def run_identity(self) -> dict:
        """The chosen run's cell type and counts, for provenance records."""
        run = self.shown_run
        return {'cell_type': (run.cell_type or None) if run is not None else None,
                'counts': self.run_counts()}

    def run_cell_type_column(self):
        """The obs column holding the chosen run's cell type: from its record,
        else the first candidate column of the loaded dataset whose values
        include it; None when there is none."""
        run, adata = self.shown_run, self.current_adata
        if run is None or not run.cell_type or adata is None:
            return None
        col = self.run_record().get('cell_type_column')
        if col and col in adata.obs.columns:
            return col
        ranked, other = self.gene_de_page._cell_type_columns(adata)
        for c in [*ranked, *other]:
            if run.cell_type in set(adata.obs[c].astype(str).unique()):
                return c
        return None

    def run_adata(self):
        """The loaded dataset's cells of the chosen run (ADR-016): every cell
        for a whole-study run, the cell type's cells for a cell type's run.
        None when no dataset is loaded or the cell type's column cannot be
        found. A view, made once per run."""
        adata, run = self.current_adata, self.shown_run
        if adata is None or run is None or not run.cell_type:
            return adata
        cached = getattr(self, '_run_cells', None)
        if cached is not None and cached[0] == (id(adata), self.run_index):
            return cached[1]
        col = self.run_cell_type_column()
        if col is None:
            self.log_message.emit(
                f"No column of the loaded dataset holds the cell type "
                f"'{run.cell_type}'; its cells cannot be selected.")
            return None
        cells = adata[(adata.obs[col].astype(str) == run.cell_type).to_numpy()]
        self._run_cells = ((id(adata), self.run_index), cells)
        return cells

    def run_donor_floors(self) -> tuple:
        """(min_cells, min_counts) the chosen run kept donors with: from its
        record, else the workspace's settings."""
        rec = self.run_record()
        min_cells = rec.get('min_cells')
        min_counts = rec.get('min_counts')
        return (int(min_cells) if min_cells is not None else self.min_cells,
                int(min_counts) if min_counts is not None else self.min_counts)

    def clear_runs(self) -> None:
        self.cell_type_runs = []
        self.run_index = -1
        self._run_cells = None
        self.run_bar.set_runs([], -1)
        self.refresh_run_bar()

    def refresh_run_bar(self) -> None:
        """Show the bar on the pages after gene DE once there is a run, and
        write its line from the chosen run."""
        from kosmic.gui.de_analysis.run_bar import detail_text
        run = self.shown_run
        show = bool(self.cell_type_runs) and self.stack.currentIndex() in self._RUN_BAR_PAGES
        self.run_bar.setVisible(show)
        if run is None:
            self.run_bar.set_detail("")
            return
        if self.significant_genes is not None and self.de_results is run.de_results:
            run.n_significant = len(self.significant_genes)
        self.run_bar.set_detail(detail_text(
            run, method=run.method or self.gene_de_page.method_token(),
            counts=run.counts or self.gene_de_page.counts_folder(),
            rule=self.significance))

    def run_label(self) -> str:
        """The chosen run's name for status lines: a cell type, 'Whole
        study', or '' when there is no run."""
        from kosmic.gui.de_analysis.run_bar import run_name
        run = self.shown_run
        return run_name(run) if run is not None else ''

    # ------------------------------------------------------------------
    # Data loading (called from main.py when scRNA pipeline provides adata)
    # ------------------------------------------------------------------

    @staticmethod
    def _isolate_from_scrna(adata):
        """Wrap 'adata' in a new AnnData for DE's own use.

        Shares the big arrays (X, layers, obsm, varm) by reference --
        no data duplication -- but gets an independent 'obs' DataFrame.
        DE mutates 'obs' (numeric-condition coercion, the "Combine
        Conditions" derived column); without this, those mutations
        would land on the exact same object the scRNA workspace still
        owns and displays, silently corrupting its state.
        """
        import anndata
        return anndata.AnnData(
            X=adata.X, obs=adata.obs.copy(), var=adata.var,
            uns=adata.uns, obsm=adata.obsm, varm=adata.varm,
            layers=adata.layers, raw=adata.raw,
        )

    def release_dataset(self) -> bool:
        """Drop the dataset this workspace holds so it can be freed.

        Called when another workspace loads a different study: one
        dataset in memory at a time. Results tables are small and stay.
        """
        import gc
        had = self.current_adata is not None
        self.current_adata = None
        self.h5ad_path = None
        self._manual_load = False
        gc.collect()
        return had

    def set_adata(self, adata):
        """Receive adata from the scRNA pipeline."""
        if adata is not None:
            adata = self._isolate_from_scrna(adata)
            self._coerce_condition_columns_to_string(adata)
        self.current_adata = adata
        # scRNA's copy can differ from the file (unsaved steps), so results
        # made from it are not kept when the file is loaded later.
        self.results_dataset = None
        if adata is not None:
            self.mark_step_complete(0)  # Step 0 = Load Data
            self.status_message.emit(
                f"Data loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes"
            )
            if hasattr(self, 'setup_page'):
                h5ad_path = getattr(self, 'h5ad_path', None) or ''
                self.setup_page.show_loaded_state(
                    label=f"Data loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes",
                    detail=str(h5ad_path) if h5ad_path else "From scRNA pipeline",
                )

    def show_scrna_not_ready(self, reason: str) -> None:
        """scRNA has data loaded but it isn't ready to hand off to DE yet
        (unnormalised, or no Control/Disease roles assigned).

        Only affects the empty state: once DE already holds valid data,
        a later scRNA readiness check that fails doesn't retroactively
        clear it -- avoids discarding a working DE session over a
        transient edit happening back in scRNA.
        """
        if self.current_adata is not None:
            return
        if hasattr(self, 'setup_page'):
            self.setup_page.show_not_ready_state(reason)

    def set_project_directory(self, directory: str):
        """Set the project directory for exports."""
        new_dir = Path(directory)
        previous_dir = self.project_dir

        # If switching to a different project, reset sidebar
        # completion state so the previous dataset's ticks don't
        # linger on the new (empty) workflow -- Known Issue #34.
        # Do this before the state-clear block so we cover both the
        # "adata loaded" and "no adata yet" switching cases.
        if previous_dir is not None and previous_dir != new_dir:
            self._step_completed = {i: False for i in self._step_completed}

        # If switching to a different project while data is loaded,
        # clear the stale adata + downstream state so we don't save
        # the wrong dataset's DE results into the new project folder.
        if (
            previous_dir is not None
            and previous_dir != new_dir
            and self.current_adata is not None
        ):
            self.current_adata = None
            self.h5ad_path = None
            self.results_dataset = None
            self._last_scrna_adata_version = None
            self._manual_load = False
            self.clear_runs()
            self.de_results = None
            self.significant_genes = None
            self.significance = None
            self.pathway_coverage = None
            self.sample_df = None
            self.pathway_de_results = None
            self.parent_stats_df = None
            self.enrichment_results = None
            self.fgsea_results = None
            self.gsea_library_results = None
            self.gsea_library_name = None
            self.gsea_library_details = None
            if hasattr(self, 'setup_page'):
                self.setup_page.hide_loaded_state()
            self.log_message.emit(
                "Project changed: previous dataset cleared. Load the dataset for the new project."
            )

        self.project_dir = new_dir
        folder_name = new_dir.name
        self.status_message.emit(f"Project: {folder_name}")
        self.project_directory_changed.emit(directory)

        # Auto-advance to the Dataset page if no data loaded yet. Only
        # activate the page when this workspace is the one on screen:
        # its on_activated() loads the study's h5ad, and doing that on
        # every study switch while the user is in scRNA held a second
        # copy of the dataset in memory (11 GB on Reichart) for nothing.
        # When DE is opened later, _activate_de calls on_activated then.
        if self.current_adata is None:
            if self.isVisible():
                self.switch_tab(0)
            else:
                page_map = self.get_page_map()
                if page_map:
                    self.stack.setCurrentIndex(page_map[0])
                    self.tab_changed.emit(0)

    def open_project_dialog(self):
        """Show a folder picker to set the project directory."""
        directory = QFileDialog.getExistingDirectory(
            self, "Open Project Directory",
            self.settings.value("last_directory", "")
        )
        if directory:
            self.settings.setValue("last_directory", directory)
            self.set_project_directory(directory)

    # ------------------------------------------------------------------
    # Inter-page wiring
    # ------------------------------------------------------------------

    def _on_mode_selected(self, mode):
        """Called when user picks scoring or discovery on the Choose Analysis page."""
        self.analysis_mode = mode
        self.mode_changed.emit(mode)
        self.mark_step_complete(1)  # step 1 = Choose Analysis
        # Auto-advance to the next step (step 2 in both modes).
        self.switch_tab(2)

    def _on_gene_de_complete(self):
        """Called when Gene DE finishes -- mode-aware step marking."""
        if self.analysis_mode == 'scoring':
            self.mark_step_complete(3)  # scoring step 3 = Gene DE
        else:
            self.mark_step_complete(2)  # discovery step 2 = Gene DE

        # Export the chosen method's DE + pseudobulk for meta-analysis (background)
        self.gene_de_page._export_for_meta()

    def _on_enrichment_complete(self):
        """Called when enrichment analysis finishes (discovery mode)."""
        if self.analysis_mode == 'discovery':
            self.mark_step_complete(3)  # discovery step 3 = Enrichment

    # ------------------------------------------------------------------
    # Workflow progress
    # ------------------------------------------------------------------

    def mark_step_complete(self, step_index: int):
        self._step_completed[step_index] = True
        self.step_completed.emit(step_index)

    def iter_completed_steps(self):
        """Yield '(step_index, done)' for every workflow step."""
        return self._step_completed.items()

    def compute_sidebar_status(self, _tab_index: int = 0) -> tuple[str, str]:
        """
        Contextual guidance for the sidebar status panel for this workspace's current page.

        Dispatches on the internal stack page rather than the tab_index argument
        (kept for signal-signature parity with ScRNAWorkspace). Returns
        '(text, state)' where state is 'info' | 'success' | 'warning' | 'error'.
        """
        stack_page = self.stack.currentIndex()
        adata = self.current_adata

        if stack_page == 0:  # Load Data
            if adata is None:
                return ("Load an h5ad file", "warning")
            n, g = adata.n_obs, adata.n_vars
            if self.sample_col and self.condition_col:
                return (f"{n:,} cells, {g:,} genes. Columns configured", "success")
            return (f"{n:,} cells, {g:,} genes. Configure columns", "warning")

        if stack_page == 1:  # Choose Analysis
            if self.analysis_mode:
                label = "Hypothesis-driven" if self.analysis_mode == 'scoring' else "Discovery"
                return (f"Mode: {label}", "success")
            return ("Choose an analysis approach", "warning")

        if stack_page == 2:  # Gene Sets (scoring)
            n_pw = len(self.pathway_gene_sets) if self.pathway_gene_sets else 0
            if n_pw > 0:
                n_genes = len({g for gs in self.pathway_gene_sets.values() for g in gs})
                return (f"{n_pw} pathways, {n_genes} genes selected", "success")
            return ("Select gene sets", "warning")

        run = self.run_label()
        prefix = f"{run}: " if run else ""

        if stack_page == 3:  # Gene DE
            if self.de_results is not None:
                n_sig = len(self.significant_genes) if self.significant_genes is not None else 0
                return (f"{prefix}DE complete: {len(self.de_results):,} genes, "
                        f"{n_sig} significant", "success")
            method = "DESeq2" if self.de_method == 'deseq2' else "t-test"
            return (f"Ready to run {method} DE", "warning")

        if stack_page == 4:  # Pathway DE (scoring)
            if self.pathway_de_results is not None:
                return (f"{prefix}Pathway DE complete: "
                        f"{len(self.pathway_de_results)} pathways", "success")
            return ("Run pathway-level DE", "warning")

        if stack_page == 5:  # Enrichment (discovery)
            if self.enrichment_results is not None and not self.enrichment_results.empty:
                return (f"{prefix}Enrichment: {len(self.enrichment_results)} terms",
                        "success")
            return ("Run enrichment analysis", "warning")

        if stack_page == 6:  # Methods
            return ("Analysis methods record", "info")

        return ("", "info")

    # ------------------------------------------------------------------
    # Theme refresh
    # ------------------------------------------------------------------

    @staticmethod
    def _coerce_condition_columns_to_string(adata):
        """
        Cast numeric obs columns used for condition labels to 'str'.

        Combo boxes in the DE Setup page populate 'ws.control_label'
        and 'ws.disease_label' as 'str(val)'. If the underlying
        'obs' column is numeric (e.g. 'HTN' stored as Int64 with
        values 0/1), then 'df['HTN'] == '0'' returns an all-False
        mask -- and every downstream plot reads as blank.

        Coerce the role-flagged condition column (preferred) plus the
        canonical 'condition' name to string. Idempotent: skips columns
        that are already object/string/categorical dtype.
        """
        candidates = [adata.uns.get('role_condition_col'), 'condition']
        seen = set()
        for col in candidates:
            if not col or col in seen:
                continue
            seen.add(col)
            if col not in adata.obs.columns:
                continue
            dtype = adata.obs[col].dtype
            if dtype.kind in 'iuf':  # int / uint / float
                adata.obs[col] = adata.obs[col].astype(str)

    def refresh_theme(self):
        """Cascade theme refresh to child pages."""
        self.run_bar.refresh_theme()
        for page in [self.setup_page, self.choose_page, self.gene_set_page,
                     self.gene_de_page, self.pathway_de_page,
                     self.enrichment_page]:
            if hasattr(page, 'refresh_theme'):
                page.refresh_theme()
