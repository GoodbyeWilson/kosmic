# DE Analysis Workspace
# Top-level workspace for differential expression analysis.
# Operates as a standalone mode in the NavigationRail, with its own
# sidebar workflow steps.
#
# Can receive adata from the scRNA pipeline or load standalone h5ad files.

from PyQt6.QtWidgets import (
    QWidget, QHBoxLayout, QStackedWidget, QFileDialog,
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
    mode_changed = pyqtSignal(str)  # "scoring" or "discovery"
    study_change_requested = pyqtSignal(str)  # accession; AppWindow switches study
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

        # Results
        self.de_results = None
        self.significant_genes = None
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
        layout = borderless(QHBoxLayout, self)

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
                self.tab_changed.emit(index)
                widget = self.stack.widget(page)
                if hasattr(widget, 'on_activated'):
                    widget.on_activated()

    def on_sidebar_step(self, index: int) -> None:
        self.switch_tab(index)

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

    def set_adata(self, adata):
        """Receive adata from the scRNA pipeline."""
        if adata is not None:
            adata = self._isolate_from_scrna(adata)
            self._coerce_condition_columns_to_string(adata)
        self.current_adata = adata
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
            self._last_scrna_adata_version = None
            self._manual_load = False
            self.de_results = None
            self.significant_genes = None
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

        if stack_page == 3:  # Gene DE
            if self.de_results is not None:
                n_sig = len(self.significant_genes) if self.significant_genes is not None else 0
                return (f"DE complete: {len(self.de_results):,} genes, "
                        f"{n_sig} significant", "success")
            method = "DESeq2" if self.de_method == 'deseq2' else "t-test"
            return (f"Ready to run {method} DE", "warning")

        if stack_page == 4:  # Pathway DE (scoring)
            if self.pathway_de_results is not None:
                return (f"Pathway DE complete: {len(self.pathway_de_results)} pathways", "success")
            return ("Run pathway-level DE", "warning")

        if stack_page == 5:  # Enrichment (discovery)
            if self.enrichment_results is not None and not self.enrichment_results.empty:
                return (f"Enrichment: {len(self.enrichment_results)} terms", "success")
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
        for page in [self.setup_page, self.choose_page, self.gene_set_page,
                     self.gene_de_page, self.pathway_de_page,
                     self.enrichment_page]:
            if hasattr(page, 'refresh_theme'):
                page.refresh_theme()
