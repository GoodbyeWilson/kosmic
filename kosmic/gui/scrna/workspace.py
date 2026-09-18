# ScRNA Workspace.
#
# 'ScRNAWorkspace' is a QWidget containing the scRNA-seq pipeline (8 tabs:
# Load Data, Gene Names, Inspect, QC, Cluster, Marker Check, Decontaminate,
# Subset). It is embedded inside 'main.AppWindow''s QStackedWidget; tab
# navigation is driven externally via 'switch_tab()'.

from typing import Optional
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QStackedWidget, QFileDialog,
)
from PyQt6.QtCore import QSettings, pyqtSignal, QTimer

from kosmic.paths import processed_data_dir
from kosmic.gui.scrna.tabs.download_tab import DownloadTab
from kosmic.gui.scrna.tabs.convert_tab import ConvertTab
from kosmic.gui.scrna.tabs.qc_tab import QCTab
from kosmic.gui.scrna.tabs.cluster_tab import ClusterTab
from kosmic.gui.scrna.tabs.filter_tab import FilterTab
from kosmic.gui.scrna.tabs.gene_names_tab import GeneNamesTab
from kosmic.gui.scrna.tabs.gene_group_tab import GeneGroupTab
from kosmic.gui.scrna.tabs.decontx_tab import DecontXTab
from kosmic.gui.help.tutorial.controller import TutorialController


# ---------------------------------------------------------------------------
# Status bar proxy — lets tabs call self.main_window.status_bar.showMessage()
# without requiring a real QStatusBar.
# ---------------------------------------------------------------------------

class _StatusBarProxy:
    """Lightweight object that forwards showMessage() to a pyqtSignal."""

    def __init__(self, signal):
        self._signal = signal

    def showMessage(self, msg, timeout=0):
        self._signal.emit(msg)


# ---------------------------------------------------------------------------
# ScRNAWorkspace — embeddable QWidget owning the full 6-tab pipeline
# ---------------------------------------------------------------------------

class ScRNAWorkspace(QWidget):
    """
    Embeddable workspace containing the scRNA-seq 6-tab pipeline.

    Owns the QStackedWidget with 6 tabs (Load, Inspect, QC, Analysis,
    Filter, DE), all inter-tab signal wiring, workflow state, welcome
    dialog, and tutorial controller.

    WORKFLOW_STEPS defines the sidebar step list rendered by the
    explorer pane when this workspace is active.

    Tabs reference *this* object as 'self.main_window' and rely on
    'set_adata()', 'set_project_directory()', 'mark_step_complete()',
    'status_bar', 'current_adata', 'inspect_tab', etc.
    """

    WORKFLOW_STEPS = (
        ("Load Data", "Fetch from GEO, or import local files."),
        ("Gene Names", "Harmonise gene symbols to current HGNC standard."),
        ("Inspect", "Explore metadata, assign conditions and samples."),
        ("Quality Control", "Filter low-quality cells, detect doublets, normalize."),
        ("Cluster", "PCA, clustering, UMAP and cell type annotation."),
        ("Marker Check", "Check the annotation against known cell-type markers."),
        ("Decontaminate", "Remove ambient RNA per sample using cell-type labels."),
        ("Subset", "Subset to specific cell types for downstream analysis."),
    )

    status_message = pyqtSignal(str)
    project_directory_changed = pyqtSignal(str)
    step_completed = pyqtSignal(int)
    steps_reset = pyqtSignal()  # emitted when all downstream steps are cleared (reload raw)
    tab_changed = pyqtSignal(int)  # emitted when active tab changes (0-based tab index)
    dataset_loaded = pyqtSignal(str)  # a study file was loaded into this workspace (path)

    def __init__(self, embedded: bool = True, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.settings = QSettings("KOSMIC", "KOSMIC")
        self.current_project_dir = None
        self.current_adata = None
        self.current_h5ad_path = None
        self.raw_h5ad_path = None
        self._adata_version = 0
        self._step_completed = {i: False for i in range(len(self.WORKFLOW_STEPS))}
        self._embedded = embedded

        self.status_bar = _StatusBarProxy(self.status_message)

        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(5, 5, 5, 5)
        main_layout.setSpacing(0)
        self.stack = QStackedWidget()
        main_layout.addWidget(self.stack)

        self.download_tab = DownloadTab(self)
        self.gene_names_tab = GeneNamesTab(self)
        self.inspect_tab = ConvertTab(self)
        self.qc_tab = QCTab(self)
        self.cluster_tab = ClusterTab(self)
        self.annotate_tab = self.cluster_tab.annotate_tab
        self.gene_group_tab = GeneGroupTab(self)
        self.decontx_tab = DecontXTab(self)
        self.filter_tab = FilterTab(self)

        for w in (self.download_tab, self.gene_names_tab, self.inspect_tab,
                  self.qc_tab, self.cluster_tab, self.gene_group_tab,
                  self.decontx_tab, self.filter_tab):
            self.stack.addWidget(w)

        self.cluster_tab.pca_complete.connect(self._on_pca_complete)
        self.cluster_tab.clustering_complete.connect(self._on_clustering_complete)
        self.filter_tab.filtered_data_ready.connect(self._on_filtered_data_ready)

        self.tutorial_controller = TutorialController(self)
        self.status_message.emit("Ready")


    def centralWidget(self):
        """Compatibility shim — TutorialController uses this for its overlay parent."""
        return self

    # ------------------------------------------------------------------
    # Tab navigation
    # ------------------------------------------------------------------

    def switch_tab(self, index: int):
        """Switch to the given tab index (0=Load, ..., 5=Marker Check, 6=Decontaminate, 7=Subset)."""
        if 0 <= index < self.stack.count():
            self.stack.setCurrentIndex(index)
            self.tab_changed.emit(index)
            self.stack.widget(index).on_tab_activated()

    def on_sidebar_step(self, index: int) -> None:
        self.switch_tab(index)

    # ------------------------------------------------------------------
    # Public API used by tabs
    # ------------------------------------------------------------------

    def set_project_directory(self, directory: str):
        # On a genuine folder change (not a re-issue of the same path),
        # reset sidebar completion ticks and drop the previous dataset's
        # AnnData -- mid-pipeline tabs read 'current_adata' directly.
        previous = self.current_project_dir
        if previous is not None and previous != directory:
            self._step_completed = {i: False for i in self._step_completed}
            if self.current_adata is not None:
                self.release_dataset("Analysis folder changed - loaded dataset cleared.")

        self.current_project_dir = directory
        self.settings.setValue("last_directory", directory)

        # Notify all tabs that have set_project_directory method
        for tab in [self.download_tab, self.gene_names_tab, self.inspect_tab,
                    self.qc_tab, self.cluster_tab, self.annotate_tab,
                    self.gene_group_tab, self.decontx_tab, self.filter_tab]:
            if hasattr(tab, 'set_project_directory'):
                tab.set_project_directory(directory)

        self.status_message.emit(f"Project: {directory}")
        self.project_directory_changed.emit(directory)

        # Auto-advance to Load Data tab (unless set_adata already navigated to Inspect)
        if self.current_adata is None:
            self.switch_tab(0)

    def set_adata(self, adata, file_path=None, **_kw):
        """
        Set the current adata object and optionally its file path.

        Central state — tabs pull from here via on_tab_activated().
        """
        self.current_adata = adata

        # Track the current working file path
        first_load = False
        if file_path is not None:
            from pathlib import Path
            self.current_h5ad_path = str(Path(file_path))
            if adata is not None:
                self.dataset_loaded.emit(self.current_h5ad_path)
            # If no raw path set yet, this is the raw file
            if self.raw_h5ad_path is None:
                self.raw_h5ad_path = self.current_h5ad_path
                first_load = True

        # Record the initial load once, so the upstream dataset size is on
        # record (best-effort; the fingerprint captures the raw dimensions).
        if first_load and adata is not None:
            self.record_provenance('load', {
                'file': Path(self.current_h5ad_path).name,
                'n_cells': int(adata.n_obs),
                'n_genes': int(adata.n_vars),
            })

        self.notify_adata_changed()

        # Broadcast to gene_names_tab so it resets for new data
        self.gene_names_tab.set_adata(adata)

        # Steps are only marked complete when the user visits them.

    def record_provenance(self, stage, params, *, data_changed=True, source=None):
        """Append a processing stage to the study's provenance sidecar.

        Writes '<study_dir>/provenance.json' next to the working h5ad (never
        rewrites the h5ad itself). No-op when no working file exists yet.
        'data_changed' stamps a fresh data fingerprint for stages that alter the
        counts / obs, so downstream staleness checks see the new state. 'source'
        is the parent h5ad a derived study (a cell-type subset) came from, stored
        so the lineage can be followed back to the upstream processing.
        """
        from pathlib import Path
        h5ad = self.current_h5ad_path
        if not h5ad or self.current_adata is None:
            return
        try:
            from kosmic import provenance
            fp = (provenance.compute_fingerprint(self.current_adata, h5ad_path=h5ad)
                  if data_changed else None)
            provenance.record_stage(Path(h5ad).parent, Path(h5ad).stem,
                                    stage, params, fingerprint=fp, source=source)
        except Exception as e:  # provenance is best-effort, never blocks work
            self.status_message.emit(f"Could not record provenance: {e}")

    def reload_raw_data(self):
        """
        Reload the original raw h5ad file and reset all downstream state.

        Called from the Download tab's "Restore imported data" button.
        """
        if self.raw_h5ad_path is None:
            self.status_message.emit("No raw data file recorded — nothing to reload")
            return

        from pathlib import Path
        raw_path = Path(self.raw_h5ad_path)
        if not raw_path.exists():
            self.status_message.emit(f"Raw file not found: {raw_path}")
            return

        # Read original file
        import anndata
        import shutil
        try:
            adata = anndata.read_h5ad(str(raw_path))
        except Exception as e:
            self.status_message.emit(f"Failed to read raw file: {e}")
            return

        # Copy raw file into processed_data/ as the new working copy
        if self.current_project_dir:
            processed_dir = processed_data_dir(self.current_project_dir)
            processed_dir.mkdir(exist_ok=True)
            working_copy = processed_dir / raw_path.name
            shutil.copy2(str(raw_path), str(working_copy))
            working_path = str(working_copy)
        else:
            working_path = str(raw_path)

        # Clear all step completion -- raw reload restarts the pipeline.
        for i in range(len(self.WORKFLOW_STEPS)):
            self._step_completed[i] = False
        self.steps_reset.emit()  # notify sidebar to uncheck all

        # Reset each tab
        for tab in [self.gene_names_tab, self.inspect_tab, self.qc_tab,
                    self.cluster_tab, self.gene_group_tab, self.decontx_tab,
                    self.filter_tab]:
            tab.reset_state()

        # Set fresh adata with the working copy path
        self.current_h5ad_path = working_path
        self.current_adata = adata
        self._adata_version += 1

        # Notify gene names tab of new data
        self.gene_names_tab.set_adata(adata)

        self.status_message.emit(f"Reloaded raw data: {raw_path.name}")

        self.switch_tab(0)
        self.mark_step_complete(0)

    def release_dataset(self, message: str = "Dataset released from memory.") -> bool:
        """Drop the loaded dataset everywhere it is held, so it can be freed.

        Setting 'current_adata' to None was not enough: every tab keeps its
        own reference (the QC, Cluster and Annotate tabs, the UMAP widget),
        so a study stayed resident until the next one overwrote those
        references. Two studies in memory was the norm, and with the atlas
        open nothing else fitted. Each tab's reset_state drops its
        references and stops its workers; the collector then runs.
        Returns True when something was released.
        """
        import gc
        had = self.current_adata is not None
        self.current_adata = None
        self.current_h5ad_path = None
        self.raw_h5ad_path = None
        self._adata_version += 1
        for tab in [self.download_tab, self.gene_names_tab, self.inspect_tab,
                    self.qc_tab, self.cluster_tab, self.annotate_tab,
                    self.gene_group_tab, self.decontx_tab, self.filter_tab]:
            reset = getattr(tab, 'reset_state', None)
            if reset is not None:
                try:
                    reset()
                except Exception:  # noqa: BLE001 -- releasing must not fail
                    pass
        gc.collect()
        if had:
            self.status_message.emit(message)
        return had

    def notify_adata_changed(self):
        """
        Bump the adata version counter.

        Tabs refresh from central state in on_tab_activated(), which is
        called by switch_tab() or when sub-tabs change.  No immediate
        broadcast — tabs pick up changes when they next become visible.
        """
        self._adata_version += 1

    def open_project_dialog(self):
        directory = QFileDialog.getExistingDirectory(
            self, "Open Project Directory",
            self.settings.value("last_directory", "")
        )
        if directory:
            self.set_project_directory(directory)

    # ------------------------------------------------------------------
    # Workflow callbacks
    # ------------------------------------------------------------------

    def _on_pca_complete(self, adata):
        """PCA finished — update central state."""
        self.current_adata = adata
        self.notify_adata_changed()
        self.status_message.emit("PCA complete — ready to cluster")

    def _on_clustering_complete(self, file_path: str):
        self.mark_step_complete(4)  # Cluster step
        self.status_message.emit(f"Clustering complete: {file_path}")
        QTimer.singleShot(3000, self._suggest_next_step)

    def _on_filtered_data_ready(self, file_path: str):
        self.mark_step_complete(7)  # Subset step
        self.status_message.emit(f"Filtered data ready: {file_path}")
        QTimer.singleShot(3000, self._suggest_next_step)

    # ------------------------------------------------------------------
    # Workflow progress
    # ------------------------------------------------------------------

    def refresh_theme(self):
        """Re-apply theme to elements that can't use global stylesheet."""
        # Project picker chrome moved to the Project workspace.
        return

    def mark_step_complete(self, step_index: int):
        """Mark a workflow step as complete and emit step_completed signal."""
        self._step_completed[step_index] = True
        self.step_completed.emit(step_index)

    def iter_completed_steps(self):
        """Yield '(step_index, done)' for every workflow step."""
        return self._step_completed.items()

    def compute_sidebar_status(self, tab_index: int) -> tuple[str, str]:
        """
        Contextual guidance for the sidebar status panel at a given tab.

        Returns '(text, state)' where state is 'info' | 'success' | 'warning' | 'error'.
        Called by AppWindow's mode-switch / tab-changed / step-completed handlers.
        """
        adata = self.current_adata
        has_data = adata is not None

        if not has_data and tab_index > 0:
            return ("Load data first (step 1)", "error")

        if tab_index == 0:  # Load Data
            if has_data:
                return ("Data loaded. Proceed to Gene Names", "success")
            return ("Import or download a dataset", "warning")

        if tab_index == 1:  # Gene Names
            harmonised = (self._step_completed.get(1, False)
                          or getattr(self.gene_names_tab, '_harmonised', False))
            if harmonised:
                return (f"Gene names harmonised ({adata.n_vars:,} genes)", "success")
            return (f"{adata.n_vars:,} genes. Check for outdated symbols", "warning")

        if tab_index == 2:  # Inspect
            has_sample = 'sample' in adata.obs.columns
            has_condition = 'condition' in adata.obs.columns
            if has_sample and has_condition:
                return ("Sample + condition assigned. Ready for QC", "success")
            return ("Assign sample and condition columns", "warning")

        if tab_index == 3:  # QC
            is_norm = hasattr(adata.X, 'max') and float(adata.X.max()) < 20
            qc_done = 'n_genes_by_counts' in adata.obs.columns
            if is_norm and qc_done:
                return (f"QC + normalization done. {adata.n_obs:,} cells", "success")
            if qc_done:
                return ("QC done. Normalize before clustering", "warning")
            return ("Filter, remove doublets, then normalize", "warning")

        if tab_index == 4:  # Cluster
            is_norm = hasattr(adata.X, 'max') and float(adata.X.max()) < 20
            has_pca = 'X_pca' in adata.obsm
            has_leiden = 'leiden' in adata.obs.columns
            has_ct = 'cell_type' in adata.obs.columns
            if not is_norm:
                return ("Data not normalized. Go to QC tab first", "error")
            if has_ct:
                n_types = adata.obs['cell_type'].nunique()
                n_clusters = adata.obs['leiden'].nunique() if has_leiden else '?'
                return (f"{n_clusters} clusters, {n_types} cell types", "success")
            if has_leiden:
                n = adata.obs['leiden'].nunique()
                return (f"{n} clusters. Annotate cell types below", "warning")
            if has_pca:
                return ("PCA done. Run clustering", "warning")
            return ("Run PCA then cluster", "warning")

        if tab_index == 5:  # Marker Check
            has_ct = 'cell_type' in adata.obs.columns
            has_leiden = 'leiden' in adata.obs.columns
            if not has_ct:
                return (("Annotate cell types first — this checks them"
                         if has_leiden else
                         "Cluster and annotate first"), "warning")
            n = adata.obs['cell_type'].nunique(dropna=True)
            return (f"{n} cell types. Check them against known markers",
                    "info")

        if tab_index == 6:  # Decontaminate
            has_ct = 'cell_type' in adata.obs.columns
            if not has_ct:
                return ("Annotate cell types first, then decontaminate", "warning")
            n_types = adata.obs['cell_type'].nunique()
            if n_types < 2:
                return ("Single cell type -- run DecontX before Subset", "warning")
            if 'decontX_counts' in adata.layers:
                return ("Decontaminated. Raw and corrected counts both kept", "success")
            return (f"{n_types} cell types. Remove ambient RNA per sample", "info")

        if tab_index == 7:  # Subset
            has_ct = 'cell_type' in adata.obs.columns
            has_leiden = 'leiden' in adata.obs.columns
            if not (has_ct or has_leiden):
                return ("No clusters or cell types. Cluster first", "error")
            if has_ct:
                n = adata.obs['cell_type'].nunique()
                return (f"{n} cell types. Select which to keep", "success")
            n = adata.obs['leiden'].nunique()
            return (f"{n} clusters. Select which to keep", "warning")

        return ("", "info")

    def _suggest_next_step(self):
        """Show the next recommended step via status_message."""
        suggestions = [
            "Next: Go to Load Data tab (Ctrl+1) to load or download your data",
            "Next: Go to Gene Names tab (Ctrl+2) to harmonise gene symbols",
            "Next: Go to Inspect tab (Ctrl+3) to explore metadata",
            "Next: Go to QC tab (Ctrl+4) to filter and normalize",
            "Next: Go to Cluster tab (Ctrl+5) for PCA, clustering, UMAP and annotation",
            "Next: Go to Marker Check tab (Ctrl+6) to sense-check the annotation",
            "Next: Go to Decontaminate tab (Ctrl+7) to remove ambient RNA",
            "Next: Go to Subset tab (Ctrl+8) to subset cell types",
        ]
        for i, msg in enumerate(suggestions):
            if not self._step_completed.get(i, False):
                self.status_message.emit(msg)
                return
        self.status_message.emit("Pipeline complete! Switch to DE Analysis mode for differential expression.")

    # ------------------------------------------------------------------
    # Tutorial
    # ------------------------------------------------------------------

    def start_tutorial(self):
        """Launch the guided tutorial."""
        self.tutorial_controller.start_tutorial()

    def resizeEvent(self, event):
        """Propagate resize to tutorial overlay."""
        super().resizeEvent(event)
        self.tutorial_controller.resize_overlay()

