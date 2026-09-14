# Decontaminate Tab (DecontX).
#
# Supervised, per-sample ambient-RNA removal on the FULL annotated
# dataset. This step sits after Cluster/Annotate and before Subset: the
# abundant source cells must be present for their ambient signal to be
# recognised and subtracted from the cells they contaminate. Three gates
# enforce that placement (annotation done, >=2 cell types, sample column
# set); when any fails the Run button is disabled with a one-line reason.
#
# Decontaminated counts are written to layers['decontX_counts']; raw
# counts are left untouched, so the step is reversible and DE can run on
# either. The layer carries through the Subset step automatically.

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import QGridLayout, QGroupBox, QLabel, QProgressBar

from kosmic.gui.shared import dialogs, run_worker
from kosmic.gui.shared.theme import NoScrollComboBox
from kosmic.gui.shared.widgets import (
    BaseWorker, PrimaryButton, SecondaryLabel, SimplePage, StatusLabel,
)

from kosmic.scrna.qc.decontx import (
    CELL_TYPE_COL,
    DECONTX_LAYER,
    check_decontx_prerequisites,
    run_decontx,
)

# Workflow step index for the Decontaminate tab (after Marker Check=5).
_STEP_INDEX = 6

# obs columns worth offering as the per-sample key: categorical-ish with
# a sensible number of levels (a sample column, not a continuous metric).
_MIN_SAMPLES = 1
_MAX_SAMPLES = 200


class DecontXWorker(BaseWorker):
    """Run supervised per-sample DecontX off the UI thread.

    Emits 'finished_ok' with '(adata, info)'.
    """

    def __init__(self, adata, sample_col: str, output_path=None):
        super().__init__()
        self.adata = adata
        self.sample_col = sample_col
        self.output_path = output_path

    def _run(self):
        self.adata, info = run_decontx(
            self.adata, self.sample_col,
            progress_callback=lambda m: self.progress.emit(m),
            progress_pct_callback=lambda p: self.progress_pct.emit(p))

        if self.output_path:
            self.progress.emit("Saving decontaminated dataset...")
            self.adata.write_h5ad(self.output_path)
        return (self.adata, info)


class DecontXTab(SimplePage):
    """Supervised ambient-RNA decontamination on the full annotated dataset."""

    help_id = "scrna/decontx"
    CONTENT_MARGINS = (16, 16, 16, 16)
    CONTENT_SPACING = 12

    log_message = pyqtSignal(str)

    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self._worker = None
        self._loading = False
        self._last_adata_version = -1
        self._setup_ui()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------
    def _setup_ui(self):
        layout = self.body_layout

        header = QLabel("Decontaminate (DecontX)")
        header.setFont(QFont("Segoe UI", 14, QFont.Weight.Bold))
        layout.addWidget(header)

        desc = SecondaryLabel(
            "Remove ambient RNA contamination using cell-type labels. "
            "DecontX subtracts transcripts that belong to other cell "
            "types present in the sample (for example cardiomyocyte "
            "counts leaking into endothelial cells). Run it here, on the "
            "full annotated dataset, before subsetting to a single cell "
            "type -- the abundant source cells must be present for their "
            "ambient signal to be identified."
        )
        desc.setWordWrap(True)
        layout.addWidget(desc)

        # --- Controls ---
        controls = QGroupBox("Settings")
        grid = QGridLayout(controls)
        grid.setSpacing(8)

        grid.addWidget(SecondaryLabel("Sample column:"), 0, 0)
        self._sample_combo = NoScrollComboBox()
        self._sample_combo.setToolTip(
            "Ambient RNA is generated during each sample's dissociation, "
            "so DecontX decontaminates one sample at a time.")
        self._sample_combo.currentIndexChanged.connect(self._on_sample_changed)
        grid.addWidget(self._sample_combo, 0, 1)
        grid.setColumnStretch(1, 1)

        layout.addWidget(controls)

        self._gate_label = StatusLabel("")
        self._gate_label.setWordWrap(True)
        layout.addWidget(self._gate_label)

        self._run_btn = PrimaryButton("Run DecontX")
        self._run_btn.clicked.connect(self._run_decontx)
        self._run_btn.setEnabled(False)
        layout.addWidget(self._run_btn)

        self._progress_bar = QProgressBar()
        self._progress_bar.setRange(0, 100)
        self._progress_bar.setValue(0)
        self._progress_bar.setVisible(False)
        layout.addWidget(self._progress_bar)

        self._result_label = SecondaryLabel("")
        self._result_label.setWordWrap(True)
        layout.addWidget(self._result_label)

        layout.addStretch(1)

        self._status_label = SecondaryLabel("")
        layout.addWidget(self._status_label)

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------
    def on_tab_activated(self):
        adata = self.main_window.current_adata
        if adata is None:
            self._gate_label.setText("No data loaded. Load, cluster and annotate first.")
            self._gate_label.set_state('info')
            self._run_btn.setEnabled(False)
            return
        version = getattr(self.main_window, '_adata_version', 0)
        if version != self._last_adata_version:
            self._last_adata_version = version
            self._populate_sample_columns(adata)
        self._refresh_gate()

    def set_adata(self, adata):
        self._last_adata_version = -1  # force refresh on next activation

    def reset_state(self):
        if self._worker is not None and self._worker.isRunning():
            self._worker.quit()
            self._worker.wait(2000)
        self._worker = None
        self._last_adata_version = -1
        self._gate_label.setText("")
        self._result_label.setText("")
        self._status_label.setText("")
        self._run_btn.setEnabled(False)

    def refresh_theme(self):
        return

    # ------------------------------------------------------------------
    # Population / gating
    # ------------------------------------------------------------------
    def _populate_sample_columns(self, adata):
        self._loading = True
        self._sample_combo.clear()
        options = []
        for col in adata.obs.columns:
            series = adata.obs[col]
            is_cat = (str(series.dtype) in ('category', 'object')
                      or hasattr(series, 'cat'))
            if not is_cat:
                continue
            n = series.nunique(dropna=True)
            if _MIN_SAMPLES <= n <= _MAX_SAMPLES:
                options.append(col)
        for col in options:
            self._sample_combo.addItem(col)
        # Default to the canonical 'sample' column assigned in Inspect.
        for pref in ('sample', 'donor', 'sample_id', 'batch'):
            idx = self._sample_combo.findText(pref)
            if idx >= 0:
                self._sample_combo.setCurrentIndex(idx)
                break
        self._loading = False

    def _selected_sample_col(self):
        col = self._sample_combo.currentText()
        return col or None

    def _on_sample_changed(self):
        if self._loading:
            return
        self._refresh_gate()

    def _refresh_gate(self):
        adata = self.main_window.current_adata
        if adata is None:
            self._run_btn.setEnabled(False)
            return

        sample_col = self._selected_sample_col()
        ok, reason = check_decontx_prerequisites(adata, sample_col)

        already = DECONTX_LAYER in adata.layers
        if not ok:
            self._gate_label.setText(reason)
            self._gate_label.set_state('warning')
            self._run_btn.setEnabled(False)
            return

        n_types = adata.obs[CELL_TYPE_COL].dropna().astype(str).nunique()
        n_samples = adata.obs[sample_col].nunique()
        if already:
            self._gate_label.setText(
                "Already decontaminated. Re-running replaces the previous "
                "corrected counts; raw counts stay untouched.")
            self._gate_label.set_state('success')
            self._run_btn.setText("Re-run DecontX")
        else:
            self._gate_label.setText(
                f"Ready: {n_types} cell types across {n_samples} sample(s). "
                f"Corrected counts are kept alongside the raw counts, so "
                f"this step is reversible.")
            self._gate_label.set_state('info')
            self._run_btn.setText("Run DecontX")
        self._run_btn.setEnabled(True)

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------
    def _run_decontx(self):
        adata = self.main_window.current_adata
        if adata is None:
            return
        sample_col = self._selected_sample_col()
        ok, reason = check_decontx_prerequisites(adata, sample_col)
        if not ok:
            dialogs.warning(self, "Cannot run DecontX", reason)
            return

        n_samples = adata.obs[sample_col].nunique()
        if not dialogs.confirm(
                self, "Run DecontX",
                "This estimates and removes ambient RNA per sample using "
                f"the cell-type labels, across {n_samples} sample(s). "
                "Raw counts are kept alongside the corrected counts, so "
                "the step is reversible.\n\nProceed?"):
            return

        self._run_btn.setEnabled(False)
        self._progress_bar.setValue(0)
        self._progress_bar.setVisible(True)
        self._set_status("Running DecontX...")
        self.log_message.emit(f"DecontX started on {n_samples} sample(s)")

        output_path = getattr(self.main_window, 'current_h5ad_path', None)
        self._worker = DecontXWorker(adata, sample_col, output_path=output_path)
        run_worker(
            self._worker,
            on_finished=self._on_finished,
            on_failed=self._on_failed,
            on_progress=self._on_progress,
            on_progress_pct=self._on_progress_pct,
        )

    def _set_status(self, msg: str):
        """Show a status message both in the tab and the app status bar."""
        self._status_label.setText(msg)
        status_bar = getattr(self.main_window, 'status_bar', None)
        if status_bar is not None:
            status_bar.showMessage(msg)

    def _on_progress(self, msg: str):
        self._set_status(msg)

    def _on_progress_pct(self, pct: int):
        self._progress_bar.setValue(int(pct))

    def _on_finished(self, payload):
        adata, info = payload
        output_path = getattr(self.main_window, 'current_h5ad_path', None)
        self.main_window.set_adata(adata, output_path, switch=False)
        if hasattr(self.main_window, 'record_provenance'):
            self.main_window.record_provenance('decontx', {
                'sample_col': self._selected_sample_col(),
                'mean_contamination': round(
                    float(info.get('mean_contamination', float('nan'))), 4),
                'n_samples': info.get('n_samples'),
                'n_samples_skipped': info.get('n_samples_skipped', 0),
            })
        self._progress_bar.setValue(100)
        self._progress_bar.setVisible(False)

        mean_contam = info.get('mean_contamination', float('nan'))
        n_skipped = info.get('n_samples_skipped', 0)
        skip_note = (f" {n_skipped} single-cell-type sample(s) kept raw."
                     if n_skipped else "")
        # Interpretation guide: typical droplet datasets sit around 5-20%.
        if mean_contam == mean_contam:  # not NaN
            if mean_contam <= 0.20:
                verdict = "typical for droplet data (5-20% is common)"
            elif mean_contam <= 0.40:
                verdict = ("on the high side; check the QC step and "
                           "sample quality")
            else:
                verdict = ("very high; check that cell-type labels and "
                           "the sample column are correct")
        else:
            verdict = ""
        self._result_label.setText(
            f"Mean estimated contamination {mean_contam:.1%} across "
            f"{info.get('n_samples', '?')} sample(s)"
            + (f" -- {verdict}." if verdict else ".")
            + f"{skip_note} Corrected counts are stored alongside the raw "
            f"counts; the DE step can use either.")
        self._set_status("Decontamination complete")
        self.main_window.mark_step_complete(_STEP_INDEX)
        self._refresh_gate()
        self.log_message.emit(info.get('message', "DecontX complete"))
        dialogs.info(self, "DecontX Complete", info.get('message', "Done."))

    def _on_failed(self, message: str):
        self._run_btn.setEnabled(True)
        self._progress_bar.setVisible(False)
        self._set_status("DecontX failed")
        self._result_label.setText("")
        self.log_message.emit(f"DecontX failed: {message}")
        dialogs.warning(self, "DecontX Failed", message)
