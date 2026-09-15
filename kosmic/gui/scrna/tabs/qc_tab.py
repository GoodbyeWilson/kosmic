# QC & Preprocessing Tab
# Three tabs: QC Filtering, Doublet Detection, Normalization.
# Each tab has controls on the left and an inline plot on the right.

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QButtonGroup,
    QCheckBox,
    QSizePolicy,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QPixmap
from pathlib import Path
import numpy as np
import pyqtgraph as pg

from kosmic.gui.shared.theme import get_color, NoScrollComboBox, NoScrollSpinBox, NoScrollDoubleSpinBox
from kosmic.gui.shared.widgets import (
    BaseWorker, SecondaryLabel, SecondaryButton, SectionHeader, SidebarPage,
    StageAccordion, StageSummaryCard, StatusLabel, show_methods_report,
)
from kosmic.gui.shared import borderless, dialogs, run_worker
from kosmic.reference.contamination import builtin_panels


# ---------------------------------------------------------------------------
# Interactive QC Histogram Widget
# ---------------------------------------------------------------------------

# Metric definitions: key -> (obs column, display name, has_min_line, has_max_line)
_METRICS = [
    ("n_genes_by_counts", "Genes / cell", True, True),
    ("total_counts", "Total counts", True, True),
    ("pct_counts_mt", "MT %", False, True),
]


class _QCHistogramWidget(QWidget):
    """
    Interactive pyqtgraph histogram with draggable threshold lines.

    Shows one metric at a time, selectable via toggle buttons.
    Emits threshold_changed(metric_key, "min"/"max", value) when a line is dragged.
    """

    threshold_changed = pyqtSignal(str, str, float)  # metric, "min"/"max", value

    def __init__(self, parent=None):
        super().__init__(parent)
        self._obs_df = None
        self._current_metric = 0  # index into _METRICS
        self._lines: dict[str, pg.InfiniteLine] = {}  # "min"/"max" -> line
        self._thresholds: dict[str, dict[str, float]] = {}  # metric -> {"min": v, "max": v}
        self._syncing = False  # guard against signal loops
        self._movable = True  # lines draggable in Fixed mode, static in MAD mode


        layout = borderless(QVBoxLayout, self)

        # --- Plot container (plot + overlaid buttons) ---
        plot_container = QWidget()
        plot_container_layout = borderless(QVBoxLayout, plot_container)

        from kosmic.gui.shared.plots import InteractivePlot
        self._plot = InteractivePlot(
            title='QC Metrics',
            left_label='Cells',
            unavailable_message='Run QC to see per-metric histograms.',
        )
        # X-axis interactive: scroll-wheel zoom + drag pan. Y stays
        # fixed (auto from histogram). Right-click exposes pyqtgraph's
        # axis-range / "View All" menu for explicit limit edits.
        self._plot.setMouseEnabled(x=True, y=False)
        self._plot.setMenuEnabled(True)
        self._plot.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        plot_container_layout.addWidget(self._plot, 1)

        # Metric selector buttons — overlaid top-right on the plot
        self._btn_bar = QWidget(self._plot)
        btn_row = QHBoxLayout(self._btn_bar)
        btn_row.setContentsMargins(0, 0, 0, 0)
        btn_row.setSpacing(2)
        self._btn_group = QButtonGroup(self)
        self._btn_group.setExclusive(True)
        self._metric_btns: list[QPushButton] = []
        for i, (_, label, _, _) in enumerate(_METRICS):
            btn = QPushButton(label)
            btn.setCheckable(True)
            btn.setChecked(i == 0)
            btn.setFixedHeight(22)
            btn.setProperty("role", "metric_toggle")
            self._btn_group.addButton(btn, i)
            btn_row.addWidget(btn)
            self._metric_btns.append(btn)
        # Log scale toggle
        self._log_btn = QPushButton("Log")
        self._log_btn.setCheckable(True)
        self._log_btn.setChecked(False)
        self._log_btn.setFixedHeight(22)
        self._log_btn.setProperty("role", "metric_toggle")
        self._log_btn.toggled.connect(self._on_log_toggled)
        btn_row.addWidget(self._log_btn)

        self._btn_bar.adjustSize()
        self._btn_group.idClicked.connect(self._on_metric_switched)

        layout.addWidget(plot_container, 1)

        # --- Info label ---
        self._info_label = SecondaryLabel("Load data to see QC histograms")
        self._info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._info_label.setContentsMargins(4, 2, 4, 2)
        layout.addWidget(self._info_label)

        # Empty-state overlay is provided by InteractivePlot; no
        # hand-rolled placeholder needed here.
        self._plot.getPlotItem().hideButtons()
        self._btn_bar.hide()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._reposition_buttons()

    def _reposition_buttons(self):
        """Keep button bar pinned to top-right of the plot widget."""
        if hasattr(self, '_btn_bar'):
            self._btn_bar.adjustSize()
            plot_rect = self._plot.rect()
            self._btn_bar.move(
                plot_rect.width() - self._btn_bar.width() - 8,
                4,
            )

    def showEvent(self, event):
        super().showEvent(event)
        # Position buttons on first show
        self._reposition_buttons()

    # -- public API --

    def set_data(self, obs_df, n_excluded: int = 0):
        """Set observation DataFrame (adata.obs) and redraw current metric.

        ``n_excluded`` is how many cells of the dataset are not in
        ``obs_df`` because their role is 'exclude'; the info line says so,
        because the pass count would otherwise look like the whole
        dataset had shrunk.
        """
        self._obs_df = obs_df
        self._n_excluded = int(n_excluded)
        # InteractivePlot's empty-state overlay is hidden in '_draw()'
        # once a metric is rendered.
        self._btn_bar.show()
        self._draw()

    def set_threshold(self, metric: str, bound: str, value: float):
        """Programmatically move a threshold line (called from spinbox changes).

        A value of 0 means the bound is off -- the filter skips it -- so
        it is removed here too, and its line goes with it. Otherwise the
        pass/removed preview would count every cell as above a cap of 0.
        """
        if self._syncing:
            return
        key, _, _, _ = _METRICS[self._current_metric]
        if value <= 0:
            self._thresholds.get(metric, {}).pop(bound, None)
            if metric == key:
                self._draw()
            else:
                self._update_info()
            return
        had_line = bound in self._lines
        self._thresholds.setdefault(metric, {})[bound] = value
        # If this metric is currently displayed, move the line
        if metric == key and had_line:
            self._syncing = True
            self._lines[bound].setValue(value)
            self._syncing = False
            self._update_info()
        elif metric == key:
            self._draw()

    def set_lines_movable(self, movable: bool):
        """Set whether threshold lines are draggable (Fixed mode) or static (MAD mode)."""
        self._movable = movable
        for line in self._lines.values():
            line.setMovable(movable)

    # -- internals --

    def _on_metric_switched(self, idx: int):
        self._current_metric = idx
        self._draw()

    def _on_log_toggled(self, checked: bool):
        self._draw()

    def _draw(self):
        """Redraw histogram for the current metric."""
        self._plot.clear_plot_items()
        self._lines.clear()

        if self._obs_df is None:
            self._info_label.setText("Load data to see QC histograms")
            self._plot.show_unavailable_message()
            return

        self._plot.hide_unavailable_message()

        key, label, has_min, has_max = _METRICS[self._current_metric]

        if key not in self._obs_df.columns:
            self._info_label.setText(f"Column '{key}' not yet computed")
            return

        values = self._obs_df[key].values.astype(float)
        use_log = self._log_btn.isChecked()
        n_total = len(values)

        # Transform values if log mode
        plot_values = np.log1p(values) if use_log else values

        # Compute histogram bins
        n_bins = min(100, max(30, n_total // 50))
        hist_y, bin_edges = np.histogram(plot_values, bins=n_bins)

        # Draw bars using BarGraphItem
        bin_widths = np.diff(bin_edges)
        bin_centers = bin_edges[:-1] + bin_widths / 2
        bar_item = pg.BarGraphItem(
            x=bin_centers, height=hist_y, width=bin_widths * 0.9,
            brush=pg.mkBrush(get_color('accent_primary') + 'cc'),
            pen=pg.mkPen(get_color('accent_primary'), width=0.5),
        )
        self._plot.addItem(bar_item)
        x_label = f"log1p({label})" if use_log else label
        self._plot.setLabel('bottom', x_label)

        # Add threshold lines (only when a threshold value has been set)
        # Transform threshold positions to match the axis scale
        thresholds = self._thresholds.get(key, {})
        if has_min and 'min' in thresholds:
            pos = np.log1p(thresholds['min']) if use_log else thresholds['min']
            line = pg.InfiniteLine(
                pos=pos, angle=90, movable=(self._movable and not use_log),
                pen=pg.mkPen('#e74c3c', style=Qt.PenStyle.DashLine, width=2),
                hoverPen=pg.mkPen('#ff6b6b', width=3),
                label='min', labelOpts={'position': 0.95, 'color': '#e74c3c'},
            )
            if not use_log:
                line.sigPositionChanged.connect(lambda ln: self._on_line_dragged('min', ln))
            self._plot.addItem(line)
            self._lines['min'] = line

        if has_max and 'max' in thresholds:
            pos = np.log1p(thresholds['max']) if use_log else thresholds['max']
            line = pg.InfiniteLine(
                pos=pos, angle=90, movable=(self._movable and not use_log),
                pen=pg.mkPen('#e74c3c', style=Qt.PenStyle.DashLine, width=2),
                hoverPen=pg.mkPen('#ff6b6b', width=3),
                label='max', labelOpts={'position': 0.95, 'color': '#e74c3c'},
            )
            if not use_log:
                line.sigPositionChanged.connect(lambda ln: self._on_line_dragged('max', ln))
            self._plot.addItem(line)
            self._lines['max'] = line

        # Fix axes: Y starts at 0, X starts at 0.
        self._plot.setYRange(0, hist_y.max() * 1.1, padding=0)
        x_min = max(0, bin_edges[0] - bin_widths[0])
        x_max = bin_edges[-1] + bin_widths[-1]

        # Widen to include the threshold lines. Ranging on the histogram
        # alone puts a marker off-screen whenever it sits outside the data
        # -- which is the normal case for an untouched upper bound, and
        # exactly when you need to see where it is.
        marker_positions = [
            np.log1p(v) if use_log else v
            for bound, v in thresholds.items()
            if (bound == 'min' and has_min) or (bound == 'max' and has_max)
        ]
        if marker_positions:
            span = max(x_max - x_min, 1.0)
            x_min = min(x_min, min(marker_positions) - span * 0.02)
            x_max = max(x_max, max(marker_positions) + span * 0.02)
            x_min = max(0, x_min)

        self._plot.setXRange(x_min, x_max, padding=0)

        self._update_info()

    def _on_line_dragged(self, bound: str, line: pg.InfiniteLine):
        """Handle a threshold line being dragged by the user."""
        if self._syncing:
            return
        key = _METRICS[self._current_metric][0]
        value = float(line.value())
        self._thresholds.setdefault(key, {})[bound] = value
        self._syncing = True
        self.threshold_changed.emit(key, bound, value)
        self._syncing = False
        self._update_info()

    def _update_info(self):
        """Update the info label with cell pass/fail counts based on all thresholds."""
        if self._obs_df is None:
            return

        n_total = len(self._obs_df)
        combined_mask = np.ones(n_total, dtype=bool)
        per_metric = []

        # Compute per-metric and combined removal
        for key, label, has_min, has_max in _METRICS:
            if key not in self._obs_df.columns:
                continue
            vals = self._obs_df[key].values.astype(float)
            th = self._thresholds.get(key, {})
            metric_mask = np.ones(n_total, dtype=bool)
            if has_min and 'min' in th:
                metric_mask &= vals >= th['min']
            if has_max and 'max' in th:
                metric_mask &= vals <= th['max']
            combined_mask &= metric_mask
            n_fail = int((~metric_mask).sum())
            if n_fail > 0:
                per_metric.append(f"{label}: -{n_fail:,}")

        n_pass = int(combined_mask.sum())
        n_removed = n_total - n_pass
        pct = (n_removed / n_total * 100) if n_total > 0 else 0

        parts = [f"{n_pass:,} / {n_total:,} pass ({n_removed:,} removed, {pct:.1f}%)"]
        if per_metric:
            parts.append("  |  ".join(per_metric))
        n_excl = getattr(self, '_n_excluded', 0)
        if n_excl:
            parts.append(f"{n_excl:,} excluded cells not shown; the same "
                         "thresholds apply to them")
        self._info_label.setText("    ".join(parts))


# ---------------------------------------------------------------------------
# Worker threads
# ---------------------------------------------------------------------------

class QCFilterWorker(BaseWorker):
    """
    Worker for running 'run_qc_pipeline' + writing the result.

    Emits 'finished_ok' with '(adata, stats)'.
    """

    def __init__(self, adata, params: dict, output_path: str):
        super().__init__()
        self.adata = adata
        self.params = params
        self.output_path = output_path

    def _run(self):
        from kosmic.scrna.qc.filter import run_qc_pipeline

        self.progress.emit("Applying QC filters...")
        adata, stats = run_qc_pipeline(self.adata, self.params)

        self.progress.emit(f"Saving to {Path(self.output_path).name}...")
        adata.write_h5ad(self.output_path)
        return (adata, stats)


class NormalizeWorker(BaseWorker):
    """
    Worker for normalizing data.

    Emits 'finished_ok' with the normalized adata.
    """

    def __init__(self, adata, target_sum: int, log_transform: bool, output_path: str):
        super().__init__()
        self.adata = adata
        self.target_sum = target_sum
        self.log_transform = log_transform
        self.output_path = output_path

    def _run(self):
        from kosmic.scrna.qc.normalize import normalize_adata

        self.progress.emit(
            f"Normalizing to {self.target_sum:,} counts per cell...")
        self.adata = normalize_adata(
            self.adata,
            target_sum=self.target_sum,
            log_transform=self.log_transform,
        )

        self.progress.emit("Saving...")
        self.adata.write_h5ad(self.output_path)
        return self.adata


class ScrubletWorker(BaseWorker):
    """
    Worker for running Scrublet doublet detection.

    Emits 'finished_ok' with a tuple '(adata, message, results_dict)'.
    """

    def __init__(self, adata, expected_doublet_rate: float, min_counts: int, output_path: str):
        super().__init__()
        self.adata = adata
        self.expected_doublet_rate = expected_doublet_rate
        self.min_counts = min_counts
        self.output_path = output_path

    def _run(self):
        try:
            from kosmic.scrna.qc.scrublet import run_scrublet
        except ImportError as e:
            raise RuntimeError(
                "Scrublet not installed. Install with: pip install scrublet") from e

        self.progress.emit("Running Scrublet doublet detection...")
        self.progress_pct.emit(10)

        self.adata, results = run_scrublet(
            self.adata,
            expected_doublet_rate=self.expected_doublet_rate,
            min_counts=self.min_counts,
        )

        self.progress.emit("Saving results...")
        self.progress_pct.emit(90)
        self.adata.write_h5ad(self.output_path)

        self.progress_pct.emit(100)
        n_doublets = results['n_doublets']
        pct = results['pct_doublets']
        n_total = results['n_total']
        msg = f"Detected {n_doublets:,} doublets ({pct:.1f}%) out of {n_total:,} cells"
        return (self.adata, msg, results)


class FilterDoubletsWorker(BaseWorker):
    """
    Worker for slicing out predicted doublets + saving the result.

    Slicing a (25k cells x 28k genes) sparse AnnData is several seconds
    of memory churn, plus another few seconds for 'write_h5ad'. Both
    must run off the main thread or the UI freezes during click->done.

    Emits 'finished_ok' with '(adata, n_removed)'.
    """

    def __init__(self, adata, output_path: str):
        super().__init__()
        self.adata = adata
        self.output_path = output_path

    def _run(self):
        n_before = self.adata.n_obs
        self.progress.emit(f"Filtering doublets from {n_before:,} cells...")
        self.progress_pct.emit(20)

        adata = self.adata[~self.adata.obs['predicted_doublet'], :].copy()
        n_after = adata.n_obs
        n_removed = n_before - n_after

        self.progress.emit(f"Saving {n_after:,} cells to disk...")
        self.progress_pct.emit(70)
        adata.write_h5ad(self.output_path)
        self.progress_pct.emit(100)

        return (adata, n_removed)


class SoupXWorker(BaseWorker):
    """
    Worker for running SoupX ambient-RNA correction.

    Emits 'finished_ok' with '(adata, message, info_dict)'. Info has
    contamination_fraction, used_filtered_as_raw, auto_estimated.
    """

    def __init__(self, adata, output_path: str, raw_path: str = None,
                 contamination_fraction: float = None,
                 tfidf_min: float = 1.0):
        super().__init__()
        self.adata = adata
        self.output_path = output_path
        self.raw_path = raw_path
        self.contamination_fraction = contamination_fraction
        self.tfidf_min = tfidf_min

    def _run(self):
        try:
            from kosmic.scrna.qc.soupx import run_soupx
        except ImportError as e:
            raise RuntimeError(
                "soupx-python not installed. Install with: pip install soupx-python") from e

        raw_adata = None
        if self.raw_path:
            self.progress.emit(f"Loading raw matrix from {self.raw_path}...")
            self.progress_pct.emit(5)
            try:
                from kosmic.scrna.load.converters import load_10x_folder
                raw_adata = load_10x_folder(self.raw_path)
            except Exception as exc:
                raise RuntimeError(
                    f"Failed to load raw 10x folder {self.raw_path!r}: {exc}") from exc

        self.progress.emit("Running SoupX...")
        self.progress_pct.emit(15)

        def _cb(msg):
            self.progress.emit(msg)

        self.adata, info = run_soupx(
            self.adata,
            raw_adata=raw_adata,
            contamination_fraction=self.contamination_fraction,
            tfidf_min=self.tfidf_min,
            progress_callback=_cb,
        )

        self.progress.emit("Saving...")
        self.progress_pct.emit(90)
        self.adata.write_h5ad(self.output_path)
        self.progress_pct.emit(100)
        return (self.adata, info['message'], info)


class QCMetricsWorker(BaseWorker):
    """
    Worker for computing QC metrics off the main thread.

    Always emits 'finished_ok' with the adata (possibly unmodified if the
    per-gene MT detection failed). Errors are logged via 'progress' rather
    than propagated, because the caller expects an adata in both cases.
    """

    def __init__(self, adata):
        super().__init__()
        self.adata = adata

    def _run(self):
        try:
            import scanpy as sc
            needs_qc = 'n_genes_by_counts' not in self.adata.obs.columns
            needs_mt = 'pct_counts_mt' not in self.adata.obs.columns

            if needs_qc or needs_mt:
                mt_genes = self.adata.var_names.str.lower().str.startswith('mt-')
                if mt_genes.sum() == 0:
                    mt_genes = self.adata.var_names.str.startswith('MT-')
                self.adata.var['mt'] = mt_genes

                if needs_qc:
                    sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], inplace=True)
                else:
                    # n_genes_by_counts already present; just compute pct_counts_mt
                    import numpy as np
                    if mt_genes.sum() > 0:
                        X = self.adata.X
                        total = np.asarray(X.sum(axis=1)).ravel()
                        mt_sum = np.asarray(X[:, mt_genes].sum(axis=1)).ravel()
                        with np.errstate(divide='ignore', invalid='ignore'):
                            pct = np.where(total > 0, mt_sum / total * 100, 0.0)
                        self.adata.obs['pct_counts_mt'] = pct
                    elif self.adata.raw is not None:
                        # MT genes filtered out during HVG selection -- use raw
                        raw_mt = self.adata.raw.var_names.str.startswith('MT-')
                        if raw_mt.sum() == 0:
                            raw_mt = self.adata.raw.var_names.str.lower().str.startswith('mt-')
                        if raw_mt.sum() > 0:
                            import numpy as np
                            X_raw = self.adata.raw.X
                            total = np.asarray(X_raw.sum(axis=1)).ravel()
                            mt_sum = np.asarray(X_raw[:, raw_mt].sum(axis=1)).ravel()
                            with np.errstate(divide='ignore', invalid='ignore'):
                                pct = np.where(total > 0, mt_sum / total * 100, 0.0)
                            self.adata.obs['pct_counts_mt'] = pct
                        else:
                            self.adata.obs['pct_counts_mt'] = 0.0
                    else:
                        self.adata.obs['pct_counts_mt'] = 0.0
        except Exception as e:
            import traceback
            self.progress.emit(
                f"QC worker error: {e}\n{traceback.format_exc()}")
        return self.adata


class ScrubletPlotWorker(BaseWorker):
    """
    Worker for generating Scrublet histogram off the main thread.

    Always emits 'finished_ok' with bytes (empty on failure) -- the caller
    interprets empty bytes as 'no plot available'.
    """

    def __init__(self, doublet_scores, threshold, dark_mode=True):
        super().__init__()
        self.doublet_scores = doublet_scores
        self.threshold = threshold
        self.dark_mode = dark_mode

    def _run(self):
        try:
            from kosmic.visualisation.scrna.qc_plots import create_scrublet_histogram
            import io
            import matplotlib.pyplot as plt
            fig = create_scrublet_histogram(
                self.doublet_scores,
                self.threshold,
                dark_mode=self.dark_mode,
            )
            if fig is None:
                return b""
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=120, bbox_inches="tight",
                        facecolor=fig.get_facecolor(), edgecolor='none')
            plt.close(fig)
            buf.seek(0)
            return buf.getvalue()
        except Exception:
            return b""


# ---------------------------------------------------------------------------
# QCTab
# ---------------------------------------------------------------------------

class QCTab(SidebarPage):
    """
    QC & Preprocessing tab: filtering, doublet detection, normalization.

    Three tabs with controls (left) and inline plots (right).
    """

    help_id = "scrna/qc"

    log_message = pyqtSignal(str)

    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.project_dir = None
        self.adata = None
        self.h5ad_path = None
        self.normalize_worker = None
        self.scrublet_worker = None
        self._qc_filter_worker = None
        self._qc_plot_worker = None
        self._scrublet_plot_worker = None
        self._soupx_worker = None
        self._filter_doublets_worker = None
        self.progress_bar = None  # wired to sidebar by main.py
        # Per-mode threshold memory, so looking at MAD does not discard
        # values typed in Fixed mode (and vice versa).
        self._mode_thresholds: dict[str, dict] = {}
        self._qc_mode: str | None = None

        self._setup_ui()

    def _setup_ui(self):
        # Sidebar lives on self.sidebar_layout (provided by SidebarPage).
        left_layout = self.sidebar_layout

        def _page(builder):
            page = QWidget()
            lay = QVBoxLayout(page)
            lay.setContentsMargins(0, 4, 0, 0)
            lay.setSpacing(4)
            builder(lay)
            lay.addStretch()
            return page

        left_layout.addWidget(SectionHeader("ANALYSIS SUMMARY"))
        self._accordion = StageAccordion(left_layout)

        self._filter_card = StageSummaryCard("Filter", icon="filter-funnel")
        self._filter_card.set_summary(["Not run yet"])
        self._filter_card.set_settings_widget(_page(self._setup_qc_panel))
        left_layout.addWidget(self._filter_card)
        self._accordion.add_card('filter', self._filter_card)

        self._doublet_card = StageSummaryCard("Doublets", icon="cells")
        self._doublet_card.set_summary(["Optional; not analysed"])
        self._doublet_card.set_settings_widget(_page(self._setup_doublet_panel))
        left_layout.addWidget(self._doublet_card)
        self._accordion.add_card('doublets', self._doublet_card)

        self._soupx_card = StageSummaryCard("Ambient (SoupX)", icon="flask-conical")
        self._soupx_card.set_summary(["Optional; not run"])
        self._soupx_card.set_settings_widget(_page(self._setup_soupx_panel))
        left_layout.addWidget(self._soupx_card)
        self._accordion.add_card('soupx', self._soupx_card)

        self._norm_card = StageSummaryCard("Normalise", icon="bar-chart-2")
        self._norm_card.set_summary(["Not run yet"])
        self._norm_card.set_settings_widget(_page(self._setup_norm_panel))
        left_layout.addWidget(self._norm_card)
        self._accordion.add_card('normalise', self._norm_card)

        left_layout.addStretch()
        self._accordion.finalize()

        # Status label
        self.status_label = QLabel("No data loaded — go to Load Data first")
        self.status_label.setProperty("role", "padded_status")
        left_layout.addWidget(self.status_label)

        report_btn = SecondaryButton("View full report")
        report_btn.setToolTip(
            "Show the methods recorded for this study (provenance).")
        report_btn.clicked.connect(
            lambda: show_methods_report(self, self.h5ad_path))
        left_layout.addWidget(report_btn)

        # ── Right: QC histogram + scrublet plot ──
        right_layout = self.content_layout
        right_layout.addWidget(self._qc_histogram, 1)

        # Scrublet plot (shown when doublet detection runs)
        self._doublet_plot_label = QLabel()
        self._doublet_plot_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._doublet_plot_label.setMinimumHeight(200)
        self._doublet_plot_label.hide()
        right_layout.addWidget(self._doublet_plot_label)

        # Lock the splitter handle so the sidebar can't be dragged wider
        # (qc_tab keeps its fixed-width sidebar invariant).
        self.splitter.handle(1).setEnabled(False)

    # ------------------------------------------------------------------
    def _get_save_path(self) -> Path:
        """Return the current h5ad path — all saves overwrite the same file."""
        return self.h5ad_path

    # Panel builders
    # ------------------------------------------------------------------

    def _setup_qc_panel(self, parent_layout):
        """QC Filtering controls."""
        left_layout = parent_layout

        self.qc_status = StatusLabel("")
        self.qc_status.setWordWrap(True)
        left_layout.addWidget(self.qc_status)

        grid = QGridLayout()
        grid.setSpacing(4)

        # Row 0: QC Mode selector
        grid.addWidget(QLabel("Mode"), 0, 0)
        self.qc_mode_combo = NoScrollComboBox()
        self.qc_mode_combo.addItems(["Fixed Thresholds", "MAD-based"])
        self.qc_mode_combo.setToolTip(
            "Fixed: type the thresholds yourself.\n"
            "MAD-based: pre-fills the threshold boxes adaptively (a number\n"
            "of median absolute deviations from the median); you can still\n"
            "adjust any value before applying.")
        self.qc_mode_combo.currentTextChanged.connect(self._toggle_qc_mode)
        # The combo does not emit for the item it starts on, so without
        # this the first switch away from Fixed has no 'previous' mode to
        # stash and the typed thresholds are lost.
        self._qc_mode = self.qc_mode_combo.currentText()
        grid.addWidget(self.qc_mode_combo, 0, 1)

        # Row 1: MAD multiplier (only enabled in MAD mode)
        self.qc_mad_label = QLabel("MAD")
        grid.addWidget(self.qc_mad_label, 1, 0)
        self.qc_mad_spin = NoScrollDoubleSpinBox()
        self.qc_mad_spin.setRange(1.0, 10.0)
        self.qc_mad_spin.setValue(5.0)
        self.qc_mad_spin.setSingleStep(0.5)
        self.qc_mad_spin.setToolTip("Number of MADs from median (log1p space). Lower = stricter. 5 is standard, 3 is aggressive.")
        self.qc_mad_spin.setEnabled(False)
        grid.addWidget(self.qc_mad_spin, 1, 1)

        # Row 2: Min genes/cell
        grid.addWidget(QLabel("Min genes"), 2, 0)
        self.qc_min_genes_spin = NoScrollSpinBox()
        self.qc_min_genes_spin.setRange(0, 50000)
        self.qc_min_genes_spin.setValue(200)
        self.qc_min_genes_spin.setToolTip("Filter cells with fewer than this many genes")
        grid.addWidget(self.qc_min_genes_spin, 2, 1)

        # Row 3: Max genes/cell
        grid.addWidget(QLabel("Max genes"), 3, 0)
        self.qc_max_genes_spin = NoScrollSpinBox()
        self.qc_max_genes_spin.setRange(0, 50000)
        self.qc_max_genes_spin.setValue(6000)
        self.qc_max_genes_spin.setSpecialValueText("off")
        self.qc_max_genes_spin.setToolTip(
            "Filter cells with more than this many genes (0 = no filter). "
            "In tissues where one cell type carries far more transcripts "
            "than the others -- cardiomyocytes in heart -- a cap removes "
            "that type, not doublets; use Scrublet for doublets.")
        grid.addWidget(self.qc_max_genes_spin, 3, 1)

        # Row 4: Min counts/cell
        grid.addWidget(QLabel("Min counts"), 4, 0)
        self.qc_min_counts_spin = NoScrollSpinBox()
        self.qc_min_counts_spin.setRange(0, 500000)
        self.qc_min_counts_spin.setValue(0)
        self.qc_min_counts_spin.setSpecialValueText("off")
        self.qc_min_counts_spin.setToolTip("Filter cells with fewer than this many total counts (0 = no filter)")
        grid.addWidget(self.qc_min_counts_spin, 4, 1)

        # Row 5: Max counts/cell
        grid.addWidget(QLabel("Max counts"), 5, 0)
        self.qc_max_counts_spin = NoScrollSpinBox()
        self.qc_max_counts_spin.setRange(0, 500000)
        self.qc_max_counts_spin.setValue(0)
        self.qc_max_counts_spin.setSpecialValueText("off")
        self.qc_max_counts_spin.setToolTip("Filter cells with more than this many total counts (0 = no filter)")
        grid.addWidget(self.qc_max_counts_spin, 5, 1)

        # Row 6: Max MT %
        grid.addWidget(QLabel("Max MT %"), 6, 0)
        self.qc_max_mt_spin = NoScrollSpinBox()
        self.qc_max_mt_spin.setRange(1, 100)
        self.qc_max_mt_spin.setValue(20)
        self.qc_max_mt_spin.setToolTip("Filter cells with more than this percentage of mitochondrial genes")
        grid.addWidget(self.qc_max_mt_spin, 6, 1)

        # Row 7: Ambient-contamination panel (optional). Scores a
        # dominant-cell-type gene panel into obs['pct_counts_<key>']
        # (percent.mito-style) on the pre-filter counts, for the DE
        # spillover diagnostic. Does NOT filter cells.
        grid.addWidget(QLabel("Contam. panel"), 7, 0)
        self.qc_contam_combo = NoScrollComboBox()
        self.qc_contam_combo.addItem("None", None)
        self._contam_panels = builtin_panels()
        for key, spec in self._contam_panels.items():
            self.qc_contam_combo.addItem(spec.get("display_name", key), key)
        self.qc_contam_combo.setToolTip(
            "Score a dominant-cell-type ambient-contamination panel per cell "
            "into pct_counts_<panel>, stored in obs so it survives subsetting "
            "and feeds the spillover diagnostic in DE. 'None' disables it.")
        grid.addWidget(self.qc_contam_combo, 7, 1)

        # Re-preview MAD thresholds when multiplier changes
        self.qc_mad_spin.valueChanged.connect(self._on_mad_spin_changed)

        left_layout.addLayout(grid)

        self.apply_qc_btn = QPushButton("Apply QC Filters")
        self.apply_qc_btn.clicked.connect(self._apply_qc_filters)
        self.apply_qc_btn.setEnabled(False)
        left_layout.addWidget(self.apply_qc_btn)

        # Histogram is created here but added to the right panel in _setup_ui
        self._qc_histogram = _QCHistogramWidget()
        self._qc_histogram.setMinimumHeight(220)

        # Connect threshold lines -> spinboxes
        self._qc_histogram.threshold_changed.connect(self._on_histogram_threshold_changed)

        # Connect spinboxes -> threshold lines
        self.qc_min_genes_spin.valueChanged.connect(
            lambda v: self._qc_histogram.set_threshold("n_genes_by_counts", "min", float(v)))
        self.qc_max_genes_spin.valueChanged.connect(
            lambda v: self._qc_histogram.set_threshold("n_genes_by_counts", "max", float(v)))
        self.qc_min_counts_spin.valueChanged.connect(
            lambda v: self._qc_histogram.set_threshold("total_counts", "min", float(v)))
        self.qc_max_counts_spin.valueChanged.connect(
            lambda v: self._qc_histogram.set_threshold("total_counts", "max", float(v)))
        self.qc_max_mt_spin.valueChanged.connect(
            lambda v: self._qc_histogram.set_threshold("pct_counts_mt", "max", float(v)))

    def _setup_doublet_panel(self, parent_layout):
        """Doublet Detection controls."""
        left_layout = parent_layout

        self.doublet_status = StatusLabel("")
        self.doublet_status.setWordWrap(True)
        left_layout.addWidget(self.doublet_status)

        grid = QGridLayout()
        grid.setSpacing(4)
        grid.addWidget(QLabel("Rate"), 0, 0)
        self.doublet_rate_spin = NoScrollSpinBox()
        self.doublet_rate_spin.setRange(1, 20)
        self.doublet_rate_spin.setValue(6)
        self.doublet_rate_spin.setSuffix("%")
        self.doublet_rate_spin.setToolTip("Expected doublet rate (typically 5-10%)")
        grid.addWidget(self.doublet_rate_spin, 0, 1)

        grid.addWidget(QLabel("Min counts"), 1, 0)
        self.doublet_min_counts_spin = NoScrollSpinBox()
        self.doublet_min_counts_spin.setRange(1, 10)
        self.doublet_min_counts_spin.setValue(2)
        self.doublet_min_counts_spin.setToolTip("Minimum counts for a gene to be used")
        grid.addWidget(self.doublet_min_counts_spin, 1, 1)
        left_layout.addLayout(grid)

        self.run_scrublet_btn = QPushButton("Run Scrublet")
        self.run_scrublet_btn.clicked.connect(self._run_scrublet)
        self.run_scrublet_btn.setEnabled(False)
        left_layout.addWidget(self.run_scrublet_btn)

        self.filter_doublets_btn = QPushButton("Filter Doublets")
        self.filter_doublets_btn.clicked.connect(self._filter_doublets)
        self.filter_doublets_btn.setEnabled(False)
        left_layout.addWidget(self.filter_doublets_btn)

    def _setup_soupx_panel(self, parent_layout):
        """SoupX ambient-RNA correction controls."""
        left_layout = parent_layout

        self.soupx_status = StatusLabel("")
        self.soupx_status.setWordWrap(True)
        left_layout.addWidget(self.soupx_status)

        # Raw 10x folder picker
        raw_row = QHBoxLayout()
        raw_row.setSpacing(4)
        self.soupx_raw_label = SecondaryLabel("(no raw 10x; will use filtered-only fallback)")
        self.soupx_raw_label.setWordWrap(True)
        raw_row.addWidget(self.soupx_raw_label, 1)
        self.soupx_browse_btn = QPushButton("Browse...")
        self.soupx_browse_btn.setToolTip(
            "CellRanger raw_feature_bc_matrix folder (with empty droplets).\n"
            "Optional but recommended -- without it, SoupX uses the filtered\n"
            "matrix as a degraded soup approximation.")
        self.soupx_browse_btn.clicked.connect(self._browse_raw_10x)
        raw_row.addWidget(self.soupx_browse_btn)
        left_layout.addLayout(raw_row)
        self._soupx_raw_path = None

        # Mode: auto-estimate vs manual
        grid = QGridLayout()
        grid.setSpacing(4)
        grid.addWidget(QLabel("Mode"), 0, 0)
        self.soupx_mode_combo = NoScrollComboBox()
        self.soupx_mode_combo.addItems(["Auto-estimate", "Manual fraction"])
        self.soupx_mode_combo.setToolTip(
            "Auto: SoupX estimates contamination from cluster markers.\n"
            "Manual: skip auto-estimation, apply a fixed fraction "
            "(useful when auto fails due to insufficient cluster separation).")
        self.soupx_mode_combo.currentTextChanged.connect(self._on_soupx_mode_changed)
        grid.addWidget(self.soupx_mode_combo, 0, 1)

        # Manual contamination fraction (only when mode = Manual)
        grid.addWidget(QLabel("Fraction"), 1, 0)
        self.soupx_frac_spin = NoScrollDoubleSpinBox()
        self.soupx_frac_spin.setRange(0.01, 0.80)
        self.soupx_frac_spin.setValue(0.10)
        self.soupx_frac_spin.setSingleStep(0.01)
        self.soupx_frac_spin.setDecimals(2)
        self.soupx_frac_spin.setToolTip(
            "Fixed contamination fraction (0.05-0.20 typical).")
        self.soupx_frac_spin.setEnabled(False)
        grid.addWidget(self.soupx_frac_spin, 1, 1)

        # tf-idf min (auto mode tuning)
        grid.addWidget(QLabel("tf-idf min"), 2, 0)
        self.soupx_tfidf_spin = NoScrollDoubleSpinBox()
        self.soupx_tfidf_spin.setRange(0.1, 5.0)
        self.soupx_tfidf_spin.setValue(1.0)
        self.soupx_tfidf_spin.setSingleStep(0.1)
        self.soupx_tfidf_spin.setDecimals(1)
        self.soupx_tfidf_spin.setToolTip(
            "Lower = more permissive marker-gene cutoff. Try 0.5 if "
            "auto-estimation reports no plausible markers.")
        grid.addWidget(self.soupx_tfidf_spin, 2, 1)

        left_layout.addLayout(grid)

        self.run_soupx_btn = QPushButton("Run SoupX")
        self.run_soupx_btn.clicked.connect(self._run_soupx)
        self.run_soupx_btn.setEnabled(False)
        left_layout.addWidget(self.run_soupx_btn)

    def _setup_norm_panel(self, parent_layout):
        """Normalization controls."""
        left_layout = parent_layout

        self.norm_status = StatusLabel("")
        self.norm_status.setWordWrap(True)
        left_layout.addWidget(self.norm_status)

        grid = QGridLayout()
        grid.setSpacing(4)
        grid.addWidget(QLabel("Target sum"), 0, 0)
        self.target_sum_spin = NoScrollSpinBox()
        self.target_sum_spin.setRange(1000, 100000)
        self.target_sum_spin.setValue(10000)
        self.target_sum_spin.setSingleStep(1000)
        self.target_sum_spin.setToolTip(
            "Every cell's counts are rescaled to sum to this value, removing\n"
            "sequencing-depth differences between cells. 10,000 is the\n"
            "scanpy convention.")
        grid.addWidget(self.target_sum_spin, 0, 1)

        self.log_transform_check = QCheckBox("Log1p transform")
        self.log_transform_check.setChecked(True)
        self.log_transform_check.setToolTip(
            "Apply log(1 + x) after rescaling so highly expressed genes\n"
            "don't dominate. Standard for clustering; leave on.")
        grid.addWidget(self.log_transform_check, 1, 0, 1, 2)
        left_layout.addLayout(grid)

        self.normalize_btn = QPushButton("Normalize & Save")
        self.normalize_btn.setToolTip("Normalize counts to target sum per cell,\nthen apply log1p transform. Saves to h5ad file.")
        self.normalize_btn.clicked.connect(self._normalize_data)
        self.normalize_btn.setEnabled(False)
        left_layout.addWidget(self.normalize_btn)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def set_project_directory(self, directory: str):
        self.project_dir = Path(directory)

    def set_data(self, adata, file_path=None):
        """Receive adata broadcast from ScRNAWorkspace.set_adata()."""
        self.adata = adata
        self.h5ad_path = Path(file_path) if file_path else None
        self._qc_plots_stale = True
        self._update_controls()

    def reset_state(self):
        """Clear all cached state (called on raw data reload)."""
        # Kill running workers
        for attr in ('normalize_worker', 'scrublet_worker', '_qc_filter_worker',
                     '_qc_plot_worker', '_scrublet_plot_worker', '_soupx_worker',
                     '_filter_doublets_worker'):
            w = getattr(self, attr, None)
            if w is not None and w.isRunning():
                w.quit()
                w.wait(2000)
        self.normalize_worker = None
        self.scrublet_worker = None
        self._qc_filter_worker = None
        self._qc_plot_worker = None
        self._scrublet_plot_worker = None
        self._soupx_worker = None
        self._filter_doublets_worker = None
        self.adata = None
        self.h5ad_path = None
        self._last_adata_version = -1
        self._last_qc_params = None
        self._doublets_removed = 0
        self._last_norm_target = None
        # Reset UI
        self.status_label.setText("No data loaded — go to Load Data first")
        self.qc_status.setText("")
        self.qc_status.set_state('info')
        self.doublet_status.setText("Not analyzed")
        self.doublet_status.set_state('info')
        self.norm_status.setText("")
        self.norm_status.set_state('info')
        self._filter_card.set_completed(False)
        self._filter_card.set_summary(["Not run yet"])
        self._doublet_card.set_completed(False)
        self._doublet_card.set_summary(["Optional; not analysed"])
        self._soupx_card.set_completed(False)
        self._soupx_card.set_summary(["Optional; not run"])
        self._norm_card.set_completed(False)
        self._norm_card.set_summary(["Not run yet"])
        self._accordion.collapse_all()
        self.apply_qc_btn.setEnabled(False)
        self.run_scrublet_btn.setEnabled(False)
        self.filter_doublets_btn.setEnabled(False)
        self.normalize_btn.setEnabled(False)
        # Clear histogram
        self._qc_histogram._obs_df = None
        self._qc_histogram._thresholds.clear()
        self._qc_histogram._plot.clear()
        self._qc_histogram._lines.clear()
        self._qc_histogram._info_label.setText("Load data to see QC histograms")
        self._qc_histogram._btn_bar.hide()

    def refresh_theme(self):
        """Re-apply theme to pyqtgraph plot."""
        from kosmic.gui.shared.theme import style_pg_plot
        style_pg_plot(self._qc_histogram._plot, title='QC Metrics', left_label='Cells')

    def on_tab_activated(self):
        """Called when tab becomes visible — pull latest adata from workspace."""
        ws = self.main_window
        if ws.current_adata is None:
            self.status_label.setText("No data loaded — go to Load Data first")
            return
        version = getattr(ws, '_adata_version', 0)
        if version != getattr(self, '_last_adata_version', -1):
            self._last_adata_version = version
            self.set_data(ws.current_adata, ws.current_h5ad_path)
        if getattr(self, '_qc_plots_stale', False):
            self._qc_plots_stale = False
            self._generate_qc_plots()

        # Auto-detect if QC + normalization already done
        if self.adata is not None:
            qc_done = 'n_genes_by_counts' in self.adata.obs.columns
            normalized = self.adata.raw is not None or (
                hasattr(self.adata.X, 'max') and float(self.adata.X.max()) < 20
            )
            if qc_done and normalized:
                ws.mark_step_complete(3)

    # ------------------------------------------------------------------
    # QC metrics + plot generation
    # ------------------------------------------------------------------

    def _ensure_qc_metrics(self):
        """
        Compute QC metrics if not already present. Called on main thread
        only when user explicitly clicks Apply (acceptable brief block).

        MT-gene detection is delegated to
        'kosmic.scrna.qc.detect_mitochondrial_genes' (auto-detects
        species). Calculating 'pct_counts_mt' only -- when n_genes /
        total_counts already exist -- is done inline because scanpy
        always recomputes the whole metric set.
        """
        if self.adata is None:
            return
        needs_qc = 'n_genes_by_counts' not in self.adata.obs.columns
        needs_mt = 'pct_counts_mt' not in self.adata.obs.columns
        if not needs_qc and not needs_mt:
            return
        try:
            import scanpy as sc
            import numpy as np
            from kosmic.scrna.qc.filter import detect_mitochondrial_genes
            detect_mitochondrial_genes(self.adata)

            if needs_qc:
                sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], inplace=True)
            else:
                # genes/counts already present, only compute pct_counts_mt
                mt_mask = self.adata.var['mt'].values
                X = self.adata.X
                total = np.asarray(X.sum(axis=1)).ravel()
                mt_sum = (np.asarray(X[:, mt_mask].sum(axis=1)).ravel()
                          if mt_mask.sum() > 0 else np.zeros(len(total)))
                with np.errstate(divide='ignore', invalid='ignore'):
                    self.adata.obs['pct_counts_mt'] = np.where(
                        total > 0, mt_sum / total * 100, 0.0)
        except Exception as e:
            self.log_message.emit(f"QC metric computation failed: {e}")

    def _generate_qc_plots(self, thresholds=None):
        """Launch QCMetricsWorker to compute metrics off main thread."""
        if self.adata is None:
            return

        # Skip the worker entirely if ALL metrics already computed
        if ('n_genes_by_counts' in self.adata.obs.columns
                and 'pct_counts_mt' in self.adata.obs.columns):
            self._on_qc_metrics_finished(self.adata, thresholds)
            return

        self.status_label.setText("Computing QC metrics...")
        if self.progress_bar:
            self.progress_bar.setRange(0, 0)  # indeterminate

        self._qc_plot_worker = QCMetricsWorker(self.adata)
        run_worker(
            self._qc_plot_worker,
            on_finished=lambda adata: self._on_qc_metrics_finished(adata, thresholds),
            on_progress=self.log_message,
        )

    def _on_qc_metrics_finished(self, adata, thresholds=None):
        """Handle QC metrics computation -- populate histogram widget."""
        self.adata = adata
        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)

        # Guaranteed fallback: compute pct_counts_mt if still missing
        if 'pct_counts_mt' not in adata.obs.columns:
            try:
                import numpy as np
                mt_mask = adata.var_names.str.startswith('MT-')
                if mt_mask.sum() == 0:
                    mt_mask = adata.var_names.str.lower().str.startswith('mt-')

                if mt_mask.sum() > 0:
                    X = adata.X
                    total = np.asarray(X.sum(axis=1)).ravel().astype(float)
                    mt_sum = np.asarray(X[:, mt_mask].sum(axis=1)).ravel().astype(float)
                    pct = np.where(total > 0, mt_sum / total * 100, 0.0)
                    source = f"{mt_mask.sum()} MT genes in var"
                elif adata.raw is not None:
                    # MT genes removed by HVG selection -- use raw
                    raw_mt = adata.raw.var_names.str.startswith('MT-')
                    if raw_mt.sum() == 0:
                        raw_mt = adata.raw.var_names.str.lower().str.startswith('mt-')
                    if raw_mt.sum() > 0:
                        X_raw = adata.raw.X
                        total = np.asarray(X_raw.sum(axis=1)).ravel().astype(float)
                        mt_sum = np.asarray(X_raw[:, raw_mt].sum(axis=1)).ravel().astype(float)
                        pct = np.where(total > 0, mt_sum / total * 100, 0.0)
                        source = f"{raw_mt.sum()} MT genes in raw"
                    else:
                        pct = np.zeros(adata.n_obs)
                        source = "no MT genes found"
                else:
                    pct = np.zeros(adata.n_obs)
                    source = "no raw, no MT genes in var"

                adata.obs['pct_counts_mt'] = pct
                self.log_message.emit(
                    f"MT%: {source} -- max={pct.max():.1f}%, mean={pct.mean():.2f}%"
                )
            except Exception as e:
                self.log_message.emit(f"MT% fallback failed: {e}")

        self._update_controls()

        # Seed threshold positions from spinbox values if no explicit thresholds
        if thresholds is None:
            thresholds = {
                "n_genes_by_counts": {
                    "min": float(self.qc_min_genes_spin.value()),
                    "max": float(self.qc_max_genes_spin.value()),
                },
                "pct_counts_mt": {"max": float(self.qc_max_mt_spin.value())},
            }

        # Push thresholds into the histogram widget
        for metric, bounds in thresholds.items():
            for bound, val in bounds.items():
                self._qc_histogram.set_threshold(metric, bound, float(val))

        # Populate histogram with the cells the thresholds come from.
        self._refresh_histogram_source()
        # status_label is updated by _update_controls() above

    def _generate_scrublet_plot(self, scores, threshold):
        """Launch ScrubletPlotWorker to create histogram."""
        self._doublet_plot_label.setText("Generating histogram...")

        self._scrublet_plot_worker = ScrubletPlotWorker(
            scores, threshold, dark_mode=True,
        )
        run_worker(
            self._scrublet_plot_worker,
            on_finished=self._on_scrublet_plot_finished,
        )

    def _on_scrublet_plot_finished(self, data: bytes):
        """Display scrublet histogram in the doublet panel."""
        if not data:
            self._doublet_plot_label.setText("Could not generate histogram")
            return
        pixmap = QPixmap()
        pixmap.loadFromData(data)
        if not pixmap.isNull():
            scaled = pixmap.scaled(
                self._doublet_plot_label.width() - 4,
                self._doublet_plot_label.height() - 4,
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation,
            )
            self._doublet_plot_label.setPixmap(scaled)

    # ------------------------------------------------------------------
    # Status / progress handlers
    # ------------------------------------------------------------------

    def _on_status(self, message: str):
        self.status_label.setText(message)
        self.log_message.emit(message)

    def _on_progress(self, value: int):
        if self.progress_bar:
            self.progress_bar.setValue(value)

    # ------------------------------------------------------------------
    # Update controls from data state
    # ------------------------------------------------------------------

    def _update_controls(self):
        """Enable/disable controls based on loaded data."""
        if self.adata is None:
            return

        n = self.adata.n_obs
        file_info = f"  [{self.h5ad_path.name}]" if self.h5ad_path else ""
        self.status_label.setText(f"{n:,} cells x {self.adata.n_vars:,} genes{file_info}")

        # QC — detect if already filtered
        qc_done = 'n_genes_by_counts' in self.adata.obs.columns
        if qc_done:
            self.qc_status.setText(f"QC metrics present — {n:,} cells")
            self.qc_status.set_state('success')
            self._filter_card.set_completed(True)
            params = getattr(self, '_last_qc_params', None)
            if params:
                self._filter_card.set_summary([
                    f"{params['min_genes']}-{params['max_genes']:,} genes · "
                    f"MT ≤ {params['max_mt']}%",
                    f"{n:,} cells kept",
                ])
            else:
                self._filter_card.set_summary([f"QC metrics present · {n:,} cells"])
        else:
            self.qc_status.setText("")
            self._filter_card.set_completed(False)
            self._filter_card.set_summary(["Not run yet"])
        self.apply_qc_btn.setEnabled(True)

        # Doublet detection
        self.run_scrublet_btn.setEnabled(True)
        if 'predicted_doublet' in self.adata.obs.columns:
            n_doublets = self.adata.obs['predicted_doublet'].sum()
            pct = (n_doublets / n) * 100
            self.doublet_status.setText(f"Detected: {n_doublets:,} doublets ({pct:.1f}%)")
            self.doublet_status.set_state('warning')
            self.filter_doublets_btn.setEnabled(True)
            self._doublet_card.set_completed(True)
            self._doublet_card.set_summary(
                [f"{n_doublets:,} doublets flagged ({pct:.1f}%)"])
        elif getattr(self, '_doublets_removed', 0):
            self.doublet_status.setText("Doublets filtered")
            self.doublet_status.set_state('success')
            self.filter_doublets_btn.setEnabled(False)
            self._doublet_card.set_completed(True)
            self._doublet_card.set_summary(
                [f"{self._doublets_removed:,} doublets removed"])
        else:
            self.doublet_status.setText("Not analyzed")
            self.doublet_status.set_state('info')
            self.filter_doublets_btn.setEnabled(False)
            self._doublet_card.set_completed(False)
            self._doublet_card.set_summary(["Optional; not analysed"])

        # SoupX -- enabled if data has integer counts (raw or in counts layer)
        has_counts_layer = 'counts' in self.adata.layers
        x_max = float(self.adata.X.max()) if hasattr(self.adata.X, 'max') else 0
        soupx_ready = has_counts_layer or x_max >= 20
        already_corrected = 'counts_pre_soupx' in self.adata.layers
        self.run_soupx_btn.setEnabled(soupx_ready and not already_corrected)
        if already_corrected:
            self.soupx_status.setText("Ambient correction already applied")
            self.soupx_status.set_state('success')
            self._soupx_card.set_completed(True)
            self._soupx_card.set_summary(["Ambient correction applied"])
        elif not soupx_ready:
            self.soupx_status.setText(
                "Needs raw counts -- run SoupX before Normalize")
            self.soupx_status.set_state('info')
            self._soupx_card.set_completed(False)
            self._soupx_card.set_summary(["Needs raw counts"])
        else:
            self.soupx_status.setText("")
            self._soupx_card.set_completed(False)
            self._soupx_card.set_summary(["Optional; not run"])

        # Normalization
        max_val = float(self.adata.X.max()) if hasattr(self.adata.X, 'max') else 0
        is_normalized = max_val < 20
        if is_normalized:
            self.norm_status.set_state('success')
            self.normalize_btn.setEnabled(False)
            self.norm_status.setText("Data is already normalized. Proceed to Cluster.")
            self._norm_card.set_completed(True)
            target = getattr(self, '_last_norm_target', None)
            self._norm_card.set_summary(
                [f"target {target:,} · log1p"] if target else ["Log-normalised"])
        else:
            self.norm_status.set_state('warning')
            self.normalize_btn.setEnabled(True)
            self.norm_status.setText("Normalize after QC filtering and doublet removal")
            self._norm_card.set_completed(False)
            self._norm_card.set_summary(["Not run yet"])


    # ------------------------------------------------------------------
    # QC filtering
    # ------------------------------------------------------------------

    _THRESHOLD_SPINS = ('qc_min_genes_spin', 'qc_max_genes_spin',
                        'qc_min_counts_spin', 'qc_max_counts_spin',
                        'qc_max_mt_spin')

    def _stash_thresholds(self, mode: str) -> None:
        """Remember the threshold values belonging to ``mode``."""
        self._mode_thresholds[mode] = {
            name: getattr(self, name).value() for name in self._THRESHOLD_SPINS}

    def _restore_thresholds(self, mode: str) -> bool:
        """Put ``mode``'s remembered values back. False if none stored."""
        stored = self._mode_thresholds.get(mode)
        if not stored:
            return False
        for name, value in stored.items():
            getattr(self, name).setValue(value)
        return True

    def _toggle_qc_mode(self, mode_text):
        """Enable/disable QC widgets based on selected mode.

        Each mode keeps its own thresholds. Switching to MAD used to
        overwrite the spin boxes with computed values and switching back
        left them there, so hand-typed limits were silently lost the
        moment you looked at what MAD suggested.
        """
        is_fixed = mode_text == "Fixed Thresholds"
        previous = getattr(self, '_qc_mode', None)
        if previous is not None and previous != mode_text:
            self._stash_thresholds(previous)
        self._qc_mode = mode_text

        self.qc_mad_spin.setEnabled(not is_fixed)
        # Threshold lines: movable in Fixed mode, static in MAD mode
        self._qc_histogram.set_lines_movable(is_fixed)

        if is_fixed:
            # Restore what was typed before, and redraw so the histogram
            # markers follow the spin boxes rather than the MAD values.
            if self._restore_thresholds(mode_text):
                self._sync_histogram_to_spins()
        else:
            self._preview_mad_thresholds()

    def _qc_obs(self):
        """The obs rows QC thresholds and histograms are built from.

        Cells marked '_role == exclude' are left out. They are still
        filtered by whatever thresholds you settle on -- the file stays
        internally consistent -- they just do not get to *define* those
        thresholds. On GSE292067 the excluded doxorubicin arm is 44% of
        the cells, so leaving it in means the median and MAD that bound
        your DCM and donor cells are set largely by a cohort you have
        already said you do not want.
        """
        from kosmic.scrna.inspect.roles import included_mask

        if self.adata is None:
            return None, 0, 0
        obs = self.adata.obs
        mask = included_mask(self.adata)
        if mask is None:
            return obs, len(obs), 0
        return obs.loc[mask], int(mask.sum()), int((~mask).sum())

    def _refresh_histogram_source(self):
        """Point the histogram at the included cells and say so."""
        obs, n_used, n_excluded = self._qc_obs()
        if obs is None:
            return
        self._qc_histogram.set_data(obs, n_excluded)
        if n_excluded:
            self.qc_status.setText(
                f"Thresholds computed from {n_used:,} cells; "
                f"{n_excluded:,} marked 'exclude' are not counted here "
                f"(they are still filtered by the thresholds you set).")
            self.qc_status.set_state('info')

    def _sync_histogram_to_spins(self):
        """Push the spin-box values onto the histogram's threshold lines."""
        pairs = (
            ("n_genes_by_counts", "min", self.qc_min_genes_spin),
            ("n_genes_by_counts", "max", self.qc_max_genes_spin),
            ("total_counts", "min", self.qc_min_counts_spin),
            ("total_counts", "max", self.qc_max_counts_spin),
            ("pct_counts_mt", "max", self.qc_max_mt_spin),
        )
        for metric, bound, spin in pairs:
            self._qc_histogram.set_threshold(metric, bound, float(spin.value()))
        self._qc_histogram._draw()

    def _on_mad_spin_changed(self, _value):
        """Re-preview MAD thresholds when the multiplier changes."""
        if self.qc_mode_combo.currentText() == "MAD-based":
            self._preview_mad_thresholds()

    def _preview_mad_thresholds(self):
        """Compute MAD thresholds and show them in spinboxes + histogram without filtering."""
        if self.adata is None:
            return
        self._ensure_qc_metrics()
        if 'n_genes_by_counts' not in self.adata.obs.columns:
            return

        from kosmic.scrna.qc.mad_thresholds import compute_mad_thresholds
        obs, _n_used, _n_excluded = self._qc_obs()
        bounds = compute_mad_thresholds(
            obs, nmads=self.qc_mad_spin.value())

        min_genes = bounds['n_genes_by_counts']['min']
        max_genes = bounds['n_genes_by_counts']['max']
        min_counts = bounds['total_counts']['min']
        max_counts = bounds['total_counts']['max']
        max_mt = bounds['pct_counts_mt']['max']

        # Update spinboxes (read-only in MAD mode, just showing values)
        self.qc_min_genes_spin.setValue(int(round(min_genes)))
        self.qc_max_genes_spin.setValue(int(round(max_genes)))
        self.qc_min_counts_spin.setValue(int(round(min_counts)))
        self.qc_max_counts_spin.setValue(int(round(max_counts)))
        self.qc_max_mt_spin.setValue(int(round(max_mt)))

        # Update histogram threshold lines
        self._qc_histogram.set_threshold("n_genes_by_counts", "min", float(min_genes))
        self._qc_histogram.set_threshold("n_genes_by_counts", "max", float(max_genes))
        self._qc_histogram.set_threshold("total_counts", "min", float(min_counts))
        self._qc_histogram.set_threshold("total_counts", "max", float(max_counts))
        self._qc_histogram.set_threshold("pct_counts_mt", "max", float(max_mt))

        # Redraw to show the new lines
        self._qc_histogram._draw()

    def _on_histogram_threshold_changed(self, metric: str, bound: str, value: float):
        """Update spinboxes when threshold lines are dragged in the histogram."""
        if metric == "n_genes_by_counts" and bound == "min":
            self.qc_min_genes_spin.setValue(int(round(value)))
        elif metric == "n_genes_by_counts" and bound == "max":
            self.qc_max_genes_spin.setValue(int(round(value)))
        elif metric == "total_counts" and bound == "min":
            self.qc_min_counts_spin.setValue(int(round(value)))
        elif metric == "total_counts" and bound == "max":
            self.qc_max_counts_spin.setValue(int(round(value)))
        elif metric == "pct_counts_mt" and bound == "max":
            self.qc_max_mt_spin.setValue(int(round(value)))

    def _apply_qc_filters(self):
        """Apply QC filters to the data (fixed thresholds or MAD-based)."""
        if self.adata is None or self.h5ad_path is None:
            return

        n_before = self.adata.n_obs
        use_mad = self.qc_mode_combo.currentText() == "MAD-based"

        if use_mad:
            self._apply_qc_mad(n_before)
        else:
            self._apply_qc_fixed(n_before)

    def _apply_qc_fixed(self, n_before):
        """Apply QC using fixed thresholds."""
        min_genes = self.qc_min_genes_spin.value()
        max_genes = self.qc_max_genes_spin.value()
        min_counts = self.qc_min_counts_spin.value()
        max_counts = self.qc_max_counts_spin.value()
        max_mt = self.qc_max_mt_spin.value()

        # Build filter description
        filters = [
            f"  - Min genes/cell: {min_genes}",
            f"  - Max genes/cell: {max_genes if max_genes > 0 else 'off'}",
        ]
        if min_counts > 0:
            filters.append(f"  - Min counts/cell: {min_counts:,}")
        if max_counts > 0:
            filters.append(f"  - Max counts/cell: {max_counts:,}")
        filters.append(f"  - Max MT %: {max_mt}%")

        contam_key = self.qc_contam_combo.currentData()
        contam_spec = self._contam_panels.get(contam_key) if contam_key else None
        if contam_spec:
            filters.append(
                f"  - Score contamination panel: {contam_spec.get('display_name', contam_key)} "
                f"(no cells removed)")

        _obs, n_included, n_excluded = self._qc_obs()
        count_line = f"Current cell count: {n_before:,}"
        if n_excluded:
            count_line += (f" ({n_included:,} included, {n_excluded:,} marked "
                           "'exclude' -- kept in the file and filtered by the "
                           "same thresholds, so a later sensitivity run can "
                           "bring them back)")
        if not dialogs.confirm(self, "Apply QC Filters", "This will filter cells with the following criteria:\n\n"
            + "\n".join(filters) + "\n\n"
            f"{count_line}\n\n"
            "Proceed?"):
            return

        # Launch the pipeline on a worker thread. Filter + metrics
        # recompute + write_h5ad are all O(nnz) on a sparse matrix and
        # block the UI on big datasets if run on the main thread.
        self.apply_qc_btn.setEnabled(False)
        if self.progress_bar:
            self.progress_bar.setRange(0, 0)  # indeterminate

        save_path = self._get_save_path()
        self._qc_n_before = n_before  # finished handler reads this back

        params = {
            'min_genes': min_genes,
            'max_genes': max_genes,
            'min_counts': min_counts,
            'max_counts': max_counts,
            'max_mt': max_mt,
        }
        if contam_spec and contam_spec.get('genes'):
            params['signature_panels'] = {contam_key: contam_spec['genes']}

        self._last_qc_params = params
        self._qc_filter_worker = QCFilterWorker(
            self.adata,
            params,
            str(save_path),
        )
        run_worker(
            self._qc_filter_worker,
            on_finished=self._on_qc_filter_finished,
            on_failed=self._on_qc_filter_failed,
            on_progress=self._on_status,
        )

    def _on_qc_filter_finished(self, payload):
        """Handle QCFilterWorker completion: rebind adata, update UI, save."""
        adata, _stats = payload
        self.apply_qc_btn.setEnabled(True)
        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)

        self.adata = adata
        save_path = self._get_save_path()
        self.main_window.set_adata(adata, str(save_path), switch=False)

        n_before = getattr(self, '_qc_n_before', adata.n_obs)
        n_after = adata.n_obs
        n_removed = n_before - n_after

        if hasattr(self.main_window, 'record_provenance'):
            p = {k: v for k, v in getattr(self, '_last_qc_params', {}).items()
                 if k != 'signature_panels'}  # gene lists too bulky for the record
            if 'signature_panels' in getattr(self, '_last_qc_params', {}):
                p['contamination_panels'] = list(self._last_qc_params['signature_panels'].keys())
            p['n_cells_before'] = int(n_before)
            p['n_cells_after'] = int(n_after)
            self.main_window.record_provenance('qc', p)
        self.main_window.mark_step_complete(3)

        self.qc_status.setText(
            f"QC Applied: {n_before:,} -> {n_after:,} cells ({n_removed:,} removed)")
        self.qc_status.set_state('success')
        self._update_controls()
        self.log_message.emit(
            f"QC: {n_before:,} -> {n_after:,} cells ({n_removed:,} removed)")

        # Refresh histogram with filtered data
        self._refresh_histogram_source()

        self.status_label.setText("Ready")
        pct = (n_removed / n_before * 100) if n_before else 0
        dialogs.info(
            self, "QC Complete",
            f"QC filtering complete!\n\n"
            f"Before: {n_before:,} cells\n"
            f"After: {n_after:,} cells\n"
            f"Removed: {n_removed:,} cells ({pct:.1f}%)\n\n"
            "Saved to file.")

    def _on_qc_filter_failed(self, message: str):
        self.apply_qc_btn.setEnabled(True)
        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)
        self.status_label.setText("Ready")
        dialogs.warning(self, "QC Failed", message)

    def _apply_qc_mad(self, n_before):
        """
        Apply QC using MAD-computed thresholds (from spinboxes, user may have overridden).

        The spinboxes already contain MAD-computed values from _preview_mad_thresholds().
        The user can adjust any value before clicking Apply. We just read the spinboxes.
        """
        # Read thresholds from spinboxes (MAD-computed, possibly user-adjusted)
        self._apply_qc_fixed(n_before)

    # ------------------------------------------------------------------
    # Scrublet doublet detection
    # ------------------------------------------------------------------

    def _run_scrublet(self):
        """Run Scrublet doublet detection."""
        if self.adata is None or self.h5ad_path is None:
            return

        expected_rate = self.doublet_rate_spin.value() / 100.0
        min_counts = self.doublet_min_counts_spin.value()

        if not dialogs.confirm(self, "Run Scrublet", f"This will run Scrublet doublet detection with:\n\n"
            f"  - Expected doublet rate: {expected_rate*100:.0f}%\n"
            f"  - Min counts: {min_counts}\n\n"
            f"Current cell count: {self.adata.n_obs:,}\n\n"
            "This may take a few minutes. Proceed?"):
            return

        self.run_scrublet_btn.setEnabled(False)
        self.filter_doublets_btn.setEnabled(False)
        if self.progress_bar:
            self.progress_bar.setValue(0)

        self.scrublet_worker = ScrubletWorker(
            self.adata,
            expected_rate,
            min_counts,
            str(self.h5ad_path)
        )
        run_worker(
            self.scrublet_worker,
            on_finished=self._on_scrublet_finished,
            on_failed=self._on_scrublet_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_progress,
        )

    def _on_scrublet_finished(self, payload):
        """Handle Scrublet completion."""
        adata, message, results = payload
        if self.progress_bar:
            self.progress_bar.setValue(0)
        self.run_scrublet_btn.setEnabled(True)

        self.adata = adata
        self.main_window.set_adata(adata, str(self.h5ad_path), switch=False)

        n_doublets = results['n_doublets']
        pct = results['pct_doublets']
        self.doublet_status.setText(f"Detected: {n_doublets:,} doublets ({pct:.1f}%)")
        self.doublet_status.set_state('warning')
        self.filter_doublets_btn.setEnabled(True)

        # Generate scrublet histogram
        if 'doublet_score' in self.adata.obs.columns:
            self._generate_scrublet_plot(
                self.adata.obs['doublet_score'].values,
                results['threshold'],
            )

        self._update_controls()
        self.log_message.emit(message)

        dialogs.info(self, "Scrublet Complete",
            f"{message}\n\n"
            f"Threshold: {results['threshold']:.3f}\n\n"
            f"New columns added:\n"
            f"  - doublet_score (continuous)\n"
            f"  - predicted_doublet (boolean)\n\n"
            f"Click 'Filter Doublets' to remove them.")

    def _on_scrublet_failed(self, message: str):
        if self.progress_bar:
            self.progress_bar.setValue(0)
        self.run_scrublet_btn.setEnabled(True)
        self.doublet_status.setText("Scrublet failed")
        self.doublet_status.set_state('error')
        dialogs.warning(self, "Scrublet Failed", message)

    def _filter_doublets(self):
        """Filter out predicted doublets."""
        if self.adata is None or self.h5ad_path is None:
            return

        if 'predicted_doublet' not in self.adata.obs.columns:
            dialogs.warning(self, "No Doublet Data",
                "Run Scrublet first to detect doublets.")
            return

        n_doublets = self.adata.obs['predicted_doublet'].sum()
        n_total = self.adata.n_obs

        if not dialogs.confirm(self, "Filter Doublets", f"This will remove {n_doublets:,} predicted doublets "
            f"({n_doublets/n_total*100:.1f}%) from the dataset.\n\n"
            f"Before: {n_total:,} cells\n"
            f"After: {n_total - n_doublets:,} cells\n\n"
            "Proceed?"):
            return

        # Slicing + write_h5ad on a 25k+ cell sparse AnnData blocks the
        # main thread for many seconds. Hand off to a worker.
        self.filter_doublets_btn.setEnabled(False)
        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)

        self._filter_doublets_worker = FilterDoubletsWorker(
            self.adata, str(self.h5ad_path))
        run_worker(
            self._filter_doublets_worker,
            on_finished=self._on_filter_doublets_finished,
            on_failed=self._on_filter_doublets_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_progress,
        )

    def _on_filter_doublets_finished(self, payload):
        adata, n_removed = payload
        if self.progress_bar:
            self.progress_bar.setValue(0)

        self._doublets_removed = int(n_removed)
        self.adata = adata
        self.main_window.set_adata(adata, str(self.h5ad_path))

        self.doublet_status.setText(
            f"Filtered: {n_removed:,} doublets removed, "
            f"{adata.n_obs:,} cells remain")
        self.doublet_status.set_state('success')
        self._update_controls()
        self.log_message.emit(
            f"Filtered {n_removed:,} doublets, {adata.n_obs:,} cells remain")
        self.status_label.setText("Ready")

        dialogs.info(
            self, "Doublets Filtered",
            f"Removed {n_removed:,} doublets.\n\n"
            f"Remaining cells: {adata.n_obs:,}\n\n"
            "Saved to file.")

    def _on_filter_doublets_failed(self, message: str):
        self.filter_doublets_btn.setEnabled(True)
        if self.progress_bar:
            self.progress_bar.setValue(0)
        self.status_label.setText("Ready")
        dialogs.warning(self, "Filter Failed", message)

    # ------------------------------------------------------------------
    # SoupX ambient correction
    # ------------------------------------------------------------------

    def _on_soupx_mode_changed(self, text: str):
        """Toggle the manual-fraction spin enable based on mode."""
        is_manual = text == "Manual fraction"
        self.soupx_frac_spin.setEnabled(is_manual)
        # tf-idf is only meaningful in auto mode
        self.soupx_tfidf_spin.setEnabled(not is_manual)

    def _browse_raw_10x(self):
        """Pick a CellRanger raw_feature_bc_matrix folder."""
        from PyQt6.QtWidgets import QFileDialog
        start = str(self.project_dir) if self.project_dir else ""
        folder = QFileDialog.getExistingDirectory(
            self, "Select CellRanger raw_feature_bc_matrix folder", start)
        if not folder:
            return
        self._soupx_raw_path = folder
        self.soupx_raw_label.setText(Path(folder).name)
        self.soupx_raw_label.setProperty("role", "")
        self.soupx_raw_label.style().unpolish(self.soupx_raw_label)
        self.soupx_raw_label.style().polish(self.soupx_raw_label)

    def _run_soupx(self):
        """Launch SoupXWorker."""
        if self.adata is None or self.h5ad_path is None:
            return

        is_manual = self.soupx_mode_combo.currentText() == "Manual fraction"
        contamination = float(self.soupx_frac_spin.value()) if is_manual else None
        tfidf_min = float(self.soupx_tfidf_spin.value())

        msg_lines = [
            f"Mode: {'manual fraction' if is_manual else 'auto-estimate'}",
        ]
        if is_manual:
            msg_lines.append(f"Contamination fraction: {contamination:.0%}")
        else:
            msg_lines.append(f"tf-idf min: {tfidf_min}")
        msg_lines.append(
            f"Raw matrix: {Path(self._soupx_raw_path).name if self._soupx_raw_path else '(filtered-only fallback)'}")

        if not dialogs.confirm(self, "Run SoupX", "This will estimate ambient RNA contamination and replace "
            "adata.X with corrected counts. Original counts are preserved "
            "in adata.layers['counts_pre_soupx'].\n\n"
            + "\n".join(msg_lines) + "\n\nProceed?"):
            return

        self.run_soupx_btn.setEnabled(False)
        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)

        save_path = self._get_save_path()
        self._soupx_worker = SoupXWorker(
            self.adata,
            str(save_path),
            raw_path=self._soupx_raw_path,
            contamination_fraction=contamination,
            tfidf_min=tfidf_min,
        )
        run_worker(
            self._soupx_worker,
            on_finished=self._on_soupx_finished,
            on_failed=self._on_soupx_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_progress,
        )

    def _on_soupx_finished(self, payload):
        adata, message, info = payload
        self.run_soupx_btn.setEnabled(True)
        if self.progress_bar:
            self.progress_bar.setValue(0)

        self.adata = adata
        save_path = self._get_save_path()
        self.main_window.set_adata(adata, str(save_path), switch=False)

        contam = info.get('contamination_fraction', 0)
        fallback_note = " (filtered-only fallback)" if info.get('used_filtered_as_raw') else ""
        self.soupx_status.setText(
            f"Applied: {contam:.1%} contamination removed{fallback_note}")
        self.soupx_status.set_state('success')
        self._update_controls()
        self.log_message.emit(message)
        self.status_label.setText("Ready")

        dialogs.info(
            self, "SoupX Complete",
            f"{message}\n\n"
            f"Original counts preserved in "
            f"adata.layers['counts_pre_soupx'].\n"
            "Saved to file.")

    def _on_soupx_failed(self, message: str):
        self.run_soupx_btn.setEnabled(True)
        if self.progress_bar:
            self.progress_bar.setValue(0)
        self.status_label.setText("Ready")
        self.soupx_status.setText("SoupX failed")
        self.soupx_status.set_state('error')
        dialogs.warning(self, "SoupX Failed", message)

    # ------------------------------------------------------------------
    # Normalization
    # ------------------------------------------------------------------

    def _normalize_data(self):
        """Normalize the data."""
        if self.adata is None or self.h5ad_path is None:
            return

        # Guard: block double-normalization
        max_val = float(self.adata.X.max()) if hasattr(self.adata.X, 'max') else 0
        if max_val < 20:
            dialogs.info(
                self, "Already Normalized",
                "Data appears already normalized (max value < 20).\n"
                "Normalizing again would corrupt the data."
            )
            self.normalize_btn.setEnabled(False)
            self.norm_status.setText("Data is already normalized. Proceed to Cluster.")
            self.norm_status.set_state('success')
            return

        # Guard: preserve raw counts before normalization
        if 'counts' not in self.adata.layers:
            self.adata.layers['counts'] = self.adata.X.copy()

        if not dialogs.confirm(self, "Normalize Data", f"This will normalize the data to {self.target_sum_spin.value():,} counts per cell"
            + (" and log-transform." if self.log_transform_check.isChecked() else ".")
            + "\n\nRaw counts will be stored in adata.raw.\n\nProceed?"):
            return

        self._last_norm_target = self.target_sum_spin.value()

        self.normalize_btn.setEnabled(False)
        save_path = self._get_save_path()
        self.normalize_worker = NormalizeWorker(
            self.adata,
            self.target_sum_spin.value(),
            self.log_transform_check.isChecked(),
            str(save_path)
        )
        self._pending_norm_path = save_path
        run_worker(
            self.normalize_worker,
            on_finished=self._on_normalize_finished,
            on_failed=self._on_normalize_failed,
            on_progress=self._on_status,
        )

    def _on_normalize_finished(self, adata):
        self.normalize_btn.setEnabled(True)
        self.adata = adata
        self.h5ad_path = self._pending_norm_path
        self.main_window.set_adata(adata, str(self.h5ad_path), switch=False)

        # Mark QC step complete after normalization
        self.main_window.mark_step_complete(3)

        self.norm_status.set_state('success')
        self.normalize_btn.setEnabled(False)
        self.norm_status.setText("Data is normalized. Proceed to Cluster.")
        self.log_message.emit("Normalization complete")

        dialogs.info(
            self, "Normalization Complete", "Normalization complete")

    def _on_normalize_failed(self, message: str):
        self.normalize_btn.setEnabled(True)
        dialogs.warning(
            self, "Normalization Failed", f"Normalization failed: {message}")
