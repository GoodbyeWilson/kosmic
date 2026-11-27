# Pathway-Level Meta-Analysis Page
#
# Pools each study's Pathway DE results across studies. In each study the
# Pathway DE page (DE workspace) scores every donor (DESeq2-normalised log2
# expression averaged over a pathway's genes) and compares disease and
# control donors with Welch's t-test, giving each pathway an effect (the
# difference in mean score, log2) and its standard error. Here those
# effects are pooled with REML (or DL) random effects, as genes are, and
# BH-corrected. Only tables made with that score and test are pooled; a
# table saved by an earlier version (which may have used one of the
# scoring methods since removed) is named, not pooled.


from PyQt6.QtWidgets import (
    QPushButton, QWidget, QVBoxLayout, QHBoxLayout,
    QSplitter, QComboBox, QSpinBox,
    QDoubleSpinBox, QFileDialog, QTabWidget, QHeaderView,
    QAbstractItemView,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QFont

import numpy as np
import pandas as pd
import pyqtgraph as pg

from kosmic.gui.shared.theme import get_color, style_pg_plot, get_font_sizes
from kosmic.gui.shared.widgets import (
    BaseWorker, Column, SettingsGroup, SidebarPage, SecondaryLabel,
    HintLabel, PrimaryButton, SecondaryButton, ResultsTable,
)
from kosmic.gui.shared import run_worker
from kosmic.numerical import neg_log10
from kosmic.paths import meta_output_dir
from kosmic import DEFAULT_FDR
from kosmic.meta_analysis.direction import count_significant_conflicts
from kosmic.de.pathway_scoring import SCORING_METHOD


_RESULT_COLUMNS = [
    Column("Pathway",   "names",            "s"),
    Column("log2FC",    "logfoldchanges",   ".3f"),
    Column("SE",        "se",               ".3f"),
    Column("CI low",    "ci_lower",         ".3f"),
    Column("CI high",   "ci_upper",         ".3f"),
    Column("Pooled P",  "pvals_pooled",     ".2e"),
    Column("FDR",       "fdr",              ".2e"),
    Column("k",         "n_studies",        "d"),
    Column("I^2",       "heterogeneity_i2", ".2f"),
]

# Shown when the per-study tables carry adjusted p-values (the Pathway DE
# results do).
_DIRECTION_COLUMNS = [
    Column("Up",       "n_up",               "d",
           tooltip="Studies with study FDR < 0.05 and log2FC > 0"),
    Column("Down",     "n_down",             "d",
           tooltip="Studies with study FDR < 0.05 and log2FC < 0"),
    Column("Conflict", "direction_conflict", "bool",
           tooltip="Significant up in one study and down in another"),
]


def usable_tables(datasets):
    """Split the per-study tables into those made with the pathway score and
    test the meta-analysis pools, and the names of the others (saved by an
    earlier version, or with another scoring method), which are not pooled."""
    usable = [d for d in datasets if d.get('scoring_method') == SCORING_METHOD]
    other = [d.get('name', '?') for d in datasets
             if d.get('scoring_method') != SCORING_METHOD]
    return usable, other


class PathwayMAWorker(BaseWorker):
    """
    Pool the per-study pathway results off the UI thread.
    Emits finished_ok with a tuple '(meta_df, method_keys)'.
    """

    def __init__(self, datasets, method_keys, min_studies, parent=None):
        super().__init__(parent)
        self.datasets = datasets
        self.method_keys = method_keys
        self.min_studies = min_studies

    def _run(self):
        from kosmic.numerical import bh_fdr
        from kosmic.meta_analysis.pooling_dispatch import get_analytical_pool_fn
        from kosmic.meta_analysis.direction import unadjusted_warning

        dfs = [d['df'] for d in self.datasets]
        warning = unadjusted_warning(dfs)
        if warning:
            self.progress.emit(f"Warning: {warning}")
        # The same vectorised pooling as the gene-level meta-analysis.
        method_key = self.method_keys[0] if self.method_keys else 'reml'
        self.progress.emit(f"Pooling ({method_key})...")
        pool_fn = get_analytical_pool_fn(method_key, self.min_studies)
        result = pool_fn(dfs, min_studies=self.min_studies)
        if result.empty:
            raise RuntimeError(f"No pathways in >= {self.min_studies} studies.")
        result['fdr'] = bh_fdr(result['pvals_pooled'])
        result = result.sort_values('pvals_pooled').reset_index(drop=True)
        return (result, self.method_keys)


class PathwayMAPage(SidebarPage):
    """Pathway-level meta-analysis (single pooling method)."""

    help_id = "meta/pathway_ma"
    SIDEBAR_SCROLLABLE = False

    log_message = pyqtSignal(str)
    analysis_complete = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._datasets = []
        self._gene_datasets = []  # for finding pseudobulk paths
        self._project_folder = None
        self._output_selection = None  # folder under meta_analysis/ (ADR-007)
        self._run_output_selection = None
        self._run_study_names = []
        self._meta_df = None  # pooled pathway-level DataFrame from the most recent run
        self._method_keys = []  # active method keys for the current run
        self._scatter_data = None
        self._forest_xlim = None
        self._pathway_gene_sets = {}
        self._setup_ui()

    @property
    def meta_df(self):
        """Pooled pathway-level DataFrame from the most recent run, or 'None'."""
        return self._meta_df

    @property
    def pathway_gene_sets(self) -> dict:
        """Gene sets used in the most recent run ({pathway_name: [genes]})."""
        return self._pathway_gene_sets

    def _setup_ui(self):
        # Sidebar lives on self.sidebar_layout (provided by SidebarPage).
        lay = self.sidebar_layout

        info_group = SettingsGroup("Loaded Studies")
        self._info_label = SecondaryLabel("No pathway DE datasets loaded.\n\n"
                                   "Run Pathway DE in the DE workspace,\n"
                                   "then scan the project folder.")
        self._info_label.setWordWrap(True)
        info_group.add_widget(self._info_label)
        lay.addWidget(info_group)

        about = HintLabel(
            "Pools each study's Pathway DE results. In each study, every donor "
            "gets one score per pathway (DESeq2-normalised log2 expression "
            "averaged over its genes) and Welch's t-test compares disease with "
            "control. Here the studies' effects are pooled with random effects, "
            "as genes are. Run Pathway DE in each study first.")
        about.setWordWrap(True)
        lay.addWidget(about)

        settings_group = SettingsGroup("Analysis Settings", collapsible=True, expanded=True)
        self._min_studies = QSpinBox()
        self._min_studies.setRange(2, 20)
        self._min_studies.setValue(2)
        settings_group.add_row("Min studies:", self._min_studies)

        self._pooling_combo = QComboBox()
        self._pooling_combo.addItem("REML", userData='reml')
        self._pooling_combo.addItem("DL", userData='dl')
        self._pooling_combo.setToolTip(
            "REML: restricted maximum likelihood (recommended).\n"
            "DL: DerSimonian-Laird method-of-moments.")
        settings_group.add_row("Pooling:", self._pooling_combo)

        lay.addWidget(settings_group)

        # Run
        self._run_btn = PrimaryButton("Run Pathway Meta-Analysis")
        self._run_btn.setMinimumHeight(32)
        self._run_btn.clicked.connect(self._run_analysis)
        lay.addWidget(self._run_btn)

        self._status_label = SecondaryLabel("")
        self._status_label.setWordWrap(True)
        lay.addWidget(self._status_label)

        thresh_group = SettingsGroup("Volcano Thresholds", collapsible=True, expanded=True)
        self._pval_spin = QDoubleSpinBox()
        self._pval_spin.setRange(0.001, 1)
        self._pval_spin.setSingleStep(0.01)
        self._pval_spin.setValue(DEFAULT_FDR)
        self._pval_spin.setDecimals(3)
        self._pval_spin.valueChanged.connect(self._refresh_volcano)
        thresh_group.add_row("P-value:", self._pval_spin)

        self._fc_spin = QDoubleSpinBox()
        self._fc_spin.setRange(0, 10)
        self._fc_spin.setSingleStep(0.1)
        self._fc_spin.setValue(0.0)
        self._fc_spin.setDecimals(2)
        self._fc_spin.valueChanged.connect(self._refresh_volcano)
        thresh_group.add_row("Min |log2FC|:", self._fc_spin)

        lay.addWidget(thresh_group)

        info2_group = SettingsGroup("Selected Pathway")
        self._pathway_info = HintLabel("Click a pathway on the volcano plot")
        self._pathway_info.setWordWrap(True)
        info2_group.add_widget(self._pathway_info)
        lay.addWidget(info2_group)

        export_group = SettingsGroup("Export", collapsible=True, expanded=True)
        self._export_btn = QPushButton("Save CSV")
        self._export_btn.setEnabled(False)
        self._export_btn.clicked.connect(self._export_csv)
        export_group.add_widget(self._export_btn)
        lay.addWidget(export_group)

        lay.addStretch()

        # ── Right panel: tabs for plots + results table ──
        right_layout = self.content_layout

        self._tabs = QTabWidget()

        # --- Tab 1: Plots (volcano + forest splitter) ---
        plots_splitter = QSplitter(Qt.Orientation.Vertical)

        from kosmic.gui.shared.plots import InteractiveVolcano
        self._volcano = InteractiveVolcano()
        self._volcano.gene_selected.connect(self._on_pathway_clicked)
        plots_splitter.addWidget(self._volcano)

        from kosmic.gui.shared.plots import InteractivePlot
        self._forest = InteractivePlot(
            title='Per-Study Forest Plot',
            left_label='Study',
            bottom_label='Log2 Fold Change',
            unavailable_message='Click a pathway on the volcano to see its per-study forest.',
        )
        self._forest.setMouseEnabled(x=True, y=False)
        self._forest.setMinimumHeight(150)
        plots_splitter.addWidget(self._forest)
        plots_splitter.setSizes([500, 200])

        self._tabs.addTab(plots_splitter, "Plots")

        # --- Tab 2: Results table ---
        results_panel = QWidget()
        results_layout = QVBoxLayout(results_panel)
        results_layout.setContentsMargins(0, 0, 0, 0)
        results_layout.setSpacing(4)

        # Toolbar above the table: explicit Plot Forest button so a row
        # click only selects, never silently switches tabs.
        toolbar = QHBoxLayout()
        toolbar.setContentsMargins(4, 4, 4, 0)
        self._plot_forest_btn = SecondaryButton("Plot Forest for Selected")
        self._plot_forest_btn.setEnabled(False)
        self._plot_forest_btn.setToolTip(
            "Render the per-study forest plot for the selected row "
            "and switch to the Plots tab."
        )
        self._plot_forest_btn.clicked.connect(self._plot_forest_for_selected)
        toolbar.addWidget(self._plot_forest_btn)
        toolbar.addStretch()
        results_layout.addLayout(toolbar)

        self._results_table = ResultsTable()
        self._results_table.set_schema(_RESULT_COLUMNS)
        self._results_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self._results_table.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection)
        self._results_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self._results_table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.Stretch)
        self._results_table.itemSelectionChanged.connect(
            self._on_results_selection_changed)
        results_layout.addWidget(self._results_table)

        self._tabs.addTab(results_panel, "Results")

        right_layout.addWidget(self._tabs)

    # ------------------------------------------------------------------
    # Data
    # ------------------------------------------------------------------

    def set_datasets(self, datasets: list, gene_datasets=None,
                     project_folder=None, pathway_gene_sets=None,
                     output_selection=None):
        self._output_selection = output_selection
        self._datasets = datasets
        self._gene_datasets = gene_datasets or []
        self._project_folder = project_folder
        if pathway_gene_sets is not None:
            self._pathway_gene_sets = pathway_gene_sets
        n = len(datasets)
        if n == 0:
            self._info_label.setText(
                "No pathway DE datasets loaded.\n\n"
                "Run Pathway DE in the DE workspace,\n"
                "then scan the project folder."
            )
        else:
            usable, other = usable_tables(datasets)
            lines = [f"{len(usable)} of {n} studies with pathway DE can be pooled:"]
            for d in usable:
                lines.append(f"  {d['name']} ({len(d['df'])} pathways)")
            if other:
                lines.append("")
                lines.append("Not pooled, made by an earlier version of Pathway "
                             "DE (re-run it there): " + ", ".join(other))
            self._info_label.setText("\n".join(lines))

    # ------------------------------------------------------------------
    # Analysis
    # ------------------------------------------------------------------

    def _run_analysis(self):
        if not self._datasets:
            self._status_label.setText("No pathway DE datasets loaded.")
            return

        pooling_key = self._pooling_combo.currentData() or 'reml'
        method_keys = [pooling_key]
        min_studies = self._min_studies.value()

        usable, other = usable_tables(self._datasets)
        if other:
            self.log_message.emit(
                "Pathway MA: not pooled, made by an earlier version of Pathway "
                "DE (re-run it there): " + ", ".join(other))
        if len(usable) < min_studies:
            self._status_label.setText(
                f"Need >= {min_studies} studies made with the current Pathway DE, "
                f"have {len(usable)}." + (
                    " Re-run Pathway DE in: " + ", ".join(other) if other else ""))
            return
        pw_gene_sets = self._pathway_gene_sets

        # Snapshot settings now so the sidecar JSON written on success
        # reflects exactly what the user picked, not whatever the
        # widgets show by the time the worker returns.
        self._last_run_settings = {
            'pooling_method':  pooling_key,
            'min_studies':     min_studies,
            'scoring_method':  SCORING_METHOD,
            'per_study_test':  "Welch's t-test",
            'n_studies_input': len(usable),
            'study_names':     [d.get('name', '?') for d in usable],
            'not_pooled':      other,
            'n_pathways_in_geneset': len(pw_gene_sets) if pw_gene_sets else 0,
        }

        # Clear stale per-pathway forest from any previous run before
        # the new pooling kicks off.
        self._forest.clear_plot_items()
        self._forest.show_unavailable_message()

        self._run_btn.setEnabled(False)
        self._status_label.setText("Running pathway meta-analysis...")
        self.log_message.emit(
            f"Pathway MA: {len(usable)} studies, pooling={pooling_key}")

        # Fixed at launch: the selection may change while it runs.
        self._run_output_selection = self._output_selection
        self._run_study_names = [d.get('name', '?') for d in usable]
        self._worker = PathwayMAWorker(usable, method_keys, min_studies)
        run_worker(
            self._worker,
            on_finished=self._on_worker_finished,
            on_failed=self._on_worker_failed,
            on_progress=self._on_worker_progress,
        )

    def _on_worker_progress(self, message):
        self._status_label.setText(message)
        self.log_message.emit(f"  {message}")

    def _on_worker_finished(self, payload):
        meta_df, method_keys = payload
        self._meta_df = meta_df
        self._method_keys = method_keys
        n_sig, n_conflict = count_significant_conflicts(
            meta_df, fdr=DEFAULT_FDR)
        conflict_text = (
            f", {n_conflict} with opposite-direction effects across studies"
            if n_conflict is not None else "")
        self._status_label.setText(
            f"{len(meta_df)} pathways tested, "
            f"{n_sig} significant (FDR < 0.05){conflict_text}.")
        self.log_message.emit(
            f"Pathway MA: {n_sig}/{len(meta_df)} significant{conflict_text}")
        self._populate_volcano(meta_df)
        self._populate_results_table(meta_df)
        self._auto_save(meta_df, method_keys)
        self._record_meta_provenance('meta_pathway')
        self._export_btn.setEnabled(True)
        self._run_btn.setEnabled(True)
        self.analysis_complete.emit()

    def _record_meta_provenance(self, stage):
        """Append the pathway meta-analysis step to the project-level sidecar so
        the combined methods view can show it. Best-effort."""
        if not self._project_folder:
            return
        try:
            from kosmic.meta_analysis.io import (
                gather_study_provenance, record_meta_stage)
            tokens = {
                label: (rec.get('fingerprint') or {}).get('token')
                for label, rec in gather_study_provenance(self._project_folder)
            }
            params = {
                'pooling_method': self._pooling_combo.currentData() or 'reml',
                'n_studies': len(self._run_study_names),
                'studies': list(self._run_study_names),
                'n_pathway_sets': len(self._pathway_gene_sets or {}),
                'study_tokens': tokens,
            }
            if self._meta_df is not None:
                n_sig, n_conflict = count_significant_conflicts(
                    self._meta_df, fdr=DEFAULT_FDR)
                params['n_significant'] = n_sig
                if n_conflict is not None:
                    params['n_significant_direction_conflict'] = n_conflict
                    params['direction_study_fdr'] = DEFAULT_FDR
            # Shared recorder: keeps the latest run, as every meta
            # stage does. A direct record_stage call appends one entry
            # per Run click.
            record_meta_stage(self._project_folder, stage, params,
                              selection=self._run_output_selection)
        except Exception as e:
            self.log_message.emit(f"Could not record pathway meta provenance: {e}")

    def _populate_results_table(self, meta_df):
        """Render the meta-analysis results in the Results tab, sorted by FDR."""
        if meta_df is None or meta_df.empty:
            self._results_table.set_data(pd.DataFrame())
            return
        schema = list(_RESULT_COLUMNS)
        if 'direction_conflict' in meta_df.columns:
            schema += _DIRECTION_COLUMNS
        self._results_table.set_schema(schema)
        cols = [c.key for c in schema if c.key in meta_df.columns]
        df = meta_df[cols].copy().sort_values('fdr').reset_index(drop=True)
        self._results_table.set_data(df)

    def _on_results_selection_changed(self):
        """Selection just enables the Plot Forest button -- nothing else."""
        has_selection = bool(
            self._results_table.selectionModel().selectedRows()
        )
        self._plot_forest_btn.setEnabled(has_selection)

    def _plot_forest_for_selected(self):
        """Render the forest for the currently selected row + switch to Plots."""
        rows = self._results_table.selectionModel().selectedRows()
        if not rows or self._meta_df is None:
            return
        row = rows[0].row()
        item = self._results_table.item(row, 0)
        if item is None:
            return
        pathway_name = item.text()
        self._on_pathway_clicked(pathway_name)
        self._tabs.setCurrentIndex(0)  # Plots tab

    def _auto_save(self, meta_df, method_keys):
        """Save the pooled DataFrame + settings sidecar to the selection's folder
        under {project}/meta_analysis/ (ADR-007).

        Both files are keyed by the pooling method, so re-running with
        the same method overwrites the pair (settings stay in sync with
        the CSV). To preserve a previous run, use the manual Save CSV
        button to copy the file elsewhere.
        """
        if not self._project_folder or meta_df is None or meta_df.empty:
            return
        try:
            output_dir = meta_output_dir(self._project_folder,
                                         self._run_output_selection)
            output_dir.mkdir(parents=True, exist_ok=True)
            method = method_keys[0] if method_keys else 'pool'
            stem = f'pathway_ma_{method}'

            csv_path = output_dir / f'{stem}.csv'
            meta_df.to_csv(csv_path, index=False)

            json_path = output_dir / f'{stem}.json'
            self._write_settings_sidecar(json_path, meta_df, method_keys)

            self.log_message.emit(
                f"Pathway MA results saved to "
                f"{output_dir.relative_to(self._project_folder).as_posix()}: "
                f"{csv_path.name} + {json_path.name}"
            )
            from kosmic.meta_analysis.io import mixed_selection_warning
            warning = mixed_selection_warning(self._run_output_selection)
            if warning:
                self.log_message.emit(warning)
        except Exception as exc:
            self.log_message.emit(f"Pathway MA auto-save failed: {exc}")

    def _write_settings_sidecar(self, path, meta_df, method_keys):
        """Write a JSON sidecar capturing settings + result summary for the run."""
        from kosmic.meta_analysis.io import write_meta_settings_sidecar

        settings = dict(getattr(self, '_last_run_settings', {}) or {})
        settings['method_keys'] = list(method_keys)

        write_meta_settings_sidecar(
            path,
            tool='KOSMIC Pathway Meta-Analysis',
            analysis_type='Pathway-Level Meta-Analysis (single pooling method)',
            parameters=settings,
            results_summary={
                'n_pathways_tested': len(meta_df),
                'n_significant':     int((meta_df['fdr'] < DEFAULT_FDR).sum()),
                'fdr_threshold':     DEFAULT_FDR,
            },
        )

    def _on_worker_failed(self, message):
        self._status_label.setText(f"Failed: {message}")
        self.log_message.emit(f"Pathway MA failed: {message}")
        self._run_btn.setEnabled(True)


    # ------------------------------------------------------------------
    # Volcano
    # ------------------------------------------------------------------

    def _populate_volcano(self, meta_df):
        """Render the pathway-level volcano using the shared widget."""
        self._meta_df = meta_df
        self._forest_xlim = None

        if meta_df is None or len(meta_df) == 0:
            self._volcano.clear_data()
            return

        # Pre-compute shared forest x-limits (scan both pooled CIs and
        # per-study effects so the per-study forest below doesn't clip).
        if 'study_effects' in meta_df.columns:
            all_lo, all_hi = [], []
            for _, r in meta_df.iterrows():
                effects = r.get('study_effects', [])
                if isinstance(effects, list):
                    for eff in effects:
                        lfc = eff.get('logfc', 0.0)
                        se = eff.get('se', float('nan'))
                        margin = 1.96 * se if np.isfinite(se) and se > 0 else 0.0
                        all_lo.append(lfc - margin)
                        all_hi.append(lfc + margin)
                if 'ci_lower' in r.index and np.isfinite(
                        r.get('ci_lower', float('nan'))):
                    all_lo.append(r['ci_lower'])
                    all_hi.append(r['ci_upper'])
            if all_lo:
                extreme = max(abs(min(all_lo)), abs(max(all_hi)))
                pad = extreme * 0.05
                self._forest_xlim = (-(extreme + pad), extreme + pad)

        x = meta_df['logfoldchanges'].values.astype(float)
        pvals = meta_df['pvals_pooled'].values.astype(float)
        pvals = np.where(np.isfinite(pvals) & (pvals > 0), pvals, 1.0)
        y = neg_log10(pvals)
        names = meta_df['names'].tolist()

        n_studies = (meta_df['n_studies'].values
                     if 'n_studies' in meta_df.columns
                     else np.ones(len(meta_df)))

        pval_thresh = self._pval_spin.value()
        fc_thresh = self._fc_spin.value()

        # Up/down/NS colours + n_studies-scaled point size
        up = get_color('plot_upregulated')
        down = get_color('plot_downregulated')
        ns_color = get_color('plot_insignificant')
        brushes = []
        sizes = []
        for i in range(len(meta_df)):
            ns = int(n_studies[i]) if np.isfinite(n_studies[i]) else 1
            base_size = max(8, min(16, 6 + ns * 2))
            if pvals[i] < pval_thresh and abs(x[i]) > fc_thresh:
                brushes.append(pg.mkBrush(up if x[i] > 0 else down))
                sizes.append(base_size)
            else:
                brushes.append(pg.mkBrush(ns_color))
                sizes.append(max(6, base_size - 2))

        fs = get_font_sizes()

        def _fmt(name, gx, gy):
            extra = ''
            row = self._meta_df[self._meta_df['names'] == name]
            if len(row) > 0 and 'n_studies' in row.columns:
                extra = f"\nStudies: {int(row.iloc[0]['n_studies'])}"
                i2 = row.iloc[0].get('heterogeneity_i2', None)
                if i2 is not None and np.isfinite(i2):
                    extra += f"\nI\u00b2: {i2:.1%}"
            p = pvals[names.index(name)] if name in names else float('nan')
            return f"{name}\nlog2FC: {gx:.3f}\np: {p:.2e}{extra}"

        self._volcano.set_scatter(
            x, y, names, brushes=brushes, sizes=sizes,
            x_label='Pooled Log2 Fold Change',
            y_label='-Log10 P-value',
            title='Pathway Meta-Analysis Volcano',
            hover_fmt=_fmt,
        )

        # Threshold lines (p-value + fold-change cutoffs)
        self._volcano.add_threshold_line(y=-np.log10(pval_thresh))
        if fc_thresh > 0:
            self._volcano.add_threshold_line(x=fc_thresh)
            self._volcano.add_threshold_line(x=-fc_thresh)

        # Pathway labels next to each point (small N so label-every is fine)
        for i, label in enumerate(names):
            short = label if len(label) <= 35 else label[:32] + '...'
            txt = pg.TextItem(text=short, color=get_color('fg_secondary'),
                              anchor=(0, 0.5))
            txt.setFont(QFont('Arial', max(7, int(fs['tick']) - 1)))
            txt.setPos(x[i] + 0.02, y[i])
            self._volcano.addItem(txt)

        self._volcano.set_symmetric_x(pad_frac=0.15)  # extra pad for labels

    def _refresh_volcano(self):
        if self._meta_df is not None and len(self._meta_df) > 0:
            self._populate_volcano(self._meta_df)

    def _on_pathway_clicked(self, pathway_name: str) -> None:
        """Volcano click -> drill down into the per-study forest plot."""
        self._show_forest(pathway_name)

    # ------------------------------------------------------------------
    # Forest plot
    # ------------------------------------------------------------------

    def _show_forest(self, pathway_name):
        self._forest.clear_plot_items()
        self._forest.hide_unavailable_message()
        self._forest.enableAutoRange(axis='x', enable=False)

        if self._meta_df is None:
            return

        row = self._meta_df[self._meta_df['names'] == pathway_name]
        if row.empty:
            return
        row = row.iloc[0]

        study_effects = row.get('study_effects', [])
        if not study_effects or not isinstance(study_effects, list):
            self._pathway_info.setText(f"{pathway_name}: no per-study data")
            return

        fg = get_color('fg_primary')
        accent = get_color('accent_primary')
        n = len(study_effects)

        style_pg_plot(self._forest, title=f'{pathway_name} -- Forest Plot',
                      left_label='', bottom_label='Log2 Fold Change')

        y_positions = list(range(n))
        x_vals = [eff['logfc'] for eff in study_effects]
        ci_lo = [eff['logfc'] - 1.96 * eff['se'] for eff in study_effects]
        ci_hi = [eff['logfc'] + 1.96 * eff['se'] for eff in study_effects]
        labels = [eff.get('dataset', f'Study {i+1}') for i, eff in enumerate(study_effects)]

        # Error bars
        for i in range(n):
            self._forest.addItem(pg.PlotDataItem(
                [ci_lo[i], ci_hi[i]], [y_positions[i], y_positions[i]],
                pen=pg.mkPen(fg, width=1.5)))

        # Study points
        scatter = pg.ScatterPlotItem(
            x=np.array(x_vals), y=np.array(y_positions, dtype=float),
            size=10, pen=pg.mkPen(fg), brush=pg.mkBrush(accent),
            symbol='o')
        self._forest.addItem(scatter)

        # Pooled diamond
        pooled = row.get('logfoldchanges', 0)
        pooled_ci_lo = row.get('ci_lower', pooled)
        pooled_ci_hi = row.get('ci_upper', pooled)
        diamond_y = -1.0
        dh = 0.3
        diamond_x = [pooled_ci_lo, pooled, pooled_ci_hi, pooled, pooled_ci_lo]
        diamond_y_pts = [diamond_y, diamond_y + dh, diamond_y,
                         diamond_y - dh, diamond_y]
        diamond = pg.PlotDataItem(diamond_x, diamond_y_pts,
                                  pen=pg.mkPen(get_color('plot_upregulated'), width=2),
                                  fillLevel=diamond_y,
                                  brush=pg.mkBrush(get_color('plot_upregulated')))
        self._forest.addItem(diamond)

        # Zero line
        self._forest.addItem(pg.InfiniteLine(
            pos=0, angle=90,
            pen=pg.mkPen(fg, width=1, style=Qt.PenStyle.DashLine)))

        # Y-axis labels
        y_ticks = [(float(i), labels[i]) for i in range(n)]
        y_ticks.append((diamond_y, 'Pooled'))
        self._forest.getAxis('left').setTicks([y_ticks])

        self._forest.setYRange(diamond_y - 0.8, n - 0.5)
        if self._forest_xlim:
            self._forest.setXRange(*self._forest_xlim, padding=0)

        # Info label with per-method p-values
        i2 = row.get('heterogeneity_i2', float('nan'))
        i2_str = f'I\u00b2={i2:.0%}' if np.isfinite(i2) else 'I\u00b2=N/A'
        pval = row.get('pvals_pooled', float('nan'))
        fdr = row.get('fdr', float('nan'))

        method_lines = []
        for k in self._method_keys:
            col = f'pvals_{k}'
            p = row.get(col, float('nan'))
            if np.isfinite(p):
                method_lines.append(f"{k}: {p:.2e}")
        method_str = "  |  ".join(method_lines) if method_lines else ""

        self._pathway_info.setText(
            f"<b>{pathway_name}</b><br>"
            f"Pooled log2FC = {pooled:.3f} "
            f"[{pooled_ci_lo:.3f}, {pooled_ci_hi:.3f}]<br>"
            f"Pooled P = {pval:.2e}  |  FDR = {fdr:.2e}<br>"
            f"{method_str}<br>"
            f"{i2_str}  |  k = {n}")

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def _export_csv(self):
        if self._meta_df is None or self._meta_df.empty:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Pathway Meta-Analysis Results", "",
            "CSV files (*.csv)"
        )
        if path:
            cols = ['names', 'logfoldchanges', 'se',
                    'ci_lower', 'ci_upper',
                    'pvals_pooled', 'fdr',
                    'n_studies', 'heterogeneity_i2']
            for k in self._method_keys:
                cols.append(f'pvals_{k}')
            cols = [c for c in cols if c in self._meta_df.columns]
            export_df = self._meta_df[cols].copy()
            export_df.to_csv(path, index=False)
            self.log_message.emit(
                f"Pathway MA results exported to {path}")

    def refresh_theme(self):
        """Re-apply theme to pyqtgraph plots."""
        style_pg_plot(self._volcano, title='Pathway Meta-Analysis Volcano',
                      left_label='-Log10 Pooled P-value',
                      bottom_label='Pooled Log2 Fold Change')
        style_pg_plot(self._forest, title='Per-Study Forest Plot',
                      left_label='Study', bottom_label='Log2 Fold Change')
