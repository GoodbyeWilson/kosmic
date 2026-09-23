# Pathway-Level Meta-Analysis Page
#
# Single-method pooling at the pathway level. The user picks one
# pooling method (DL, HKSJ, REML, Stouffer's Z, sign test); the page
# runs it across studies and BH-corrects. Multi-method "consensus"
# (max-FDR across methods) is only on the gene-level page.


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
from kosmic import DEFAULT_FDR, MIN_STUDIES, PATHWAY_VIF_RHO
from kosmic.meta_analysis.direction import (
    add_wald_pvalues, count_significant_conflicts,
)


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

# Shown only when the per-study tables carry adjusted p-values: the
# pathway-DE results, or the Wald p-values added to the DESeq2+VIF
# summaries.
_DIRECTION_COLUMNS = [
    Column("Up",       "n_up",               "d",
           tooltip="Studies with study FDR < 0.05 and log2FC > 0"),
    Column("Down",     "n_down",             "d",
           tooltip="Studies with study FDR < 0.05 and log2FC < 0"),
    Column("Conflict", "direction_conflict", "bool",
           tooltip="Significant up in one study and down in another"),
]


# Available pooling methods. (key, display, tooltip)
_PATHWAY_METHODS = [
    ('dl', 'DL', 'DerSimonian-Laird random effects (effect-size pooling)'),
    ('hksj', 'HKSJ', 'DL with Hartung-Knapp t-CI (more conservative for small k)'),
    ('reml', 'REML', 'Restricted Maximum Likelihood random effects'),
    ('stouffer', "Stouffer's Z",
     "Precision-weighted combination of per-study p-values"),
    ('sign', 'Sign Test',
     "Counts studies showing same effect direction (binomial null)"),
]


# ---------------------------------------------------------------------------
# Pooling helpers (module-level)
# ---------------------------------------------------------------------------


def _pool_advanced(base_df, adv_method):
    """
    Run DL+HKSJ, REML, or Stouffer across all pathways at once.

    Returns a Series of p-values indexed by 'base_df.index' (pathway name).
    Falls back to the base pooled p-value if a pathway has fewer than 2
    valid per-study effects.
    """
    out = pd.Series(index=base_df.index, dtype=float)

    # Each advanced method goes through a single vectorised pool_fn call
    # across all pathways. Same input shape, different pool function.
    if adv_method == 'hartung_knapp':
        from kosmic.meta_analysis.pooling.dl import dl_fast
        def pool_fn(per_study):
            return dl_fast(per_study, min_studies=MIN_STUDIES, hksj=True)
    elif adv_method == 'reml':
        from kosmic.meta_analysis.pooling.reml import reml_fast
        def pool_fn(per_study):
            return reml_fast(per_study, min_studies=MIN_STUDIES)
    elif adv_method == 'stouffer':
        from kosmic.meta_analysis.pooling.stouffer import stouffer_fast
        def pool_fn(per_study):
            return stouffer_fast(per_study, min_studies=MIN_STUDIES, weighted=True)
    else:
        return base_df['pvals_pooled'].copy()

    per_study = _build_per_study_dfs(base_df)
    if not per_study:
        return base_df['pvals_pooled'].copy()
    result = pool_fn(per_study)
    pval_map = dict(zip(result['names'], result['pvals_pooled']))
    for name in base_df.index:
        out[name] = float(np.clip(
            pval_map.get(name, base_df.loc[name, 'pvals_pooled']),
            1e-300, 1.0))
    return out


def _build_per_study_dfs(base_df):
    """
    Transpose base_df's per-pathway study_effects lists into per-study
    DataFrames suitable for dl_fast / reml_fast / stouffer_fast. One
    DataFrame per dataset, each holding a (pathway → logfc, se, pval) row
    for every pathway that dataset contributed to. DL/REML ignore pvals
    at proportion_top=1.0; Stouffer uses them.
    """
    study_records = {}  # dataset -> list of {names, logfc, se, pval}
    for name, row in base_df.iterrows():
        effects_data = row.get('study_effects', [])
        if not isinstance(effects_data, list):
            continue
        for d in effects_data:
            logfc = d.get('logfc', np.nan)
            se = d.get('se', np.nan)
            if not (np.isfinite(logfc) and np.isfinite(se) and se > 0):
                continue
            pval = d.get('pval', np.nan)
            if not np.isfinite(pval) or pval <= 0 or pval >= 1:
                pval = 0.5  # neutral fallback for downstream pool fns
            ds = d.get('dataset', 'unknown')
            study_records.setdefault(ds, []).append({
                'names': name, 'logfc': logfc, 'se': se, 'pval': pval,
            })

    per_study = []
    for ds, records in study_records.items():
        per_study.append(pd.DataFrame({
            'names': [r['names'] for r in records],
            'logfoldchanges': [r['logfc'] for r in records],
            'se': [r['se'] for r in records],
            'pvals': [r['pval'] for r in records],
            'dataset': [ds] * len(records),
        }))
    return per_study


class PathwayMAWorker(BaseWorker):
    """
    Run pathway-level pooling off the UI thread.

    Emits finished_ok with a tuple '(meta_df, method_keys)'.
    """

    def __init__(self, datasets, gene_datasets, project_folder,
                 method_keys, min_studies,
                 scoring_method='cpm_welch', rho=PATHWAY_VIF_RHO,
                 pathway_gene_sets=None, parent=None):
        super().__init__(parent)
        self.datasets = datasets
        self.gene_datasets = gene_datasets
        self.project_folder = project_folder
        self.method_keys = method_keys
        self.min_studies = min_studies
        self.scoring_method = scoring_method  # 'cpm_welch' or 'deseq2_vif'
        self.rho = rho
        self.pathway_gene_sets = pathway_gene_sets or {}

    def _build_deseq2_vif_datasets(self):
        """
        Build per-study pathway-level DataFrames from gene-level
        DESeq2 results with CAMERA-style VIF adjustment.

        For each study and each pathway:
          logFC = IV-weighted mean of per-gene logFCs
          SE = sqrt(sum(w_i^2 * se_i^2) / (sum(w_i))^2) * sqrt(VIF)
          VIF = 1 + (m - 1) * rho
        where w_i = 1/se_i^2 and rho is the inter-gene correlation.
        """
        if not self.gene_datasets or not self.pathway_gene_sets:
            return None

        rho = self.rho
        results_per_study = []

        for ds in self.gene_datasets:
            df = ds if isinstance(ds, pd.DataFrame) else ds.get('df')
            if df is None:
                continue

            gene_col = 'names' if 'names' in df.columns else df.columns[0]
            if 'logfoldchanges' not in df.columns or 'se' not in df.columns:
                continue

            gene_df = df.set_index(gene_col)
            rows = []

            for pw_name, pw_genes in self.pathway_gene_sets.items():
                # Find genes present in this study
                available = [g for g in pw_genes
                             if g.upper() in gene_df.index or g in gene_df.index]
                if len(available) < 2:
                    continue

                lfcs, ses = [], []
                for g in available:
                    key = g.upper() if g.upper() in gene_df.index else g
                    if key not in gene_df.index:
                        continue
                    r = gene_df.loc[key]
                    if isinstance(r, pd.DataFrame):
                        r = r.iloc[0]
                    lfc = r.get('logfoldchanges', np.nan)
                    se = r.get('se', np.nan)
                    if (np.isfinite(lfc) and np.isfinite(se)
                            and se > 0):
                        lfcs.append(float(lfc))
                        ses.append(float(se))

                m = len(lfcs)
                if m < 2:
                    continue

                lfcs = np.array(lfcs)
                ses = np.array(ses)
                w = 1.0 / (ses ** 2)
                w_sum = w.sum()

                # IV-weighted mean logFC
                pw_lfc = float((w * lfcs).sum() / w_sum)

                # SE with VIF adjustment
                vif = 1.0 + (m - 1) * rho
                pw_se_raw = float(np.sqrt((w ** 2 * ses ** 2).sum())
                                  / w_sum)
                pw_se = pw_se_raw * np.sqrt(vif)

                rows.append({
                    'names': pw_name,
                    'logfoldchanges': pw_lfc,
                    'se': pw_se,
                    'n_genes': m,
                })

            if rows:
                # So the direction columns can be counted.
                results_per_study.append(add_wald_pvalues(pd.DataFrame(rows)))

        return results_per_study if results_per_study else None

    def _run(self):
        from kosmic.numerical import bh_fdr
        from kosmic.meta_analysis.pooling_dispatch import get_analytical_pool_fn

        # Step 1: Build per-study pathway DataFrames
        if self.scoring_method == 'deseq2_vif':
            self.progress.emit(
                "Building pathway summaries from DESeq2 "
                f"(VIF, rho={self.rho})...")
            vif_dfs = self._build_deseq2_vif_datasets()
            if vif_dfs:
                dfs = vif_dfs
            else:
                self.progress.emit(
                    "DESeq2+VIF failed; falling back to CPM+Welch")
                dfs = [d['df'] for d in self.datasets]
        else:
            dfs = [d['df'] for d in self.datasets]

        from kosmic.meta_analysis.direction import unadjusted_warning
        warning = unadjusted_warning(dfs)
        if warning:
            self.progress.emit(f"Warning: {warning}")

        # Step 2: Pool across studies using the same vectorised
        # code as gene-level MA (no reimplementation).
        method_key = self.method_keys[0] if self.method_keys else 'dl'
        self.progress.emit(f"Pooling ({method_key})...")
        pool_fn = get_analytical_pool_fn(method_key, self.min_studies)
        result = pool_fn(dfs, min_studies=self.min_studies)

        if result.empty:
            raise RuntimeError(
                f"No pathways in >= {self.min_studies} studies.")

        # Step 3: BH FDR
        result['fdr'] = bh_fdr(result['pvals_pooled'])

        result = result.sort_values('pvals_pooled').reset_index(
            drop=True)
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
        self._meta_df = None  # pooled pathway-level DataFrame from the most recent run
        self._method_keys = []  # active method keys for the current run
        self._scatter_data = None
        self._forest_xlim = None
        self._method_checkboxes = {}
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

        settings_group = SettingsGroup("Analysis Settings", collapsible=True, expanded=True)
        self._min_studies = QSpinBox()
        self._min_studies.setRange(2, 20)
        self._min_studies.setValue(2)
        settings_group.add_row("Min studies:", self._min_studies)

        self._scoring_combo = QComboBox()
        self._scoring_combo.addItems([
            "DESeq2 + VIF (recommended)",
            "CPM + Welch (legacy)"])
        self._scoring_combo.setToolTip(
            "DESeq2 + VIF: pathway summaries from per-gene DESeq2 "
            "results with CAMERA-style VIF correction.\n\n"
            "CPM + Welch: mean log2-CPM, Welch t-test. Legacy.")
        settings_group.add_row("Scoring:", self._scoring_combo)

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
                     project_folder=None, pathway_gene_sets=None):
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
            lines = [f"{n} studies with pathway DE:"]
            for d in datasets:
                lines.append(f"  {d['name']} ({len(d['df'])} pathways)")
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

        if len(self._datasets) < min_studies:
            self._status_label.setText(
                f"Need >= {min_studies} studies, have "
                f"{len(self._datasets)}.")
            return

        scoring = ('deseq2_vif'
                   if self._scoring_combo.currentIndex() == 0
                   else 'cpm_welch')
        pw_gene_sets = self._pathway_gene_sets

        # Snapshot settings now so the sidecar JSON written on success
        # reflects exactly what the user picked, not whatever the
        # widgets show by the time the worker returns.
        self._last_run_settings = {
            'pooling_method':  pooling_key,
            'min_studies':     min_studies,
            'scoring_method':  scoring,
            'vif_rho':         PATHWAY_VIF_RHO if scoring == 'deseq2_vif' else None,
            'n_studies_input': len(self._datasets),
            'study_names':     [d.get('name', '?') for d in self._datasets],
            'n_pathways_in_geneset': len(pw_gene_sets) if pw_gene_sets else 0,
        }

        # Clear stale per-pathway forest from any previous run before
        # the new pooling kicks off.
        self._forest.clear_plot_items()
        self._forest.show_unavailable_message()

        self._run_btn.setEnabled(False)
        self._status_label.setText("Running pathway meta-analysis...")
        self.log_message.emit(
            f"Pathway MA: {len(self._datasets)} studies, "
            f"pooling={pooling_key}, scoring={scoring}")

        self._worker = PathwayMAWorker(
            self._datasets, self._gene_datasets, self._project_folder,
            method_keys, min_studies,
            scoring_method=scoring, rho=PATHWAY_VIF_RHO,
            pathway_gene_sets=pw_gene_sets,
        )
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
                'n_studies': len(self._datasets or []),
                'studies': [d.get('name', '?') for d in (self._datasets or [])],
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
            record_meta_stage(self._project_folder, stage, params)
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
        """Save the pooled DataFrame + settings sidecar to {project}/meta_analysis/.

        Both files are keyed by the pooling method, so re-running with
        the same method overwrites the pair (settings stay in sync with
        the CSV). To preserve a previous run, use the manual Save CSV
        button to copy the file elsewhere.
        """
        if not self._project_folder or meta_df is None or meta_df.empty:
            return
        try:
            output_dir = meta_output_dir(self._project_folder)
            output_dir.mkdir(parents=True, exist_ok=True)
            method = method_keys[0] if method_keys else 'pool'
            stem = f'pathway_ma_{method}'

            csv_path = output_dir / f'{stem}.csv'
            meta_df.to_csv(csv_path, index=False)

            json_path = output_dir / f'{stem}.json'
            self._write_settings_sidecar(json_path, meta_df, method_keys)

            self.log_message.emit(
                f"Pathway MA results saved: {csv_path.name} + {json_path.name}"
            )
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
