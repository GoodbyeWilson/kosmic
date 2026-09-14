# Individual Methods Page (Methods Comparison, step 1 of 2).
#
# Evaluates each selected pooling method standalone (no intersection).
# Produces a per-method summary table (n_called at FDR<0.05,
# n_meta_discovered, n_obvious, direction stats) plus a per-method
# Results table with an All / Obvious / Meta-discovered partition filter.
#
# Reuses the consensus pipeline's shared-DE infrastructure but skips
# the merge step so per-method results stay separate. Results are
# cached on the page for the downstream Consensus Evaluation page.

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor, QKeySequence, QShortcut
from PyQt6.QtWidgets import (
    QAbstractItemView, QApplication, QButtonGroup, QCheckBox, QComboBox, QHBoxLayout, QHeaderView, QLabel, QRadioButton, QSpinBox, QTabWidget, QTableWidget,
    QVBoxLayout, QWidget,
)

from kosmic.gui.shared.theme import get_color
from kosmic.gui.shared.widgets import BaseWorker, Column, ResultsTable, SettingsGroup, SidebarPage, SecondaryLabel, PrimaryButton


# Reuse the single source of truth for method labels.
from kosmic.gui.meta_analysis.pages.gene_ma_page import _METHODS as METHODS_REGISTRY
from kosmic import DEFAULT_FDR
from kosmic.gui.shared import dialogs, run_worker


def _copy_table_to_clipboard(table: QTableWidget, select_all_if_empty=True):
    """
    Copy selected cells from a QTableWidget to the clipboard as TSV.

    If nothing is selected and 'select_all_if_empty' is True, the whole
    table is copied. Headers are always included.
    """
    ranges = table.selectedRanges()
    if not ranges and select_all_if_empty:
        rows = list(range(table.rowCount()))
        cols = list(range(table.columnCount()))
    elif not ranges:
        return
    else:
        row_set = set()
        col_set = set()
        for r in ranges:
            for i in range(r.topRow(), r.bottomRow() + 1):
                row_set.add(i)
            for j in range(r.leftColumn(), r.rightColumn() + 1):
                col_set.add(j)
        rows = sorted(row_set)
        cols = sorted(col_set)

    if not rows or not cols:
        return

    headers = []
    for j in cols:
        h = table.horizontalHeaderItem(j)
        headers.append(h.text() if h is not None else '')
    lines = ['\t'.join(headers)]
    for i in rows:
        cells = []
        for j in cols:
            it = table.item(i, j)
            cells.append(it.text() if it is not None else '')
        lines.append('\t'.join(cells))

    QApplication.clipboard().setText('\n'.join(lines))


def _wire_table_copy(table: QTableWidget):
    """Enable Ctrl+C copy (TSV) and Ctrl+A select-all on a QTableWidget."""
    table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
    copy_sc = QShortcut(QKeySequence.StandardKey.Copy, table)
    copy_sc.setContext(Qt.ShortcutContext.WidgetShortcut)
    copy_sc.activated.connect(
        lambda t=table: _copy_table_to_clipboard(t, select_all_if_empty=False))
    select_all_sc = QShortcut(QKeySequence.StandardKey.SelectAll, table)
    select_all_sc.setContext(Qt.ShortcutContext.WidgetShortcut)
    select_all_sc.activated.connect(table.selectAll)


# Worker: per-method pooling (no merge)

class IndividualMethodsWorker(BaseWorker):
    """Run N pooling methods with shared CC permutation. No intersection."""

    def __init__(self, datasets, labels, params, parent=None):
        super().__init__(parent)
        self.datasets = datasets
        self.labels = labels
        self.params = params

    def _run(self):
        from kosmic.meta_analysis.pooling_dispatch import (
            get_analytical_pool_fn)

        method_keys = self.params['method_keys']
        calibration = self.params.get('calibration', 'analytical')
        min_st = self.params.get('min_studies', 3)
        n_cc_perms = self.params.get('cc_n_perms', 1000)
        n_methods = len(method_keys)

        # Step 1: Analytical pooling per method
        analytical_dfs = {}
        for i, key in enumerate(method_keys):
            self.progress.emit(
                f"Analytical pooling {i+1}/{n_methods}: {key}...")
            pool_fn = get_analytical_pool_fn(key, min_st)
            meta_df = pool_fn(self.datasets, min_studies=min_st)
            analytical_dfs[key] = meta_df

        # Step 2: CC perm calibration (optional)
        if calibration == 'cc_perm':
            pb_paths = self.params.get('pseudobulk_paths', [])
            self.progress.emit(
                f"CC perm: {len(pb_paths)} pseudobulk path(s) received")
            if not pb_paths:
                raise RuntimeError(
                    "No pseudobulk data. Run DE with 'Save pseudobulk' "
                    "first, or switch to Analytical calibration.")

            pb_data = self._load_pseudobulk(pb_paths)
            self.progress.emit(
                f"CC perm: {len(pb_data)}/{len(pb_paths)} pseudobulk "
                f"file(s) loaded")
            if not pb_data:
                raise RuntimeError("Failed to load pseudobulk data.")

            from kosmic.meta_analysis.cc_permutation import (
                consensus_cc_permutation)

            cc_de_method = self.params.get('de_method', 'deseq2_fast')
            self.progress.emit(
                f"CC permutation: {n_cc_perms} perms, {cc_de_method}, "
                f"{len(pb_data)} studies...")
            method_dfs = consensus_cc_permutation(
                pb_data, analytical_dfs,
                n_perms=n_cc_perms, seed=42,
                de_method=cc_de_method,
                min_studies=min_st,
                progress_cb=lambda cur, tot: self.progress.emit(
                    f"CC permutation {cur}/{tot}"),
                status_cb=lambda msg: self.progress.emit(msg),
            )
            # Post-CC-perm sanity check: if no method got a
            # pvals_analytical column it means calibration didn't
            # run (empty null distributions). Tell the user.
            calibrated_any = any(
                'pvals_analytical' in df.columns
                for df in method_dfs.values()
                if df is not None and len(df) > 0)
            if not calibrated_any:
                raise RuntimeError(
                    "CC permutation completed but produced no "
                    "calibrated p-values (empty null distribution for "
                    "every method). Likely causes: per-permutation DE "
                    "failed silently, or disease/control labels don't "
                    "match the pseudobulk condition values. Check the "
                    "log for DE errors.")
            self.progress.emit(
                f"CC permutation complete: calibrated p-values on "
                f"{sum(1 for df in method_dfs.values() if df is not None and 'pvals_analytical' in df.columns)} "
                f"/{len(method_dfs)} methods")
        else:
            method_dfs = analytical_dfs

        # Step 3: Per-method FDR + contributing-study counts.
        # DO NOT merge. Keep per-method outputs separate.
        from kosmic.numerical import bh_fdr

        # Count per-study FDR-sig calls per gene (used for the
        # obvious-vs-meta-discovered partition). Uses whichever
        # gene-name column the raw DE input has.
        n_sig_per_gene = {}
        for ds in self.datasets:
            if 'names' in ds.columns:
                gene_col = 'names'
            elif 'gene' in ds.columns:
                gene_col = 'gene'
            elif 'gene_name' in ds.columns:
                gene_col = 'gene_name'
            else:
                self.progress.emit(
                    f"Skipping dataset: no gene-name column "
                    f"(cols: {list(ds.columns)[:8]})")
                continue
            pcol = 'pvals_adj' if 'pvals_adj' in ds.columns else (
                'fdr' if 'fdr' in ds.columns else None)
            if pcol is None:
                continue
            sig_names = ds.loc[ds[pcol] < DEFAULT_FDR, gene_col]
            for g in sig_names:
                n_sig_per_gene[g] = n_sig_per_gene.get(g, 0) + 1

        # Build per-gene directional-agreement lookup: for each gene,
        # (n_agreeing_studies, n_valid_studies) given the gene's pooled
        # sign. We compute this OFFLINE from raw per-study logFCs,
        # keyed by gene name, then look up per method because pooled
        # sign is always the sign of the method's pooled logFC.
        # Stored as {gene: [+/-1 signs across studies]} then per-
        # method agreement is derived.
        per_study_signs: Dict[str, List[int]] = {}
        for ds in self.datasets:
            gene_col = None
            if 'names' in ds.columns:
                gene_col = 'names'
            elif 'gene' in ds.columns:
                gene_col = 'gene'
            elif 'gene_name' in ds.columns:
                gene_col = 'gene_name'
            if gene_col is None or 'logfoldchanges' not in ds.columns:
                continue
            for _, row in ds[[gene_col, 'logfoldchanges']].iterrows():
                g = row[gene_col]
                lfc = row['logfoldchanges']
                if not np.isfinite(lfc) or lfc == 0:
                    continue
                per_study_signs.setdefault(g, []).append(
                    1 if lfc > 0 else -1)

        # Add 'fdr', 'n_contributing_sig', 'n_agreeing', 'n_valid',
        # 'frac_agreeing' to each per-method df.
        for key, df in method_dfs.items():
            if df is None or len(df) == 0:
                continue
            if 'pvals_pooled' not in df.columns:
                self.progress.emit(
                    f"Method {key} returned no 'pvals_pooled' column "
                    f"(cols: {list(df.columns)[:8]}) -- skipping FDR")
                continue
            if 'names' not in df.columns:
                self.progress.emit(
                    f"Method {key} returned no 'names' column "
                    f"(cols: {list(df.columns)[:8]}) -- skipping partition")
                continue
            df['fdr'] = bh_fdr(df['pvals_pooled'])
            df['n_contributing_sig'] = df['names'].map(
                lambda g: n_sig_per_gene.get(g, 0)).astype(int)

            # Directional concordance per gene: count valid studies
            # whose per-study logFC sign matches the pooled sign.
            if 'logfoldchanges' in df.columns:
                pooled = df['logfoldchanges'].to_numpy(dtype=float)
                names = df['names'].tolist()
                n_agreeing = np.zeros(len(df), dtype=int)
                n_valid = np.zeros(len(df), dtype=int)
                for i, g in enumerate(names):
                    signs = per_study_signs.get(g, [])
                    p_sign = 1 if pooled[i] > 0 else (
                        -1 if pooled[i] < 0 else 0)
                    n_valid[i] = len(signs)
                    if p_sign != 0 and signs:
                        n_agreeing[i] = sum(1 for s in signs if s == p_sign)
                df['n_agreeing'] = n_agreeing
                df['n_valid'] = n_valid
                with np.errstate(divide='ignore', invalid='ignore'):
                    frac = np.where(n_valid > 0,
                                    n_agreeing / n_valid.astype(float),
                                    np.nan)
                df['frac_agreeing'] = frac

        return {
            'method_dfs': method_dfs,
            'method_keys': method_keys,
            'n_sig_per_gene': n_sig_per_gene,
            'params': self.params,
        }

    def _load_pseudobulk(self, pb_paths):
        """
        Skip per-path failures so a single unloadable study does not
        kill the whole individual-methods run.
        """
        from kosmic.meta_analysis.io import load_pseudobulk_set
        return load_pseudobulk_set(
            pb_paths, labels=self.labels,
            progress_cb=self.progress.emit, skip_errors=True)


# Page UI

class MethodsComparisonPage(QWidget):
    """Methods Comparison -- step 1: evaluate each method standalone."""

    help_id = "meta/methods_comparison"

    log_message = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._datasets = None
        self._dataset_labels = None
        self._project_folder = None
        self._method_results = None  # dict set on Discovery Run finish
        self._loo_results = None     # dict set on Cross-Study Validation finish
        self._worker = None
        self._pseudobulk_paths: List[str] = []  # provided by workspace
        self.progress_bar = None
        self._setup_ui()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _setup_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(12, 12, 12, 12)
        outer.setSpacing(10)

        # Header
        title = QLabel("Individual Methods")
        title.setProperty("role", "page_title")
        outer.addWidget(title)

        subtitle = SecondaryLabel(
            "Evaluate each selected pooling method standalone (no "
            "intersection). Per-method outputs let you decide which methods "
            "are worth trusting and for which purpose. Results feed into "
            "the next step (Consensus Evaluation) where combinations can be "
            "built interactively on cached data.")
        subtitle.setWordWrap(True)
        outer.addWidget(subtitle)

        # Sidebar + tabs body (composition: page hosts a SidebarPage).
        self._body = SidebarPage()
        outer.addWidget(self._body, 1)
        ctrl_lay = self._body.sidebar_layout
        split = self._body.splitter

        methods_box = SettingsGroup("Methods", collapsible=True, expanded=True)
        self._method_checks: Dict[str, QCheckBox] = {}
        default_on = {'dl', 'reml', 'stouffer', 'fisher',
                      'sumrank', 'gwop'}
        for key, label in METHODS_REGISTRY:
            cb = QCheckBox(label)
            cb.setChecked(key in default_on)
            self._method_checks[key] = cb
            methods_box.add_widget(cb)
        ctrl_lay.addWidget(methods_box)

        settings_box = SettingsGroup("Settings", collapsible=True, expanded=True)
        self._cal_combo = QComboBox()
        self._cal_combo.addItems(["Analytical", "CC permutation"])
        settings_box.add_row("Calibration:", self._cal_combo)

        self._de_combo = QComboBox()
        self._de_combo.addItems(["DESeq2", "Welch CPM", "Welch Raw"])
        settings_box.add_row("DE method:", self._de_combo)

        self._min_studies = QSpinBox()
        self._min_studies.setRange(2, 20)
        self._min_studies.setValue(3)
        settings_box.add_row("Min studies:", self._min_studies)

        self._fdr_spin = QComboBox()
        self._fdr_spin.addItems(["0.01", "0.05", "0.1", "0.25"])
        self._fdr_spin.setCurrentText("0.05")
        self._fdr_spin.currentIndexChanged.connect(self._on_fdr_changed)
        settings_box.add_row("FDR threshold:", self._fdr_spin)

        self._cc_perms = QSpinBox()
        self._cc_perms.setRange(100, 10000)
        self._cc_perms.setValue(1000)
        self._cc_perms.setSingleStep(100)
        settings_box.add_row("# CC perms:", self._cc_perms)

        ctrl_lay.addWidget(settings_box)

        self._run_btn = PrimaryButton("Run Individual Methods")
        self._run_btn.setMinimumHeight(32)
        self._run_btn.clicked.connect(self._run)
        ctrl_lay.addWidget(self._run_btn)

        self._status = SecondaryLabel("")
        self._status.setWordWrap(True)
        ctrl_lay.addWidget(self._status)

        ctrl_lay.addStretch()

        # --- RIGHT: results (two tabs: Discovery + Cross-Study Validation) ---
        right_lay = self._body.content_layout
        right_lay.setSpacing(8)

        self._right_tabs = QTabWidget()
        right_lay.addWidget(self._right_tabs)

        # Tab 1: Discovery (per-method pooling output)
        discovery_tab = QWidget()
        disc_lay = QVBoxLayout(discovery_tab)
        disc_lay.setContentsMargins(8, 8, 8, 8)
        disc_lay.setSpacing(8)

        comp_box = SettingsGroup(
            "List composition (counts at FDR threshold + effect size)")
        self._summary_table = ResultsTable()
        self._summary_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch)
        self._summary_table.verticalHeader().setVisible(False)
        self._summary_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self._summary_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self._summary_table.set_schema([
            Column("Method",                  "method",       "s"),
            Column("n called",                "n_called",     "d"),
            Column("n obvious (n>=1)",        "n_obvious",    "d"),
            Column("n meta-discovered (n=0)", "n_meta",       "d",
                   "Meta-discovered = FDR-significant in this method but NOT "
                   "individually FDR-significant in any contributing study. "
                   "Requires effect-size pooling to surface."),
            Column("% meta",                  "pct_meta",     ".1%"),
            Column("mean |logFC|",            "mean_abs_lfc", ".3f",
                   "Mean absolute log2 fold change across the method's entire "
                   "FDR<0.05 call set. Smaller means the method reaches into "
                   "subthreshold effects; larger means it requires bigger signal."),
        ])
        comp_box.add_widget(self._summary_table)

        disc_lay.addWidget(comp_box, 1)

        concord_box = SettingsGroup(
            "Directional concordance (whole call set vs meta-discovered "
            "subset)")
        self._concord_table = ResultsTable()
        self._concord_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch)
        self._concord_table.verticalHeader().setVisible(False)
        self._concord_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self._concord_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self._concord_table.set_schema([
            Column("Method",                "method",              "s"),
            Column("mean concord % (all)",  "mean_concord",        ".1%",
                   "Mean directional concordance across the method's ENTIRE "
                   "FDR<0.05 call set: for each gene, fraction of valid studies "
                   "whose per-study logFC sign matches the pooled sign. "
                   "Averaged across genes."),
            Column("% unanimous (all)",     "pct_unanimous",       ".1%",
                   "Fraction of the method's ENTIRE call set with unanimous "
                   "directional agreement (all valid studies match pooled sign)."),
            Column("mean concord % (meta)", "meta_mean_concord",   ".1%", na_text="-",
                   tooltip=("Mean directional concordance computed ONLY on the "
                            "meta-discovered (n=0) subset. Tests whether the method's "
                            "meta-discoveries are directionally robust.")),
            Column("% unanimous (meta)",    "meta_pct_unanimous",  ".1%", na_text="-",
                   tooltip=("Fraction of meta-discoveries with unanimous directional "
                            "agreement. The strongest defensibility metric for the "
                            "paper's meta-discovery claim.")),
        ])
        concord_box.add_widget(self._concord_table)

        disc_lay.addWidget(concord_box, 1)

        self._right_tabs.addTab(discovery_tab, "Discovery")

        # Tab 2: Gene Lists (per-method drill-down with partition filter).
        # Moved out of Discovery because three stacked group-boxes was
        # too crowded for one tab.
        gene_tab = QWidget()
        gene_tab_lay = QVBoxLayout(gene_tab)
        gene_tab_lay.setContentsMargins(8, 8, 8, 8)
        gene_tab_lay.setSpacing(8)

        results_box = SettingsGroup("Per-method gene list")
        rb_lay = QVBoxLayout()
        rb_lay.setContentsMargins(0, 0, 0, 0)
        rb_lay.setSpacing(6)
        results_box.add_layout(rb_lay)

        # Method picker
        method_row = QHBoxLayout()
        method_row.addWidget(QLabel("Method:"))
        self._method_picker = QComboBox()
        self._method_picker.currentIndexChanged.connect(
            self._on_method_changed)
        method_row.addWidget(self._method_picker, 1)
        rb_lay.addLayout(method_row)

        # Per-method distribution label (shows counts at each n_sig bin)
        self._distribution_label = SecondaryLabel("")
        self._distribution_label.setWordWrap(True)
        rb_lay.addWidget(self._distribution_label)

        # Granular partition filter by n contributing studies FDR-sig.
        # Mirrors the LOO 'By n contributing studies FDR-sig' breakdown
        # so Individual Methods output is directly comparable.
        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("Show:"))
        self._partition_group = QButtonGroup(self)
        self._rb_all = QRadioButton("All")
        self._rb_n0 = QRadioButton("n=0 (meta-discovered)")
        self._rb_n1 = QRadioButton("n=1")
        self._rb_n2 = QRadioButton("n=2")
        self._rb_n3plus = QRadioButton("n>=3 (obvious)")
        self._rb_all.setChecked(True)
        self._partition_buttons = [
            self._rb_all, self._rb_n0, self._rb_n1,
            self._rb_n2, self._rb_n3plus,
        ]
        for rb in self._partition_buttons:
            self._partition_group.addButton(rb)
            filter_row.addWidget(rb)
            rb.toggled.connect(self._populate_results_table)
        rb_lay.addLayout(filter_row)

        self._results_table = ResultsTable()
        self._results_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self._results_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self._results_table.setSortingEnabled(True)
        _wire_table_copy(self._results_table)
        rb_lay.addWidget(self._results_table, 1)

        self._gene_count_label = SecondaryLabel("")
        rb_lay.addWidget(self._gene_count_label)

        gene_tab_lay.addWidget(results_box, 1)
        self._right_tabs.addTab(gene_tab, "Gene Lists")

        # Tab 3: Cross-Study Validation (per-method LOO + held-out AUC)
        validation_tab = self._build_validation_tab()
        self._right_tabs.addTab(validation_tab, "Cross-Study Validation")

        # methods_comparison wants a 1:3 split favoring content (results
        # tabs are dense). Override the SidebarPage default 0:1.
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 3)

    # ------------------------------------------------------------------
    # Validation tab: per-method LOO + AUC
    # ------------------------------------------------------------------

    def _build_validation_tab(self) -> QWidget:
        tab = QWidget()
        lay = QVBoxLayout(tab)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(10)

        info = SecondaryLabel(
            "Held-out K-fold LOO for each method. Per fold: build the "
            "method's K-1 call set at the chosen FDR, test whether each "
            "gene replicates in direction on the held-out study, and "
            "compute held-out case/control AUC from the top-1000 gene "
            "list (requires pseudobulk)."
        )
        info.setWordWrap(True)
        lay.addWidget(info)

        ctrl_row = QHBoxLayout()
        self._val_run_btn = PrimaryButton("Run LOO Validation")
        self._val_run_btn.setMinimumHeight(28)
        self._val_run_btn.setEnabled(False)
        self._val_run_btn.clicked.connect(self._run_validation)
        ctrl_row.addWidget(self._val_run_btn)
        ctrl_row.addStretch()
        lay.addLayout(ctrl_row)

        self._val_status = SecondaryLabel(
            "Run the Discovery step first (left panel) to select methods. "
            "Then click Run LOO Validation.")
        self._val_status.setWordWrap(True)
        lay.addWidget(self._val_status)

        # Table 1: Direction replication by partition of the method's
        # own FDR<0.05 call set. Each dir % is on that specific partition;
        # "overall" is on the full call set (all partitions combined).
        dir_box = SettingsGroup(
            "Direction replication by partition (method's own FDR<0.05 "
            "call set, per n contributing studies FDR-sig)")
        self._val_table = ResultsTable()
        self._val_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch)
        self._val_table.verticalHeader().setVisible(False)
        self._val_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self._val_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self._val_table.set_schema([
            Column("Method",           "method",      "s"),
            Column("n=0 count",        "n_meta",      ".0f",
                   "Average count of meta-discovered genes (FDR-sig in 0 "
                   "contributing studies) across K LOO folds."),
            Column("n=0 dir % (meta)", "dir_meta",    ".1%", na_text="-",
                   tooltip=("Direction-replication rate on the n=0 (meta-discovered) "
                            "subset only. The paper's key validation metric.")),
            Column("n=1 dir %",        "dir_1",       ".1%", na_text="-",
                   tooltip=("Direction-replication rate on genes called by exactly 1 "
                            "contributing study (independent of the n=0 result).")),
            Column("n=2 dir %",        "dir_2",       ".1%", na_text="-"),
            Column("n>=3 dir %",       "dir_3plus",   ".1%", na_text="-"),
            Column("overall dir %",    "overall_dir", ".1%", na_text="-",
                   tooltip=("Direction-replication rate across the method's entire "
                            "FDR<0.05 call set (all partitions combined).")),
        ])
        dir_box.add_widget(self._val_table)

        lay.addWidget(dir_box, 1)

        # Table 2: Fixed-size top-1000 comparison. Completely different
        # gene list from the solo FDR<0.05 call set above -- this is
        # each method's top 1000 genes by p-value, so gene-list size is
        # constant across methods and AUC differences reflect method
        # quality, not list size.
        auc_box = SettingsGroup(
            "Top-1000 fixed-size comparison (held-out case/control "
            "prediction + wrong-direction rate)")
        self._auc_table = ResultsTable()
        self._auc_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch)
        self._auc_table.verticalHeader().setVisible(False)
        self._auc_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self._auc_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self._auc_table.set_schema([
            Column("Method",                "method", "s"),
            Column("AUC (sum-of-signed-z)", "auc_z",  ".3f", na_text="-",
                   tooltip=("Held-out case/control AUC using sum-of-signed-z module "
                            "score on each method's top-1000 genes (by its own p-value). "
                            "Averaged across K folds. Requires pseudobulk. Higher is "
                            "better; 0.5 = null.")),
            Column("AUC (UCell-rank)",      "auc_u",  ".3f", na_text="-",
                   tooltip=("Same as AUC z but using UCell-style rank-based module "
                            "score -- strict parity with the SumRank paper's validation.")),
            Column("wrong-direction %",     "wd_pct", ".1%", na_text="-",
                   tooltip=("Percentage of the top-1000 genes whose direction in the "
                            "held-out study opposes the K-1 pooled direction.")),
        ])
        auc_box.add_widget(self._auc_table)

        lay.addWidget(auc_box, 1)

        # Paper-ready headline summary
        self._val_headline = QLabel("")
        self._val_headline.setWordWrap(True)
        self._val_headline.setProperty("role", "validation_headline")
        lay.addWidget(self._val_headline)

        return tab

    # ------------------------------------------------------------------
    # External activation hooks
    # ------------------------------------------------------------------

    def set_datasets(self, datasets, labels, project_folder):
        self._datasets = datasets
        self._dataset_labels = labels
        self._project_folder = project_folder

    def set_pseudobulk_paths(self, paths):
        """
        Receive pseudobulk file paths discovered by the workspace
        (delegated to GeneMAPage._find_pseudobulk_paths to avoid
        duplicating the multi-location discovery + caching logic).
        """
        self._pseudobulk_paths = list(paths) if paths else []

    def get_method_results(self) -> Optional[dict]:
        return self._method_results

    def get_loo_results(self) -> Optional[dict]:
        """Per-method LOO results from Cross-Study Validation, if run."""
        return self._loo_results

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def _on_fdr_changed(self):
        if self._method_results is not None:
            self._populate_summary_table()
            self._update_distribution_label()
            self._populate_results_table()

    def _run(self):
        if not self._datasets:
            dialogs.warning(self, "No data",
                                "Import studies first.")
            return
        selected = [k for k, cb in self._method_checks.items() if cb.isChecked()]
        if not selected:
            dialogs.warning(self, "No methods",
                                "Select at least one method.")
            return

        cal = ['analytical', 'cc_perm'][self._cal_combo.currentIndex()]
        de_methods = ['deseq2', 'welch_cpm', 'welch_raw']
        de = de_methods[self._de_combo.currentIndex()]

        params = {
            'method_keys': selected,
            'calibration': cal,
            'min_studies': self._min_studies.value(),
            'cc_n_perms': self._cc_perms.value(),
            'de_method': de,
        }

        # Pseudobulk paths if CC perm. Prefer workspace-provided paths
        # (these come from GeneMAPage._find_pseudobulk_paths which
        # does the heavy lifting: multi-location candidates, SMB re-auth,
        # persistent cache). Fall back to our simple glob only if the
        # workspace didn't provide any.
        if cal == 'cc_perm':
            pb_paths = self._pseudobulk_paths or self._discover_pseudobulk_paths()
            if not pb_paths:
                dialogs.warning(
                    self, "No pseudobulk",
                    "CC permutation requires saved pseudobulk files. "
                    "Re-run per-study DE with 'Save pseudobulk' enabled, "
                    "or switch to Analytical calibration.\n\nTip: try "
                    "running CC permutation on the Consensus page first -- "
                    "that warms the pseudobulk path cache, and the paths "
                    "are then shared with this page.")
                return
            params['pseudobulk_paths'] = pb_paths

        # Unpack datasets: the workspace passes a list of dicts
        # ({'path', 'name', 'df'}); workers expect a list of DataFrames
        # with the usual 'names'/'logfoldchanges'/'se'/'pvals' columns.
        # Mirrors gene_ma_page.py's _run_consensus pattern.
        if self._datasets and isinstance(self._datasets[0], dict):
            de_dfs = [d['df'] for d in self._datasets]
            ds_labels = [d['name'] for d in self._datasets]
        else:
            de_dfs = list(self._datasets)
            ds_labels = list(self._dataset_labels or [])

        # Disable UI + clear displayed results so stale data doesn't
        # mask a failed run + kick off worker
        self._run_btn.setEnabled(False)
        self._status.setText(f"Running ({cal})...")
        self._method_results = None
        self._summary_table.setRowCount(0)
        self._concord_table.setRowCount(0)
        self._results_table.setRowCount(0)
        self._results_table.setColumnCount(0)
        self._method_picker.blockSignals(True)
        self._method_picker.clear()
        self._method_picker.blockSignals(False)
        self._gene_count_label.setText("")
        if self.progress_bar is not None:
            self.progress_bar.setRange(0, 0)

        self._worker = IndividualMethodsWorker(
            de_dfs, ds_labels, params, parent=self)
        run_worker(
            self._worker,
            on_finished=self._on_finished,
            on_failed=self._on_failed,
            on_progress=self._on_progress,
        )

    def _discover_pseudobulk_paths(self) -> List[str]:
        """Find per-study pseudobulk.pkl files under the project folder."""
        if not self._project_folder or not self._dataset_labels:
            return []
        project = Path(self._project_folder)
        paths = []
        for label in self._dataset_labels:
            # Match the naming convention used elsewhere (study_folder/
            # analysis/pseudobulk.pkl or similar). Try a few common patterns.
            candidates = list(project.glob(f"{label}/*pseudobulk*.pkl"))
            candidates += list(project.glob(f"{label}/**/pseudobulk*.pkl"))
            if candidates:
                paths.append(str(candidates[0]))
        return paths

    def _on_progress(self, msg):
        self._status.setText(msg)
        self.log_message.emit(f"  individual methods: {msg}")

    def _on_failed(self, msg):
        self._run_btn.setEnabled(True)
        if self.progress_bar is not None:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)
        self._status.setText(f"Failed: {msg.splitlines()[0]}")
        dialogs.warning(self, "Run failed", msg)

    def _on_finished(self, result):
        self._run_btn.setEnabled(True)
        if self.progress_bar is not None:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(100)
        self._method_results = result
        n_methods = len(result['method_keys'])
        self._status.setText(
            f"Done -- {n_methods} methods, "
            f"{result['params'].get('calibration', 'analytical')} calibration")

        # Populate method picker
        self._method_picker.blockSignals(True)
        self._method_picker.clear()
        key_to_label = dict(METHODS_REGISTRY)
        for k in result['method_keys']:
            self._method_picker.addItem(key_to_label.get(k, k), userData=k)
        self._method_picker.blockSignals(False)

        self._populate_summary_table()
        self._update_distribution_label()
        self._populate_results_table()

        # Validation tab now runnable: methods are known
        self._val_run_btn.setEnabled(True)
        self._val_status.setText(
            f"Ready to LOO-validate {n_methods} methods. "
            f"Click Run LOO Validation.")

    # ------------------------------------------------------------------
    # Validation-tab handlers
    # ------------------------------------------------------------------

    def _run_validation(self):
        if self._method_results is None:
            dialogs.warning(
                self, "No methods",
                "Run the Discovery step first to pick methods.")
            return
        if not self._datasets:
            dialogs.warning(self, "No data", "Import studies first.")
            return

        method_keys = self._method_results.get('method_keys', [])
        if not method_keys:
            return

        # Unpack datasets if list of dicts
        if self._datasets and isinstance(self._datasets[0], dict):
            de_dfs = [d['df'] for d in self._datasets]
            ds_labels = [d['name'] for d in self._datasets]
        else:
            de_dfs = list(self._datasets)
            ds_labels = list(self._dataset_labels or [])

        try:
            fdr_thr = float(self._fdr_spin.currentText())
        except ValueError:
            fdr_thr = DEFAULT_FDR

        params = {
            'calibration': 'analytical',
            'min_studies': self._min_studies.value(),
            'consensus_fdr_threshold': fdr_thr,
            'replication_threshold': fdr_thr,
            'top_n': 1000,
        }

        self._val_run_btn.setEnabled(False)
        self._val_status.setText("Running per-method LOO...")
        self._val_table.setRowCount(0)
        self._val_headline.setText("")
        if self.progress_bar is not None:
            self.progress_bar.setRange(0, 0)

        from kosmic.gui.meta_analysis.per_method_loo_worker import PerMethodLOOWorker
        self._val_worker = PerMethodLOOWorker(
            de_dfs, ds_labels, method_keys, params,
            pseudobulk_paths=self._pseudobulk_paths or None,
            parent=self)
        run_worker(
            self._val_worker,
            on_finished=self._on_val_finished,
            on_failed=self._on_val_failed,
            on_progress=self._on_val_progress,
        )

    def _on_val_progress(self, msg):
        self._val_status.setText(msg)
        self.log_message.emit(f"  val: {msg}")

    def _on_val_failed(self, msg):
        self._val_run_btn.setEnabled(True)
        if self.progress_bar is not None:
            self.progress_bar.setRange(0, 100)
        self._val_status.setText(f"Failed: {msg.splitlines()[0]}")
        dialogs.warning(self, "LOO validation failed", msg)

    def _on_val_finished(self, result):
        self._val_run_btn.setEnabled(True)
        if self.progress_bar is not None:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(100)
        # Cache so the Consensus Evaluation page can pull LOO metrics
        self._loo_results = result
        self._populate_val_table(result)
        auc_note = ("" if result.get('pseudobulk_available')
                    else "  (pseudobulk unavailable; AUC columns blank)")
        self._val_status.setText(
            f"Done -- {len(result['method_keys'])} methods LOO'd"
            + auc_note)

    def _populate_val_table(self, result):
        per_method = result['per_method']
        method_keys = result['method_keys']
        key_to_label = dict(METHODS_REGISTRY)

        # Percent columns stored in 0-1 range (widget's '.1%' adds the
        # trailing '%' and numeric sort works).  NaN → '-' via na_text.
        val_rows = []
        auc_rows = []
        paper_lines = []
        def _partition_n(p, n_sig):
            return float(p.get(n_sig, {}).get('mean_n', 0.0))

        def _partition_dir_pct(p, n_sig):
            return float(p.get(n_sig, {}).get(
                'mean_replication_pct_direction_only', 0.0)) / 100.0

        for key in method_keys:
            res = per_method[key]
            partitions = res.get('partition_by_n_contributing', {})

            n_meta = _partition_n(partitions, 0)
            dir_meta = (_partition_dir_pct(partitions, 0)
                        if n_meta > 0.5 else float('nan'))
            dir_1 = _partition_dir_pct(partitions, 1)
            dir_2 = _partition_dir_pct(partitions, 2)
            n_3plus = 0.0
            weighted_dir_pct = 0.0
            for n_sig, p in partitions.items():
                if n_sig >= 3:
                    n_i = float(p.get('mean_n', 0.0))
                    n_3plus += n_i
                    weighted_dir_pct += n_i * float(p.get(
                        'mean_replication_pct_direction_only', 0.0)) / 100.0
            dir_3plus = (weighted_dir_pct / n_3plus
                         if n_3plus > 0 else float('nan'))
            overall_dir = float(res.get(
                'avg_replication_pct_direction', 0.0)) / 100.0

            tn = res.get('top_n_per_method', {}).get(key, {})
            auc_z = tn.get('mean_auc_sum_of_z', float('nan'))
            auc_u = tn.get('mean_auc_ucell', float('nan'))
            wd_raw = tn.get('mean_wrong_dir_pct', float('nan'))
            wd_frac = float('nan') if wd_raw != wd_raw else wd_raw / 100.0

            val_rows.append({
                'method':      key_to_label.get(key, key),
                'n_meta':      n_meta,
                'dir_meta':    dir_meta,
                'dir_1':       dir_1,
                'dir_2':       dir_2,
                'dir_3plus':   dir_3plus,
                'overall_dir': overall_dir,
            })
            auc_rows.append({
                'method': key_to_label.get(key, key),
                'auc_z':  auc_z,
                'auc_u':  auc_u,
                'wd_pct': wd_frac,
            })

            if n_meta >= 0.5:
                label = key_to_label.get(key, key)
                auc_z_s = "-" if auc_z != auc_z else f"{auc_z:.3f}"
                auc_u_s = "-" if auc_u != auc_u else f"{auc_u:.3f}"
                paper_lines.append(
                    f"  {label}: {n_meta:.0f} meta-discoveries, "
                    f"{dir_meta * 100.0:.1f}% direction rep on that subset, "
                    f"top-1000 AUC z={auc_z_s} UCell={auc_u_s}")

        accent = QColor(get_color('accent_primary'))

        def _val_decorate(item, row, col_i):
            # Accent the n=0 dir % cell (col 2) when a real value is shown.
            if col_i == 2 and item.text() != "-":
                item.setForeground(accent)

        self._val_table.set_data(pd.DataFrame(val_rows), decorate=_val_decorate)
        self._auc_table.set_data(pd.DataFrame(auc_rows))

        if paper_lines:
            self._val_headline.setText(
                "Paper-ready per-method validation (meta-discovery "
                "direction rep is computed on the n=0 subset; AUCs are "
                "on each method's top-1000 for fair comparison):\n"
                + "\n".join(paper_lines))
        else:
            self._val_headline.setText(
                "No meta-discoveries from any method at this FDR.")

    # ------------------------------------------------------------------
    # Tables
    # ------------------------------------------------------------------

    def _current_fdr_threshold(self) -> float:
        try:
            return float(self._fdr_spin.currentText())
        except ValueError:
            return DEFAULT_FDR

    def _populate_summary_table(self):
        if self._method_results is None:
            return
        method_dfs = self._method_results['method_dfs']
        method_keys = self._method_results['method_keys']
        fdr_thr = self._current_fdr_threshold()
        key_to_label = dict(METHODS_REGISTRY)

        # Percent columns are stored in 0-1 range so the widget's '.1%'
        # format renders them with a trailing '%' and sorts numerically.
        rows = []
        for key in method_keys:
            df = method_dfs.get(key)
            stats = {
                'method':       key_to_label.get(key, key),
                'n_called':     0, 'n_obvious': 0, 'n_meta': 0,
                'pct_meta':          0.0,
                'mean_abs_lfc':      0.0,
                'mean_concord':      0.0, 'pct_unanimous': 0.0,
                'meta_mean_concord': float('nan'),
                'meta_pct_unanimous': float('nan'),
            }
            if df is not None and len(df) > 0:
                sig = df[df['fdr'] < fdr_thr]
                n_called = len(sig)
                n_obvious = int((sig['n_contributing_sig'] >= 1).sum())
                n_meta = int((sig['n_contributing_sig'] == 0).sum())
                stats['n_called'] = n_called
                stats['n_obvious'] = n_obvious
                stats['n_meta'] = n_meta
                stats['pct_meta'] = (n_meta / n_called) if n_called > 0 else 0.0
                if 'logfoldchanges' in sig.columns and n_called > 0:
                    stats['mean_abs_lfc'] = float(
                        np.abs(sig['logfoldchanges']).mean())
                if 'frac_agreeing' in sig.columns and n_called > 0:
                    frac = sig['frac_agreeing'].to_numpy(dtype=float)
                    finite = np.isfinite(frac)
                    if finite.sum() > 0:
                        stats['mean_concord'] = float(frac[finite].mean())
                        stats['pct_unanimous'] = (
                            float(np.isclose(frac[finite], 1.0).sum())
                            / n_called)

                meta_sig = sig[sig['n_contributing_sig'] == 0]
                if ('frac_agreeing' in meta_sig.columns
                        and len(meta_sig) > 0):
                    mfrac = meta_sig['frac_agreeing'].to_numpy(dtype=float)
                    mfinite = np.isfinite(mfrac)
                    if mfinite.sum() > 0:
                        stats['meta_mean_concord'] = float(mfrac[mfinite].mean())
                        stats['meta_pct_unanimous'] = (
                            float(np.isclose(mfrac[mfinite], 1.0).sum())
                            / len(meta_sig))
            rows.append(stats)

        summary_df = pd.DataFrame(rows)

        # Accent-colour the n_meta cell when > 0 (col index 3 in summary).
        accent = QColor(get_color('accent_primary'))

        def _summary_decorate(item, row, col_i):
            if col_i == 3 and int(row['n_meta']) > 0:
                item.setForeground(accent)

        self._summary_table.set_data(summary_df, decorate=_summary_decorate)
        self._concord_table.set_data(summary_df)

    def _on_method_changed(self):
        """
        Refresh the distribution label AND the Results table when the
        method dropdown changes.
        """
        self._update_distribution_label()
        self._populate_results_table()

    def _update_distribution_label(self):
        """
        Show the selected method's gene-count distribution across
        n_contributing_sig bins, mirroring the LOO 'By n contributing
        studies FDR-sig' headline.
        """
        if self._method_results is None:
            self._distribution_label.setText("")
            return
        idx = self._method_picker.currentIndex()
        if idx < 0:
            self._distribution_label.setText("")
            return
        key = self._method_picker.itemData(idx)
        df = self._method_results['method_dfs'].get(key)
        if df is None or len(df) == 0:
            self._distribution_label.setText("")
            return
        fdr_thr = self._current_fdr_threshold()
        sig = df[df['fdr'] < fdr_thr]
        if len(sig) == 0:
            self._distribution_label.setText(
                "No genes pass the FDR threshold.")
            return
        n_by = sig['n_contributing_sig'].value_counts().to_dict()
        n0 = int(n_by.get(0, 0))
        n1 = int(n_by.get(1, 0))
        n2 = int(n_by.get(2, 0))
        n3p = int(sum(v for k, v in n_by.items() if k >= 3))
        total = n0 + n1 + n2 + n3p
        self._distribution_label.setText(
            f"Distribution by n contributing studies FDR-sig: "
            f"n=0 {n0}  |  n=1 {n1}  |  n=2 {n2}  |  n>=3 {n3p}  "
            f"(total {total})")

    def _populate_results_table(self):
        if self._method_results is None:
            return
        method_dfs = self._method_results['method_dfs']

        # Which method is currently picked
        idx = self._method_picker.currentIndex()
        if idx < 0:
            return
        key = self._method_picker.itemData(idx)
        df = method_dfs.get(key)
        if df is None or len(df) == 0:
            self._results_table.setRowCount(0)
            self._results_table.setColumnCount(0)
            self._gene_count_label.setText("No results for this method.")
            return

        # Apply FDR threshold and partition filter
        fdr_thr = self._current_fdr_threshold()
        sig = df[df['fdr'] < fdr_thr].copy()
        if self._rb_n0.isChecked():
            sig = sig[sig['n_contributing_sig'] == 0]
            partition_label = "n=0 (meta-discovered)"
        elif self._rb_n1.isChecked():
            sig = sig[sig['n_contributing_sig'] == 1]
            partition_label = "n=1"
        elif self._rb_n2.isChecked():
            sig = sig[sig['n_contributing_sig'] == 2]
            partition_label = "n=2"
        elif self._rb_n3plus.isChecked():
            sig = sig[sig['n_contributing_sig'] >= 3]
            partition_label = "n>=3 (obvious)"
        else:
            partition_label = "all"

        sig = sig.sort_values('fdr').reset_index(drop=True)

        # Schema: order matches the original column order; include only
        # columns that actually exist in 'sig'.
        _COL_SPEC = [
            Column('Gene',              'names',              's'),
            Column('log2FC',            'logfoldchanges',     '.3f'),
            Column('FDR',               'fdr',                '.2e'),
            Column('pvals_pooled',      'pvals_pooled',       '.2e'),
            Column('n studies FDR-sig', 'n_contributing_sig', 'd'),
            Column('n agreeing',        'n_agreeing',         'd'),
            Column('n measuring',       'n_valid',            'd'),
            Column('concord',           'frac_agreeing',      '.0%'),
            Column('I²',                'heterogeneity_i2',   '.3f'),
            Column('n studies',         'n_studies',          'd'),
        ]
        cols = [c for c in _COL_SPEC if c.key in sig.columns]
        self._results_table.set_schema(cols)
        self._results_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents)
        self._results_table.set_data(sig)

        self._gene_count_label.setText(
            f"{len(sig)} genes ({partition_label}) at FDR < {fdr_thr}")
