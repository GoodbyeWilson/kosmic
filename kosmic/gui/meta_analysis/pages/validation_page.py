# Discovery step 5: how much should the pooled gene list be believed?
#
# Four questions, four tabs, deliberately under one heading because in
# practice they are asked together:
#
#   - Cross-Dataset Reproducibility: do the studies agree gene by gene?
#     Needs only the loaded per-study DE, so it runs before any pooling.
#   - LOO Validation: hold each study out, refit, does the call survive?
#     Reuses the last pooling run's exact settings so the held-out test
#     is apples-to-apples, and so needs a pooling run first.
#   - Olink Panels: can these proteins be measured on a routine assay?
#     A clinical-measurability flag, NOT evidence the protein changed.
#   - GWAS Overlap: not built yet.
#
# The first two moved here from GeneMAPage, where they were two buttons
# in a sidebar card and two tabs among nine.
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QAbstractItemView, QDialog, QLabel, QPushButton, QTabWidget,
    QVBoxLayout, QWidget,
)

from kosmic import DEFAULT_FDR
from kosmic.gui.shared import dialogs, run_worker
from kosmic.gui.shared.widgets import (
    Column, PrimaryButton, ResultsTable, SecondaryLabel, SectionHeader,
    SidebarTabbedPage, StageAccordion, StageSummaryCard,
)
from kosmic.gui.meta_analysis.pages.gene_ma.workers import (
    ConsensusLOOWorker,
    LOOBenchmarkWorker,
    ReproducibilityWorker,
)


class ValidationPage(SidebarTabbedPage):
    """Discovery step 5: internal stability plus external evidence."""

    help_id = "meta/validation"

    log_message = pyqtSignal(str)
    # GeneMAPage owns the per-mode result cache, so the LOO result goes
    # back to it rather than this page reaching into its state.
    loo_complete = pyqtSignal(dict)

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._datasets = []
        self._meta_df = None
        self._run_context = {}
        self._pw_worker = None
        self._loo_worker = None
        self._loo_benchmark_worker = None
        self._loo_result = None
        self._project_folder = None
        self._selection = None
        # Assigned by the workspace, as on the other pages.
        self.progress_bar = None
        self._setup_ui()

    # -- UI ------------------------------------------------------------

    def _setup_ui(self):
        lay = self.sidebar_layout
        lay.addWidget(SectionHeader("VALIDATION"))
        self._accordion = StageAccordion(lay)

        self._studies_card = StageSummaryCard("Studies", icon="dataset")
        self._studies_card.set_read_only(True)
        self._studies_card.set_summary(["No studies loaded yet", ""])
        lay.addWidget(self._studies_card)

        self._pooling_card = StageSummaryCard("Pooled result", icon="bar-chart-2")
        self._pooling_card.set_read_only(True)
        self._pooling_card.set_summary(["No pooling run yet", ""])
        lay.addWidget(self._pooling_card)

        lay.addStretch()
        self._accordion.finalize()

        # Actions sit below the cards, always reachable.
        self._repro_btn = PrimaryButton("Run Cross-Dataset Reproducibility")
        self._repro_btn.setToolTip(
            "Gene-by-gene agreement between the loaded studies. Needs no "
            "pooling run -- it reads the per-study DE directly.")
        self._repro_btn.clicked.connect(self._run_repro_analysis)
        lay.addWidget(self._repro_btn)

        self._loo_btn = PrimaryButton("Run LOO Validation")
        self._loo_btn.setToolTip(
            "Hold each study out in turn, refit, and check the call "
            "survives. Reuses the last pooling run's exact settings.")
        self._loo_btn.clicked.connect(self._run_loo_validation)
        lay.addWidget(self._loo_btn)

        self._status_label = SecondaryLabel("")
        self._status_label.setWordWrap(True)
        lay.addWidget(self._status_label)

        # -- Tabs --
        from kosmic.gui.meta_analysis.pages.gene_ma.reproducibility import (
            ReproducibilityTab)
        self._repro_tab = ReproducibilityTab()
        self.add_tab(self._repro_tab, "Cross-Dataset Reproducibility")

        from kosmic.gui.meta_analysis.pages.gene_ma.loo import LOOTab
        self._loo_tab = LOOTab(
            n_methods_getter=lambda: len(self._ctx('method_keys') or []),
            mode_getter=lambda: self._ctx('mode') or '',
        )
        self._loo_tab.characterise_clicked.connect(
            self._show_replicator_characterisation)
        self._loo_tab.benchmark_clicked.connect(self._run_loo_benchmark)
        self.add_tab(self._loo_tab, "LOO Validation")

        self._olink_tab = _OlinkTab()
        self.add_tab(self._olink_tab, "Olink Panels")

        self._gwas_tab = _ComingSoonTab(
            "GWAS Overlap",
            "Does a GWAS hit for this trait sit near the gene?\n\n"
            "The overlap calculation exists in 'dev/bench/gwas_overlap.py' "
            "but has not been brought into the app yet, so there is "
            "nothing to show here.")
        self.add_tab(self._gwas_tab, "GWAS Overlap")

        self._sync_buttons()

    def set_project_folder(self, folder, selection=None):
        """Where the run's provenance sidecar lives.

        'selection' is the output folder of the pooling run being
        validated (ADR-007).
        """
        self._project_folder = folder
        self._selection = selection

    def _ctx(self, key):
        """One read point for the pooling run's captured settings."""
        return self._run_context.get(key)

    # -- Inputs from the workspace -------------------------------------

    def set_studies(self, datasets):
        """Per-study DE, available as soon as studies are imported."""
        self._datasets = list(datasets or [])
        n = len(self._datasets)
        names = ", ".join(d['name'] for d in self._datasets[:3])
        if n > 3:
            names += f", +{n - 3} more"
        self._studies_card.set_summary(
            [f"{n} stud{'y' if n == 1 else 'ies'} loaded", names or ""]
            if n else ["No studies loaded yet", ""])
        self._sync_buttons()

    def set_run_context(self, context):
        """Take the last pooling run's results and exact settings."""
        self._run_context = dict(context or {})
        self._meta_df = self._run_context.get('meta_df')
        n_genes = 0 if self._meta_df is None else len(self._meta_df)
        methods = self._run_context.get('method_keys') or []
        self._pooling_card.set_summary(
            [f"{n_genes:,} genes pooled",
             f"method: {'+'.join(methods)}" if methods else ""]
            if n_genes else ["No pooling run yet", ""])
        self._olink_tab.populate(self._meta_df)
        self._sync_buttons()

    def _sync_buttons(self):
        """Say why a button is off rather than just grey it out."""
        n = len(self._datasets)
        self._repro_btn.setEnabled(n >= 2)

        # LOO folds over the studies the pooling run actually used, not
        # over everything loaded, so its guard reads the run context --
        # matching the check inside _run_loo_validation.
        has_run = bool(self._run_context.get('params'))
        n_folds = len(self._run_context.get('de_dfs') or [])
        self._loo_btn.setEnabled(has_run and n_folds >= 3)

        notes = []
        if n < 2:
            notes.append("Reproducibility needs at least 2 studies.")
        if not has_run:
            notes.append("LOO needs a pooling run first.")
        elif n_folds < 3:
            notes.append(
                f"LOO needs at least 3 pooled studies; the last run "
                f"used {n_folds}.")
        self._status_label.setText("  ".join(notes))


    def _run_repro_analysis(self):
        """Run reproducibility analysis on a worker thread (no consensus needed)."""
        if not self._datasets or len(self._datasets) < 2:
            dialogs.warning(
                self, "Too Few Studies",
                f"Reproducibility analysis needs at least 2 studies; "
                f"have {len(self._datasets) if self._datasets else 0}. "
                "Load datasets on the Import page first.")
            return

        de_dfs = [d['df'] for d in self._datasets]
        labels = [d['name'] for d in self._datasets]
        K = len(de_dfs)

        self._repro_btn.setEnabled(False)
        self._repro_tab.set_status(f"Computing reproducibility (K={K})...")
        self.log_message.emit(f"Reproducibility: K={K} studies")

        self._pw_worker = ReproducibilityWorker(de_dfs, labels)
        run_worker(
            self._pw_worker,
            on_finished=self._on_repro_finished,
            on_failed=self._on_repro_failed,
            on_progress=self._on_repro_progress,
        )

    def _on_repro_progress(self, msg):
        self._repro_tab.set_status(msg)
        self.log_message.emit(f"  repro: {msg}")

    def _on_repro_finished(self, result):
        self._repro_btn.setEnabled(True)
        K = result['n_studies']
        summary = result['summary']
        self._repro_tab.populate(result)
        self._repro_tab.set_status(
            f"Reproducibility complete (K={K}, "
            f"{summary['modes']['strict']['n_unique_DEGs']} "
            f"unique strict DEGs).")
        self.log_message.emit(
            f"Reproducibility done: strict mean orphan rate "
            f"{summary['modes']['strict']['mean_orphan_pct']:.1f}%")

        from kosmic.meta_analysis.io import record_meta_stage
        strict = summary['modes']['strict']
        record_meta_stage(self._project_folder, 'meta_reproducibility', {
            'n_studies': K,
            'studies': [d.get('name', '?') for d in (self._datasets or [])],
            'n_unique_strict_degs': strict['n_unique_DEGs'],
            'mean_orphan_pct': round(float(strict['mean_orphan_pct']), 2),
        }, selection=self._selection)

    def _on_repro_failed(self, msg):
        self._repro_btn.setEnabled(True)
        self._repro_tab.set_status(f"Reproducibility failed: {msg}")
        self.log_message.emit(f"Reproducibility failed: {msg}")
        dialogs.warning(self, "Reproducibility Failed", msg)

    # --- LOO Validation tab ---

    def _run_loo_validation(self):
        """Kick off the LOO worker, reusing the most recent consensus settings."""
        if self._ctx('params') is None:
            dialogs.warning(
                self, "No Consensus Run Yet",
                "Run a normal consensus analysis first. LOO validation "
                "reuses those exact settings (methods, calibration, "
                "min studies) to keep the held-out test apples-to-apples.")
            return

        de_dfs = self._ctx('de_dfs')
        labels = self._ctx('labels')
        if not de_dfs or len(de_dfs) < 3:
            dialogs.warning(
                self, "Too Few Studies",
                f"LOO validation requires at least 3 studies; "
                f"have {len(de_dfs) if de_dfs else 0}.")
            return

        params = dict(self._ctx('params'))
        if params.get('calibration') == 'cc_perm':
            est_min = len(de_dfs) * 4  # ~4 min/fold
            if not dialogs.confirm(self, "Confirm Long Run", f"LOO with CC permutation will run the consensus "
                f"pipeline {len(de_dfs)} times (one fold per study). "
                f"Estimated time: ~{est_min} minutes for K={len(de_dfs)}. "
                f"Continue?"):
                return

        self._loo_btn.setEnabled(False)
        self._loo_tab.clear()
        self._loo_tab.set_status("Starting LOO worker...")

        self.log_message.emit(
            f"LOO: starting K={len(de_dfs)} folds, "
            f"calibration={params.get('calibration', 'analytical')}")

        self._loo_worker = ConsensusLOOWorker(
            de_dfs, labels, params,
            replication_mode='direction_only',
            pb_paths=self._ctx('pb_paths') or [])
        self._loo_worker.fold_progress.connect(self._on_loo_fold_progress)
        run_worker(
            self._loo_worker,
            on_finished=self._on_loo_finished,
            on_failed=self._on_loo_failed,
            on_progress=self._on_loo_progress,
        )

    def _on_loo_progress(self, message):
        self._loo_tab.set_status(message)
        self.log_message.emit(f"  loo: {message}")

    def _on_loo_fold_progress(self, fold_one_based, n_folds, name):
        if self.progress_bar is not None:
            self.progress_bar.setRange(0, n_folds)
            self.progress_bar.setValue(fold_one_based)

    def _on_loo_finished(self, result):
        self._loo_btn.setEnabled(True)
        if self.progress_bar is not None:
            self.progress_bar.setRange(0, 0)
        self._loo_result = result

        self.loo_complete.emit(result)

        self._loo_tab.populate(result)
        self._loo_tab.set_status("LOO complete.")
        self.log_message.emit(
            f"LOO complete: avg replication "
            f"{result['avg_replication_pct']:.1f}% "
            f"(range {result['min_replication_pct']:.1f}-"
            f"{result['max_replication_pct']:.1f}%)")

        from kosmic.meta_analysis.io import record_meta_stage
        params = self._ctx('params') or {}
        record_meta_stage(self._project_folder, 'meta_loo', {
            'n_folds': len(self._ctx('de_dfs') or []),
            'studies': list(self._ctx('labels') or []),
            'replication_mode': 'direction_only',
            'pooling_methods': list(params.get('method_keys') or []),
            'calibration': params.get('calibration', 'analytical'),
            'avg_replication_pct': round(
                float(result['avg_replication_pct']), 2),
            'min_replication_pct': round(
                float(result['min_replication_pct']), 2),
            'max_replication_pct': round(
                float(result['max_replication_pct']), 2),
        }, selection=self._selection)

    def _on_loo_failed(self, message):
        self._loo_btn.setEnabled(True)
        if self.progress_bar is not None:
            self.progress_bar.setRange(0, 0)
        self._loo_tab.set_status(f"LOO failed: {message.splitlines()[0]}")
        self.log_message.emit(f"LOO failed: {message}")
        dialogs.warning(self, "LOO Failed", message)

    # --- LOO method benchmark ---
    def _run_loo_benchmark(self, level='current'):
        """
        Run LOO across method configurations.

        level:
            'current'    -- each chosen method alone + the consensus of
                            all chosen methods (fast)
            'quick'      -- 18 configs (singles x 2 Top50 + multi-subsets)
            'full'       -- 53 configs (3 states per method)
            'exhaustive' -- 127 configs (each variant independent)
        """
        from kosmic.meta_analysis.consensus_loo import (
            enumerate_method_configurations)

        if self._ctx('params') is None:
            dialogs.warning(
                self, "No Consensus Run Yet",
                "Run a normal consensus analysis first. The benchmark "
                "re-uses its method list and settings.")
            return

        de_dfs = self._ctx('de_dfs')
        labels = self._ctx('labels')
        method_keys = list(self._ctx('params').get('method_keys', []))
        if not de_dfs or len(de_dfs) < 3:
            dialogs.warning(
                self, "Too Few Studies",
                f"Benchmark needs at least 3 studies; "
                f"have {len(de_dfs) if de_dfs else 0}.")
            return

        # 'quick' fixes the Top50 dimension before its 31-subset sweep.
        top50 = bool(self._ctx('top50'))
        per_method_top50 = {k: top50 for k in (self._ctx('method_keys') or [])}

        if level == 'current':
            n_configs = (len(method_keys)
                         + (1 if len(method_keys) > 1 else 0))
        else:
            try:
                n_configs = len(enumerate_method_configurations(
                    level, user_top50_per_method=per_method_top50))
            except Exception:
                n_configs = -1

        if level in ('full', 'exhaustive') or (
                level == 'quick' and n_configs > 30):
            est_min = max(1, int(n_configs * 1.0))
            if not dialogs.confirm(self, "Confirm Combinatorial Benchmark", f"This will run {n_configs} LOO analyses. "
                f"Estimated time at ~1 minute per config: "
                f"~{est_min} minutes.\n\n"
                f"All configs use analytical calibration (CC perm "
                f"would take far longer). Continue?"):
                return

        if level == 'current' and not method_keys:
            dialogs.warning(
                self, "No Methods",
                "The most recent consensus had no methods selected.")
            return

        self._loo_tab.benchmark_btn.setEnabled(False)
        self._loo_tab.set_benchmark_status(
            f"Benchmark starting (level={level}): K={len(de_dfs)} "
            f"studies, {n_configs} configs (analytical calibration).")
        self._loo_tab.benchmark_table.setRowCount(0)
        self._loo_tab.benchmark_table.setColumnCount(0)

        self.log_message.emit(
            f"LOO benchmark ({level}): methods={method_keys}, "
            f"K={len(de_dfs)}, configs={n_configs}")

        self._loo_benchmark_worker = LOOBenchmarkWorker(
            de_dfs, labels, method_keys,
            dict(self._ctx('params')),
            replication_mode='direction_only',
            level=level,
            per_method_top50=per_method_top50)
        self._loo_benchmark_worker.config_progress.connect(
            self._on_benchmark_config_progress)
        run_worker(
            self._loo_benchmark_worker,
            on_finished=self._on_benchmark_finished,
            on_failed=self._on_benchmark_failed,
            on_progress=self._on_benchmark_progress,
        )

    # --- Replicator characterisation ---
    def _show_replicator_characterisation(self):
        """
        Show feature comparison between replicated and failed consensus genes,
        pooled across LOO folds.
        """
        if self._loo_result is None:
            dialogs.warning(
                self, "No LOO Result",
                "Run an LOO validation first. The characterisation "
                "pools per-fold per_gene tables to compare features "
                "of replicated vs failed consensus genes.")
            return

        try:
            from kosmic.meta_analysis.consensus_loo import (
                characterise_per_gene,
                characterise_per_study,
                characterise_replicators,
            )
            rep = characterise_replicators(self._loo_result)
            per_study = characterise_per_study(self._loo_result)
            per_gene = characterise_per_gene(self._loo_result)
        except Exception as e:
            dialogs.warning(
                self, "Characterisation Failed", str(e))
            return

        dlg = QDialog(self)
        dlg.setWindowTitle("Replicator Characterisation")
        dlg.resize(1000, 700)
        lay = QVBoxLayout(dlg)

        n_total = rep['n_genes_total']
        counts = rep['n_per_reason']
        n_rep = int(counts.get('replicated', 0))
        n_wd  = int(counts.get('wrong_direction', 0))
        n_ns  = int(counts.get('not_significant', 0))
        n_ms  = int(counts.get('missing', 0))
        header = QLabel(
            f"Pooled across {rep['n_folds']} LOO folds: {n_total} "
            f"gene-fold rows (one per (gene, fold) pair -- a gene "
            f"that's consensus in every fold contributes "
            f"{rep['n_folds']} rows).  "
            f"replicated {n_rep}, wrong_direction {n_wd}, "
            f"not_significant {n_ns}, missing {n_ms}.")
        header.setWordWrap(True)
        header.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse
            | Qt.TextInteractionFlag.TextSelectableByKeyboard)
        lay.addWidget(header)

        sub_tabs = QTabWidget()
        lay.addWidget(sub_tabs, stretch=1)

        # ----- Tab 1: feature comparison -----
        t1 = QWidget()
        t1_lay = QVBoxLayout(t1)
        t1_lay.setContentsMargins(4, 4, 4, 4)

        explain = SecondaryLabel(
            "Compare the mean / median of each feature between "
            "replicated and wrong_direction groups. Features that "
            "differ substantially are candidates for a post-hoc "
            "filter that would improve strict replication.\n"
            "Typical pattern when a filter is useful: replicated "
            "genes have larger |log2FC|, lower heterogeneity_i2, "
            "more studies, and smaller SE than wrong_direction "
            "failures.")
        explain.setWordWrap(True)
        t1_lay.addWidget(explain)

        summary = rep['feature_summary']
        feat_table = ResultsTable()
        feat_table.setAlternatingRowColors(True)
        feat_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)

        features = ('abs_consensus_log2FC', 'heterogeneity_i2',
                    'n_studies', 'se_approx', 'consensus_fdr')
        _NICE = {
            'abs_consensus_log2FC': '|log2FC|',
            'heterogeneity_i2':     'I2',
            'n_studies':            'n_studies',
            'se_approx':            'SE',
            'consensus_fdr':        'Cons. FDR',
        }
        cols = [Column('Reason', 'reason', 's'),
                Column('N', 'n', 'd')]
        for feat in features:
            for stat in ('mean', 'median'):
                cols.append(Column(
                    f'{_NICE[feat]} {stat}', f'{feat}_{stat}', '.3f'))

        def _decorate(item, row, col_i):
            # Swap to scientific for very small feature values.
            if col_i < 2:
                return
            key = cols[col_i].key
            v = row.get(key, float('nan'))
            try:
                f = float(v)
            except (TypeError, ValueError):
                return
            if np.isfinite(f) and abs(f) < 1e-3 and f != 0:
                item.setText(f'{f:.2e}')

        feat_table.set_schema(cols)
        feat_table.set_data(summary, decorate=_decorate)
        t1_lay.addWidget(feat_table)

        hint = self._interpret_replicator_summary(summary, rep)
        hint_label = QLabel(hint)
        hint_label.setWordWrap(True)
        hint_label.setProperty("role", "headline")
        hint_label.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse
            | Qt.TextInteractionFlag.TextSelectableByKeyboard)
        t1_lay.addWidget(hint_label)
        sub_tabs.addTab(t1, "Feature comparison")

        # ----- Tab 2: per held-out study -----
        t2 = QWidget()
        t2_lay = QVBoxLayout(t2)
        t2_lay.setContentsMargins(4, 4, 4, 4)
        t2_explain = SecondaryLabel(
            "One row per LOO fold (= one per held-out study). High "
            "wrong_direction % in a single study suggests a "
            "study-specific artifact (e.g. DESeq2 normalisation flip, "
            "batch effect) rather than a gene-level problem.")
        t2_explain.setWordWrap(True)
        t2_lay.addWidget(t2_explain)

        ps_table = ResultsTable()
        ps_table.setAlternatingRowColors(True)
        ps_table.set_schema([
            Column('Held out',     'study_held_out',      's'),
            Column('Consensus',    'n_consensus',         'd'),
            Column('Replicated',   'n_replicated',        'd'),
            Column('Wrong dir',    'n_wrong_direction',   'd'),
            Column('Not sig',      'n_not_significant',   'd'),
            Column('Missing',      'n_missing',           'd'),
            Column('Replicated %', 'pct_replicated',      '.1f'),
            Column('Wrong dir %',  'pct_wrong_direction', '.1f'),
            Column('Not sig %',    'pct_not_significant', '.1f'),
        ])
        ps_table.set_data(per_study)
        t2_lay.addWidget(ps_table)
        sub_tabs.addTab(t2, "Per held-out study")

        # ----- Tab 3: per-gene breakdown -----
        t3 = QWidget()
        t3_lay = QVBoxLayout(t3)
        t3_lay.setContentsMargins(4, 4, 4, 4)

        n_unique_genes = len(per_gene)
        n_wd_gt0 = int((per_gene['n_wrong_direction'] > 0).sum()) \
            if len(per_gene) else 0
        t3_summary = SecondaryLabel(
            f"{n_unique_genes} unique consensus genes across all folds. "
            f"{n_wd_gt0} had at least one wrong-direction fold. "
            f"Top of the list = genes most often contradicted by "
            f"their held-out study. Sorted by wrong-direction fold "
            f"count descending. If the same ~5% of genes account for "
            f"~50% of wrong-direction failures, those specific genes "
            f"are the concrete filter targets. Table capped at top "
            f"500 to keep the UI responsive.")
        t3_summary.setWordWrap(True)
        t3_lay.addWidget(t3_summary)

        pg_shown = per_gene.head(500) if len(per_gene) else per_gene
        pg_table = ResultsTable()
        pg_table.setAlternatingRowColors(True)
        pg_table.set_schema([
            Column('Gene',          'gene',                 's'),
            Column('N folds',       'n_folds',              'd'),
            Column('Replicated',    'n_replicated',         'd'),
            Column('Wrong dir',     'n_wrong_direction',    'd'),
            Column('Not sig',       'n_not_significant',    'd'),
            Column('Missing',       'n_missing',            'd'),
            Column('Wrong dir %',   'pct_wrong_direction',  '.1f'),
            Column('Replicated %',  'pct_replicated',       '.1f'),
            Column('|log2FC| mean', 'abs_log2FC_mean',      '.3f'),
            Column('|log2FC| max',  'abs_log2FC_max',       '.3f'),
            Column('I2 mean',       'i2_mean',              '.2f'),
            Column('I2 max',        'i2_max',               '.2f'),
        ])
        pg_table.set_data(pg_shown)
        t3_lay.addWidget(pg_table)
        sub_tabs.addTab(t3, "Per-gene breakdown")

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dlg.accept)
        lay.addWidget(close_btn)

        dlg.exec()

    def _interpret_replicator_summary(self, summary, rep):
        """
        Plain-English read on whether replicated vs wrong-direction differ
        on the key features. Heuristic, not a formal test.
        """
        if len(summary) == 0:
            return "(no rows in the summary)"
        df = summary.set_index('reason')
        if ('replicated' not in df.index
                or 'wrong_direction' not in df.index):
            return "Not enough data to compare groups."
        if df.loc['replicated', 'n'] == 0:
            return "No replicated genes in this run."
        if df.loc['wrong_direction', 'n'] == 0:
            return ("No wrong-direction failures in this run. "
                    "A filter on the replicator features would "
                    "have nothing to remove.")

        lines = []
        lines.append(
            f"n replicated = {int(df.loc['replicated','n'])},  "
            f"n wrong_direction = {int(df.loc['wrong_direction','n'])}.")

        # |log2FC|
        r_lfc = df.loc['replicated', 'abs_consensus_log2FC_median']
        w_lfc = df.loc['wrong_direction', 'abs_consensus_log2FC_median']
        if pd.notna(r_lfc) and pd.notna(w_lfc):
            ratio = r_lfc / w_lfc if w_lfc > 0 else float('inf')
            arrow = ">" if r_lfc > w_lfc else "<"
            lines.append(
                f"|log2FC|: replicated {arrow} wrong "
                f"({r_lfc:.2f} vs {w_lfc:.2f}; "
                f"ratio {ratio:.1f}x)" +
                ("  <-- filter candidate" if ratio > 1.3 else ""))

        # I^2
        r_i2 = df.loc['replicated', 'heterogeneity_i2_median']
        w_i2 = df.loc['wrong_direction', 'heterogeneity_i2_median']
        if pd.notna(r_i2) and pd.notna(w_i2):
            arrow = "<" if r_i2 < w_i2 else ">"
            lines.append(
                f"Heterogeneity I2: replicated {arrow} wrong "
                f"({r_i2:.2f} vs {w_i2:.2f})" +
                ("  <-- filter candidate"
                 if (w_i2 - r_i2) > 0.1 else ""))

        # n_studies
        r_n = df.loc['replicated', 'n_studies_mean']
        w_n = df.loc['wrong_direction', 'n_studies_mean']
        if pd.notna(r_n) and pd.notna(w_n):
            arrow = ">" if r_n > w_n else "<"
            lines.append(
                f"n_studies: replicated {arrow} wrong "
                f"({r_n:.1f} vs {w_n:.1f})")

        # SE
        r_se = df.loc['replicated', 'se_approx_median']
        w_se = df.loc['wrong_direction', 'se_approx_median']
        if pd.notna(r_se) and pd.notna(w_se):
            arrow = "<" if r_se < w_se else ">"
            lines.append(
                f"SE: replicated {arrow} wrong "
                f"({r_se:.3f} vs {w_se:.3f})")

        return "\n".join(lines)

    def _on_benchmark_progress(self, message):
        self._loo_tab.set_benchmark_status(message)
        self.log_message.emit(f"  bench: {message}")

    def _on_benchmark_config_progress(self, i, n, label):
        if self.progress_bar is not None:
            self.progress_bar.setRange(0, n)
            self.progress_bar.setValue(i)

    def _on_benchmark_finished(self, results):
        self._loo_tab.benchmark_btn.setEnabled(True)
        if self.progress_bar is not None:
            self.progress_bar.setRange(0, 0)
        self._loo_tab.populate_benchmark(results)
        self._loo_tab.set_benchmark_status(
            f"Benchmark complete: {len(results)} configs compared "
            f"across the same K folds.")

    def _on_benchmark_failed(self, message):
        self._loo_tab.benchmark_btn.setEnabled(True)
        if self.progress_bar is not None:
            self.progress_bar.setRange(0, 0)
        self._loo_tab.set_benchmark_status(
            f"Benchmark failed: {message.splitlines()[0]}")
        dialogs.warning(self, "Benchmark Failed", message)


class _OlinkTab(QWidget):
    """Which pooled genes are measurable on a routine Olink assay."""

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 12, 12, 12)
        lay.setSpacing(8)

        self._caption = QLabel(
            "Clinical-measurability flag, NOT a validation: panel "
            "membership says nothing about whether the protein is "
            "actually changed in disease.")
        self._caption.setWordWrap(True)
        self._caption.setProperty("role", "padded_status")
        lay.addWidget(self._caption)

        self._body = SecondaryLabel("Run a meta-analysis to see panel coverage.")
        self._body.setWordWrap(True)
        lay.addWidget(self._body)
        lay.addStretch()

    def populate(self, meta_df):
        if meta_df is None or meta_df.empty:
            self._body.setText("Run a meta-analysis to see panel coverage.")
            return
        from kosmic.meta_analysis.olink_panels import panel_summary
        try:
            summary = panel_summary(meta_df, fdr_threshold=DEFAULT_FDR)
        except (KeyError, ValueError) as exc:
            self._body.setText(f"Panel coverage unavailable: {exc}")
            return
        n = summary['n_consensus']
        lines = [f"{n:,} significant genes at FDR < {DEFAULT_FDR}.", ""]
        for info in summary['panels'].values():
            lines.append(
                f"{info['label']}: {info['n_on_panel']:,} "
                f"({info['pct']:.1f}%)")
        self._body.setText("\n".join(lines))


class _ComingSoonTab(QWidget):
    """Placeholder for a validation route that is not built yet."""

    def __init__(self, title: str, detail: str,
                 parent: Optional[QWidget] = None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 12, 12, 12)
        lay.setSpacing(8)
        heading = QLabel(f"{title} — coming soon")
        heading.setProperty("role", "stage_card_title")
        lay.addWidget(heading)
        body = SecondaryLabel(detail)
        body.setWordWrap(True)
        lay.addWidget(body)
        lay.addStretch()
