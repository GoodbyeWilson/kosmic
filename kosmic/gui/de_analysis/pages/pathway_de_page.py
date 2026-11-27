# DE Analysis -- Pathway DE page (hypothesis mode).
#
# Works on the gene DE run chosen in the run bar (ADR-016): a cell type's
# run scores that cell type's cells, with the donors its gene DE kept, and
# its results are saved in that run's folder under pathway_scoring/.
#
# One score per donor for each chosen pathway (DESeq2-normalised log2
# expression averaged over the pathway's genes), compared between disease
# and control donors with Welch's t-test. The sidebar holds the run's
# settings; above the results, the run bar and a pathway picker choose what
# is shown. Tabs: the picked pathway's per-donor scores; GSEA of the chosen
# pathways; the ambient-RNA check; which pathways had enough genes measured.

from pathlib import Path

import numpy as np
import pandas as pd
from PyQt6.QtWidgets import QFileDialog, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget

from kosmic import PATHWAY_GENE_DETECTION_PCT, PATHWAY_MIN_COVERAGE
from kosmic.de.pathway_scoring import SCORING_DESCRIPTION, SCORING_METHOD
from kosmic.gui.de_analysis.workers import PathwayScoringWorker
from kosmic.gui.shared import dialogs, run_worker
from kosmic.gui.shared.plots import InteractivePlot
from kosmic.gui.shared.theme import NoScrollComboBox, NoScrollDoubleSpinBox, get_color
from kosmic.gui.shared.widgets import (
    Column, HintLabel, PrimaryButton, ResultsTableView, SecondaryLabel, SectionHeader,
    SettingsGroup, SidebarTabbedPage,
)
from kosmic.scrna.counts import COUNT_SOURCE_TEXT

#: Prefix of a parent category's entry in the pathway picker.
_CATEGORY = "Category: "


class PathwayDEPage(SidebarTabbedPage):
    """Pathway DE: per-donor pathway scores, tested with Welch's t-test."""

    help_id = "de/pathway_de"

    def __init__(self, workspace):
        super().__init__()
        self.ws = workspace
        self.scoring_worker = None
        self.parent_scoring_worker = None
        self._pw_cache = {}
        self._fgsea_details = {}
        self._setup_ui()
        # The run bar's choice (ADR-016).
        self.ws.run_changed.connect(self._on_run_changed)

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------

    def _setup_ui(self):
        left = self.sidebar_layout
        left.addWidget(SectionHeader("PATHWAY SCORES"))
        about = HintLabel(SCORING_DESCRIPTION)
        about.setWordWrap(True)
        left.addWidget(about)

        settings = SettingsGroup("Settings", collapsible=True, expanded=True)
        count_row = QHBoxLayout()
        count_row.addWidget(QLabel("Counts:"))
        self.count_source_label = SecondaryLabel(COUNT_SOURCE_TEXT['raw'])
        self.count_source_label.setToolTip(
            "The counts the scores are computed from, chosen under Counts on "
            "the Load Data page.")
        count_row.addWidget(self.count_source_label, 1)
        settings.add_layout(count_row)

        detect_row = QHBoxLayout()
        detect_row.addWidget(QLabel("Gene detection ≥:"))
        self.gene_detection_spin = NoScrollDoubleSpinBox()
        self.gene_detection_spin.setRange(0, 100)
        self.gene_detection_spin.setValue(PATHWAY_GENE_DETECTION_PCT * 100)
        self.gene_detection_spin.setSingleStep(1)
        self.gene_detection_spin.setSuffix(" %")
        self.gene_detection_spin.setToolTip(
            "A gene counts towards a pathway's score only if it is detected "
            "(non-zero count) in at least this percentage of cells in disease "
            "or in control.")
        detect_row.addWidget(self.gene_detection_spin, 1)
        settings.add_layout(detect_row)

        coverage_row = QHBoxLayout()
        coverage_row.addWidget(QLabel("Min pathway coverage ≥:"))
        self.min_coverage_spin = NoScrollDoubleSpinBox()
        self.min_coverage_spin.setRange(0, 100)
        self.min_coverage_spin.setValue(PATHWAY_MIN_COVERAGE * 100)
        self.min_coverage_spin.setSingleStep(5)
        self.min_coverage_spin.setSuffix(" %")
        self.min_coverage_spin.setToolTip(
            "A pathway is scored only if at least this percentage of its "
            "genes pass the gene detection filter (see Genes measured).")
        coverage_row.addWidget(self.min_coverage_spin, 1)
        settings.add_layout(coverage_row)
        left.addWidget(settings)
        left.addStretch()

        self.run_pathway_btn = PrimaryButton("Run Pathway Analysis")
        self.run_pathway_btn.setEnabled(False)
        self.run_pathway_btn.clicked.connect(self._on_run_pathway_clicked)
        left.addWidget(self.run_pathway_btn)

        self.export_btn = QPushButton("Export CSV...")
        self.export_btn.setEnabled(False)
        self.export_btn.setToolTip(
            "Save every pathway's result (difference, p-value, FDR) as a table. "
            "Each run is already saved in the project for Figure Export.")
        self.export_btn.clicked.connect(self._export_csv)
        left.addWidget(self.export_btn)

        self.status_label = SecondaryLabel("Run Gene DE first, then Run Pathway Analysis.")
        self.status_label.setWordWrap(True)
        left.addWidget(self.status_label)

        # -- Above the results: which pathway the tabs show ---------------
        # (the run bar is put above this by the workspace).
        picker_row = QHBoxLayout()
        picker_row.setContentsMargins(12, 6, 12, 6)
        title = QLabel("Pathway")
        title.setStyleSheet("font-weight: 600;")
        picker_row.addWidget(title)
        self.pathway_combo = NoScrollComboBox()
        self.pathway_combo.setMinimumWidth(320)
        self.pathway_combo.setToolTip(
            "The pathway the per-donor scores, the GSEA curve and the ambient "
            "RNA check show. Each entry gives the pathway's direction and FDR.")
        self.pathway_combo.currentIndexChanged.connect(self._on_pathway_picked)
        picker_row.addWidget(self.pathway_combo, 1)
        picker_row.addStretch(1)
        self._picker = QWidget()
        self._picker.setLayout(picker_row)
        self.content_layout.insertWidget(0, self._picker)

        self.plot_tabs = self.tabs

        # -- Tab: per-donor scores -----------------------------------------
        score_tab = QWidget()
        score_layout = QVBoxLayout(score_tab)
        self.score_stats_label = SecondaryLabel("")
        self.score_stats_label.setWordWrap(True)
        score_layout.addWidget(self.score_stats_label)
        self._score_plot = InteractivePlot(
            title='Pathway score',
            left_label='Pathway score (per donor)',
            unavailable_message='Run Pathway Analysis to see the per-donor scores.',
        )
        score_layout.addWidget(self._score_plot, 1)
        self.plot_tabs.addTab(score_tab, "Pathway scores")

        # -- Tab: GSEA of the chosen pathways ---------------------------------
        fgsea_tab = QWidget()
        fgsea_layout = QVBoxLayout(fgsea_tab)
        fgsea_controls = QHBoxLayout()
        self.fgsea_run_btn = PrimaryButton("Run GSEA")
        self.fgsea_run_btn.clicked.connect(self._on_run_fgsea)
        fgsea_controls.addWidget(self.fgsea_run_btn)
        fgsea_controls.addStretch()
        fgsea_layout.addLayout(fgsea_controls)
        self.fgsea_status_label = SecondaryLabel(
            "Preranked GSEA of the chosen pathways: every tested gene ranked by "
            "its DESeq2 Wald statistic; a pathway is tested when 15-500 of its "
            "genes are in the ranking. It asks whether a pathway's genes stand "
            "out from the rest of the transcriptome.")
        self.fgsea_status_label.setWordWrap(True)
        fgsea_layout.addWidget(self.fgsea_status_label)
        self._fgsea_nes_plot = InteractivePlot(
            title='Pathway NES', left_label='',
            unavailable_message='Run GSEA to see the pathways\' NES.')
        fgsea_layout.addWidget(self._fgsea_nes_plot, 2)
        self._fgsea_plot = InteractivePlot(
            title='GSEA enrichment', left_label='Running enrichment score',
            unavailable_message='Run GSEA; the picked pathway\'s curve is shown here.')
        fgsea_layout.addWidget(self._fgsea_plot, 1)
        self.plot_tabs.addTab(fgsea_tab, "GSEA of these pathways")

        # -- Tab: ambient RNA check -------------------------------------------
        contam_tab = QWidget()
        contam_layout = QVBoxLayout(contam_tab)
        contam_controls = QHBoxLayout()
        contam_controls.addWidget(QLabel("Contamination:"))
        self.spill_metric_combo = NoScrollComboBox()
        self.spill_metric_combo.setToolTip(
            "A per-cell contamination metric from obs (pct_counts_<panel>, "
            "scored in the QC step, or DecontX's estimate).")
        contam_controls.addWidget(self.spill_metric_combo)
        self.spill_run_btn = PrimaryButton("Run check")
        self.spill_run_btn.clicked.connect(self._on_run_spillover)
        contam_controls.addWidget(self.spill_run_btn)
        contam_controls.addStretch()
        contam_layout.addLayout(contam_controls)
        self.spill_result_label = SecondaryLabel(
            "Scores the picked pathway per donor and asks how much of its "
            "disease effect is explained by ambient RNA from other cell types. "
            "Pick a contamination metric, then Run check.")
        self.spill_result_label.setWordWrap(True)
        contam_layout.addWidget(self.spill_result_label)
        self._spill_plot = InteractivePlot(
            title='Contamination vs pathway score',
            left_label='Pathway score (per donor)',
            unavailable_message=(
                'Score a contamination panel in the QC step, run the pathway '
                'analysis, then Run check.'),
        )
        contam_layout.addWidget(self._spill_plot, 1)
        self.plot_tabs.addTab(contam_tab, "Ambient RNA check")

        # -- Tab: genes measured ----------------------------------------------
        coverage_tab = QWidget()
        coverage_layout = QVBoxLayout(coverage_tab)
        self.coverage_status_label = SecondaryLabel(
            "How many of each pathway's genes were measured. A pathway is scored "
            "only if enough of its genes pass the gene detection filter. Run the "
            "pathway analysis to fill this in.")
        self.coverage_status_label.setWordWrap(True)
        coverage_layout.addWidget(self.coverage_status_label)
        self.coverage_table = ResultsTableView()
        self.coverage_table.set_schema([
            Column('Pathway', 'pathway', 's'),
            Column('Annotated', 'annotated', 'd'),
            Column('Present', 'present', 'd'),
            Column('Detected', 'detected', 'd'),
            Column('Coverage', 'coverage', '.0%'),
            Column('Status', 'status', 's'),
        ])
        coverage_layout.addWidget(self.coverage_table, 1)
        self.plot_tabs.addTab(coverage_tab, "Genes measured")

    # ------------------------------------------------------------------
    # Activation and the run chosen
    # ------------------------------------------------------------------

    def on_activated(self):
        """Called when page becomes visible."""
        self._populate_spill_metrics()
        self._fgsea_ready_state()
        self._refresh_count_source()
        self.run_pathway_btn.setEnabled(
            self.ws.current_adata is not None and bool(self.ws.pathway_gene_sets))
        if self.ws.pathway_de_results is None:
            self._try_load_previous_results()
        if self.ws.pathway_de_results is not None:
            self.status_label.setText("Results shown. Click Run to re-analyse.")
        self._populate_pathway_combo()

    def _try_load_previous_results(self):
        """Show the saved pathway results of the run chosen in the run bar
        (ADR-016), with its per-donor scores when they were saved."""
        folder = self.ws.run_saved_folder()
        if folder is None:
            return
        csv_path = folder / "pathway_de_results.csv"
        if not csv_path.exists():
            return
        try:
            self.ws.pathway_de_results = pd.read_csv(csv_path)
            from kosmic.de.pathway_scoring import PER_DONOR_SCORES_FILE, load_per_donor_scores
            self.ws.pathway_score_data = load_per_donor_scores(folder / PER_DONOR_SCORES_FILE)
            parent = folder / "parent_pathway_de_results.csv"
            if parent.exists():
                self.ws.parent_stats_df = pd.read_csv(parent)
            self._pw_cache = {}
            self.export_btn.setEnabled(True)
            self.ws.log_message.emit(f"Loaded previous pathway DE results from {csv_path.name}")
            if self.ws.analysis_mode == 'scoring':
                self.ws.mark_step_complete(4)  # scoring step 4 = Pathway DE
        except Exception as e:
            self.ws.log_message.emit(f"Could not load previous pathway results: {e}")

    def _on_run_changed(self):
        """Another gene DE run was chosen (ADR-016): clear the results shown,
        then show that run's saved ones."""
        self._pw_cache = {}
        self.ws.pathway_de_results = None
        self.ws.pathway_score_data = None
        self.ws.parent_stats_df = None
        self.export_btn.setEnabled(False)
        self._fgsea_details = {}
        for plot in (self._score_plot, self._fgsea_nes_plot, self._fgsea_plot):
            plot.clear_plot_items()
            plot.show_unavailable_message()
        self.score_stats_label.setText("")
        self._try_load_previous_results()
        folder = self.ws.run_saved_folder()
        if folder is not None:
            from kosmic.de.fgsea import FGSEA_RESULTS_FILE, load_fgsea_results
            res = load_fgsea_results(folder / FGSEA_RESULTS_FILE)
            if res is not None:
                self.ws.fgsea_results = res
                self._draw_fgsea_nes(res)
        if self.ws.pathway_de_results is None:
            self.status_label.setText(
                f"No saved pathway results for {self.ws.run_label() or 'this run'}. "
                "Click Run Pathway Analysis.")
        else:
            self.status_label.setText("Results shown. Click Run to re-analyse.")
        self._populate_pathway_combo()

    # ------------------------------------------------------------------
    # The pathway picker
    # ------------------------------------------------------------------

    @staticmethod
    def _result_text(row) -> str:
        """'down, FDR 0.002' for a pathway's result row, or 'not tested'."""
        lfc, fdr = row.get('logfoldchanges'), row.get('pvals_adj')
        if lfc is None or fdr is None or not np.isfinite(fdr):
            return "not tested"
        return f"{'up' if lfc > 0 else 'down'}, FDR {fdr:.2g}"

    def _populate_pathway_combo(self):
        """List the chosen pathways, most significant first, each with its
        direction and FDR; then the parent categories, when there are any."""
        current = self.pathway_combo.currentData()
        self.pathway_combo.blockSignals(True)
        self.pathway_combo.clear()
        names = list((self.ws.pathway_gene_sets or {}).keys())
        res = self.ws.pathway_de_results
        rows = {}
        if res is not None and 'names' in res.columns:
            rows = {str(r['names']): r for _, r in res.iterrows()}
        order = sorted(names, key=lambda n: (rows[n]['pvals_adj'] if n in rows
                                             and np.isfinite(rows[n]['pvals_adj'])
                                             else np.inf, n.lower()))
        for name in order:
            detail = self._result_text(rows[name]) if name in rows else (
                "not scored" if res is not None else "")
            self.pathway_combo.addItem(f"{name}  ·  {detail}" if detail else name, name)
        parents = getattr(self.ws, 'parent_gene_sets', None) or {}
        pstats = getattr(self.ws, 'parent_stats_df', None)
        prow = {}
        if pstats is not None and 'names' in getattr(pstats, 'columns', []):
            prow = {str(r['names']): r for _, r in pstats.iterrows()}
        if parents:
            self.pathway_combo.insertSeparator(self.pathway_combo.count())
            for name in parents:
                detail = self._result_text(prow[name]) if name in prow else ""
                text = f"{_CATEGORY}{name}" + (f"  ·  {detail}" if detail else "")
                self.pathway_combo.addItem(text, _CATEGORY + name)
        idx = self.pathway_combo.findData(current) if current else -1
        self.pathway_combo.setCurrentIndex(idx if idx >= 0 else 0)
        self.pathway_combo.blockSignals(False)
        self._picker.setVisible(self.pathway_combo.count() > 0)
        self._on_pathway_picked()

    def _picked(self):
        """(pathway name, its genes, is a parent category) of the picker's
        entry, or (None, [], False)."""
        data = self.pathway_combo.currentData()
        if not data:
            return None, [], False
        if data.startswith(_CATEGORY):
            name = data[len(_CATEGORY):]
            return name, list((self.ws.parent_gene_sets or {}).get(name, [])), True
        return data, list((self.ws.pathway_gene_sets or {}).get(data, [])), False

    def _result_row(self, name, is_category):
        res = (getattr(self.ws, 'parent_stats_df', None) if is_category
               else self.ws.pathway_de_results)
        if res is None or 'names' not in getattr(res, 'columns', []):
            return None
        hit = res[res['names'].astype(str) == name]
        return None if hit.empty else hit.iloc[0]

    def _on_pathway_picked(self, *_):
        name, genes, is_category = self._picked()
        if not name:
            return
        self._show_score_stats(name, is_category)
        self._draw_fgsea_enrichment(name)
        if self.ws.current_adata is None or not genes:
            return
        if name in self._pw_cache:
            self._show_pathway_plots(self._pw_cache[name])
            return
        from kosmic.de.de_analysis import tested_gene_set
        from kosmic.gui.de_analysis.workers import PathwayDataWorker
        cells = self.ws.run_adata()
        if cells is None:
            return
        self._pw_plot_worker = PathwayDataWorker(
            cells, genes, self.ws.sample_col, self.ws.condition_col,
            self.ws.disease_label, self.ws.control_label, name,
            method=SCORING_METHOD,
            precomputed_scores=None if is_category else self._stored_pathway_score(name),
            allowed_genes=tested_gene_set(getattr(self.ws, 'de_results', None)),
        )
        run_worker(self._pw_plot_worker, on_finished=self._on_pw_plot_ready)

    def _show_score_stats(self, name, is_category):
        row = self._result_row(name, is_category)
        if row is None:
            self.score_stats_label.setText(
                f"{name}: not scored (too few of its genes measured; see Genes "
                "measured)." if self.ws.pathway_de_results is not None else "")
            return
        parts = [f"{name}: difference (disease - control) "
                 f"{row['logfoldchanges']:+.3f} log2"]
        if np.isfinite(row.get('pvals', np.nan)):
            parts.append(f"p = {row['pvals']:.2g}")
        if np.isfinite(row.get('pvals_adj', np.nan)):
            parts.append(f"FDR = {row['pvals_adj']:.2g}")
        n_d, n_c = row.get('disease_samples'), row.get('control_samples')
        if n_d is not None and n_c is not None and np.isfinite(n_d) and np.isfinite(n_c):
            parts.append(f"{int(n_d)} disease, {int(n_c)} control donors")
        self.score_stats_label.setText("; ".join(parts) + " (Welch's t-test).")

    def _stored_pathway_score(self, pathway_name):
        """Per-donor (scores, sample_df) for a pathway from the last run, or None."""
        data = getattr(self.ws, 'pathway_score_data', None)
        if not data or pathway_name not in data['names']:
            return None
        col = data['names'].index(pathway_name)
        return data['matrix'][:, col], data['sample_df']

    def _on_pw_plot_ready(self, results):
        if results is None:
            return
        self._pw_cache[results['pathway_name']] = results
        if results['pathway_name'] == self._picked()[0]:
            self._show_pathway_plots(results)

    def _show_pathway_plots(self, results):
        """Per-donor scores of the picked pathway: a box per arm, a dot per donor."""
        import pyqtgraph as pg

        self._score_plot.clear_plot_items()
        scores = results.get('pathway_scores')
        if scores is None:
            self._score_plot.show_unavailable_message(
                f"{results['pathway_name']} was not scored: too few of its genes "
                "were measured (see Genes measured).")
            return
        self._score_plot.hide_unavailable_message()
        d_scores = scores[results['disease_mask']]
        c_scores = scores[results['control_mask']]
        if len(d_scores) == 0 or len(c_scores) == 0:
            return
        fg = get_color('fg_primary')
        dot = pg.mkBrush(get_color('fg_secondary'))
        box_w = 0.5

        def _box(x_pos, vals, color):
            """Tukey box plot at x_pos: IQR box, median, 1.5-IQR whiskers."""
            vals = np.asarray(vals, dtype=float)
            q1, med, q3 = np.percentile(vals, [25, 50, 75])
            iqr = q3 - q1
            if iqr > 0:
                lo = float(vals[vals >= q1 - 1.5 * iqr].min())
                hi = float(vals[vals <= q3 + 1.5 * iqr].max())
            else:
                lo, hi = float(vals.min()), float(vals.max())
            xl, xr = x_pos - box_w / 2, x_pos + box_w / 2
            pen = pg.mkPen(color, width=2)
            cap = box_w / 4
            self._score_plot.plot([xl, xr, xr, xl, xl], [q1, q1, q3, q3, q1], pen=pen)
            self._score_plot.plot([xl, xr], [med, med], pen=pg.mkPen(color, width=3))
            self._score_plot.plot([x_pos, x_pos], [q3, hi], pen=pen)
            self._score_plot.plot([x_pos, x_pos], [q1, lo], pen=pen)
            self._score_plot.plot([x_pos - cap, x_pos + cap], [hi, hi], pen=pen)
            self._score_plot.plot([x_pos - cap, x_pos + cap], [lo, lo], pen=pen)

        _box(0, c_scores, get_color('plot_control'))
        _box(1, d_scores, get_color('plot_disease'))
        rng = np.random.default_rng(42)
        for x0, vals in ((0.0, c_scores), (1.0, d_scores)):
            self._score_plot.plot(
                x0 + rng.uniform(-0.12, 0.12, len(vals)), vals, pen=None,
                symbolBrush=dot, symbolPen=pg.mkPen(fg, width=0.5),
                symbolSize=8, symbol='o')
        self._score_plot.getPlotItem().getAxis('bottom').setTicks(
            [[(0, results['control_label']), (1, results['disease_label'])]])
        self._score_plot.setXRange(-0.5, 1.5)
        self._score_plot.setTitle(
            f"{results['pathway_name']} -- per-donor pathway score", color=fg)

    # ------------------------------------------------------------------
    # Scoring
    # ------------------------------------------------------------------

    def _on_run_pathway_clicked(self):
        if self.ws.current_adata is None:
            return
        self.run_pathway_btn.setEnabled(False)
        self.status_label.setText("Running pathway scoring...")
        self._start_scoring()

    def _get_output_dir(self):
        """Where the scoring results are saved: the chosen run's folder (ADR-016)."""
        folder = self.ws.run_folder()
        if folder is not None:
            return folder
        return Path.cwd() / "de_analysis_results" / "pathway_scoring"

    def _scoring_cells(self):
        """The chosen run's cells (ADR-016), or None after saying why they
        cannot be selected."""
        adata = self.ws.run_adata()
        if adata is None and self.ws.current_adata is not None:
            dialogs.warning(
                self, "Cell type not found",
                f"No column of the loaded dataset holds the cell type "
                f"'{self.ws.run_label()}', so its cells cannot be scored. Load the "
                "dataset the per-cell-type run was made on.")
        return adata

    def _refresh_count_source(self):
        self.count_source_label.setText(
            COUNT_SOURCE_TEXT.get(self._count_source(), self._count_source()))

    def _count_source(self):
        """Which counts the loaded dataset holds: 'raw' or a layer name."""
        from kosmic.scrna.counts import loaded_count_source
        return loaded_count_source(self.ws.current_adata)

    def _coverage_kwargs(self):
        return {
            'gene_detection_pct': self.gene_detection_spin.value() / 100.0,
            'min_coverage': self.min_coverage_spin.value() / 100.0,
        }

    def _start_scoring(self):
        gene_de = self.ws.de_results
        if gene_de is None or len(gene_de) == 0:
            dialogs.warning(
                self, "Gene DE Required",
                "Run Gene DE first: the pathway scores use the same DESeq2 "
                "normalisation and donors as the gene DE run.")
            self.run_pathway_btn.setEnabled(True)
            return
        cells = self._scoring_cells()
        if cells is None:
            self.run_pathway_btn.setEnabled(True)
            return
        min_cells, min_counts = self.ws.run_donor_floors()
        self.scoring_worker = PathwayScoringWorker(
            cells, self.ws.sample_col, self.ws.condition_col,
            self.ws.control_label, self.ws.disease_label,
            self.ws.pathway_gene_sets, self._get_output_dir(),
            method=SCORING_METHOD, gene_de_results=gene_de,
            min_cells=min_cells, min_counts=min_counts,
            **self._coverage_kwargs(),
        )
        if self.ws.progress_bar:
            self.ws.progress_bar.setRange(0, 0)  # indeterminate
        run_worker(self.scoring_worker, on_finished=self._on_scoring_finished,
                   on_failed=self._on_scoring_failed, on_progress=self._on_progress)

    def _on_scoring_finished(self, payload):
        pathway_de_df, scores, coverage_report, message = payload
        self.ws.log_message.emit(message)
        self.ws.pathway_de_results = pathway_de_df
        self.ws.pathway_score_data = scores
        self._pw_cache = {}
        self._populate_coverage_table(coverage_report)
        self.export_btn.setEnabled(True)
        if self.ws.current_hierarchy is not None and self.ws.parent_gene_sets:
            self._start_parent_scoring()
        else:
            self._finalize_analysis()

    def _populate_coverage_table(self, coverage_report):
        """Fill Genes measured from the worker's per-pathway report."""
        self.ws.pathway_coverage_report = coverage_report
        if coverage_report is None or coverage_report.empty:
            self.coverage_status_label.setText("No coverage report.")
            self.coverage_table.set_data(pd.DataFrame())
            return
        n_scored = int((coverage_report['status'] == 'scored').sum())
        n_dropped = len(coverage_report) - n_scored
        self.coverage_status_label.setText(
            f"{n_scored} pathways scored, {n_dropped} not scored "
            f"(gene detected in ≥ {self.gene_detection_spin.value():.0f}% of cells; "
            f"≥ {self.min_coverage_spin.value():.0f}% of a pathway's genes detected).")
        ordered = (coverage_report
                   .assign(_scored=(coverage_report['status'] == 'scored'))
                   .sort_values(['_scored', 'coverage'], ascending=[True, True])
                   .drop(columns='_scored')
                   .reset_index(drop=True))
        self.coverage_table.set_data(ordered)

    def _on_scoring_failed(self, message: str):
        self.ws.log_message.emit(f"Pathway scoring failed: {message}")
        self.status_label.setText(f"Failed: {message}")
        self._finalize_analysis(ok=False)

    def _start_parent_scoring(self):
        cells = self.ws.run_adata()
        if cells is None:
            self._finalize_analysis()
            return
        min_cells, min_counts = self.ws.run_donor_floors()
        self.parent_scoring_worker = PathwayScoringWorker(
            cells, self.ws.sample_col, self.ws.condition_col,
            self.ws.control_label, self.ws.disease_label,
            self.ws.parent_gene_sets, self._get_output_dir(),
            level_label='parent_', method=SCORING_METHOD,
            min_cells=min_cells, min_counts=min_counts,
            **self._coverage_kwargs(),
        )
        run_worker(self.parent_scoring_worker,
                   on_finished=self._on_parent_scoring_finished,
                   on_failed=self._on_parent_scoring_failed,
                   on_progress=self._on_progress)

    def _on_parent_scoring_finished(self, payload):
        stats_df, _scores, _coverage, message = payload
        self.ws.log_message.emit(message)
        self.ws.parent_stats_df = stats_df
        self._finalize_analysis()

    def _on_parent_scoring_failed(self, message: str):
        self.ws.log_message.emit(f"Parent-level scoring failed: {message}")
        self._finalize_analysis()

    def _finalize_analysis(self, ok=True):
        self.run_pathway_btn.setEnabled(True)
        self.ws.status_message.emit("Ready")
        if self.ws.progress_bar:
            self.ws.progress_bar.setRange(0, 100)
            self.ws.progress_bar.setValue(0)
        if not ok:
            return
        self.ws.log_message.emit("Pathway analysis complete!")
        res = self.ws.pathway_de_results
        n_sig = 0
        if res is not None and 'pvals_adj' in res.columns:
            n_sig = int((pd.to_numeric(res['pvals_adj'], errors='coerce') < 0.05).sum())
        n = 0 if res is None else len(res)
        self.status_label.setText(f"Complete: {n} pathways scored, {n_sig} at FDR < 0.05.")
        if self.ws.analysis_mode == 'scoring':
            self.ws.mark_step_complete(4)
        self._record_pathway_provenance()
        self._populate_pathway_combo()

    def _record_pathway_provenance(self):
        """Append a pathway_de stage to the study's provenance, tagged with the
        upstream data token it ran against. Best-effort; never blocks."""
        h5ad = getattr(self.ws, 'h5ad_path', None)
        if not h5ad or self.ws.current_adata is None:
            return
        try:
            from kosmic import provenance
            study_dir = Path(h5ad).parent
            upstream = provenance.compute_fingerprint(
                self.ws.current_adata, h5ad_path=h5ad)['token']
            cov = self._coverage_kwargs()
            params = {
                'scoring_method': SCORING_METHOD,
                'test': "Welch's t-test",
                'count_source': self._count_source(),
                'gene_detection_pct': round(cov['gene_detection_pct'], 4),
                'min_coverage': round(cov['min_coverage'], 4),
                'n_gene_sets': len(self.ws.pathway_gene_sets or {}),
                'mode': getattr(self.ws, 'analysis_mode', None),
                # The gene DE run it was made on (ADR-016).
                **self.ws.run_identity(),
            }
            # One entry per gene DE run (ADR-016): a rerun replaces only the
            # entry of the same cell type and counts.
            provenance.record_stage(
                study_dir, Path(h5ad).stem, 'pathway_de', params,
                consumed_token=upstream,
                replace_matching={'cell_type': params.get('cell_type'),
                                  'counts': params.get('counts')})
        except Exception as e:
            self.ws.log_message.emit(f"Could not record pathway provenance: {e}")

    def _on_progress(self, message):
        self.status_label.setText(message)
        self.ws.log_message.emit(message)

    def _export_csv(self):
        res = self.ws.pathway_de_results
        if res is None or res.empty:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export pathway results", "", "CSV files (*.csv)")
        if path:
            res.to_csv(path, index=False)
            self.ws.log_message.emit(f"Pathway results exported to {path}")

    # ------------------------------------------------------------------
    # Ambient RNA check
    # ------------------------------------------------------------------

    def _populate_spill_metrics(self):
        """Fill the contamination-metric selector; Run needs a metric."""
        adata = self.ws.current_adata
        self.spill_metric_combo.blockSignals(True)
        self.spill_metric_combo.clear()
        metrics = []
        if adata is not None:
            metrics = [c for c in adata.obs.columns if str(c).startswith('pct_counts_')]
            if 'decontX_contamination' in adata.obs.columns:
                metrics.append('decontX_contamination')
        for col in metrics:
            label = ('DecontX (per-cell)' if col == 'decontX_contamination'
                     else str(col)[len('pct_counts_'):] or str(col))
            self.spill_metric_combo.addItem(label, col)
        self.spill_metric_combo.blockSignals(False)
        self.spill_run_btn.setEnabled(bool(metrics))
        if not metrics:
            self.spill_result_label.setText(
                "No contamination metric found. Re-run the QC step with a "
                "'Contam. panel' selected, then come back here.")

    def _on_run_spillover(self):
        adata = self._scoring_cells()   # the chosen run's cells (ADR-016)
        if adata is None:
            return
        contam_key = self.spill_metric_combo.currentData()
        pathway, genes, _ = self._picked()
        if not contam_key or not pathway or not genes:
            return
        from kosmic.gui.de_analysis.workers import SpilloverWorker
        self.spill_run_btn.setEnabled(False)
        self.spill_result_label.setText("Running contamination check...")
        self._spill_worker = SpilloverWorker(
            adata, pathway, genes,
            self.ws.sample_col, self.ws.condition_col,
            self.ws.disease_label, self.ws.control_label, contam_key,
            method=SCORING_METHOD,
        )
        run_worker(self._spill_worker, on_finished=self._on_spillover_ready,
                   on_failed=self._on_spillover_failed)

    def _on_spillover_failed(self, message):
        self.spill_run_btn.setEnabled(True)
        self.spill_result_label.setText(f"Could not run check: {message}")

    def _on_spillover_ready(self, res):
        self.spill_run_btn.setEnabled(True)
        if not res:
            self.spill_result_label.setText("No result (too few samples?).")
            return
        raw, adj = res['raw_effect'], res['adjusted_effect']
        pct = res['pct_attributable']
        p = res['adjusted_p']
        rho = res['spearman_score_contam']
        nd, nc = res['n_disease'], res['n_control']
        metric = str(res['contam_key'])[len('pct_counts_'):]
        pct_txt = "n/a" if pct != pct else f"{pct:.0f}%"          # NaN-safe
        if raw == 0:
            verdict = "No raw disease effect to attribute."
        elif pct == pct and pct >= 60:
            verdict = (f"Most of the effect (~{pct_txt}) is explained by "
                       f"{metric} contamination -- treat with caution.")
        elif p < 0.05:
            verdict = (f"Effect survives adjustment for {metric} "
                       f"(adjusted p={p:.2g}) -- not primarily contamination.")
        else:
            verdict = (f"Effect is no longer significant after adjusting for "
                       f"{metric} (adjusted p={p:.2g}).")
        self.spill_result_label.setText(
            f"{res['pathway']}: disease effect {raw:+.3f} -> {adj:+.3f} after "
            f"adjusting for {metric}; {pct_txt} attributable; "
            f"Spearman({metric}, score)={rho:+.2f}; n={nd} vs {nc}.  {verdict}")
        self._render_spill_plot(res)

    def _render_spill_plot(self, res):
        """Scatter of per-donor contamination vs pathway score, by condition."""
        import pyqtgraph as pg
        x = np.asarray(res['contamination'], dtype=float)
        y = np.asarray(res['score'], dtype=float)
        is_dis = np.asarray(res['is_disease'], dtype=float).astype(bool)
        fg = get_color('fg_primary')
        self._spill_plot.clear_plot_items()
        self._spill_plot.hide_unavailable_message()
        self._spill_plot.plot(
            x[~is_dis], y[~is_dis], pen=None, symbol='o', symbolSize=9,
            symbolBrush=pg.mkBrush(get_color('plot_control')),
            symbolPen=pg.mkPen(fg, width=0.5))
        self._spill_plot.plot(
            x[is_dis], y[is_dis], pen=None, symbol='o', symbolSize=9,
            symbolBrush=pg.mkBrush(get_color('plot_disease')),
            symbolPen=pg.mkPen(fg, width=0.5))
        metric = str(res['contam_key'])[len('pct_counts_'):]
        self._spill_plot.getPlotItem().setLabel('bottom', f"{metric} contamination (%)")
        self._spill_plot.setTitle(
            f"{res['pathway']}: contamination vs score (disease vs control)", color=fg)

    # ------------------------------------------------------------------
    # GSEA of the chosen pathways
    # ------------------------------------------------------------------

    def _fgsea_ready_state(self):
        gene_de = self.ws.de_results
        have_de = gene_de is not None and len(gene_de) > 0
        have_sets = bool(self.ws.pathway_gene_sets)
        self.fgsea_run_btn.setEnabled(have_de and have_sets)
        if not have_de:
            self.fgsea_status_label.setText("Run Gene DE first: GSEA ranks its genes.")
        elif not have_sets:
            self.fgsea_status_label.setText("Select gene sets first.")

    def _on_run_fgsea(self):
        # GSEA's null is the whole transcriptome, so it ranks the genome-wide
        # gene DE the run kept, not the hypothesis-restricted table.
        gene_de = getattr(self.ws, 'gene_de_genome_wide', None)
        if gene_de is None:
            gene_de = self.ws.de_results
        if gene_de is None or len(gene_de) == 0 or not self.ws.pathway_gene_sets:
            return
        if 'statistic' not in gene_de.columns:
            self.fgsea_status_label.setText(
                "This DE result has no 'statistic' column to rank by (it was saved "
                "by an earlier version). Run Gene DE again.")
            return
        from kosmic.gui.de_analysis.workers import FgseaWorker
        self.fgsea_run_btn.setEnabled(False)
        self.fgsea_status_label.setText("Running preranked GSEA...")
        self._fgsea_worker = FgseaWorker(gene_de, dict(self.ws.pathway_gene_sets))
        run_worker(self._fgsea_worker, on_finished=self._on_fgsea_ready,
                   on_failed=self._on_fgsea_failed)

    def _on_fgsea_failed(self, message):
        self.fgsea_run_btn.setEnabled(True)
        self.fgsea_status_label.setText(f"GSEA failed: {message}")

    def _on_fgsea_ready(self, payload):
        from kosmic.de.fgsea import GSEA_MIN_SIZE, GSEA_MAX_SIZE
        self.fgsea_run_btn.setEnabled(True)
        res, details = payload if isinstance(payload, tuple) else (payload, {})
        self._fgsea_details = details or {}
        chosen = list(self.ws.pathway_gene_sets or {})
        tested = set() if res is None else set(res['names'].astype(str))
        untested = [n for n in chosen if n not in tested]
        note = ""
        if untested:
            note = (f" Not tested, fewer than {GSEA_MIN_SIZE} (or more than "
                    f"{GSEA_MAX_SIZE}) of their genes in the ranking: "
                    + ", ".join(untested) + ".")
        if res is None or len(res) == 0:
            self.fgsea_status_label.setText("No pathway could be tested." + note)
            for p in (self._fgsea_nes_plot, self._fgsea_plot):
                p.clear_plot_items()
                p.show_unavailable_message()
            return
        n_sig = int((pd.to_numeric(res['pvals_adj'], errors='coerce') < 0.05).sum())
        self.fgsea_status_label.setText(
            f"GSEA complete: {len(res)} of {len(chosen)} pathways tested, {n_sig} at "
            "FDR < 0.05. NES > 0: towards the disease-up end of the ranking; "
            "NES < 0: towards the control-up end." + note)
        # Kept for Figure Export and saved, so it survives a reload.
        self.ws.fgsea_results = res
        if getattr(self.ws, 'project_dir', None):
            try:
                from kosmic.de.fgsea import FGSEA_RESULTS_FILE, save_fgsea_results
                save_fgsea_results(res, self.ws.run_folder() / FGSEA_RESULTS_FILE)
            except Exception as e:
                self.ws.log_message.emit(f"Could not save GSEA results: {e}")
        self._draw_fgsea_nes(res)
        self._draw_fgsea_enrichment(self._picked()[0])

    def _draw_fgsea_nes(self, res):
        """NES bar chart across pathways, sorted, coloured by FDR significance."""
        import pyqtgraph as pg
        res = res.copy()
        res['nes'] = pd.to_numeric(res['nes'], errors='coerce')
        res['pvals_adj'] = pd.to_numeric(res['pvals_adj'], errors='coerce')
        res = res.dropna(subset=['nes']).sort_values('nes', ascending=True)
        self._fgsea_nes_plot.clear_plot_items()
        self._fgsea_nes_plot.hide_unavailable_message()
        fg = get_color('fg_primary')
        sig_color = (60, 170, 160, 220)      # teal  = FDR < 0.05
        ns_color = (215, 120, 140, 200)      # pink  = not significant
        names = res['names'].astype(str).tolist()
        for i, (name, v, p) in enumerate(zip(names, res['nes'], res['pvals_adj'])):
            color = sig_color if (np.isfinite(p) and p < 0.05) else ns_color
            self._fgsea_nes_plot.addItem(pg.BarGraphItem(
                x0=[0], x1=[float(v)], y=[float(i)], height=0.7,
                brush=color, pen=pg.mkPen('k', width=0.5)))
            label = pg.TextItem(text=name.replace('_', ' '), color=fg,
                                anchor=(1, 0.5) if v >= 0 else (0, 0.5))
            label.setPos(0, i)
            self._fgsea_nes_plot.addItem(label)
        self._fgsea_nes_plot.addItem(pg.InfiniteLine(pos=0, angle=90, pen=pg.mkPen(fg, width=1)))
        self._fgsea_nes_plot.setYRange(-0.5, len(names) - 0.5)
        self._fgsea_nes_plot.getPlotItem().setLabel('bottom', 'Normalized Enrichment Score')
        self._fgsea_nes_plot.setTitle('Pathway NES', color=fg)

    def _draw_fgsea_enrichment(self, term):
        """Running-enrichment curve + hit rug for the picked pathway."""
        import pyqtgraph as pg
        details = self._fgsea_details
        if not term or not details or term not in details:
            self._fgsea_plot.clear_plot_items()
            self._fgsea_plot.show_unavailable_message(
                "Run GSEA; the picked pathway's curve is shown here."
                if not details else f"{term} was not tested by GSEA.")
            return
        d = details[term]
        res = np.asarray(d['RES'], dtype=float)
        hits = np.asarray(d['hits'], dtype=int)
        n = len(res)
        if n == 0:
            return
        x = np.arange(n, dtype=float)
        fg = get_color('fg_primary')
        es_line = (60, 160, 90)
        self._fgsea_plot.clear_plot_items()
        self._fgsea_plot.hide_unavailable_message()
        self._fgsea_plot.plot(x, res, pen=pg.mkPen(es_line, width=2))
        self._fgsea_plot.addItem(pg.InfiniteLine(pos=0, angle=0, pen=pg.mkPen(fg, width=1)))
        peak = int(np.argmax(np.abs(res)))
        self._fgsea_plot.plot(
            [float(x[peak])], [float(res[peak])], pen=None, symbol='o',
            symbolSize=9, symbolBrush=pg.mkBrush(*es_line), symbolPen=pg.mkPen(fg, width=0.5))
        span = float(res.max() - res.min()) or 1.0
        base = float(res.min()) - 0.06 * span
        tick = 0.04 * span
        for h in hits:
            self._fgsea_plot.plot([float(h), float(h)], [base, base - tick],
                                  pen=pg.mkPen(fg, width=1))
        self._fgsea_plot.setXRange(0, n)
        self._fgsea_plot.setTitle(f"{term} -- GSEA (NES = {d['nes']:+.2f})", color=fg)
        self._fgsea_plot.getPlotItem().setLabel(
            'bottom', 'Gene rank (disease-up  ->  control-up)')

    def refresh_theme(self):
        """Re-apply theme to pyqtgraph plots."""
        from kosmic.gui.shared.theme import style_pg_plot
        style_pg_plot(self._score_plot, title='Pathway score',
                      left_label='Pathway score (per donor)', bottom_label='')
        style_pg_plot(self._spill_plot, title='Contamination vs pathway score',
                      left_label='Pathway score (per donor)', bottom_label='')
