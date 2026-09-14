# DE Analysis — Pathway DE Page (Step 3)
# Pathway scoring + three interactive plot tabs:
#   1. Pathway Activity  (bar chart from pathway_de_results)
#   2. Patient Dotplot   (per-sample scores, pathway + gene combos)
#   3. Top DE Genes      (drill-down from pathway → gene level)

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
)
from pathlib import Path
import numpy as np
import pandas as pd

from kosmic.de.pathway_scoring import SCORING_METHODS, SCORING_METHOD_HELP, DEFAULT_METHOD
from kosmic.paths import de_analysis_dir, de_pathway_scoring_dir
from kosmic import PATHWAY_GENE_DETECTION_PCT, PATHWAY_MIN_COVERAGE
from kosmic.gui.shared.theme import get_color, NoScrollComboBox, NoScrollDoubleSpinBox
from kosmic.gui.shared.plots import InteractivePathwayPlot
from kosmic.gui.de_analysis.workers import PathwayScoringWorker
from kosmic.gui.shared import dialogs, run_worker
from kosmic.gui.shared.widgets import SettingsGroup, SidebarTabbedPage, SecondaryLabel, HintLabel, PrimaryButton, ResultsTableView, Column


class PathwayDEPage(SidebarTabbedPage):
    """Step 3: pathway scoring + interactive plot tabs."""

    help_id = "de/pathway_de"

    def __init__(self, workspace):
        super().__init__()
        self.ws = workspace

        # Worker references (scoring only)
        self.scoring_worker = None
        self.parent_scoring_worker = None

        self._setup_ui()

    def _setup_ui(self):
        # Sidebar lives on self.sidebar_layout (provided by SidebarTabbedPage).
        left_layout = self.sidebar_layout

        status_group = SettingsGroup("Status")
        self.status_label = SecondaryLabel("Waiting for Gene DE analysis to complete...")
        self.status_label.setWordWrap(True)
        status_group.add_widget(self.status_label)
        left_layout.addWidget(status_group)

        method_group = SettingsGroup("Scoring Method", collapsible=True, expanded=True)
        method_row = QHBoxLayout()
        method_row.addWidget(QLabel("Method:"))
        self.scoring_method_combo = NoScrollComboBox()
        for key in SCORING_METHODS:
            self.scoring_method_combo.addItem(SCORING_METHODS[key], key)
        default_idx = list(SCORING_METHODS.keys()).index(DEFAULT_METHOD)
        self.scoring_method_combo.setCurrentIndex(default_idx)
        self.scoring_method_combo.currentIndexChanged.connect(self._on_scoring_method_changed)
        method_row.addWidget(self.scoring_method_combo)
        method_group.add_layout(method_row)

        self.scoring_desc_label = HintLabel(SCORING_METHOD_HELP.get(DEFAULT_METHOD, ''))
        self.scoring_desc_label.setWordWrap(True)
        method_group.add_widget(self.scoring_desc_label)

        count_row = QHBoxLayout()
        count_row.addWidget(QLabel("Count source:"))
        self.count_source_combo = NoScrollComboBox()
        self.count_source_combo.addItems(["Raw counts", "Decontaminated (DecontX)"])
        self.count_source_combo.setToolTip(
            "Raw counts, or DecontX-corrected counts (layers['decontX_counts']).\n"
            "Decontaminated is available once the Decontaminate step has run."
        )
        count_row.addWidget(self.count_source_combo)
        method_group.add_layout(count_row)

        # Gene detection + pathway coverage thresholds. Genes below the
        # detection floor are dropped from a pathway's score; pathways below
        # the coverage floor (fraction of annotated genes still detected) are
        # not scored at all.
        detect_row = QHBoxLayout()
        detect_row.addWidget(QLabel("Gene detection ≥:"))
        self.gene_detection_spin = NoScrollDoubleSpinBox()
        self.gene_detection_spin.setRange(0, 100)
        self.gene_detection_spin.setValue(PATHWAY_GENE_DETECTION_PCT * 100)
        self.gene_detection_spin.setSingleStep(1)
        self.gene_detection_spin.setSuffix(" %")
        self.gene_detection_spin.setToolTip(
            "A gene contributes to a pathway score only if detected (non-zero "
            "count) in at least this fraction of cells in disease OR control. "
            "Same donor-blind '% of cells' metric as the gene DE detection "
            "filter, but set separately here for scoring.")
        detect_row.addWidget(self.gene_detection_spin, 1)
        method_group.add_layout(detect_row)

        coverage_row = QHBoxLayout()
        coverage_row.addWidget(QLabel("Min pathway coverage ≥:"))
        self.min_coverage_spin = NoScrollDoubleSpinBox()
        self.min_coverage_spin.setRange(0, 100)
        self.min_coverage_spin.setValue(PATHWAY_MIN_COVERAGE * 100)
        self.min_coverage_spin.setSingleStep(5)
        self.min_coverage_spin.setSuffix(" %")
        self.min_coverage_spin.setToolTip(
            "A pathway is scored only if at least this fraction of the "
            "pathway's annotated genes survive the gene detection filter "
            "above. Keeps scores reflecting the biology measurable in this "
            "cell type rather than the full database annotation.")
        coverage_row.addWidget(self.min_coverage_spin, 1)
        method_group.add_layout(coverage_row)

        self.run_pathway_btn = PrimaryButton("Run Pathway Analysis")
        self.run_pathway_btn.setEnabled(False)
        self.run_pathway_btn.clicked.connect(self._on_run_pathway_clicked)
        method_group.add_widget(self.run_pathway_btn)

        left_layout.addWidget(method_group)

        pathway_group = SettingsGroup("Pathway", collapsible=True, expanded=True)
        self.dotplot_pathway_combo = NoScrollComboBox()
        self.dotplot_pathway_combo.currentTextChanged.connect(self._on_pathway_combo_changed)
        pathway_group.add_widget(self.dotplot_pathway_combo)

        self._pw_info_label = SecondaryLabel("")
        self._pw_info_label.setWordWrap(True)
        pathway_group.add_widget(self._pw_info_label)

        left_layout.addWidget(pathway_group)

        left_layout.addStretch()

        # ── Right panel — Tabbed interactive plots (provided by SidebarTabbedPage) ──
        self.plot_tabs = self.tabs

        # ── Tab 1: Pathway Activity ──
        pathway_tab = QWidget()
        pathway_layout = QVBoxLayout(pathway_tab)

        self.interactive_pathway = InteractivePathwayPlot()
        pathway_layout.addWidget(self.interactive_pathway, 1)

        pathway_controls = QHBoxLayout()

        # Level switcher (only visible when hierarchy exists)
        self.pathway_level_label = QLabel("Level:")
        pathway_controls.addWidget(self.pathway_level_label)
        self.pathway_level_combo = NoScrollComboBox()
        self.pathway_level_combo.addItems(["Subfamilies", "Parent Categories"])
        self.pathway_level_combo.currentTextChanged.connect(self._on_pathway_level_changed)
        pathway_controls.addWidget(self.pathway_level_combo)
        self.pathway_level_label.setVisible(False)
        self.pathway_level_combo.setVisible(False)

        reset_pathway_btn = QPushButton("Reset View")
        reset_pathway_btn.clicked.connect(self.interactive_pathway.reset_view)
        pathway_controls.addWidget(reset_pathway_btn)

        pathway_controls.addStretch()
        info1 = SecondaryLabel("Scroll to zoom, drag to pan, hover for details")
        pathway_controls.addWidget(info1)

        pathway_layout.addLayout(pathway_controls)
        self.plot_tabs.addTab(pathway_tab, "Pathway Activity")

        # ── Tab 2: Score Bar Chart (replaces Patient Dotplot) ──
        score_tab = QWidget()
        score_layout = QVBoxLayout(score_tab)

        from kosmic.gui.shared.plots import InteractivePlot
        self._score_plot = InteractivePlot(
            title='Pathway Score',
            left_label='Pathway score (per donor)',
            unavailable_message='Run pathway scoring to see per-pathway scores.',
        )
        score_layout.addWidget(self._score_plot, 1)

        self.plot_tabs.addTab(score_tab, "Pathway Score")

        # ── Tab 3: Contamination (ambient-spillover diagnostic) ──
        contam_tab = QWidget()
        contam_layout = QVBoxLayout(contam_tab)

        contam_controls = QHBoxLayout()
        contam_controls.addWidget(QLabel("Contamination:"))
        self.spill_metric_combo = NoScrollComboBox()
        self.spill_metric_combo.setToolTip(
            "A per-cell contamination metric from obs (pct_counts_<panel>), "
            "scored in the QC step.")
        contam_controls.addWidget(self.spill_metric_combo)
        contam_controls.addWidget(QLabel("Pathway:"))
        self.spill_pathway_combo = NoScrollComboBox()
        contam_controls.addWidget(self.spill_pathway_combo)
        self.spill_run_btn = PrimaryButton("Run check")
        self.spill_run_btn.clicked.connect(self._on_run_spillover)
        contam_controls.addWidget(self.spill_run_btn)
        contam_controls.addStretch()
        contam_layout.addLayout(contam_controls)

        self.spill_result_label = SecondaryLabel(
            "Scores this pathway per donor and asks how much of the "
            "disease effect is explained by contamination. Pick a metric "
            "and pathway, then Run check.")
        self.spill_result_label.setWordWrap(True)
        contam_layout.addWidget(self.spill_result_label)

        self._spill_plot = InteractivePlot(
            title='Contamination vs pathway score',
            left_label='Pathway score (per donor)',
            unavailable_message=(
                'Score a contamination panel in the QC step, run pathway '
                'scoring, then Run check.'),
        )
        contam_layout.addWidget(self._spill_plot, 1)

        self.plot_tabs.addTab(contam_tab, "Contamination")

        # ── Tab 4: Enrichment (fgsea) — independent gene-set test ──
        fgsea_tab = QWidget()
        fgsea_layout = QVBoxLayout(fgsea_tab)

        fgsea_controls = QHBoxLayout()
        self.fgsea_run_btn = PrimaryButton("Run fgsea")
        self.fgsea_run_btn.clicked.connect(self._on_run_fgsea)
        fgsea_controls.addWidget(self.fgsea_run_btn)
        fgsea_controls.addWidget(QLabel("Pathway:"))
        self.fgsea_pathway_combo = NoScrollComboBox()
        self.fgsea_pathway_combo.currentTextChanged.connect(
            self._on_fgsea_pathway_changed)
        fgsea_controls.addWidget(self.fgsea_pathway_combo)
        fgsea_controls.addStretch()
        fgsea_layout.addLayout(fgsea_controls)

        self.fgsea_status_label = SecondaryLabel(
            "Preranked GSEA of your pathways against the DESeq2 gene ranks — "
            "the independent test that corroborates the per-donor scores. "
            "Needs Gene DE run first.")
        self.fgsea_status_label.setWordWrap(True)
        fgsea_layout.addWidget(self.fgsea_status_label)

        # Primary figure: NES bar chart across all pathways, coloured by
        # significance (the standard fgsea summary plot).
        self._fgsea_nes_plot = InteractivePlot(
            title='Pathway NES (fgsea)', left_label='',
            unavailable_message='Run fgsea to see pathway NES.')
        fgsea_layout.addWidget(self._fgsea_nes_plot, 2)

        # Secondary: running-enrichment curve for the pathway picked above.
        self._fgsea_plot = InteractivePlot(
            title='GSEA enrichment', left_label='Running enrichment score',
            unavailable_message='Pick a pathway to see its enrichment curve.')
        fgsea_layout.addWidget(self._fgsea_plot, 1)

        self.plot_tabs.addTab(fgsea_tab, "Enrichment (fgsea)")

        # ── Tab 5: Coverage (which pathways scored, and why) ──
        coverage_tab = QWidget()
        coverage_layout = QVBoxLayout(coverage_tab)
        self.coverage_status_label = SecondaryLabel(
            "Per-pathway gene detection and coverage. A pathway is scored only "
            "if enough of its annotated genes are detected in this cell type. "
            "Run pathway scoring to populate.")
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
        self.plot_tabs.addTab(coverage_tab, "Coverage")

        self.plot_tabs.currentChanged.connect(self._on_tab_changed)

    # ------------------------------------------------------------------
    # Activation
    # ------------------------------------------------------------------

    def on_activated(self):
        """Called when page becomes visible."""
        has_hierarchy = self.ws.current_hierarchy is not None
        self.pathway_level_label.setVisible(has_hierarchy)
        self.pathway_level_combo.setVisible(has_hierarchy)

        self._populate_spill_controls()
        self._fgsea_ready_state()
        self._refresh_count_source()

        # Enable Run if data + gene sets are available
        if self.ws.current_adata is not None and self.ws.pathway_gene_sets:
            self.run_pathway_btn.setEnabled(True)
            if not self.dotplot_pathway_combo.count():
                self._populate_pathway_combo()

        # If results already in memory, just show status
        if self.ws.pathway_de_results is not None:
            self.status_label.setText("Pathway analysis complete. Click Run to re-analyse.")
            self.interactive_pathway.set_data(self.ws.pathway_de_results)
            return

        # Try to load previous results from disk
        self._try_load_previous_results()

        if self.ws.pathway_de_results is None:
            self.status_label.setText("Ready -- select a scoring method and click Run.")

    # ------------------------------------------------------------------
    # Helper: populate pathway combo
    # ------------------------------------------------------------------

    def _try_load_previous_results(self):
        """Check disk for previously saved pathway DE results."""
        import pandas as pd
        if not self.ws.project_dir:
            return

        # Pathway DE results CSV (in pathway_scoring/ subfolder)
        csv_path = de_pathway_scoring_dir(self.ws.project_dir) / "pathway_de_results.csv"
        if not csv_path.exists():
            return
        try:
            pathway_de = pd.read_csv(csv_path)
            self.ws.pathway_de_results = pathway_de
            # The CSV holds only summary stats; the per-donor matrix isn't on
            # disk, so the plot recomputes until the next run.
            self.ws.pathway_score_data = None
            self._pw_cache = {}

            # Plot directly from standard format (InteractivePathwayPlot accepts both)
            self.interactive_pathway.set_data(pathway_de)

            self.status_label.setText("Previous results loaded. Click Run to re-analyse.")
            self.ws.log_message.emit(f"Loaded previous pathway DE results from {csv_path.name}")
            if self.ws.analysis_mode == 'scoring':
                self.ws.mark_step_complete(4)  # scoring step 4 = Pathway DE
        except Exception as e:
            self.ws.log_message.emit(f"Could not load previous pathway results: {e}")

    def _populate_pathway_combo(self):
        """Fill the pathway selector from workspace gene sets."""
        self.dotplot_pathway_combo.blockSignals(True)
        self.dotplot_pathway_combo.clear()
        for pathway_name in self.ws.pathway_gene_sets.keys():
            self.dotplot_pathway_combo.addItem(pathway_name)
        self.dotplot_pathway_combo.blockSignals(False)
        if self.dotplot_pathway_combo.count() > 0:
            self._on_pathway_combo_changed(self.dotplot_pathway_combo.currentText())

    def _on_run_pathway_clicked(self):
        """Run pathway scoring with the selected method."""
        if self.ws.current_adata is None:
            return
        self.run_pathway_btn.setEnabled(False)
        method_desc = SCORING_METHODS.get(self._get_scoring_method(), "")
        self.status_label.setText(f"Running pathway scoring ({method_desc})...")
        self._start_scoring()

    def _get_output_dir(self):
        if self.ws.project_dir:
            return de_analysis_dir(self.ws.project_dir)
        return Path.cwd() / "de_analysis_results"

    # ------------------------------------------------------------------
    # Scoring workers
    # ------------------------------------------------------------------

    def _on_scoring_method_changed(self, index):
        key = self.scoring_method_combo.currentData()
        if key:
            self.scoring_desc_label.setText(SCORING_METHOD_HELP.get(key, SCORING_METHODS.get(key, "")))
        # The per-donor Pathway Score plot depends on the method -- re-render
        # the currently-selected pathway so the figure tracks the choice.
        current_pw = self.dotplot_pathway_combo.currentText()
        if current_pw and self.ws.current_adata is not None:
            self._on_pathway_combo_changed(current_pw)

    def _get_scoring_method(self):
        return self.scoring_method_combo.currentData() or DEFAULT_METHOD

    def _refresh_count_source(self):
        """Enable the DecontX count source only when the layer is present."""
        adata = self.ws.current_adata
        has_decontx = (adata is not None
                       and 'decontX_counts' in getattr(adata, 'layers', {}))
        self.count_source_combo.model().item(1).setEnabled(has_decontx)
        if not has_decontx and self.count_source_combo.currentIndex() == 1:
            self.count_source_combo.setCurrentIndex(0)

    def _selected_counts_layer(self):
        """Return the counts layer to score, or None for raw counts."""
        if self.count_source_combo.currentIndex() == 1:
            return 'decontX_counts'
        return None

    def _coverage_kwargs(self):
        """Gene-detection + pathway-coverage thresholds from the spinboxes."""
        return {
            'gene_detection_pct': self.gene_detection_spin.value() / 100.0,
            'min_coverage': self.min_coverage_spin.value() / 100.0,
        }

    def _stored_pathway_score(self, pathway_name):
        """Per-donor (scores, sample_df) for a pathway from the last run, or None.

        None when no run has retained the matrix (e.g. results loaded from
        disk), in which case the plot recomputes the score.
        """
        data = getattr(self.ws, 'pathway_score_data', None)
        if not data or pathway_name not in data['names']:
            return None
        col = data['names'].index(pathway_name)
        return data['matrix'][:, col], data['sample_df']

    def _start_scoring(self):
        output_dir = self._get_output_dir()
        method = self._get_scoring_method()
        gene_de = self.ws.de_results
        if method == 'deseq2' and (gene_de is None or len(gene_de) == 0):
            dialogs.warning(
                self, "Gene DE Required",
                "The DESeq2 pathway score is part of the DESeq2 framework — "
                "run Gene DE (DESeq2) first, then come back here so the pathway "
                "scores and the fgsea ranks come off the same DESeq2 run.")
            self.run_pathway_btn.setEnabled(True)
            if self.ws.progress_bar:
                self.ws.progress_bar.setRange(0, 100)
                self.ws.progress_bar.setValue(0)
            return

        self.scoring_worker = PathwayScoringWorker(
            self.ws.current_adata, self.ws.sample_col, self.ws.condition_col,
            self.ws.control_label, self.ws.disease_label,
            self.ws.pathway_gene_sets, output_dir,
            method=method,
            gene_de_results=gene_de,
            counts_layer=self._selected_counts_layer(),
            **self._coverage_kwargs(),
        )
        if self.ws.progress_bar:
            self.ws.progress_bar.setRange(0, 0)  # indeterminate
        run_worker(
            self.scoring_worker,
            on_finished=self._on_scoring_finished,
            on_failed=self._on_scoring_failed,
            on_progress=self._on_progress,
        )

    def _on_scoring_finished(self, payload):
        pathway_de_df, scores, coverage_report, message = payload
        self.ws.log_message.emit(message)
        self.ws.pathway_de_results = pathway_de_df
        # Retain the per-donor score matrix so the Pathway Score plot reads it
        # instead of recomputing; drop any per-pathway plot cache from a prior run.
        self.ws.pathway_score_data = scores
        self._pw_cache = {}
        self._populate_coverage_table(coverage_report)

        # Update interactive pathway plot (accepts either column-name format)
        self.interactive_pathway.set_data(pathway_de_df)

        # Run parent-level scoring if hierarchy exists
        if self.ws.current_hierarchy is not None and self.ws.parent_gene_sets:
            self._start_parent_scoring()
        else:
            self._finalize_analysis()

    def _populate_coverage_table(self, coverage_report):
        """Fill the Coverage tab from the worker's per-pathway report."""
        self.ws.pathway_coverage_report = coverage_report
        if coverage_report is None or coverage_report.empty:
            self.coverage_status_label.setText(
                "No coverage report (per-cell scoring methods are not gated by "
                "gene detection).")
            self.coverage_table.set_data(pd.DataFrame())
            return
        n_scored = int((coverage_report['status'] == 'scored').sum())
        n_dropped = len(coverage_report) - n_scored
        self.coverage_status_label.setText(
            f"{n_scored} pathways scored, {n_dropped} dropped "
            f"(detection ≥ {self.gene_detection_spin.value():.0f}% of cells, "
            f"coverage ≥ {self.min_coverage_spin.value():.0f}% of annotated genes).")
        # Show dropped pathways first so the exclusions are immediately visible.
        ordered = (coverage_report
                   .assign(_scored=(coverage_report['status'] == 'scored'))
                   .sort_values(['_scored', 'coverage'], ascending=[True, True])
                   .drop(columns='_scored')
                   .reset_index(drop=True))
        self.coverage_table.set_data(ordered)

    def _on_scoring_failed(self, message: str):
        self.ws.log_message.emit(f"Pathway scoring failed: {message}")
        self._finalize_analysis()

    def _start_parent_scoring(self):
        output_dir = self._get_output_dir()
        method = self._get_scoring_method()
        self.parent_scoring_worker = PathwayScoringWorker(
            self.ws.current_adata, self.ws.sample_col, self.ws.condition_col,
            self.ws.control_label, self.ws.disease_label,
            self.ws.parent_gene_sets, output_dir,
            level_label='parent_', method=method,
            counts_layer=self._selected_counts_layer(),
            **self._coverage_kwargs(),
        )
        run_worker(
            self.parent_scoring_worker,
            on_finished=self._on_parent_scoring_finished,
            on_failed=self._on_parent_scoring_failed,
            on_progress=self._on_progress,
        )

    def _on_parent_scoring_finished(self, payload):
        stats_df, _scores, _coverage, message = payload
        self.ws.log_message.emit(message)
        self.ws.parent_stats_df = stats_df
        self._finalize_analysis()

    def _on_parent_scoring_failed(self, message: str):
        self.ws.log_message.emit(f"Parent-level scoring failed: {message}")
        self._finalize_analysis()

    def _finalize_analysis(self):
        self.ws.log_message.emit("Pathway analysis complete!")
        method_desc = SCORING_METHODS.get(self._get_scoring_method(), "")
        self.status_label.setText(f"Complete -- {method_desc}")
        self.run_pathway_btn.setEnabled(True)
        self.ws.status_message.emit("Ready")
        if self.ws.progress_bar:
            self.ws.progress_bar.setRange(0, 100)
            self.ws.progress_bar.setValue(0)
        # Mark step complete (scoring step 4 = Pathway DE)
        if self.ws.analysis_mode == 'scoring':
            self.ws.mark_step_complete(4)
        self._record_pathway_provenance()

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
                'scoring_method': self._get_scoring_method(),
                'count_source': self._selected_counts_layer() or 'raw',
                'gene_detection_pct': round(cov['gene_detection_pct'], 4),
                'min_coverage': round(cov['min_coverage'], 4),
                'n_gene_sets': len(self.ws.pathway_gene_sets or {}),
                'mode': getattr(self.ws, 'analysis_mode', None),
            }
            provenance.record_stage(
                study_dir, Path(h5ad).stem, 'pathway_de', params,
                consumed_token=upstream)
        except Exception as e:
            self.ws.log_message.emit(f"Could not record pathway provenance: {e}")

    # ------------------------------------------------------------------
    # Progress
    # ------------------------------------------------------------------

    def _on_progress(self, message):
        self.status_label.setText(message)
        self.ws.log_message.emit(message)

    # ------------------------------------------------------------------
    # Pathway Activity controls
    # ------------------------------------------------------------------

    def _on_pathway_level_changed(self, level_text):
        if level_text == "Parent Categories" and self.ws.parent_stats_df is not None:
            self.interactive_pathway.set_data(self.ws.parent_stats_df)
        elif self.ws.pathway_de_results is not None:
            self.interactive_pathway.set_data(self.ws.pathway_de_results)

    # ------------------------------------------------------------------
    # Tab switching — enable/disable pathway selector
    # ------------------------------------------------------------------

    def _on_tab_changed(self, index):
        """Enable pathway combo only for pathway-specific tabs."""
        needs_pathway = index == 1  # Pathway Score tab
        self.dotplot_pathway_combo.setEnabled(needs_pathway)
        self._pw_info_label.setVisible(needs_pathway)

    # ------------------------------------------------------------------
    # Pathway expression plots (Score Bar, Gene Bars, Heatmap)
    # ------------------------------------------------------------------

    def _on_pathway_combo_changed(self, pathway_name):
        """When user selects a pathway, update Score Bar, Gene Bars, and Heatmap."""
        if not pathway_name or self.ws.current_adata is None:
            return

        # Cache is keyed by (pathway, method) -- the per-donor score depends
        # on the selected scoring method, so a method switch must recompute.
        if not hasattr(self, '_pw_cache'):
            self._pw_cache = {}
        method = self._get_scoring_method()
        if (pathway_name, method) in self._pw_cache:
            self._show_pathway_plots(self._pw_cache[(pathway_name, method)])
            return

        # Compute on the fly using Explorer's worker
        from kosmic.gui.de_analysis.workers import PathwayDataWorker as _PathwayDataWorker

        genes = list(self.ws.pathway_gene_sets.get(pathway_name, []))
        if not genes:
            return

        from kosmic.de.de_analysis import tested_gene_set
        self._pw_plot_worker = _PathwayDataWorker(
            self.ws.current_adata, genes,
            self.ws.sample_col, self.ws.condition_col,
            self.ws.disease_label, self.ws.control_label,
            pathway_name,
            method=method,
            counts_layer=self._selected_counts_layer(),
            precomputed_scores=self._stored_pathway_score(pathway_name),
            allowed_genes=tested_gene_set(getattr(self.ws, 'de_results', None)),
        )
        run_worker(
            self._pw_plot_worker,
            on_finished=self._on_pw_plot_ready,
        )

    def _on_pw_plot_ready(self, results):
        if results is None:
            return
        # Cache locally, keyed by (pathway, method).
        if not hasattr(self, '_pw_cache'):
            self._pw_cache = {}
        self._pw_cache[(results['pathway_name'], results.get('method'))] = results
        self._show_pathway_plots(results)

    def _show_pathway_plots(self, results):
        """Render score bar chart for the selected pathway."""
        import pyqtgraph as pg

        pw_name = results['pathway_name']
        disease = results['disease_label']
        control = results['control_label']
        pathway_scores = results['pathway_scores']
        disease_mask = results['disease_mask']
        control_mask = results['control_mask']

        fg = get_color('fg_primary')
        ctrl_color = get_color('plot_control')
        dis_color = get_color('plot_disease')

        # Parse themed dot color (rgba string -> QColor)
        import re as _re
        _dot_str = get_color('plot_dot')
        _dot_match = _re.match(r'rgba?\((\d+),\s*(\d+),\s*(\d+),?\s*(\d*)\)', _dot_str)
        if _dot_match:
            _dr, _dg, _db = int(_dot_match.group(1)), int(_dot_match.group(2)), int(_dot_match.group(3))
            _da = int(_dot_match.group(4)) if _dot_match.group(4) else 255
        else:
            _dr, _dg, _db, _da = 255, 255, 255, 180

        self._score_plot.clear_plot_items()
        self._score_plot.hide_unavailable_message()

        d_scores = pathway_scores[disease_mask]
        c_scores = pathway_scores[control_mask]

        if len(d_scores) > 0 and len(c_scores) > 0:
            box_w = 0.5

            def _box(x_pos, vals, color):
                """Tukey box plot at x_pos: IQR box, median, 1.5-IQR whiskers."""
                vals = np.asarray(vals, dtype=float)
                if len(vals) == 0:
                    return
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

            _box(0, c_scores, ctrl_color)
            _box(1, d_scores, dis_color)

            # Jittered per-donor points on top of the boxes.
            rng = np.random.default_rng(42)
            self._score_plot.plot(
                rng.uniform(-0.12, 0.12, len(c_scores)), c_scores, pen=None,
                symbolBrush=pg.mkBrush(_dr, _dg, _db, _da),
                symbolPen=pg.mkPen(fg, width=0.5), symbolSize=8, symbol='o',
            )
            self._score_plot.plot(
                1.0 + rng.uniform(-0.12, 0.12, len(d_scores)), d_scores, pen=None,
                symbolBrush=pg.mkBrush(_dr, _dg, _db, _da),
                symbolPen=pg.mkPen(fg, width=0.5), symbolSize=8, symbol='o',
            )

            self._score_plot.getPlotItem().getAxis('bottom').setTicks([[(0, control), (1, disease)]])
            self._score_plot.setXRange(-0.5, 1.5)
            self._score_plot.setTitle(f"{pw_name} -- per-donor pathway score", color=fg)

    # ------------------------------------------------------------------
    # Contamination tab — ambient-spillover diagnostic
    # ------------------------------------------------------------------

    def _populate_spill_controls(self):
        """Fill the contamination-metric and pathway selectors and enable Run
        only when both a metric and pathways are available."""
        adata = self.ws.current_adata

        self.spill_metric_combo.blockSignals(True)
        self.spill_metric_combo.clear()
        metrics = []
        if adata is not None:
            metrics = [c for c in adata.obs.columns
                       if str(c).startswith('pct_counts_')]
            if 'decontX_contamination' in adata.obs.columns:
                metrics.append('decontX_contamination')
        for col in metrics:
            if col == 'decontX_contamination':
                label = 'decontX (per-cell)'
            else:
                label = str(col)[len('pct_counts_'):] or str(col)
            self.spill_metric_combo.addItem(label, col)
        self.spill_metric_combo.blockSignals(False)

        self.spill_pathway_combo.blockSignals(True)
        self.spill_pathway_combo.clear()
        for name in (self.ws.pathway_gene_sets or {}).keys():
            self.spill_pathway_combo.addItem(name)
        self.spill_pathway_combo.blockSignals(False)

        ready = bool(metrics) and self.spill_pathway_combo.count() > 0
        self.spill_run_btn.setEnabled(ready)
        if not metrics:
            self.spill_result_label.setText(
                "No contamination metric found. Re-run the QC step with a "
                "'Contam. panel' selected, then come back here.")

    def _on_run_spillover(self):
        adata = self.ws.current_adata
        if adata is None:
            return
        contam_key = self.spill_metric_combo.currentData()
        pathway = self.spill_pathway_combo.currentText()
        genes = list((self.ws.pathway_gene_sets or {}).get(pathway, []))
        if not contam_key or not pathway or not genes:
            return

        from kosmic.gui.de_analysis.workers import SpilloverWorker
        self.spill_run_btn.setEnabled(False)
        self.spill_result_label.setText("Running contamination check...")
        self._spill_worker = SpilloverWorker(
            adata, pathway, genes,
            self.ws.sample_col, self.ws.condition_col,
            self.ws.disease_label, self.ws.control_label, contam_key,
            method=self._get_scoring_method(),
            counts_layer=self._selected_counts_layer(),
        )
        run_worker(
            self._spill_worker,
            on_finished=self._on_spillover_ready,
            on_failed=self._on_spillover_failed,
        )

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
        # Plain-language verdict.
        if raw == 0:
            verdict = "No raw disease effect to attribute."
        elif pct == pct and pct >= 60:
            verdict = (f"Most of the effect (~{pct_txt}) is explained by "
                       f"{metric} contamination — treat with caution.")
        elif p < 0.05:
            verdict = (f"Effect survives adjustment for {metric} "
                       f"(adjusted p={p:.2g}) — not primarily contamination.")
        else:
            verdict = (f"Effect is no longer significant after adjusting for "
                       f"{metric} (adjusted p={p:.2g}).")

        self.spill_result_label.setText(
            f"{res['pathway']} — disease effect {raw:+.3f} → {adj:+.3f} after "
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
        dis_color = get_color('plot_disease')
        ctrl_color = get_color('plot_control')

        self._spill_plot.clear_plot_items()
        self._spill_plot.hide_unavailable_message()

        self._spill_plot.plot(
            x[~is_dis], y[~is_dis], pen=None, symbol='o', symbolSize=9,
            symbolBrush=pg.mkBrush(ctrl_color), symbolPen=pg.mkPen(fg, width=0.5))
        self._spill_plot.plot(
            x[is_dis], y[is_dis], pen=None, symbol='o', symbolSize=9,
            symbolBrush=pg.mkBrush(dis_color), symbolPen=pg.mkPen(fg, width=0.5))

        metric = str(res['contam_key'])[len('pct_counts_'):]
        self._spill_plot.getPlotItem().setLabel('bottom', f"{metric} contamination (%)")
        self._spill_plot.setTitle(
            f"{res['pathway']}: contamination vs score "
            f"(disease vs control)", color=fg)

    # ------------------------------------------------------------------
    # Enrichment (fgsea) tab — independent gene-set test on DESeq2 ranks
    # ------------------------------------------------------------------

    def _fgsea_ready_state(self):
        """Gate Run on gene DE + gene sets being available."""
        gene_de = self.ws.de_results
        have_de = gene_de is not None and len(gene_de) > 0
        have_sets = bool(self.ws.pathway_gene_sets)
        self.fgsea_run_btn.setEnabled(have_de and have_sets)
        if not have_de:
            self.fgsea_status_label.setText(
                "Run Gene DE first — fgsea ranks genes by the DESeq2 result.")
        elif not have_sets:
            self.fgsea_status_label.setText("Select gene sets first.")

    def _on_run_fgsea(self):
        # Preranked GSEA's null is the whole transcriptome, so it must rank the
        # genome-wide gene DE -- not the hypothesis-restricted table (which
        # holds only the committed pathway genes). The genome-wide ranking is
        # preserved by the gene DE run; fall back to de_results in discovery
        # mode, where de_results is already genome-wide.
        gene_de = getattr(self.ws, 'gene_de_genome_wide', None)
        if gene_de is None:
            gene_de = self.ws.de_results
        if gene_de is None or len(gene_de) == 0 or not self.ws.pathway_gene_sets:
            return
        from kosmic.gui.de_analysis.workers import FgseaWorker
        self.fgsea_run_btn.setEnabled(False)
        self.fgsea_status_label.setText("Running preranked GSEA (fgsea)...")
        self._fgsea_worker = FgseaWorker(gene_de, dict(self.ws.pathway_gene_sets))
        run_worker(self._fgsea_worker,
                   on_finished=self._on_fgsea_ready,
                   on_failed=self._on_fgsea_failed)

    def _on_fgsea_failed(self, message):
        self.fgsea_run_btn.setEnabled(True)
        self.fgsea_status_label.setText(f"fgsea failed: {message}")

    def _on_fgsea_ready(self, payload):
        self.fgsea_run_btn.setEnabled(True)
        res, details = payload if isinstance(payload, tuple) else (payload, {})
        self._fgsea_details = details or {}
        if res is None or len(res) == 0:
            self.fgsea_status_label.setText("fgsea returned no results (coverage?).")
            self.fgsea_pathway_combo.clear()
            for p in (self._fgsea_nes_plot, self._fgsea_plot):
                p.clear_plot_items()
                p.show_unavailable_message()
            return
        n_sig = int((pd.to_numeric(res['pvals_adj'], errors='coerce') < 0.05).sum())
        self.fgsea_status_label.setText(
            f"fgsea complete — {len(res)} pathways, {n_sig} at FDR < 0.05. "
            "NES > 0 = enriched toward the disease-up end of the ranking, "
            "NES < 0 toward control-up.")

        # Retain for the figure export (read live) and persist so it survives a
        # reload -- the same reload-not-recompute contract as the pathway scores.
        self.ws.fgsea_results = res
        if getattr(self.ws, 'project_dir', None):
            try:
                from kosmic.de.fgsea import save_fgsea_results, FGSEA_RESULTS_FILE
                save_fgsea_results(
                    res, de_pathway_scoring_dir(self.ws.project_dir) / FGSEA_RESULTS_FILE)
            except Exception as e:
                self.ws.log_message.emit(f"Could not save fgsea results: {e}")

        # Primary: NES bar chart (all pathways, coloured by significance).
        self._draw_fgsea_nes(res)

        # Pathway selector (most significant first) + draw the top curve.
        res_sig = res.sort_values('pvals_adj', kind='stable')
        self.fgsea_pathway_combo.blockSignals(True)
        self.fgsea_pathway_combo.clear()
        for name in res_sig['names'].astype(str):
            if name in self._fgsea_details:
                self.fgsea_pathway_combo.addItem(name)
        self.fgsea_pathway_combo.blockSignals(False)
        if self.fgsea_pathway_combo.count() > 0:
            self._draw_fgsea_enrichment(self.fgsea_pathway_combo.currentText())

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
        nes_vals = res['nes'].to_numpy()
        padj = res['pvals_adj'].to_numpy()

        for i, (name, v, p) in enumerate(zip(names, nes_vals, padj)):
            color = sig_color if (np.isfinite(p) and p < 0.05) else ns_color
            self._fgsea_nes_plot.addItem(pg.BarGraphItem(
                x0=[0], x1=[float(v)], y=[float(i)], height=0.7,
                brush=color, pen=pg.mkPen('k', width=0.5)))
            label = pg.TextItem(text=name.replace('_', ' '), color=fg,
                                anchor=(1, 0.5) if v >= 0 else (0, 0.5))
            label.setPos(0, i)
            self._fgsea_nes_plot.addItem(label)

        self._fgsea_nes_plot.addItem(pg.InfiniteLine(
            pos=0, angle=90, pen=pg.mkPen(fg, width=1)))
        self._fgsea_nes_plot.setYRange(-0.5, len(names) - 0.5)
        self._fgsea_nes_plot.getPlotItem().setLabel(
            'bottom', 'Normalized Enrichment Score')
        self._fgsea_nes_plot.setTitle('Pathway NES (fgsea)', color=fg)

    def _on_fgsea_pathway_changed(self, term):
        if term:
            self._draw_fgsea_enrichment(term)

    def _draw_fgsea_enrichment(self, term):
        """Native running-enrichment (GSEA) curve + hit rug for one pathway."""
        import pyqtgraph as pg
        details = getattr(self, '_fgsea_details', None)
        if not details or term not in details:
            return
        d = details[term]
        res = np.asarray(d['RES'], dtype=float)
        hits = np.asarray(d['hits'], dtype=int)
        n = len(res)
        if n == 0:
            return
        x = np.arange(n, dtype=float)
        fg = get_color('fg_primary')
        es_line = (60, 160, 90)                          # GSEA-classic green

        self._fgsea_plot.clear_plot_items()
        self._fgsea_plot.hide_unavailable_message()

        # Running enrichment score curve.
        self._fgsea_plot.plot(x, res, pen=pg.mkPen(es_line, width=2))
        # Zero baseline.
        self._fgsea_plot.addItem(pg.InfiniteLine(
            pos=0, angle=0, pen=pg.mkPen(fg, width=1)))
        # Peak ES marker.
        peak = int(np.argmax(np.abs(res)))
        self._fgsea_plot.plot(
            [float(x[peak])], [float(res[peak])], pen=None, symbol='o',
            symbolSize=9, symbolBrush=pg.mkBrush(*es_line),
            symbolPen=pg.mkPen(fg, width=0.5))
        # Hit rug beneath the curve.
        span = float(res.max() - res.min()) or 1.0
        base = float(res.min()) - 0.06 * span
        tick = 0.04 * span
        for h in hits:
            self._fgsea_plot.plot([float(h), float(h)], [base, base - tick],
                                  pen=pg.mkPen(fg, width=1))

        self._fgsea_plot.setXRange(0, n)
        self._fgsea_plot.setTitle(
            f"{term} — GSEA (NES = {d['nes']:+.2f})", color=fg)
        self._fgsea_plot.getPlotItem().setLabel(
            'bottom', 'Gene rank (disease-up  →  control-up)')

    def refresh_theme(self):
        """Re-apply theme to pyqtgraph plots."""
        from kosmic.gui.shared.theme import style_pg_plot
        style_pg_plot(self._score_plot, title='Pathway Score',
                      left_label='Pathway score (per donor)', bottom_label='')
        if hasattr(self, '_spill_plot'):
            style_pg_plot(self._spill_plot, title='Contamination vs pathway score',
                          left_label='Pathway score (per donor)', bottom_label='')
        if hasattr(self, 'interactive_pathway'):
            self.interactive_pathway.refresh_theme()

