# Enrichment Page: preranked GSEA of every tested gene against a gene-set
# library (the bundled MSigDB Hallmark sets, or one downloaded from Enrichr,
# which includes KEGG, Reactome and the GO libraries).
#
# In the DE workspace the page works on the gene DE run chosen in the run
# bar (ADR-016): genes are ranked by the run's DESeq2 Wald statistic, the
# results are saved in that run's folder and, when another run is chosen,
# the page clears what it shows and shows that run's saved results. The
# Meta-Analysis workspace hosts the same page on the pooled gene list,
# ranked by the pooled z-statistic, and has no run bar.

from PyQt6.QtWidgets import (
    QLabel, QPushButton, QWidget, QVBoxLayout, QTableWidget, QFileDialog,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QCursor, QColor, QBrush

from kosmic.gui.shared.theme import (
    get_color, NoScrollComboBox, NoScrollDoubleSpinBox, NoScrollSpinBox,
)
from kosmic.gui.shared.widgets import (
    Column, ResultsTable, SidebarTabbedPage, SecondaryLabel, HintLabel,
    SectionHeader, PrimaryButton, StageAccordion, StageSummaryCard,
)
from kosmic.gui.shared.plots.interactive_enrichment_bar import InteractiveEnrichmentBar
from kosmic.gui.shared import run_worker
from kosmic import DEFAULT_FDR

# Published libraries bundled with KOSMIC, used exactly as released
# (kosmic.reference.pathways.libraries). Listed before the Enrichr ones.
_BUNDLED_LIBRARIES = [
    ("Hallmark 2026.2 (human)", "MSigDB_Hallmark_2026.2_Hs"),
    ("Hallmark 2026.2 (mouse)", "MSigDB_Hallmark_2026.2_Mm"),
]

_ENRICHR_LIBRARIES = [
    ("KEGG 2026", "KEGG_2026"),
    ("Reactome 2024", "Reactome_Pathways_2024"),
    ("MSigDB Hallmark 2020", "MSigDB_Hallmark_2020"),
    ("WikiPathways 2024 Human", "WikiPathways_2024_Human"),
    ("BioPlanet 2019", "BioPlanet_2019"),
    ("GO Biological Process 2025", "GO_Biological_Process_2025"),
    ("GO Cellular Component 2025", "GO_Cellular_Component_2025"),
    ("GO Molecular Function 2025", "GO_Molecular_Function_2025"),
]

#: The method this page runs, as recorded with each run.
METHOD = 'gsea'


class EnrichmentPage(SidebarTabbedPage):
    """Preranked GSEA of the DE result shown, against a chosen library."""

    help_id = "de/enrichment"

    # The DE table column the genes are ranked by: DESeq2's Wald statistic
    # (log2FC / SE), the usual ranking for preranked GSEA of DESeq2 results,
    # so a noisy fold change of a lowly expressed gene does not reach the
    # ends of the ranking. The meta-analysis host ranks by the pooled z.
    rank_column = 'statistic'
    rank_description = "the DESeq2 Wald statistic (log2 fold change / its SE)"

    enrichment_complete = pyqtSignal()

    def __init__(self, workspace):
        super().__init__()
        self.ws = workspace
        self._worker = None
        self._results_df = None
        self._gsea_library_id = None
        self._setup_ui()
        run_changed = getattr(self.ws, 'run_changed', None)
        if run_changed is not None:
            run_changed.connect(self._on_run_changed)

    def _setup_ui(self):
        left_layout = self.sidebar_layout
        left_layout.addWidget(SectionHeader("PRERANKED GSEA"))

        self._about = HintLabel(
            "Every tested gene is ranked by " + self.rank_description + ". "
            "For each gene set, GSEA asks whether its genes gather towards the "
            "top (up in disease) or the bottom (down) of the ranking. No "
            "significance cut-off is needed.")
        self._about.setWordWrap(True)
        left_layout.addWidget(self._about)

        self._accordion = StageAccordion(left_layout)
        self._card_pages = {}
        for key in ('library', 'sets'):
            page = QWidget()
            page_lay = QVBoxLayout(page)
            page_lay.setContentsMargins(0, 4, 0, 0)
            page_lay.setSpacing(4)
            self._card_pages[key] = (page, page_lay)

        # -- Library ---------------------------------------------------
        lib_lay = self._card_pages['library'][1]
        lib_lay.addWidget(QLabel("Library:"))
        self._library_combo = NoScrollComboBox()
        lib_lay.addWidget(self._library_combo)
        self._populate_libraries()
        lib_lay.addStretch()

        # -- Gene sets: which are tested, and how strictly they are called
        from kosmic.de.fgsea import GSEA_MIN_SIZE, GSEA_MAX_SIZE
        sets_lay = self._card_pages['sets'][1]
        size_tip = (
            "A gene set is tested only when this many of its genes are in the\n"
            f"ranking. {GSEA_MIN_SIZE} and {GSEA_MAX_SIZE} are the defaults of the "
            "Broad Institute's GSEA\nand of GSEApy: smaller sets give unstable "
            "scores, larger ones are\ntoo broad to interpret.")
        sets_lay.addWidget(QLabel("Min genes per set (in the ranking):"))
        self._gsea_min_size = NoScrollSpinBox()
        self._gsea_min_size.setRange(1, 10000)
        self._gsea_min_size.setValue(GSEA_MIN_SIZE)
        self._gsea_min_size.setToolTip(size_tip)
        sets_lay.addWidget(self._gsea_min_size)
        sets_lay.addWidget(QLabel("Max genes per set (in the ranking):"))
        self._gsea_max_size = NoScrollSpinBox()
        self._gsea_max_size.setRange(1, 100000)
        self._gsea_max_size.setValue(GSEA_MAX_SIZE)
        self._gsea_max_size.setToolTip(size_tip)
        sets_lay.addWidget(self._gsea_max_size)

        sets_lay.addWidget(QLabel("FDR threshold:"))
        self._p_thresh = NoScrollDoubleSpinBox()
        self._p_thresh.setRange(0.0001, 1.0)
        self._p_thresh.setValue(DEFAULT_FDR)
        self._p_thresh.setSingleStep(0.01)
        self._p_thresh.setDecimals(4)
        self._p_thresh.setToolTip("FDR below which a gene set is called significant.")
        sets_lay.addWidget(self._p_thresh)
        sets_lay.addStretch()

        for key, title, icon in (('library', "Library", "dataset"),
                                 ('sets', "Gene sets", "filter-funnel")):
            page, _ = self._card_pages[key]
            card = StageSummaryCard(title, icon=icon)
            card.set_settings_widget(page)
            left_layout.addWidget(card)
            self._accordion.add_card(key, card)
            setattr(self, f"_{key}_card", card)

        left_layout.addStretch()
        self._accordion.finalize()

        # Actions below the accordion, so collapsing a card never hides Run.
        self._run_btn = PrimaryButton("Run GSEA")
        self._run_btn.setMinimumHeight(32)
        self._run_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._run_btn.clicked.connect(self._run_enrichment)
        left_layout.addWidget(self._run_btn)

        self._export_btn = QPushButton("Export CSV...")
        self._export_btn.setEnabled(False)
        self._export_btn.setToolTip(
            "Save a copy of the results table wherever you like. Each run is "
            "already written into the project for Figure Export; this is for "
            "taking a copy elsewhere.")
        self._export_btn.clicked.connect(self._export_csv)
        left_layout.addWidget(self._export_btn)

        self._ranking_label = SecondaryLabel("")
        self._ranking_label.setWordWrap(True)
        left_layout.addWidget(self._ranking_label)

        self._status_label = SecondaryLabel("")
        self._status_label.setWordWrap(True)
        left_layout.addWidget(self._status_label)

        for _sig in (self._library_combo.currentIndexChanged,
                     self._p_thresh.valueChanged,
                     self._gsea_min_size.valueChanged,
                     self._gsea_max_size.valueChanged):
            _sig.connect(self._refresh_card_summaries)
        self._refresh_card_summaries()

        # -- Right content area: tabs (provided by SidebarTabbedPage) --
        self._tabs = self.tabs
        self._table = ResultsTable()
        self._table.setAlternatingRowColors(True)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSortingEnabled(True)
        self._tabs.addTab(self._table, "Results Table")

        self._plot_widget = InteractiveEnrichmentBar()
        self._tabs.addTab(self._plot_widget, "Bar Plot")

    def _populate_libraries(self):
        self._library_combo.blockSignals(True)
        self._library_combo.clear()
        for display, lib_id in _BUNDLED_LIBRARIES:
            self._library_combo.addItem(f"MSigDB: {display}", lib_id)
        for display, lib_id in _ENRICHR_LIBRARIES:
            self._library_combo.addItem(f"Enrichr: {display}", lib_id)
        self._library_combo.blockSignals(False)

    def _rank_by(self) -> str:
        """The DE table column the genes are ranked by."""
        return self.rank_column

    def _refresh_card_summaries(self, *_):
        """Two short lines per collapsed card."""
        if not hasattr(self, '_sets_card'):
            return  # still building the sidebar
        library = self._library_combo.currentText() or "no library"
        self._library_card.set_summary([library.split(": ")[-1], ""])
        self._sets_card.set_summary([
            f"{self._gsea_min_size.value()}-{self._gsea_max_size.value()} genes",
            f"FDR<{self._p_thresh.value():g}"])

    # ------------------------------------------------------------------
    # The ranking
    # ------------------------------------------------------------------

    def _ranking_table(self):
        """The gene table GSEA ranks: the genome-wide ranking the DE run
        kept (present even in hypothesis mode, so the null is the whole
        transcriptome), else the DE table shown. None when there is none."""
        gene_de = getattr(self.ws, 'gene_de_genome_wide', None)
        if gene_de is None:
            gene_de = getattr(self.ws, 'de_results', None)
        if gene_de is None or len(gene_de) == 0 or 'names' not in gene_de.columns:
            return None
        return gene_de

    def _refresh_ranking_note(self):
        """'n genes will be ranked', or why nothing can be."""
        gene_de = self._ranking_table()
        if gene_de is None:
            text = "No DE result yet: run Gene DE first."
        elif self.rank_column not in gene_de.columns:
            text = (f"This DE result has no '{self.rank_column}' column to rank by "
                    "(it was saved by an earlier version). Run the DE again.")
        else:
            n = int(gene_de[self.rank_column].notna().sum())
            text = f"{n:,} genes will be ranked."
        self._ranking_label.setText(text)

    def on_activated(self):
        """Called when page becomes visible."""
        self._refresh_ranking_note()

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def _run_enrichment(self):
        if self._worker is not None and self._worker.isRunning():
            return
        self._run_gsea()

    def _run_gsea(self):
        gene_de = self._ranking_table()
        if gene_de is None:
            self._status_label.setText("No DE result to rank. Run Gene DE first.")
            return
        rank_col = self._rank_by()
        if rank_col not in gene_de.columns:
            self._status_label.setText(
                f"This DE result has no '{rank_col}' column to rank by "
                "(it was saved by an earlier version). Run the DE again.")
            return
        lib_id = self._library_combo.currentData()
        if not lib_id:
            return
        lo, hi = self._gsea_min_size.value(), self._gsea_max_size.value()
        if lo > hi:
            self._status_label.setText(
                f"The minimum set size ({lo}) is above the maximum ({hi}).")
            return

        self._gsea_library_id = lib_id
        self.ws.gsea_library_name = self._library_combo.currentText()
        self._run_btn.setEnabled(False)
        self._status_label.setText(f"Running preranked GSEA ({lib_id})...")

        # The worker ranks by the column it is given under the name
        # 'logfoldchanges'.
        rank_df = gene_de[['names', rank_col]].copy()
        rank_df.columns = ['names', 'logfoldchanges']
        self._n_ranked = int(rank_df['logfoldchanges'].notna().sum())

        from kosmic.gui.de_analysis.workers import LibraryGseaWorker
        self._worker = LibraryGseaWorker(rank_df, lib_id, min_size=lo, max_size=hi)
        if self.ws.progress_bar:
            self.ws.progress_bar.setRange(0, 0)  # indeterminate
        run_worker(
            self._worker,
            on_finished=self._on_gsea_finished,
            on_failed=self._on_enrichment_failed,
            on_progress=self._on_worker_progress,
        )

    def _on_worker_progress(self, message):
        self._status_label.setText(message)
        self.ws.log_message.emit(message)

    def _on_enrichment_failed(self, message: str):
        self._run_btn.setEnabled(True)
        if self.ws.progress_bar:
            self.ws.progress_bar.setRange(0, 100)
            self.ws.progress_bar.setValue(0)
        self._status_label.setText(f"Failed: {message}")
        self.ws.log_message.emit(f"Enrichment failed: {message}")

    def _on_gsea_finished(self, payload):
        res, details, message = payload
        self.ws.gsea_library_details = details or {}
        self._run_btn.setEnabled(True)
        if self.ws.progress_bar:
            self.ws.progress_bar.setRange(0, 100)
            self.ws.progress_bar.setValue(0)
        self._status_label.setText(message)
        self.ws.log_message.emit(message)
        self.ws.status_message.emit(message)
        self._show_gsea_results(res)
        self._library_source = getattr(self._worker, 'library_source', None)
        self._record_study_enrichment(METHOD)

    # ------------------------------------------------------------------
    # Provenance
    # ------------------------------------------------------------------

    def _enrichment_params(self, method=METHOD) -> dict:
        """What one GSEA run used and found, for the provenance."""
        df = self._results_df
        n_terms = 0 if df is None else len(df)
        n_sig = 0
        if df is not None and 'pvals_adj' in getattr(df, 'columns', []):
            n_sig = int((df['pvals_adj'] < self._p_thresh.value()).sum())
        return {
            'method':      METHOD,
            'library':     self._library_combo.currentText(),
            # Where the library came from: bundled file or Enrichr download,
            # its version or address, a SHA-1 of its content, set and gene
            # counts and the time it was loaded.
            'library_source': getattr(self, '_library_source', None),
            'rank_by':     self._rank_by(),
            'fdr':         self._p_thresh.value(),
            'gsea_min_size': self._gsea_min_size.value(),
            'gsea_max_size': self._gsea_max_size.value(),
            'n_genes_ranked': getattr(self, '_n_ranked', None),
            'n_terms':     n_terms,
            'n_terms_sig': n_sig,
            # The gene DE run it was made on (ADR-016); the DE workspace only.
            **(self.ws.run_identity() if hasattr(self.ws, 'run_identity') else {}),
        }

    def _record_study_enrichment(self, method=METHOD):
        """Record a DE-workspace GSEA run in the study's provenance.

        One entry per library: rerunning a library replaces its entry,
        running another library adds one. The meta-analysis host has no
        study file and records into its own sidecar instead.
        """
        h5ad = getattr(self.ws, 'h5ad_path', None)
        if not h5ad:
            return
        try:
            from pathlib import Path

            from kosmic import provenance
            study_dir = Path(h5ad).parent
            params = self._enrichment_params()
            de_stage = provenance.last_stage(study_dir, 'gene_de')
            provenance.record_stage(
                study_dir, Path(h5ad).stem, 'enrichment', params,
                consumed_token=(de_stage or {}).get('consumed_token'),
                replace_matching={'method': METHOD, 'library': params['library'],
                                  'cell_type': params.get('cell_type'),
                                  'counts': params.get('counts')})
        except Exception as e:  # provenance is best-effort
            self.ws.log_message.emit(f"Could not record enrichment provenance: {e}")

    def _results_folder(self):
        """The folder of the gene DE run chosen in the run bar (ADR-016), or
        None to not save. The meta-analysis host overrides this."""
        if not getattr(self.ws, 'project_dir', None):
            return None
        run_folder = getattr(self.ws, 'run_folder', None)
        if run_folder is not None:
            return run_folder()     # the chosen run's folder (ADR-016)
        from kosmic.paths import de_pathway_scoring_dir
        return de_pathway_scoring_dir(self.ws.project_dir)

    def _gsea_results_path(self):
        """Where a GSEA run is saved, or None to not save it: the chosen
        run's folder (ADR-016), under the one file name the Figure Export
        library GSEA page reloads; each run replaces the last. The
        meta-analysis host overrides this."""
        folder = self._results_folder()
        if folder is None:
            return None
        from kosmic.de.fgsea import GSEA_LIBRARY_RESULTS_FILE
        return folder / GSEA_LIBRARY_RESULTS_FILE

    def _on_run_changed(self):
        """Another gene DE run was chosen (ADR-016): clear what the page
        shows, then show that run's saved GSEA, if it has one."""
        self._results_df = None
        self._table.setRowCount(0)
        self._export_btn.setEnabled(False)
        self._plot_widget.clear_plot_items()
        self._plot_widget.show_unavailable_message()
        self._refresh_ranking_note()
        name = self.ws.run_label() or "this run"
        folder = self.ws.run_saved_folder()
        if folder is None:
            self._status_label.setText("")
            return
        from kosmic.de import saved
        from kosmic.de.fgsea import GSEA_LIBRARY_RESULTS_FILE, load_fgsea_results
        res = load_fgsea_results(folder / GSEA_LIBRARY_RESULTS_FILE)
        if res is not None:
            details, library = saved.read_gsea_details(folder / GSEA_LIBRARY_RESULTS_FILE)
            self.ws.gsea_library_details = details or {}
            self.ws.gsea_library_name = library or None
            self._show_gsea_results(res, save=False)
            self._status_label.setText(f"Saved GSEA of {name} shown.")
            return
        self._status_label.setText(f"No saved GSEA for {name}.")

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def _show_gsea_results(self, res, save=True):
        """Populate the table + bar plot with the GSEA (NES) results; with
        'save', also save them for Figure Export."""
        import pandas as pd

        if res is None or res.empty:
            self._status_label.setText("GSEA returned no results.")
            self._results_df = pd.DataFrame()
            return

        self._results_df = res
        self.ws.gsea_library_results = res

        # Persist so Figure Export reloads the exact run instead of recomputing.
        path = self._gsea_results_path() if save else None
        if path is not None:
            from kosmic.de.fgsea import save_fgsea_results
            from kosmic.de.saved import save_gsea_details
            save_fgsea_results(res, path)
            # The running-enrichment curves the mountain figure draws.
            save_gsea_details(path, self.ws.gsea_library_details or {},
                              self.ws.gsea_library_name or '')

        alpha = self._p_thresh.value()
        down = QBrush(QColor(get_color('plot_downregulated')))
        up = QBrush(QColor(get_color('plot_upregulated')))

        display = res.copy()
        display['_le'] = display['leading_edge'].astype(str)
        display['_le_n'] = display['_le'].str.split(';').apply(
            lambda gs: len([g for g in gs if g.strip()]))
        display['_le_disp'] = display['_le'].str.slice(0, 200)
        from kosmic.reference.pathways.libraries import term_label
        display['_label'] = display['names'].map(term_label)

        from kosmic.de.fgsea import PERMUTATIONS, format_permutation_value
        floor = (f"Estimated from {PERMUTATIONS:,} gene-set permutations; 0 means "
                 f"none reached the observed score, so the value is below "
                 f"{1 / PERMUTATIONS:g} and is shown as < {1 / PERMUTATIONS:g}.")
        cols = [
            Column('Gene set',           '_label',    's'),
            Column('NES',                'nes',       '.2f'),
            Column('ES',                 'es',        '.2f'),
            Column('P-value',            'pvals',     '.2e', tooltip=floor),
            Column('FDR',                'pvals_adj', '.2e', tooltip=floor),
            Column('Leading edge',       '_le_n',     'd'),
            Column('Leading-edge genes', '_le_disp',  's'),
        ]

        def _decorate(item, row, col_i):
            key = cols[col_i].key
            if key in ('pvals', 'pvals_adj') and row[key] == 0:
                # Sorting still uses the stored 0.
                item.setText(format_permutation_value(0))
            if key == 'pvals_adj' and pd.notna(row['pvals_adj']) and row['pvals_adj'] < alpha:
                f = item.font()
                f.setBold(True)
                item.setFont(f)
            elif key == 'nes' and pd.notna(row['nes']):
                item.setForeground(down if row['nes'] < 0 else up)
            elif key == '_le_disp':
                item.setToolTip(row['_le'])
            elif key == '_label':
                item.setToolTip(str(row['names']))      # the library's own name

        self._table.set_schema(cols)
        self._table.set_data(display, decorate=_decorate)
        self._export_btn.setEnabled(True)
        self.enrichment_complete.emit()
        try:
            self._plot_widget.set_gsea(res, alpha)
        except Exception as e:
            self._plot_widget.show_unavailable_message(f"Plot rendering failed: {e}")

    def _export_csv(self):
        if self._results_df is None or self._results_df.empty:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Enrichment Results", "", "CSV files (*.csv)"
        )
        if path:
            self._results_df.to_csv(path, index=False)
            self.ws.log_message.emit(f"Enrichment results exported to {path}")
