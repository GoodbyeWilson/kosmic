# Enrichment Page (DE Discovery mode): over-representation analysis of
# significant DE genes via GO (local elim) or Enrichr-based ORA, run
# separately for upregulated and downregulated genes.

from PyQt6.QtWidgets import (
    QLabel, QPushButton,
    QCheckBox, QWidget, QVBoxLayout, QHBoxLayout,
    QTableWidget, QHeaderView,
    QFileDialog,
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
from kosmic import (
    DEFAULT_FDR, DEFAULT_LFC_THRESHOLD,
    DE_ENRICH_IF_DEFAULT_PCT, DE_ENRICH_IF_ON_BY_DEFAULT,
)

# Libraries grouped by source
_GO_LIBRARIES = [
    ("GO Biological Process", "GO_BP"),
    ("GO Molecular Function", "GO_MF"),
    ("GO Cellular Component", "GO_CC"),
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

# Enrichment methods. ORA tests the significant-gene list against a background;
# preranked GSEA (fgsea) is threshold-free, ranking every gene -- the same
# rank-based method used for the a-priori panel, run here against a library.
_METHODS = [
    ("Over-representation (ORA)", "ora"),
    ("Preranked GSEA (fgsea)", "gsea"),
]


class EnrichmentPage(SidebarTabbedPage):
    """
    Discovery mode endpoint: pathway enrichment of significant DE genes.

    Runs a single ORA pass over all significant genes, then annotates
    each term with the majority direction of its overlapping DE genes
    (the Direction column). It does NOT run up and down separately.
    """

    help_id = "de/enrichment"

    # Column used to rank genes for preranked GSEA. Overridden by the
    # meta-analysis host; DE leaves it at the log2FC default.
    rank_column = 'logfoldchanges'

    enrichment_complete = pyqtSignal()

    def __init__(self, workspace):
        super().__init__()
        self.ws = workspace
        self._worker = None
        self._results_df = None
        self._lfc_map = {}  # gene -> log2FC for direction annotation
        self._setup_ui()

    def _setup_ui(self):
        # Sidebar lives on self.sidebar_layout (provided by SidebarTabbedPage).
        #
        # Three cards, because there are three decisions: which genes go
        # in, how they are tested, and how strictly terms are called. It
        # was five SettingsGroups, with Method and Library split apart
        # even though picking GSEA rewrites the library list, and with
        # Run buried inside the Settings group so collapsing it hid the
        # button.
        left_layout = self.sidebar_layout
        left_layout.addWidget(SectionHeader("ENRICHMENT SETUP"))
        self._accordion = StageAccordion(left_layout)

        self._card_pages = {}
        for key in ('genes', 'method', 'terms'):
            page = QWidget()
            page_lay = QVBoxLayout(page)
            page_lay.setContentsMargins(0, 4, 0, 0)
            page_lay.setSpacing(4)
            self._card_pages[key] = (page, page_lay)

        genes_lay = self._card_pages['genes'][1]
        genes_lay.addWidget(HintLabel(
            "Which genes are handed to the test. ORA uses this list; "
            "preranked GSEA ignores the cutoffs and ranks every gene, "
            "so they grey out."))

        # pct cells filter -- mirrors the gene DE page. Critical for
        # snRNA where the genome-wide FDR is dominated by dropout-noise
        # genes that drag every real signal's adjusted p toward 1.0.
        # Checkbox above its spinbox, not beside it. Side by side the
        # two minimum widths added up to 362px inside a 274px card.
        # Every other control here is label-above-control anyway.
        self._pct_filter_check = QCheckBox("Detection filter")
        self._pct_filter_check.setChecked(DE_ENRICH_IF_ON_BY_DEFAULT)
        self._pct_filter_check.setToolTip(
            "Sensitivity analysis, off by default.\n\n"
            "Drops genes detected in fewer than this % of cells in "
            "EITHER arm, then recomputes BH FDR on what is left -- so "
            "the gene list tested here stops matching the one the DE "
            "page exported.\n\n"
            "Leave off to enrich exactly the significant genes your DE "
            "run produced. The ORA background is the tested-gene "
            "universe either way.")

        self._pct_filter = NoScrollDoubleSpinBox()
        self._pct_filter.setRange(0, 100)
        self._pct_filter.setValue(DE_ENRICH_IF_DEFAULT_PCT * 100)
        self._pct_filter.setSingleStep(5)
        self._pct_filter.setSuffix(" %")
        self._pct_filter.setEnabled(DE_ENRICH_IF_ON_BY_DEFAULT)
        self._pct_filter.setToolTip(self._pct_filter_check.toolTip())
        # Capped so the pair always fits one sidebar row: a spinbox
        # asks for far more than it needs to show two digits.
        self._pct_filter.setMinimumWidth(64)
        self._pct_filter.setMaximumWidth(96)

        # One row, because it is one setting: the box switches the
        # filter on and the number is its threshold.
        pct_row = QHBoxLayout()
        pct_row.setContentsMargins(0, 0, 0, 0)
        pct_row.setSpacing(6)
        pct_row.addWidget(self._pct_filter_check, 1)
        pct_row.addWidget(self._pct_filter)
        genes_lay.addLayout(pct_row)

        self._pct_filter_check.toggled.connect(self._pct_filter.setEnabled)
        self._pct_filter_check.toggled.connect(lambda _c: self._update_gene_count())
        self._pct_filter.valueChanged.connect(lambda _v: self._update_gene_count())

        genes_lay.addWidget(QLabel("FDR cutoff:"))
        self._fdr_spin = NoScrollDoubleSpinBox()
        self._fdr_spin.setRange(0.0001, 1.0)
        self._fdr_spin.setValue(DEFAULT_FDR)
        self._fdr_spin.setSingleStep(0.01)
        self._fdr_spin.setDecimals(4)
        self._fdr_spin.valueChanged.connect(self._update_gene_count)
        genes_lay.addWidget(self._fdr_spin)

        genes_lay.addWidget(QLabel("|log2FC| cutoff:"))
        self._lfc_spin = NoScrollDoubleSpinBox()
        self._lfc_spin.setRange(0.0, 5.0)
        self._lfc_spin.setValue(DEFAULT_LFC_THRESHOLD)
        self._lfc_spin.setSingleStep(0.1)
        self._lfc_spin.setDecimals(2)
        self._lfc_spin.valueChanged.connect(self._update_gene_count)
        genes_lay.addWidget(self._lfc_spin)

        # Restricting to one direction before the test asks a different
        # question from reading the Direction column afterwards: up-only
        # ORA tests "what is upregulated", where the combined run tests
        # "what changed" and then reports which way most of it went.
        genes_lay.addWidget(QLabel("Direction:"))
        self._direction_combo = NoScrollComboBox()
        self._direction_combo.addItem("Up and down together", "both")
        self._direction_combo.addItem("Upregulated only", "up")
        self._direction_combo.addItem("Downregulated only", "down")
        self._direction_combo.currentIndexChanged.connect(
            self._update_gene_count)
        genes_lay.addWidget(self._direction_combo)

        self._gene_count_label = SecondaryLabel("")
        self._gene_count_label.setWordWrap(True)
        genes_lay.addWidget(self._gene_count_label)
        genes_lay.addStretch()

        # -- Method: one decision chain, so one card ------------------
        method_lay = self._card_pages['method'][1]
        method_lay.addWidget(HintLabel(
            "The method decides which libraries are available, and the "
            "library decides whether an algorithm applies -- so they "
            "belong together."))

        method_lay.addWidget(QLabel("Method:"))
        self._method_combo = NoScrollComboBox()
        for display, method_id in _METHODS:
            self._method_combo.addItem(display, method_id)
        self._method_combo.setToolTip(
            "ORA: over-representation of the significant-gene list.\n"
            "Preranked GSEA: threshold-free, ranks all genes (unbiased).")
        self._method_combo.currentIndexChanged.connect(self._on_method_changed)
        method_lay.addWidget(self._method_combo)

        method_lay.addWidget(QLabel("Library:"))
        self._library_combo = NoScrollComboBox()
        self._library_combo.currentIndexChanged.connect(self._on_library_changed)
        method_lay.addWidget(self._library_combo)
        self._populate_libraries(gsea_only=False)

        self._algo_label = QLabel("Algorithm:")
        method_lay.addWidget(self._algo_label)
        self._algo_combo = NoScrollComboBox()
        self._algo_combo.addItem("elim (DAG-aware)", "elim")
        self._algo_combo.addItem("Fisher + BH", "fisher")
        self._algo_combo.setToolTip(
            "elim: DAG-aware, removes redundant parent terms (recommended).\n"
            "Fisher + BH: standard over-representation with FDR correction.")
        method_lay.addWidget(self._algo_combo)
        method_lay.addStretch()

        # -- Term filters ---------------------------------------------
        terms_lay = self._card_pages['terms'][1]
        terms_lay.addWidget(HintLabel(
            "Applied to the pathway terms that come back, not to the "
            "genes going in."))

        terms_lay.addWidget(QLabel("P threshold:"))
        self._p_thresh = NoScrollDoubleSpinBox()
        self._p_thresh.setRange(0.0001, 1.0)
        self._p_thresh.setValue(DEFAULT_FDR)
        self._p_thresh.setSingleStep(0.01)
        self._p_thresh.setDecimals(4)
        self._p_thresh.setToolTip(
            "FDR cutoff for calling a pathway term significant.")
        terms_lay.addWidget(self._p_thresh)

        terms_lay.addWidget(QLabel("Min genes per term:"))
        self._min_genes = NoScrollSpinBox()
        self._min_genes.setRange(2, 500)
        self._min_genes.setValue(5)
        self._min_genes.setToolTip(
            "Skip terms with fewer than this many genes in the library. "
            "Very small terms give unstable p-values.")
        terms_lay.addWidget(self._min_genes)
        terms_lay.addStretch()

        self._on_library_changed()

        for key, title, icon in (
                ('genes', "Genes tested", "dataset"),
                ('method', "Method", "bar-chart-2"),
                ('terms', "Term filters", "filter-funnel")):
            page, _ = self._card_pages[key]
            card = StageSummaryCard(title, icon=icon)
            card.set_settings_widget(page)
            left_layout.addWidget(card)
            self._accordion.add_card(key, card)
            setattr(self, f"_{key}_card", card)

        left_layout.addStretch()
        self._accordion.finalize()

        # Actions below the accordion: collapsing a card must never hide
        # Run, which is exactly what the old Settings group did.
        self._run_btn = PrimaryButton("Run Enrichment")
        self._run_btn.setMinimumHeight(32)
        self._run_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._run_btn.clicked.connect(self._run_enrichment)
        left_layout.addWidget(self._run_btn)

        self._export_btn = QPushButton("Export CSV...")
        self._export_btn.setEnabled(False)
        self._export_btn.setToolTip(
            "Save a copy of the results table wherever you like. A GSEA "
            "run is already written into the project for Figure Export; "
            "this is for taking a copy elsewhere.")
        self._export_btn.clicked.connect(self._export_csv)
        left_layout.addWidget(self._export_btn)

        self._status_label = SecondaryLabel("")
        self._status_label.setWordWrap(True)
        left_layout.addWidget(self._status_label)

        for _sig in (self._method_combo.currentIndexChanged,
                     self._library_combo.currentIndexChanged,
                     self._direction_combo.currentIndexChanged,
                     self._algo_combo.currentIndexChanged,
                     self._fdr_spin.valueChanged,
                     self._lfc_spin.valueChanged,
                     self._p_thresh.valueChanged,
                     self._min_genes.valueChanged,
                     self._pct_filter.valueChanged):
            _sig.connect(self._refresh_card_summaries)
        self._pct_filter_check.toggled.connect(self._refresh_card_summaries)
        self._refresh_card_summaries()

        # -- Right content area: tabs (provided by SidebarTabbedPage) --
        self._tabs = self.tabs

        # Results table (with Direction column)
        self._table = ResultsTable()
        self._table.setAlternatingRowColors(True)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSortingEnabled(True)
        hdr = self._table.horizontalHeader()
        hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for col in range(1, 8):
            hdr.setSectionResizeMode(col, QHeaderView.ResizeMode.ResizeToContents)
        hdr.setSectionResizeMode(8, QHeaderView.ResizeMode.Stretch)
        self._tabs.addTab(self._table, "Results Table")

        # Bar plot tab (native pyqtgraph; ORA and GSEA share the widget)
        self._plot_widget = InteractiveEnrichmentBar()
        self._tabs.addTab(self._plot_widget, "Bar Plot")

        # Footer note below the tab widget.
        dir_note = HintLabel(
            "Direction shows the majority direction of overlapping DE genes, "
            "not a statistical test of pathway-level activity. "
            "Use Hypothesis-driven mode for pathway-level DE testing."
        )
        dir_note.setWordWrap(True)
        self.content_layout.addWidget(dir_note)

    def _on_library_changed(self):
        """Show/hide GO-specific algorithm selector."""
        lib_id = self._library_combo.currentData()
        is_go = lib_id is not None and lib_id in ("GO_BP", "GO_MF", "GO_CC")
        self._algo_label.setVisible(is_go)
        self._algo_combo.setVisible(is_go)

    def _populate_libraries(self, gsea_only):
        """Fill the library combo. Preranked GSEA needs fetchable gene-set
        libraries, so the local-GO methods are offered for ORA only."""
        self._library_combo.blockSignals(True)
        self._library_combo.clear()
        if not gsea_only:
            for display, lib_id in _GO_LIBRARIES:
                self._library_combo.addItem(f"GO: {display}", lib_id)
        for display, lib_id in _ENRICHR_LIBRARIES:
            self._library_combo.addItem(f"Enrichr: {display}", lib_id)
        self._library_combo.blockSignals(False)

    def _on_method_changed(self):
        """Switch between ORA and preranked GSEA."""
        is_gsea = self._method_combo.currentData() == "gsea"
        self._populate_libraries(gsea_only=is_gsea)
        # GSEA ranks every gene, so the significant-gene cutoffs don't apply
        # (the pct-cells filter stays -- it scopes the ranking universe).
        self._fdr_spin.setEnabled(not is_gsea)
        self._lfc_spin.setEnabled(not is_gsea)
        self._run_btn.setText("Run GSEA" if is_gsea else "Run Enrichment")
        self._on_library_changed()

    def _scoped_de_results(self):
        """
        Return '(scoped_df, total_tested, fdr_scope_size)'.

        Applies the pct-cells pre-filter when active, then recomputes
        BH-FDR within the surviving subset (Independent
        Filtering). The returned dataframe has its 'pvals_adj'
        column scoped to whatever the user has filtered to. Used by
        both the live gene-count label and '_get_sig_genes'.
        """
        de = self.ws.de_results
        if de is None:
            return None, 0, 0

        full_total = len(de)
        df = de.copy()

        pct_active = self._pct_filter_check.isChecked()
        pct_threshold = self._pct_filter.value() / 100.0 if pct_active else 0.0

        if (pct_threshold > 0
                and 'pct_disease' in df.columns
                and 'pct_control' in df.columns):
            keep = (
                (df['pct_disease'].fillna(0) >= pct_threshold)
                | (df['pct_control'].fillna(0) >= pct_threshold)
            )
            df = df[keep].copy()
            if 'pvals' in df.columns and not df.empty:
                from kosmic.numerical import bh_fdr
                df['pvals_adj'] = bh_fdr(df['pvals'].astype(float).values)

        return df, full_total, len(df)

    def _refresh_card_summaries(self, *_):
        """Two short lines per collapsed card.

        Collapsed cards are a fixed height, so a line that wraps is a
        line that gets cut off. Everything here stays terse, and
        anything at its default is left out rather than spelled out.
        """
        if not hasattr(self, '_terms_card'):
            return  # still building the sidebar
        is_gsea = self._method_combo.currentData() == "gsea"

        # Genes tested
        if is_gsea:
            headline, detail = "all genes ranked", "threshold-free"
        else:
            count = (self._gene_count_label.text() or "").splitlines()
            headline = count[0] if count and count[0] else "no genes yet"
            detail = (f"FDR<{self._fdr_spin.value():g}, "
                      f"|LFC|>{self._lfc_spin.value():g}")
            direction = self._direction_combo.currentData()
            if direction != "both":
                detail += f", {direction} only"
        self._genes_card.set_summary([headline, detail])

        # Method: the short name, then the library.
        label = "Preranked GSEA" if is_gsea else "ORA"
        library = self._library_combo.currentText() or "no library"
        self._method_card.set_summary([label, library.split(": ")[-1]])

        # Term filters
        terms = f"FDR<{self._p_thresh.value():g}, min {self._min_genes.value()} genes"
        # Not isVisible(): that is False whenever the page itself is not
        # shown, so the algorithm silently vanished from the summary.
        # The library is what decides whether an algorithm applies.
        algo = ("GO " + self._algo_combo.currentData()
                if self._library_combo.currentData() in ("GO_BP", "GO_MF", "GO_CC")
                else "")
        self._terms_card.set_summary([terms, algo])

    def _update_gene_count(self):
        """Update the live count of significant genes with direction breakdown."""
        df, full_total, fdr_scope = self._scoped_de_results()
        if df is None:
            self._gene_count_label.setText("No DE results available.")
            return

        fdr = self._fdr_spin.value()
        lfc = self._lfc_spin.value()
        sig = (df['pvals_adj'] < fdr) & (df['logfoldchanges'].abs() > lfc)
        n_up = (sig & (df['logfoldchanges'] > 0)).sum()
        n_down = (sig & (df['logfoldchanges'] < 0)).sum()

        scope_note = ""
        if fdr_scope < full_total:
            scope_note = (f"\nFDR scope: {fdr_scope:,} genes "
                          f"(filtered from {full_total:,})")
        self._gene_count_label.setText(
            f"{n_up + n_down:,} significant"
            f"{scope_note}\n"
            f"{n_up:,} up, {n_down:,} down"
        )
        self._refresh_card_summaries()

    def on_activated(self):
        """Called when page becomes visible."""
        self._update_gene_count()

    def _get_sig_genes(self):
        """
        Return (study_genes, background_genes, lfc_map) or (None, None, None).

        The background is the *post-pct-filter* gene set so the ORA's
        denominator matches the FDR scope. Without that, we'd be
        testing against a 30k-gene background, which inflates pathway p-values.
        """
        df, _, _ = self._scoped_de_results()
        if df is None or df.empty:
            return None, None, None

        fdr = self._fdr_spin.value()
        lfc = self._lfc_spin.value()
        sig = (df['pvals_adj'] < fdr) & (df['logfoldchanges'].abs() > lfc)

        direction = self._direction_combo.currentData()
        if direction == 'up':
            sig &= df['logfoldchanges'] > 0
        elif direction == 'down':
            sig &= df['logfoldchanges'] < 0

        study = df.loc[sig, 'names'].tolist()
        background = df['names'].tolist()

        # Build gene -> log2FC map for post-hoc direction annotation
        lfc_map = dict(zip(df['names'], df['logfoldchanges']))

        return study, background, lfc_map

    # ------------------------------------------------------------------
    # Run enrichment (single pass on all significant genes)
    # ------------------------------------------------------------------

    def _run_enrichment(self):
        if self._worker is not None and self._worker.isRunning():
            return

        if self._method_combo.currentData() == "gsea":
            self._run_gsea()
            return

        study, background, lfc_map = self._get_sig_genes()
        if not study:
            self._status_label.setText(
                "No significant genes at current thresholds. "
                "Adjust cutoffs or run Gene DE first."
            )
            return

        self._lfc_map = lfc_map
        self._run_btn.setEnabled(False)
        self._status_label.setText(
            f"Running enrichment on {len(study):,} significant genes..."
        )

        lib_id = self._library_combo.currentData()
        if not lib_id:
            return

        if lib_id in ("GO_BP", "GO_MF", "GO_CC"):
            ontology = lib_id.split("_")[1]
            from kosmic.gui.de_analysis.workers import GOEnrichmentWorker
            self._worker = GOEnrichmentWorker(
                study_genes=study,
                background_genes=background,
                ontology=ontology,
                method=self._algo_combo.currentData(),
                alpha=self._p_thresh.value(),
                min_genes=self._min_genes.value(),
            )
            on_finished_cb = self._on_go_finished
        else:
            from kosmic.gui.de_analysis.workers import EnrichrORAWorker
            self._worker = EnrichrORAWorker(
                study_genes=study,
                background_genes=background,
                library_name=lib_id,
                alpha=self._p_thresh.value(),
                min_genes=self._min_genes.value(),
            )
            on_finished_cb = self._on_enrichr_finished

        if self.ws.progress_bar:
            self.ws.progress_bar.setRange(0, 100)
            self.ws.progress_bar.setValue(0)
        run_worker(
            self._worker,
            on_finished=on_finished_cb,
            on_failed=self._on_enrichment_failed,
            on_progress=self._on_worker_progress,
            on_progress_pct=self._on_worker_pct,
        )

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------

    def _on_worker_progress(self, message):
        self._status_label.setText(message)
        self.ws.log_message.emit(message)

    def _on_worker_pct(self, pct):
        if self.ws.progress_bar:
            self.ws.progress_bar.setValue(pct)

    def _on_go_finished(self, payload):
        results_df, message = payload
        self._run_btn.setEnabled(True)
        if self.ws.progress_bar:
            self.ws.progress_bar.setValue(0)
        self._status_label.setText(message)
        self.ws.log_message.emit(message)
        self.ws.status_message.emit(message)
        norm = self._normalise_go(results_df)
        self._annotate_direction_and_display(norm)

    def _on_enrichr_finished(self, payload):
        results_df, message = payload
        self._run_btn.setEnabled(True)
        if self.ws.progress_bar:
            self.ws.progress_bar.setValue(0)
        self._status_label.setText(message)
        self.ws.log_message.emit(message)
        self.ws.status_message.emit(message)
        norm = self._normalise_enrichr(results_df)
        self._annotate_direction_and_display(norm)

    def _on_enrichment_failed(self, message: str):
        self._run_btn.setEnabled(True)
        if self.ws.progress_bar:
            self.ws.progress_bar.setValue(0)
        self._status_label.setText(f"Failed: {message}")
        self.ws.log_message.emit(f"Enrichment failed: {message}")

    # ------------------------------------------------------------------
    # Preranked GSEA (fgsea) against a library
    # ------------------------------------------------------------------

    def _run_gsea(self):
        # Rank the genome-wide DE. Prefer the retained genome-wide ranking
        # (present even in hypothesis mode) so the null is the whole
        # transcriptome; fall back to the current table in discovery mode.
        gene_de = getattr(self.ws, 'gene_de_genome_wide', None)
        if gene_de is None:
            gene_de, _, _ = self._scoped_de_results()
        if gene_de is None or len(gene_de) == 0 or 'logfoldchanges' not in gene_de.columns:
            self._status_label.setText(
                "No usable DE ranking. Run Gene DE first (discovery mode gives "
                "the unbiased genome-wide ranking).")
            return

        lib_id = self._library_combo.currentData()
        if not lib_id:
            return

        self._lfc_map = dict(zip(gene_de['names'], gene_de['logfoldchanges']))
        self.ws.gsea_library_name = self._library_combo.currentText()
        self._run_btn.setEnabled(False)
        self._status_label.setText(f"Running preranked GSEA ({lib_id})...")

        # Which column ranks the genes. DE ranks by log2FC; the
        # meta-analysis host points this at the pooled z-statistic,
        # which divides the effect by its pooled SE and so already
        # accounts for between-study heterogeneity. The worker always
        # receives the column named 'logfoldchanges'.
        rank_col = getattr(self, 'rank_column', 'logfoldchanges')
        if rank_col not in gene_de.columns:
            rank_col = 'logfoldchanges'
        rank_df = gene_de[['names', rank_col]].copy()
        rank_df.columns = ['names', 'logfoldchanges']

        from kosmic.gui.de_analysis.workers import LibraryGseaWorker
        self._worker = LibraryGseaWorker(rank_df, lib_id)
        if self.ws.progress_bar:
            self.ws.progress_bar.setRange(0, 0)  # indeterminate
        run_worker(
            self._worker,
            on_finished=self._on_gsea_finished,
            on_failed=self._on_enrichment_failed,
            on_progress=self._on_worker_progress,
        )

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

    def _show_gsea_results(self, res):
        """Populate the shared table + bar plot with the fgsea (NES) results."""
        import pandas as pd

        if res is None or res.empty:
            self._status_label.setText("GSEA returned no results.")
            self._results_df = pd.DataFrame()
            return

        self._results_df = res
        self.ws.gsea_library_results = res

        # Persist so Figure Export reloads the exact run instead of recomputing.
        if getattr(self.ws, 'project_dir', None):
            from kosmic.paths import de_pathway_scoring_dir
            from kosmic.de.fgsea import save_fgsea_results, GSEA_LIBRARY_RESULTS_FILE
            save_fgsea_results(
                res, de_pathway_scoring_dir(self.ws.project_dir)
                / GSEA_LIBRARY_RESULTS_FILE)

        alpha = self._p_thresh.value()
        down = QBrush(QColor(get_color('plot_downregulated')))
        up = QBrush(QColor(get_color('plot_upregulated')))

        display = res.copy()
        display['_le'] = display['leading_edge'].astype(str)
        display['_le_n'] = display['_le'].str.split(';').apply(
            lambda gs: len([g for g in gs if g.strip()]))
        display['_le_disp'] = display['_le'].str.slice(0, 200)

        cols = [
            Column('Gene set',           'names',     's'),
            Column('NES',                'nes',       '.2f'),
            Column('ES',                 'es',        '.2f'),
            Column('P-value',            'pvals',     '.2e'),
            Column('FDR',                'pvals_adj', '.2e'),
            Column('Leading edge',       '_le_n',     'd'),
            Column('Leading-edge genes', '_le_disp',  's'),
        ]

        def _decorate(item, row, col_i):
            key = cols[col_i].key
            if key == 'pvals_adj' and pd.notna(row['pvals_adj']) and row['pvals_adj'] < alpha:
                f = item.font()
                f.setBold(True)
                item.setFont(f)
            elif key == 'nes' and pd.notna(row['nes']):
                item.setForeground(down if row['nes'] < 0 else up)
            elif key == '_le_disp':
                item.setToolTip(row['_le'])

        self._table.set_schema(cols)
        self._table.set_data(display, decorate=_decorate)
        self._export_btn.setEnabled(True)
        self.enrichment_complete.emit()
        self._render_gsea_bar(res, alpha)

    def _render_gsea_bar(self, res, alpha):
        """Diverging NES bar of the top up/down gene sets (native renderer)."""
        try:
            self._plot_widget.set_gsea(res, alpha)
        except Exception as e:
            self._plot_widget.show_unavailable_message(f"Plot rendering failed: {e}")

    def _annotate_direction_and_display(self, df):
        """Add direction based on majority vote of overlapping gene log2FCs."""
        import pandas as pd

        if df is None or df.empty:
            self._status_label.setText("No enriched terms found.")
            self._results_df = pd.DataFrame()
            return

        lfc_map = getattr(self, '_lfc_map', {})

        directions = []
        pct_up_list = []
        for _, row in df.iterrows():
            genes_str = row.get('Genes', '')
            if isinstance(genes_str, list):
                genes = genes_str
            else:
                genes = [g.strip() for g in str(genes_str).split(',') if g.strip()]

            if not genes or not lfc_map:
                directions.append('')
                pct_up_list.append(50)
                continue

            n_up = sum(1 for g in genes if lfc_map.get(g, 0) > 0)
            n_down = sum(1 for g in genes if lfc_map.get(g, 0) < 0)
            total = n_up + n_down

            if total == 0:
                directions.append('')
                pct_up_list.append(50)
            elif n_up > n_down:
                directions.append(f'{n_up}/{total} up')
                pct_up_list.append(n_up / total * 100)
            elif n_down > n_up:
                directions.append(f'{n_down}/{total} down')
                pct_up_list.append((1 - n_down / total) * 100)
            else:
                directions.append(f'{n_up}/{total} mixed')
                pct_up_list.append(50)

        df = df.copy()
        df.insert(0, 'Direction', directions)

        alpha = self._p_thresh.value()
        n_sig = (df['FDR'] < alpha).sum()
        self._status_label.setText(
            f"{len(df)} terms, {n_sig} significant (FDR < {alpha})"
        )

        self._show_results(df)

    # ------------------------------------------------------------------
    # Normalise results to common format
    # ------------------------------------------------------------------

    def _normalise_go(self, df):
        import pandas as pd
        if df is None or df.empty:
            return pd.DataFrame()

        norm = pd.DataFrame()
        norm['Term'] = df.get('GO_ID', '').astype(str) + ' ' + df.get('Term', '').astype(str)
        norm['P_value'] = df['P_value']
        norm['FDR'] = df.get('P_adjusted', df['P_value'])
        norm['Fold_Enrichment'] = df.get('Fold_Enrichment', 1.0)
        norm['Gene_Count'] = df.get('Genes', df.get('Gene_Count', '')).apply(
            lambda x: len(x) if isinstance(x, list) else 0
        )
        norm['Term_Size'] = df.get('Term_Size', norm['Gene_Count'])
        norm['Genes'] = df.get('Genes', '').apply(
            lambda x: ', '.join(x) if isinstance(x, list) else str(x)
        )
        return norm

    def _normalise_enrichr(self, df):
        import pandas as pd
        if df is None or df.empty:
            return pd.DataFrame()

        norm = pd.DataFrame()
        norm['Term'] = df['Term']
        norm['P_value'] = df['P_value']
        norm['FDR'] = df['FDR']
        norm['Fold_Enrichment'] = df['Fold_Enrichment']
        norm['Gene_Count'] = df['Gene_Count']
        norm['Term_Size'] = df.get('Term_Size', df['Gene_Count'])
        norm['Genes'] = df['Genes'].apply(
            lambda x: ', '.join(x) if isinstance(x, list) else str(x)
        )
        return norm

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def _show_results(self, df):
        """Populate the results table and bar plot from combined DataFrame."""
        self._results_df = df
        self.ws.enrichment_results = df

        alpha = self._p_thresh.value()
        lfc_map = getattr(self, '_lfc_map', {})
        up_color = QBrush(QColor(get_color('plot_upregulated')))
        down_color = QBrush(QColor(get_color('plot_downregulated')))

        # Derive the per-row display columns once, as DataFrame columns,
        # so the widget schema stays declarative.
        display = df.copy()
        display['_genes_str'] = display.get('Genes', '').astype(str)
        display['_gene_list'] = display['_genes_str'].str.split(',').apply(
            lambda gs: [g.strip() for g in gs if g.strip()])
        display['_n_overlap'] = display['_gene_list'].str.len()
        display['_n_up'] = display['_gene_list'].apply(
            lambda gs: sum(1 for g in gs if lfc_map.get(g, 0) > 0))
        display['_n_down'] = display['_gene_list'].apply(
            lambda gs: sum(1 for g in gs if lfc_map.get(g, 0) < 0))
        display['_genes_display'] = display['_genes_str'].str.slice(0, 200)
        if 'Term_Size' not in display.columns:
            display['Term_Size'] = display['_n_overlap']

        cols = [
            Column('Term',         'Term',            's'),
            Column('P-value',      'P_value',         '.2e'),
            Column('FDR',          'FDR',             '.2e'),
            Column('Fold Enrich.', 'Fold_Enrichment', '.2f'),
            Column('DEG',          '_n_overlap',      'd'),
            Column('Up',           '_n_up',           'd'),
            Column('Down',         '_n_down',         'd'),
            Column('Term Size',    'Term_Size',       'd'),
            Column('Genes',        '_genes_display',  's'),
        ]

        def _decorate(item, row, col_i):
            key = cols[col_i].key
            if key == 'FDR' and row['FDR'] < alpha:
                f = item.font()
                f.setBold(True)
                item.setFont(f)
            elif key == '_n_up' and row['_n_up'] > row['_n_down']:
                item.setForeground(up_color)
            elif key == '_n_down' and row['_n_down'] > row['_n_up']:
                item.setForeground(down_color)
            elif key == '_genes_display':
                item.setToolTip(row['_genes_str'])

        self._table.set_schema(cols)
        self._table.set_data(display, decorate=_decorate)
        self._export_btn.setEnabled(True)
        self.enrichment_complete.emit()

        self._render_bar_plot(df, alpha)

    def _render_bar_plot(self, df, alpha):
        """Render -log10(P) bars per term, coloured by direction (native renderer)."""
        try:
            self._plot_widget.set_ora(df, alpha)
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
