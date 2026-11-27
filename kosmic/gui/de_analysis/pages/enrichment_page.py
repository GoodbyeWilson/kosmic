# Enrichment Page (DE Discovery mode): over-representation analysis (ORA)
# of the significant DE genes, or preranked GSEA of every gene, against a
# gene-set library (bundled MSigDB Hallmark, or one downloaded from Enrichr,
# which includes the GO libraries).
#
# In the DE workspace the page works on the gene DE run chosen in the run
# bar (ADR-016): it saves its results in that run's folder and, when
# another run is chosen, clears what it shows and shows that run's saved
# results. The Meta-Analysis workspace hosts the same page on the pooled
# gene list and has no run bar.

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

# Enrichment methods. ORA tests the significant-gene list against a background;
# preranked GSEA (fgsea) is threshold-free, ranking every gene -- the same
# rank-based method used for the a-priori panel, run here against a library.
_METHODS = [
    ("Over-representation (ORA)", "ora"),
    ("Preranked GSEA (fgsea)", "gsea"),
]


def _now() -> str:
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


class EnrichmentPage(SidebarTabbedPage):
    """
    Discovery mode endpoint: pathway enrichment of significant DE genes.

    Runs a single ORA pass over all significant genes, then annotates
    each term with the majority direction of its overlapping DE genes
    (the Direction column). It does NOT run up and down separately.
    """

    help_id = "de/enrichment"

    # Column used to rank genes for preranked GSEA when the page offers no
    # choice. Overridden by the meta-analysis host (the pooled z).
    rank_column = 'logfoldchanges'

    # The rankings a preranked GSEA can use: (label, DE table column). The
    # first is the default. The Wald statistic (log2FC / SE; Welch's t for
    # the t-test methods) is what preranked GSEA of DESeq2 results usually
    # ranks by: a noisy fold change of a lowly expressed gene does not
    # reach the extremes of the ranking. Empty: no choice, 'rank_column'.
    RANK_CHOICES = (
        ("Test statistic (DESeq2: Wald, log2FC / SE)", 'statistic'),
        ("Log2 fold change", 'logfoldchanges'),
    )

    enrichment_complete = pyqtSignal()

    def __init__(self, workspace):
        super().__init__()
        self.ws = workspace
        self._worker = None
        self._results_df = None
        self._lfc_map = {}  # gene -> log2FC for direction annotation
        self._setup_ui()
        run_changed = getattr(self.ws, 'run_changed', None)
        if run_changed is not None:
            run_changed.connect(self._on_run_changed)

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

        # Preranked GSEA only: what orders the genes.
        self._rank_label = QLabel("Rank genes by:")
        self._rank_combo = NoScrollComboBox()
        for display, column in self.RANK_CHOICES:
            self._rank_combo.addItem(display, column)
        self._rank_combo.setToolTip(
            "Log2 fold change: the effect size alone.\n"
            "Test statistic: the fold change divided by its standard error, so\n"
            "genes measured with little precision (few counts) do not reach the\n"
            "ends of the ranking. Preranked GSEA of DESeq2 results usually\n"
            "ranks by it.")
        method_lay.addWidget(self._rank_label)
        method_lay.addWidget(self._rank_combo)
        self._rank_label.setVisible(False)
        self._rank_combo.setVisible(False)

        method_lay.addWidget(QLabel("Library:"))
        self._library_combo = NoScrollComboBox()
        method_lay.addWidget(self._library_combo)
        self._populate_libraries()

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

        self._min_genes_label = QLabel("Min genes per term:")
        terms_lay.addWidget(self._min_genes_label)
        self._min_genes = NoScrollSpinBox()
        self._min_genes.setRange(2, 500)
        self._min_genes.setValue(5)
        self._min_genes.setToolTip(
            "Skip terms with fewer than this many genes in the library. "
            "Very small terms give unstable p-values.")
        terms_lay.addWidget(self._min_genes)

        # Preranked GSEA only: which gene sets are tested, by how many of
        # their genes are in the ranking. The defaults are GSEA's own.
        from kosmic.de.fgsea import GSEA_MIN_SIZE, GSEA_MAX_SIZE
        size_tip = (
            "Preranked GSEA tests a gene set only when this many of its genes\n"
            f"are in the ranking. {GSEA_MIN_SIZE} and {GSEA_MAX_SIZE} are the "
            "defaults of the Broad\nInstitute's GSEA and of GSEApy: smaller sets "
            "give unstable scores,\nlarger ones are too broad to interpret.")
        self._gsea_min_label = QLabel("Min genes per set (in the ranking):")
        self._gsea_min_size = NoScrollSpinBox()
        self._gsea_min_size.setRange(1, 10000)
        self._gsea_min_size.setValue(GSEA_MIN_SIZE)
        self._gsea_min_size.setToolTip(size_tip)
        self._gsea_max_label = QLabel("Max genes per set (in the ranking):")
        self._gsea_max_size = NoScrollSpinBox()
        self._gsea_max_size.setRange(1, 100000)
        self._gsea_max_size.setValue(GSEA_MAX_SIZE)
        self._gsea_max_size.setToolTip(size_tip)
        for w in (self._gsea_min_label, self._gsea_min_size,
                  self._gsea_max_label, self._gsea_max_size):
            terms_lay.addWidget(w)
            w.setVisible(False)
        terms_lay.addStretch()

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
                     self._fdr_spin.valueChanged,
                     self._lfc_spin.valueChanged,
                     self._p_thresh.valueChanged,
                     self._min_genes.valueChanged,
                     self._gsea_min_size.valueChanged,
                     self._gsea_max_size.valueChanged,
                     self._pct_filter.valueChanged):
            _sig.connect(self._refresh_card_summaries)
        self._pct_filter_check.toggled.connect(self._refresh_card_summaries)

        # Defaults: preranked GSEA of the Hallmark sets (the first bundled
        # library), which ranks every gene and needs no significance cutoff.
        self._method_combo.setCurrentIndex(self._method_combo.findData("gsea"))
        self._on_method_changed()
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

    def _populate_libraries(self):
        """Fill the library combo: the same libraries serve ORA and GSEA."""
        self._library_combo.blockSignals(True)
        self._library_combo.clear()
        for display, lib_id in _BUNDLED_LIBRARIES:
            self._library_combo.addItem(f"MSigDB: {display}", lib_id)
        for display, lib_id in _ENRICHR_LIBRARIES:
            self._library_combo.addItem(f"Enrichr: {display}", lib_id)
        self._library_combo.blockSignals(False)

    def _rank_by(self) -> str:
        """The DE table column a preranked GSEA ranks the genes by."""
        if self.RANK_CHOICES:
            return self._rank_combo.currentData() or self.RANK_CHOICES[0][1]
        return self.rank_column

    def _on_method_changed(self):
        """Switch between ORA and preranked GSEA."""
        is_gsea = self._method_combo.currentData() == "gsea"
        offer_rank = is_gsea and bool(self.RANK_CHOICES)
        self._rank_label.setVisible(offer_rank)
        self._rank_combo.setVisible(offer_rank)
        # The term-size filter is ORA's; GSEA has its own set-size bounds.
        self._min_genes_label.setVisible(not is_gsea)
        self._min_genes.setVisible(not is_gsea)
        for w in (self._gsea_min_label, self._gsea_min_size,
                  self._gsea_max_label, self._gsea_max_size):
            w.setVisible(is_gsea)
        # GSEA ranks every gene, so the significant-gene cutoffs don't apply
        # (the pct-cells filter stays -- it scopes the ranking universe).
        self._fdr_spin.setEnabled(not is_gsea)
        self._lfc_spin.setEnabled(not is_gsea)
        self._run_btn.setText("Run GSEA" if is_gsea else "Run Enrichment")

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
        if is_gsea:
            terms = (f"FDR<{self._p_thresh.value():g}, sets of "
                     f"{self._gsea_min_size.value()}-{self._gsea_max_size.value()} genes")
        else:
            terms = f"FDR<{self._p_thresh.value():g}, min {self._min_genes.value()} genes"
        self._terms_card.set_summary([terms, ""])

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
        self._sync_significance()
        self._update_gene_count()

    def _sync_significance(self):
        """Start the FDR and |log2FC| boxes from the rule of the DE run shown
        (ADR-015), each time a different run is shown; values typed here
        stay until then."""
        rule = getattr(self.ws, 'significance', None)
        if not rule or rule == getattr(self, '_synced_rule', None):
            return
        self._synced_rule = rule
        fdr, lfc, _assumed = rule
        for spin, value in ((self._fdr_spin, fdr), (self._lfc_spin, lfc)):
            spin.blockSignals(True)
            spin.setValue(value)
            spin.blockSignals(False)

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
        self._library_source = getattr(self._worker, 'library_source', None)
        self._record_study_enrichment('ora')

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
        self._gsea_library_id = lib_id
        self.ws.gsea_library_name = self._library_combo.currentText()
        self._run_btn.setEnabled(False)
        self._status_label.setText(f"Running preranked GSEA ({lib_id})...")

        # Which column ranks the genes: the Rank genes by choice on the DE
        # page; the meta-analysis host's pooled z-statistic, which divides
        # the effect by its pooled SE and so already accounts for
        # between-study heterogeneity. The worker always receives the
        # column named 'logfoldchanges'.
        rank_col = self._rank_by()
        if rank_col not in gene_de.columns:
            if self.RANK_CHOICES and rank_col != 'logfoldchanges':
                self._status_label.setText(
                    f"This DE result has no '{rank_col}' column to rank by "
                    "(it was saved by an earlier version). Run the DE again, "
                    "or rank by log2 fold change.")
                self._run_btn.setEnabled(True)
                return
            rank_col = 'logfoldchanges'
        rank_df = gene_de[['names', rank_col]].copy()
        rank_df.columns = ['names', 'logfoldchanges']

        from kosmic.gui.de_analysis.workers import LibraryGseaWorker
        lo, hi = self._gsea_min_size.value(), self._gsea_max_size.value()
        if lo > hi:
            self._status_label.setText(
                f"The minimum set size ({lo}) is above the maximum ({hi}).")
            self._run_btn.setEnabled(True)
            return
        self._worker = LibraryGseaWorker(rank_df, lib_id, min_size=lo, max_size=hi)
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
        self._library_source = getattr(self._worker, 'library_source', None)
        self._record_study_enrichment('gsea')

    # ------------------------------------------------------------------
    # Provenance
    # ------------------------------------------------------------------

    # The adjusted p-value column of each method's results table: ORA's
    # is normalised to 'FDR', the fgsea summary keeps 'pvals_adj'.
    _ADJ_P_COLUMN = {'ora': 'FDR', 'gsea': 'pvals_adj'}

    def _enrichment_params(self, method) -> dict:
        """What one enrichment run used and found, for the provenance."""
        df = self._results_df
        n_terms = 0 if df is None else len(df)
        col = self._ADJ_P_COLUMN[method]
        n_sig = 0
        if df is not None and col in getattr(df, 'columns', []):
            n_sig = int((df[col] < self._p_thresh.value()).sum())
        study, background, _ = self._get_sig_genes()
        return {
            'method':      method,
            'library':     self._library_combo.currentText(),
            # Where the library came from: bundled file or Enrichr download,
            # its version or address, a SHA-1 of its content, set and gene
            # counts and the time it was loaded.
            'library_source': getattr(self, '_library_source', None),
            'direction':   self._direction_combo.currentData(),
            'rank_by':     self._rank_by() if method == 'gsea' else None,
            'fdr':         self._p_thresh.value(),
            'gene_fdr':    self._fdr_spin.value(),
            'gene_lfc':    self._lfc_spin.value(),
            'min_genes':   self._min_genes.value() if method != 'gsea' else None,
            'gsea_min_size': self._gsea_min_size.value() if method == 'gsea' else None,
            'gsea_max_size': self._gsea_max_size.value() if method == 'gsea' else None,
            'n_genes_in':  len(study or []),
            'n_background': len(background or []),
            'n_terms':     n_terms,
            'n_terms_sig': n_sig,
            # The gene DE run it was made on (ADR-016); the DE workspace only.
            **(self.ws.run_identity() if hasattr(self.ws, 'run_identity') else {}),
        }

    def _record_study_enrichment(self, method):
        """Record a DE-workspace enrichment run in the study's provenance.

        One entry per method and library: rerunning a library replaces its
        entry, running another library adds one. The meta-analysis host
        has no study file and records into its own sidecar instead.
        """
        h5ad = getattr(self.ws, 'h5ad_path', None)
        if not h5ad:
            return
        try:
            from pathlib import Path

            from kosmic import provenance
            study_dir = Path(h5ad).parent
            params = self._enrichment_params(method)
            de_stage = provenance.last_stage(study_dir, 'gene_de')
            provenance.record_stage(
                study_dir, Path(h5ad).stem, 'enrichment', params,
                consumed_token=(de_stage or {}).get('consumed_token'),
                replace_matching={'method': method, 'library': params['library'],
                                  'cell_type': params.get('cell_type'),
                                  'counts': params.get('counts')})
        except Exception as e:  # provenance is best-effort
            self.ws.log_message.emit(f"Could not record enrichment provenance: {e}")

    def _results_folder(self):
        """Where this page saves an over-representation run: the folder of
        the gene DE run chosen in the run bar (ADR-016), or None to not save
        it. The meta-analysis host overrides this."""
        if not getattr(self.ws, 'project_dir', None):
            return None
        run_folder = getattr(self.ws, 'run_folder', None)
        if run_folder is not None:
            return run_folder()     # the chosen run's folder (ADR-016)
        from kosmic.paths import de_pathway_scoring_dir
        return de_pathway_scoring_dir(self.ws.project_dir)

    def _gsea_results_path(self):
        """Where a library GSEA run is saved, or None to not save it.

        The chosen run's folder (ADR-016), under the one file name the
        Figure Export library GSEA page reloads; each run replaces the
        last. The meta-analysis host overrides this.
        """
        folder = self._results_folder()
        if folder is None:
            return None
        from kosmic.de.fgsea import GSEA_LIBRARY_RESULTS_FILE
        return folder / GSEA_LIBRARY_RESULTS_FILE

    def _on_run_changed(self):
        """Another gene DE run was chosen (ADR-016): clear what the page
        shows, then show that run's saved results of the chosen method."""
        self._results_df = None
        self._table.setRowCount(0)
        self._export_btn.setEnabled(False)
        self._plot_widget.clear_plot_items()
        self._plot_widget.show_unavailable_message()
        self._sync_significance()
        self._update_gene_count()
        de = self.ws.de_results
        self._lfc_map = (dict(zip(de['names'], de['logfoldchanges']))
                         if de is not None and 'logfoldchanges' in de.columns else {})
        name = self.ws.run_label() or "this run"
        folder = self.ws.run_saved_folder()
        if folder is None:
            self._status_label.setText("")
            return
        from kosmic.de import saved
        from kosmic.de.fgsea import GSEA_LIBRARY_RESULTS_FILE, load_fgsea_results
        if self._method_combo.currentData() == "gsea":
            res = load_fgsea_results(folder / GSEA_LIBRARY_RESULTS_FILE)
            if res is not None:
                details, library = saved.read_gsea_details(folder / GSEA_LIBRARY_RESULTS_FILE)
                self.ws.gsea_library_details = details or {}
                self.ws.gsea_library_name = library or None
                self._show_gsea_results(res, save=False)
                self._status_label.setText(f"Saved library GSEA of {name} shown.")
                return
        else:
            ora = saved.read_enrichment(folder)
            if ora is not None:
                self._show_results(ora)
                self._status_label.setText(f"Saved over-representation run of {name} shown.")
                return
        self._status_label.setText(f"No saved run of this method for {name}.")

    def _show_gsea_results(self, res, save=True):
        """Populate the shared table + bar plot with the fgsea (NES) results;
        with 'save', also save them for Figure Export."""
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
        # Saved so Figure Export can draw the run when this workspace no
        # longer holds it (ADR-014).
        folder = self._results_folder()
        if folder is not None:
            from kosmic.de.saved import save_enrichment
            try:
                save_enrichment(folder, df)
            except OSError as exc:
                self.ws.log_message.emit(f"Could not save the enrichment table: {exc}")

    # ------------------------------------------------------------------
    # Normalise results to common format
    # ------------------------------------------------------------------

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
        from kosmic.reference.pathways.libraries import term_label
        display['_term'] = display['Term'].map(term_label)

        cols = [
            Column('Term',         '_term',           's'),
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
            elif key == '_term':
                item.setToolTip(str(row['Term']))       # the library's own name

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
