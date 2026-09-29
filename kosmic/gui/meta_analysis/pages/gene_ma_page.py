# Consensus Meta-Analysis Page.
#
# Runs the user-selected pooling methods with optional CC-permutation
# calibration and reports genes significant in all chosen methods.

from __future__ import annotations

from pathlib import Path

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QSplitter, QGridLayout, QSpinBox, QComboBox,
    QCheckBox, QLineEdit,
    QRadioButton, QButtonGroup,
)
from PyQt6.QtCore import QSettings, Qt, pyqtSignal

import numpy as np
import pyqtgraph as pg

from kosmic.gui.shared.theme import get_color, style_pg_plot
from kosmic.gui.shared.widgets import (
    SidebarTabbedPage,
    SecondaryLabel, SectionHeader, PrimaryButton,
    StageAccordion, StageSummaryCard,
)

from kosmic.paths import meta_output_dir
from kosmic.numerical import neg_log10
from kosmic.meta_analysis.olink_panels import (
    DEFAULT_PANELS,
    annotate_consensus,
    panel_summary,
)
from kosmic import DEFAULT_FDR
from kosmic.gui.shared import dialogs, run_worker
from kosmic.gui.meta_analysis.pseudobulk_discovery import (
    clear_memory_cache as clear_pb_path_cache,
    find_pseudobulk_paths,
)
from kosmic.gui.meta_analysis.pages.gene_ma.workers import (
    GeneMAWorker,
)


# --- Method registry ---

_METHODS = [
    ('dl', 'DL'),
    ('reml', 'REML'),
    ('stouffer', 'Stouffer'),
    ('fisher', 'Fisher'),
    ('sumrank', 'SumRank'),
    ('gwop', 'gwOP'),
]

_METHOD_KEYS_BASE = {name: key for key, name in _METHODS}
_METHOD_KEYS_TOP50 = {
    'dl': 'dl_top50',
    'reml': 'reml_top50',
    'sumrank': 'sumrank_top50',
}

# Effect-size methods composed with HKSJ small-sample correction.
_HKSJ_MAP = {
    'dl': 'dl_hksj',
    'reml': 'reml_hksj',
    'dl_top50': 'dl_top50_hksj',
    'reml_top50': 'reml_top50_hksj',
}

# DE-method combo index → canonical method key used in DE filenames.
# Combo order: "DESeq2", "Welch CPM", "Welch Raw".
_DE_METHOD_COMBO_MAP = ['deseq2', 'welch_cpm', 'welch_raw']

# Allowable filename-key fallbacks per DE method. Tight aliases only --
# never silently substitute a different DE algorithm.
_DE_METHOD_ALIASES = {
    'deseq2':     ['deseq2'],
    'welch_cpm':  ['welch_cpm', 'welch_cpm_eb'],
    'welch_raw':  ['welch_raw'],
}


# Fixed left-axis width so the volcano and forest align, and so
# the floating legend can sit clear of the axis.
_AXIS_W = 50


class GeneMAPage(SidebarTabbedPage):
    """Discovery meta-analysis: genome-wide pooling, tabbed results.

    One method or several. Pooling several and taking their overlap is
    the 'consensus' option within this page -- it is not what the page
    itself is.
    """

    help_id = "meta/gene_ma"

    log_message = pyqtSignal(str)
    analysis_complete = pyqtSignal()
    # (text, color) -- forwards to workspace.data_status so the sidebar
    # shows live progress, not stale "Done: X/Y" from the prior run.
    sidebar_status = pyqtSignal(str, str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._datasets = []
        self._labels = []
        self._project_folder = None
        # Output folder of the selection (ADR-007): set with the datasets,
        # and frozen per run so Enrichment and Validation record beside it.
        self._output_selection = None
        self._run_output_selection = None
        self._pending_output_selection = None
        self._meta_df = None
        self._method_dfs = None
        self._method_keys = None
        self._worker = None
        self._last_logged_perm_decile = -1
        self._scatter_data = None
        self._highlight_scatter = None
        self._forest_xlim = None
        self._mode = 'exploratory'
        self._pathway_ma_results = None
        self._pathway_gene_sets = {}
        self._method_cache = {}
        # LOO worker reuses these so the user doesn't have to re-pick methods.
        self._last_consensus_params = None
        self._last_consensus_pb_paths = None
        self._last_consensus_de_dfs = None
        self._last_consensus_labels = None
        self._loo_result = None
        self._loo_worker = None
        self._last_selected_gene = None
        self._last_single_family = None
        # Per-mode result cache: switching modes preserves prior work.
        # Stores raw data; widgets re-render on restore.
        # (Reproducibility cache lives on the ReproducibilityTab.)
        self._results_by_mode = {
            'exploratory': None, 'hypothesis': None, 'comparison': None,
        }
        self.progress_bar = None
        self._setup_ui()

    def set_datasets(self, datasets, labels, project_folder,
                     pathway_ma_results=None, pathway_gene_sets=None,
                     output_selection=None):
        """
        Called when page becomes visible.

        Parameters
        ----------
        pathway_ma_results : DataFrame, optional
            Results from the Pathway-Level MA page (used in hypothesis mode
            to identify significant pathways).
        pathway_gene_sets : dict, optional
            {pathway_name: [gene_list]} -- the gene sets that were tested
            at the pathway level.  Used to map significant pathways back
            to their constituent genes.
        output_selection : str, optional
            Folder under meta_analysis/ for this selection's results and
            methods record (ADR-007).
        """
        if (list(labels or []) != list(self._labels)
                or project_folder != self._project_folder):
            # Another selection: results pooled from the previous one (and
            # the per-method and per-mode caches of them) no longer describe
            # it, and would otherwise stay on screen under the new cell type.
            if self._meta_df is not None:
                self._clear_results()
            self._method_cache = {}
            self._results_by_mode = dict.fromkeys(self._results_by_mode)
            self._run_output_selection = None
        self._output_selection = output_selection
        self._datasets = datasets or []
        self._labels = labels or []
        self._project_folder = project_folder
        self._sync_min_studies()
        self._refresh_ma_cards()
        self._pathway_ma_results = pathway_ma_results
        self._pathway_gene_sets = pathway_gene_sets or {}
        n = len(self._datasets)
        self._info_label.setText(
            f"{n} studies loaded" if n > 0 else "No studies loaded yet.")

        # Warm the pseudobulk path cache so first CC-perm run doesn't pay for the scan.
        if self._datasets:
            self.log_message.emit(
                f"Pseudobulk: warming cache for {len(self._datasets)} datasets")
            try:
                paths = find_pseudobulk_paths(
                    self._datasets, self._project_folder,
                    log_cb=self.log_message.emit)
                self.log_message.emit(
                    f"Pseudobulk: cache warm-up found "
                    f"{len(paths)}/{len(self._datasets)} paths")
            except Exception as e:
                self.log_message.emit(
                    f"Pseudobulk cache warm-up failed: {e}")

    # --- UI setup ---
    def _setup_ui(self):
        lay = self.sidebar_layout

        # Same pattern as the DE page: a summary card per concern that
        # expands in place. Six SettingsGroups stacked in one scroll made
        # it impossible to see at a glance what a run would do.
        lay.addWidget(SectionHeader("ANALYSIS SETUP"))
        self._accordion = StageAccordion(lay)

        # Built here, mounted into the cards at the end of this method.
        self._card_pages = {}
        for key in ('method', 'settings'):
            page = QWidget()
            page_lay = QVBoxLayout(page)
            page_lay.setContentsMargins(0, 4, 0, 0)
            page_lay.setSpacing(4)
            self._card_pages[key] = (page, page_lay)

        self._info_label = SecondaryLabel("No studies loaded yet.")
        self._mode_label = SecondaryLabel("Mode: Exploratory")

        # --- Method Mode ---
        from kosmic.gui.meta_analysis.dialogs.ma_settings_dialog import (
            METHOD_FAMILIES)

        # The cards already draw a frame and carry a title, so the
        # contents go in bare. A SettingsGroup here nested a titled box
        # inside a titled box and repeated the heading.
        mb_lay = QHBoxLayout()
        mb_lay.setContentsMargins(0, 0, 0, 0)
        self._single_radio = QRadioButton("Single")
        self._consensus_radio = QRadioButton("Consensus")
        self._single_radio.setChecked(True)
        mode_group = QButtonGroup(self)
        mode_group.addButton(self._single_radio)
        mode_group.addButton(self._consensus_radio)
        self._single_radio.toggled.connect(self._on_method_mode_toggled)
        self._single_radio.toggled.connect(self._sync_run_button_text)
        mb_lay.addWidget(self._single_radio)
        mb_lay.addWidget(self._consensus_radio)
        mb_lay.addStretch()
        self._card_pages['method'][1].addLayout(mb_lay)

        method_grid = QGridLayout()
        method_grid.setContentsMargins(0, 0, 0, 0)
        method_grid.setColumnStretch(0, 0)
        method_grid.setColumnStretch(1, 1)

        self._family_checks: dict = {}
        self._family_combos: dict = {}
        self._family_btn_group = QButtonGroup(self)
        self._family_btn_group.setExclusive(True)  # single mode default
        self._family_btn_group.buttonToggled.connect(
            self._on_family_selection_changed)

        for i, (fam_key, fam_label, _methods) in enumerate(METHOD_FAMILIES):
            cb = QCheckBox(f"{fam_label}:")
            cb.setChecked(i == 0)
            self._family_checks[fam_key] = cb
            self._family_btn_group.addButton(cb)
            combo = QComboBox()
            combo.setMinimumContentsLength(12)
            self._family_combos[fam_key] = combo
            method_grid.addWidget(cb, i, 0)
            method_grid.addWidget(combo, i, 1)

        self._single_combo = list(self._family_combos.values())[0]

        self._populate_method_dropdowns()
        self._card_pages['method'][1].addLayout(method_grid)

        sg = QVBoxLayout()
        sg.setContentsMargins(0, 0, 0, 0)
        sg.setSpacing(6)

        self._top50_cb = QCheckBox("Top 50% selection")
        self._top50_cb.setToolTip(
            "Pre-filter to the top half of genes by per-study p-value "
            "before pooling. Helps SumRank substantially, leaves DL "
            "essentially unchanged.")
        self._top50_cb.setChecked(False)
        sg.addWidget(self._top50_cb)

        # Kept short: the card viewport is ~274px and the qualifier ran
        # to 542px, so the text was simply cut off. The detail is in the
        # tooltip, which is where it was already spelled out anyway.
        self._hksj_cb = QCheckBox("HKSJ correction")
        self._hksj_cb.setToolTip(
            "Hartung-Knapp-Sidik-Jonkman adjustment: t-test on k-1 df "
            "+ adjusted SE. Applied to DL/REML only. Conservative; "
            "recommended by Cochrane for small k (<20) but reduces power.")
        self._hksj_cb.setChecked(False)
        sg.addWidget(self._hksj_cb)

        # --- Calibration, DE method, Min studies ---
        grid = QGridLayout()
        grid.setColumnStretch(0, 0)
        grid.setColumnStretch(1, 1)
        r = 0
        grid.addWidget(QLabel("Calibration:"), r, 0)
        self._cal_combo = QComboBox()
        self._cal_combo.addItems(["Analytical", "CC permutation"])
        grid.addWidget(self._cal_combo, r, 1)
        r += 1

        grid.addWidget(QLabel("DE method:"), r, 0)
        self._de_combo = QComboBox()
        self._de_combo.addItems(["DESeq2", "Welch CPM", "Welch Raw"])
        grid.addWidget(self._de_combo, r, 1)
        r += 1

        grid.addWidget(QLabel("Min studies:"), r, 0)
        self._min_studies_spin = QSpinBox()
        self._min_studies_spin.setRange(2, 20)
        self._min_studies_spin.setValue(3)
        self._min_studies_spin.setToolTip(
            "How many studies must report a gene for it to be pooled."
            + "\n\n" +
            "A gene needs a finite effect size AND standard error from"
            " this many studies; below that it is dropped entirely."
            + "\n\n" +
            "Capped at the number of studies loaded, because a threshold"
            " above that would silently discard every gene.")
        grid.addWidget(self._min_studies_spin, r, 1)

        sg.addLayout(grid)

        # --- Independent filter (Bourgon-style cross-study pct filter) ---
        from kosmic import META_IF_DEFAULT_PCT, META_IF_ON_BY_DEFAULT
        # The threshold belongs to the checkbox, so it is indented
        # under it and greys out with it. On separate rows with its
        # own label it read as a second, independent setting.
        self._if_check = QCheckBox("Independent filter")
        self._if_check.setToolTip(
            "Before pooling, drop genes whose detection rate averaged "
            "across studies falls below the threshold, then recompute "
            "BH FDR on what is left.\n\n"
            "Off by default. It shrinks the pooled gene universe, so "
            "the meta arm stops testing the same genes as a pooled "
            "~study + condition run -- on the DCM EC data, 8,691 "
            "genes against 11,949.")
        self._if_check.setChecked(bool(META_IF_ON_BY_DEFAULT))

        self._if_pct_spin = QSpinBox()
        self._if_pct_spin.setRange(0, 50)
        self._if_pct_spin.setSuffix(" %")
        self._if_pct_spin.setToolTip(self._if_check.toolTip())
        self._if_pct_spin.setValue(int(round(float(META_IF_DEFAULT_PCT) * 100)))
        self._if_pct_spin.setMinimumWidth(64)
        self._if_pct_spin.setMaximumWidth(96)
        # The threshold only means anything when the filter is on.
        self._if_pct_spin.setEnabled(self._if_check.isChecked())
        self._if_check.toggled.connect(self._if_pct_spin.setEnabled)

        # One row, because it is one setting. On separate rows with
        # its own label the threshold read as a second control.
        if_row = QHBoxLayout()
        if_row.setContentsMargins(0, 0, 0, 0)
        if_row.setSpacing(6)
        if_row.addWidget(self._if_check, 1)
        if_row.addWidget(self._if_pct_spin)
        sg.addLayout(if_row)

        self._card_pages['settings'][1].addLayout(sg)

        # --- Run ---
        self._run_btn = PrimaryButton("Run Meta-Analysis")
        self._run_btn.setMinimumHeight(32)
        self._run_btn.clicked.connect(self._run_consensus)

        self._status_label = SecondaryLabel("")
        self._status_label.setWordWrap(True)
        lay.addWidget(self._status_label)


        # --- Post-Consensus Analyses ---
        # Buttons for analyses that run after a consensus; results render in their own tabs.
        # These are graph controls, so they live in the strip under the
        # plots (built with the tabs, below). The group that used to hold
        # them claimed they had "moved to the View menu"; they had not --
        # nothing added them anywhere, so they were unreachable, and once
        # the group lost its last parent Qt deleted the combo underneath
        # self._yaxis_combo and every volcano draw died on it.
        self._yaxis_combo = QComboBox()
        self._yaxis_combo.addItems(["-log10(FDR)", "Z-statistic"])
        self._yaxis_combo.setCurrentIndex(1)
        self._yaxis_combo.currentIndexChanged.connect(self._redraw_volcano)

        # Significance toggles populate _filter_menu; _sig_toggles maps level -> action.
        self._sig_toggles = []

        # Olink panel rings overlay (don't replace) the significance shading.
        self._olink_cvd_check = QCheckBox("Olink CVD III")
        self._olink_cvd_check.setChecked(False)
        self._olink_cvd_check.setToolTip(
            "Highlight consensus genes that are on the Olink "
            "Cardiovascular III panel (91 proteins, routine clinical "
            "cardiac plasma assay).")
        self._olink_cvd_check.toggled.connect(self._redraw_volcano)

        self._olink_explore_check = QCheckBox("Olink Explore 3072")
        self._olink_explore_check.setChecked(False)
        self._olink_explore_check.setToolTip(
            "Highlight consensus genes that are on any Olink Explore "
            "sub-panel (~2925 unique proteins, research-grade plasma "
            "assay). The panel is large -- expect many highlights; "
            "untick to declutter.")
        self._olink_explore_check.toggled.connect(self._redraw_volcano)


        self._gene_search = QLineEdit()
        self._gene_search.setPlaceholderText("Type gene name...")
        self._gene_search.textChanged.connect(self._on_gene_search)
        self._gene_search.returnPressed.connect(self._on_gene_search_entered)

        from PyQt6.QtWidgets import QCompleter
        from PyQt6.QtCore import QStringListModel
        self._gene_completer_model = QStringListModel()
        self._gene_completer = QCompleter(self._gene_completer_model, self)
        self._gene_completer.setCaseSensitivity(
            Qt.CaseSensitivity.CaseInsensitive)
        self._gene_completer.setFilterMode(
            Qt.MatchFlag.MatchContains)
        self._gene_completer.setMaxVisibleItems(12)
        self._gene_search.setCompleter(self._gene_completer)

        # Gene search and the selected-gene readout are not specific to
        # one view -- they matter as much while reading the table as
        # while reading the volcano -- so they sit in a strip above the
        # tabs rather than in the sidebar or inside one tab.
        self._gene_info_label = SecondaryLabel("")
        self._gene_info_label.setWordWrap(False)
        self._gene_info_label.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse)

        # Hypothesis mode: pathways are picked on the PathwayExplorer page
        # and passed in via set_hypothesis_pathways().
        self._hyp_pathways = {}  # {pathway: [genes]}

        # --- the cards, in the order a run is set up ---
        self._dataset_card = StageSummaryCard("Studies", icon="dataset")
        self._dataset_card.set_read_only(True)
        self._dataset_card.set_summary(["No studies loaded yet", ""])
        lay.addWidget(self._dataset_card)

        for key, title, icon in (
                ('method', "Pooling method", "bar-chart-2"),
                ('settings', "Settings", "filter-funnel")):
            page, page_lay = self._card_pages[key]
            page_lay.addStretch()
            card = StageSummaryCard(title, icon=icon)
            card.set_settings_widget(page)
            lay.addWidget(card)
            self._accordion.add_card(key, card)
            setattr(self, f"_{key}_card", card)

        lay.addStretch()
        self._accordion.finalize()

        # The action stays out of the accordion -- collapsing a card must
        # never hide Run.
        lay.addWidget(self._run_btn)
        self._sync_run_button_text()

        # === Right panel (tabs provided by SidebarTabbedPage) ===
        self._tabs = self.tabs

        # Gene search and the selected-gene readout belong to whichever
        # view you are reading: under the volcano on Plots, under the
        # table on Results. One widget, moved between the two, so there
        # is only ever one box and one readout to keep in step.
        self._search_strip = QWidget()
        strip_lay = QHBoxLayout(self._search_strip)
        strip_lay.setContentsMargins(0, 0, 0, 0)
        strip_lay.setSpacing(8)
        strip_lay.addWidget(QLabel("Find gene:"))
        # Without a floor the readout takes every spare pixel and the
        # box collapses to a sliver.
        self._gene_search.setMinimumWidth(160)
        self._gene_search.setMaximumWidth(240)
        strip_lay.addWidget(self._gene_search)
        strip_lay.addStretch(1)

        # -- Tab 0: Plots --
        plots_widget = QWidget()
        plots_layout = QVBoxLayout(plots_widget)
        plots_layout.setContentsMargins(0, 0, 0, 0)

        plot_splitter = QSplitter(Qt.Orientation.Vertical)


        from kosmic.gui.shared.plots import InteractiveVolcano
        self._volcano_plot = InteractiveVolcano()
        self._volcano_plot.setMinimumHeight(100)
        self._volcano_plot.getAxis('left').setWidth(_AXIS_W)
        self._volcano_plot.gene_selected.connect(
            self._on_volcano_gene_selected)

        # Floating top-right filter button on the volcano. Menu items per agreement level
        # are populated after each analysis run.
        from PyQt6.QtWidgets import QMenu
        self._filter_btn = QPushButton("Filter \u25be")
        self._filter_btn.setParent(self._volcano_plot)
        self._filter_btn.setObjectName("volcano_filter_btn")
        self._filter_btn.setVisible(False)
        self._filter_menu = QMenu(self._filter_btn)
        self._filter_btn.setMenu(self._filter_menu)
        # The selected-gene readout belongs to the plot it describes, so
        # it floats over the volcano rather than sitting in the strip
        # underneath. Top-left, because top-right is the Filter button.
        self._gene_info_label.setParent(self._volcano_plot)
        self._gene_info_label.setObjectName("volcano_gene_legend")
        self._gene_info_label.setWordWrap(True)
        self._gene_info_label.setVisible(False)

        self._volcano_plot.installEventFilter(self)  # for resize repositioning

        plot_splitter.addWidget(self._volcano_plot)

        from kosmic.gui.shared.plots import InteractivePlot
        self._forest_plot = InteractivePlot(
            title='Per-Gene Forest Plot',
            bottom_label='Log2 Fold Change',
            unavailable_message='Click a gene on the volcano to see its per-study forest.',
        )
        self._forest_plot.setMinimumHeight(80)
        self._forest_plot.setMouseEnabled(x=True, y=False)
        plot_splitter.addWidget(self._forest_plot)

        plot_splitter.setSizes([400, 200])
        plots_layout.addWidget(plot_splitter)

        # Controls that act on the plots, under the plots. Gene search
        # highlights on the volcano and the readout describes whatever
        # was clicked -- neither is a setting for the run, so neither
        # belongs in the sidebar with the ones that are.
        plot_controls = QHBoxLayout()
        plot_controls.setContentsMargins(6, 0, 6, 0)
        plot_controls.setSpacing(8)
        plot_controls.addWidget(self._search_strip, 1)
        plot_controls.addSpacing(12)
        plot_controls.addWidget(QLabel("Y-axis:"))
        plot_controls.addWidget(self._yaxis_combo)
        plot_controls.addSpacing(8)
        plot_controls.addWidget(QLabel("Mark:"))
        plot_controls.addWidget(self._olink_cvd_check)
        plot_controls.addWidget(self._olink_explore_check)
        plots_layout.addLayout(plot_controls)
        self._plot_controls = plot_controls

        self._tabs.addTab(plots_widget, "Plots")

        # -- Tab 1b: Pathway Forest (one panel per sub-pathway) --
        from kosmic.gui.meta_analysis.pages.gene_ma.pathway_forest import (
            PathwayForestTab)
        self._pathway_forest_tab = PathwayForestTab()
        self._tabs.addTab(self._pathway_forest_tab, "Pathway Forest")

        # -- Tab 1: Results --
        from kosmic.gui.meta_analysis.pages.gene_ma.results import ResultsTab
        self._results_tab = ResultsTab()
        self._results_tab.gene_selected.connect(self._on_results_gene_selected)
        self._tabs.addTab(self._results_tab, "Results")

        # -- Tab 2: Venn (area-proportional, QPainter circles) --
        from kosmic.gui.meta_analysis.pages.gene_ma.venn import VennTab
        self._venn_tab = VennTab()
        self._venn_tab.log_message.connect(self.log_message)
        self._tabs.addTab(self._venn_tab, "Venn")

        self._tabs.currentChanged.connect(self._place_search_strip)

        # Make headline labels selectable so users can copy numbers into papers.
        _SELECTABLE_LABELS = ('_status_label',)
        flags = (Qt.TextInteractionFlag.TextSelectableByMouse
                 | Qt.TextInteractionFlag.TextSelectableByKeyboard)
        for attr in _SELECTABLE_LABELS:
            lab = getattr(self, attr, None)
            if lab is not None:
                lab.setTextInteractionFlags(flags)

    # --- Mode switching ---
    def set_mode(self, mode):
        """
        Set 'hypothesis' or 'exploratory' (Discovery) mode.

        Called by the workspace based on the user's mode-chooser selection.
        Methods Comparison is a separate workflow with its own pages, not
        a mode handled here.
        """
        previous = self._mode

        # Snapshot the current display under the previous mode before swapping.
        if self._meta_df is not None:
            if previous not in self._results_by_mode:
                self._results_by_mode[previous] = {}
            self._results_by_mode[previous] = {
                'meta_df': self._meta_df,
                'method_dfs': self._method_dfs,
                'method_keys': self._method_keys,
            }

        self._mode = mode

        if mode == 'hypothesis':
            self._mode_label.setText("Mode: Hypothesis")
            self._run_btn.setText("Run Focused Analysis")
            for i in range(self._tabs.count()):
                label = self._tabs.tabText(i)
                if label in ('Venn', 'Enrichment', 'Validation',
                             'Cross-Dataset Reproducibility'):
                    self._tabs.setTabVisible(i, False)
                elif label == 'Pathway Forest':
                    self._tabs.setTabVisible(i, True)
        else:
            self._mode_label.setText("Mode: Discovery")
            self._run_btn.setText("Run Discovery Analysis")
            for i in range(self._tabs.count()):
                label = self._tabs.tabText(i)
                if label == 'Pathway Forest':
                    self._tabs.setTabVisible(i, False)
                else:
                    self._tabs.setTabVisible(i, True)

        if previous is None or previous == mode:
            return

        self._clear_results()

        cached = self._results_by_mode.get(mode)
        if cached is not None:
            self._meta_df = cached['meta_df']
            self._method_dfs = cached['method_dfs']
            self._method_keys = cached['method_keys']
            self._apply_translational_annotation()
            self._populate_volcano()
            self._results_tab.populate(self._meta_df, self._method_keys)
            self._venn_tab.populate(self._meta_df, self._method_keys)
            self._update_panel_summary_label()
            self._update_gene_completer()
            cached_loo = cached.get('loo_result') if isinstance(cached, dict) else None
            if cached_loo is not None:
                self._loo_result = cached_loo

    def _clear_results(self):
        """Clear all result widgets and cached state."""
        self._meta_df = None
        self._method_dfs = None
        self._method_keys = None
        self._scatter_data = None
        self._highlight_scatter = None
        self._forest_xlim = None

        self._volcano_plot.clear_data()
        self._forest_plot.clear_plot_items()
        self._forest_plot.show_unavailable_message()
        style_pg_plot(self._volcano_plot, title='Volcano',
                      left_label='|Z-statistic|',
                      bottom_label='Pooled Log2 Fold Change')
        self._update_forest_labels()

        self._results_tab.clear()
        self._set_gene_legend("")

        self._venn_tab.clear()

        if hasattr(self, '_pathway_forest_tab'):
            self._pathway_forest_tab.clear()

        self._loo_result = None

        self._gene_search.blockSignals(True)
        self._gene_search.clear()
        self._gene_search.blockSignals(False)

        self._tabs.setCurrentIndex(0)
        self._status_label.setText("")

    # --- Method selection helpers (Single / Consensus radio design) ---

    def _update_forest_labels(self):
        """Set forest plot axis labels and x-range based on the current method selection."""
        keys = self._get_method_keys()
        if len(keys) == 1 and keys[0] in self._RANK_FOREST_KEYS:
            title = 'Per-Study Ranks'
            bottom = 'Signed Rank'
            self._forest_plot.setXRange(-1.1, 1.1)
        elif len(keys) == 1 and keys[0] in self._PVALUE_FOREST_KEYS:
            title = 'Per-Study Evidence'
            bottom = '-log10(p-value)'
            self._forest_plot.enableAutoRange(axis='x')
        else:
            title = 'Per-Gene Forest Plot'
            bottom = 'Log2 Fold Change'
            if self._forest_xlim:
                self._forest_plot.setXRange(
                    *self._forest_xlim, padding=0)
            else:
                self._forest_plot.enableAutoRange(axis='x')
        style_pg_plot(self._forest_plot, title=title,
                      left_label='', bottom_label=bottom)

    def _on_family_selection_changed(self, button, checked):
        """Cache current results and restore cached results for the new selection."""
        current_keys = self._get_method_keys()
        run_keys = self._method_keys or []

        if set(current_keys) == set(run_keys):
            return

        if run_keys and self._meta_df is not None:
            cache_key = tuple(sorted(run_keys))
            self._method_cache[cache_key] = {
                'meta_df': self._meta_df,
                'method_dfs': self._method_dfs,
                'method_keys': list(run_keys),
            }

        new_key = tuple(sorted(current_keys))
        cached = self._method_cache.get(new_key)
        if cached is not None:
            self._meta_df = cached['meta_df']
            self._method_dfs = cached['method_dfs']
            self._method_keys = cached['method_keys']
            self._populate_volcano()
            self._results_tab.populate(self._meta_df, self._method_keys)
            self._venn_tab.populate(self._meta_df, self._method_keys)
            self._update_panel_summary_label()
            self._update_gene_completer()
            self._update_forest_labels()
            # Re-show the forest for the last selected gene
            last = self._last_selected_gene
            if last:
                self._show_gene_forest(last)
                self._results_tab.select_gene(last)
            self._status_label.setText(
                f"Restored cached results for "
                f"{', '.join(current_keys)}.")
        else:
            self._clear_results()
            self._update_forest_labels()
            if current_keys:
                self._status_label.setText(
                    "Method changed -- click Run to update results.")

    def _sync_run_button_text(self, _=None):
        """Label the button with what it will actually do."""
        if hasattr(self, '_run_btn'):
            self._run_btn.setText(f"Run {self._run_noun()}")
        self._refresh_ma_cards()

    def _run_noun(self):
        """'Consensus' only when several methods are actually pooled.

        Calling a single-method run a consensus is simply wrong -- there
        is nothing to reach consensus between -- and it made the log
        unreadable when comparing one method against another.
        """
        consensus = (hasattr(self, '_consensus_radio')
                     and self._consensus_radio.isChecked())
        return "Consensus" if consensus else "Meta-analysis"

    def _refresh_ma_cards(self, _=None):
        """Two lines per card saying what a run would actually do."""
        if not hasattr(self, '_method_card'):
            return

        n = len(getattr(self, '_datasets', []) or [])
        if n:
            names = [str(d.get('name', '?')) for d in self._datasets]
            first = f"{n} studies pooled"
            second = ", ".join(names[:2]) + (" ..." if n > 2 else "")
        else:
            first, second = "No studies loaded yet", "Select them on Select Studies"
        self._dataset_card.set_summary([first, second])

        consensus = self._consensus_radio.isChecked()
        keys = self._get_method_keys() if hasattr(self, '_get_method_keys') else []
        self._method_card.set_summary([
            "Consensus of several methods" if consensus else "Single method",
            (f"{len(keys)} method(s) selected" if keys else "no method selected"),
        ])

        de_method = (self._de_combo.currentText()
                     if hasattr(self, '_de_combo') else "?")
        min_studies = (self._min_studies_spin.value()
                       if hasattr(self, '_min_studies_spin') else "?")
        self._settings_card.set_summary([
            f"DE input: {de_method}",
            f"a gene needs {min_studies} of {n or '?'} studies",
        ])

    def _on_method_mode_toggled(self, checked):
        """
        Toggle between Single Method and Consensus UI.

        Single: checkboxes are mutually exclusive (one family at a time).
        Consensus: checkboxes are independent (multi-select).
        Remembers the last single-mode selection when switching back.
        """
        is_single = self._single_radio.isChecked()
        if not is_single:
            for fam_key, cb in self._family_checks.items():
                if cb.isChecked():
                    self._last_single_family = fam_key
                    break
            self._family_btn_group.setExclusive(False)
            for cb in self._family_checks.values():
                cb.setChecked(True)
        else:
            self._family_btn_group.setExclusive(False)
            restore = self._last_single_family
            for fam_key, cb in self._family_checks.items():
                cb.setChecked(fam_key == restore)
            if not any(cb.isChecked()
                       for cb in self._family_checks.values()):
                first = list(self._family_checks.values())[0]
                first.setChecked(True)
            self._family_btn_group.setExclusive(True)

    def _populate_method_dropdowns(self):
        """Fill Single and Consensus dropdowns from the enabled methods in QSettings."""
        from kosmic.gui.meta_analysis.dialogs.ma_settings_dialog import (
            METHOD_FAMILIES, load_enabled_methods)

        enabled = load_enabled_methods(QSettings("KOSMIC", "KOSMIC"))

        self._single_combo.clear()
        key_to_label = dict(_METHODS)
        default_set = False
        for _fk, _fl, methods in METHOD_FAMILIES:
            for m_key, m_label in methods:
                if m_key not in enabled:
                    continue
                self._single_combo.addItem(
                    key_to_label.get(m_key, m_label), userData=m_key)
                if m_key == 'reml' and not default_set:
                    self._single_combo.setCurrentIndex(
                        self._single_combo.count() - 1)
                    default_set = True
        if not default_set and self._single_combo.count() > 0:
            self._single_combo.setCurrentIndex(0)

        for fam_key, _fam_label, methods in METHOD_FAMILIES:
            combo = self._family_combos.get(fam_key)
            cb = self._family_checks.get(fam_key)
            if combo is None:
                continue
            enabled_in_family = [(m_key, m_label) for m_key, m_label
                                 in methods if m_key in enabled]
            combo.clear()
            if len(enabled_in_family) == 0:
                combo.addItem("(none)")
                combo.setEnabled(False)
                combo.setVisible(True)
                if cb:
                    cb.setChecked(False)
                    cb.setEnabled(False)
            elif len(enabled_in_family) == 1:
                m_key, m_label = enabled_in_family[0]
                combo.addItem(
                    key_to_label.get(m_key, m_label), userData=m_key)
                combo.setEnabled(False)
                combo.setVisible(True)
                if cb:
                    cb.setEnabled(True)
            else:
                for m_key, m_label in enabled_in_family:
                    combo.addItem(
                        key_to_label.get(m_key, m_label), userData=m_key)
                combo.setEnabled(True)
                combo.setVisible(True)
                if cb:
                    cb.setEnabled(True)
                combo.setCurrentIndex(0)

    def refresh_method_controls(self):
        """Re-populate dropdowns after settings change."""
        self._populate_method_dropdowns()

    # --- Run ---
    def _get_method_keys(self):
        """Get selected method keys, applying Top 50% and HKSJ modifiers."""
        raw_keys = []
        for fam_key, combo in self._family_combos.items():
            cb = self._family_checks.get(fam_key)
            if cb and not cb.isChecked():
                continue
            key = combo.currentData()
            if key is not None:
                raw_keys.append(key)

        top50 = self._top50_cb.isChecked()
        hksj = self._hksj_cb.isChecked()

        keys = []
        for k in raw_keys:
            # Top 50% applied first so HKSJ composes on top.
            if top50:
                k = _METHOD_KEYS_TOP50.get(k, k)
            if hksj:
                k = _HKSJ_MAP.get(k, k)
            if k not in keys:
                keys.append(k)
        return keys

    def _reload_datasets_for_de_method(self, method_key: str
                                        ) -> tuple[int, int, list[str]]:
        """
        Re-discover DE result files in the project folder and
        reload self._datasets/_labels with the files matching
        'method_key' (e.g. 'deseq2', 'welch_cpm', 'welch_raw').

        This swaps the DE *method* for the studies already selected on
        the Select Studies page -- it does not re-select studies.
        Discovery walks the whole project folder, so without the filter
        below an unticked study (or a per-cell-type result sitting beside
        the whole-dataset one) would be silently pooled back in.
        Externally imported datasets have no discoverable file, so they
        are carried through unchanged.

        If a selected study only has files from a different method, it's
        loaded with that method and a warning is logged.

        Returns
        -------
        (n_loaded, n_total_discovered, missing_accessions)
            n_loaded: number of studies now in self._datasets
            n_total_discovered: number of de_results entries found
            missing_accessions: accessions that fell back to a
                non-preferred DE method (informational)
        """
        try:
            from kosmic.meta_analysis.io import (
                discover_de_results, load_de_results, select_de_entries,
            )
        except ImportError as e:
            self.log_message.emit(f"  DE method reload skipped: {e}")
            return (len(self._datasets), 0, [])

        preferred_methods = _DE_METHOD_ALIASES.get(
            method_key, [method_key])
        discovered = discover_de_results(Path(self._project_folder))
        de_entries = [e for e in discovered
                      if e.get('file_type') == 'de_results']

        # Externals were imported from outside the project, so discovery
        # never sees them; carry them through untouched.
        externals = [d for d in self._datasets if d.get('external')]
        external_names = {d['name'] for d in externals}
        # The studies ticked on the Select Studies page, plus the
        # externals (which are selected by virtue of being imported).
        selected = {str(name) for name in self._labels}

        chosen, fallback, skipped = select_de_entries(
            de_entries, preferred_methods, allowed=selected or None)

        new_datasets = list(externals)
        new_labels = [d['name'] for d in externals]
        for entry in chosen:
            acc = entry['accession']
            if acc in external_names:
                continue
            try:
                df = load_de_results(str(entry['path']), dataset_name=acc)
                new_datasets.append({
                    'path': str(entry['path']),
                    'name': acc,
                    'df':   df,
                })
                new_labels.append(acc)
            except Exception as exc:
                self.log_message.emit(
                    f"  failed to load {acc}: {exc}")

        if skipped:
            self.log_message.emit(
                f"  not selected, skipped: {', '.join(skipped)}")

        self._datasets = new_datasets
        self._labels = new_labels
        self._sync_min_studies()
        # Different DE method may live in a different analysis folder.
        clear_pb_path_cache()
        return (len(new_datasets), len(de_entries), fallback)

    def _sync_min_studies(self):
        """Cap 'Min studies' at the number of studies actually loaded.

        The default of 3 predates two-study designs. With two studies it
        drops every gene -- pooling requires a finite effect and SE from
        at least 'min_studies' of them, so nothing ever qualifies -- and
        the run returns an empty result with no explanation.
        """
        spin = getattr(self, '_min_studies_spin', None)
        if spin is None:
            return
        n_loaded = len(self._datasets)
        if n_loaded < 2:
            return
        spin.setMaximum(n_loaded)
        if spin.value() > n_loaded:
            spin.setValue(n_loaded)
            self.log_message.emit(
                f"  'Min studies' lowered to {n_loaded} -- only {n_loaded} "
                f"studies are loaded, so a higher threshold would drop "
                f"every gene.")

    def _current_de_method_key(self) -> str:
        """Canonical key for the currently-selected DE method combo."""
        idx = self._de_combo.currentIndex()
        if 0 <= idx < len(_DE_METHOD_COMBO_MAP):
            return _DE_METHOD_COMBO_MAP[idx]
        return 'deseq2'

    def _run_consensus(self):
        self.log_message.emit(f"Run clicked (mode={getattr(self, '_mode', '?')})")
        # Dispatch based on mode
        if getattr(self, '_mode', 'exploratory') == 'hypothesis':
            self._run_hypothesis()
            return

        method_keys = self._get_method_keys()
        self.log_message.emit(f"Exploratory: method_keys = {method_keys}")
        if len(method_keys) < 1:
            self.log_message.emit("Aborting -- no methods selected")
            dialogs.warning(
                self, "No Methods Selected",
                "Select at least one pooling method.")
            return

        # Reload DE files for the selected method (Import page only loads once).
        de_method_key = self._current_de_method_key()
        n_loaded, n_total, fallback = self._reload_datasets_for_de_method(
            de_method_key)
        self.log_message.emit(
            f"Exploratory: DE method = {de_method_key}; reloaded "
            f"{n_loaded} studies from {n_total} discovered entries.")

        if not self._datasets:
            self.log_message.emit("Aborting -- no datasets")
            dialogs.warning(
                self, "No Data",
                f"No DE result files matching DE method "
                f"'{de_method_key}' found in the project folder. "
                f"Check that DE has been run and saved for each "
                f"study, or pick a different DE method.")
            return

        # Modal warning on silent fallback: without it, mixed-method runs
        # surface only as identical numbers across DE-method changes.
        if fallback:
            self.log_message.emit(
                f"  studies falling back to non-preferred DE method: "
                f"{', '.join(fallback)}")
            n_fb = len(fallback)
            n_match = n_loaded - n_fb
            msg = (
                f"DE method picked: {de_method_key}.\n\n"
                f"Of {n_loaded} loaded studies:\n"
                f"  {n_match} have files matching '{de_method_key}'\n"
                f"  {n_fb} fell back to a different DE method:\n     "
                + "\n     ".join(fallback) + "\n\n"
                f"Continuing would mix DE methods across studies, "
                f"which is usually not what you want. To run "
                f"'{de_method_key}' for the missing studies, go to "
                f"the DE workspace and re-run DE for them with "
                f"that method.\n\n"
                f"Continue anyway with mixed methods?")
            if not dialogs.confirm(self, "Mixed DE methods detected", msg):
                self._status_label.setText(
                    f"Aborted: {n_fb} studies missing "
                    f"{de_method_key} DE files.")
                return

        cal_key = ['analytical', 'cc_perm'][self._cal_combo.currentIndex()]
        self.log_message.emit(f"Exploratory: calibration = {cal_key}")

        pb_paths = []
        if cal_key == 'cc_perm':
            pb_paths = find_pseudobulk_paths(
                self._datasets, self._project_folder,
                log_cb=self.log_message.emit)
            self.log_message.emit(
                f"Exploratory: found {len(pb_paths)}/{len(self._datasets)} pseudobulk files")
            if not pb_paths:
                self.log_message.emit("Aborting -- no pseudobulk")
                dialogs.warning(
                    self, "No Pseudobulk",
                    "CC permutation requires pseudobulk CSV files.")
                return

        # Independent filter: shrink the gene universe before pooling so
        # BH FDR runs on the survivors. Per-dataset DataFrames are
        # filtered to that universe before being handed to the worker.
        from kosmic import META_IF_FILTER_STAT
        from kosmic.meta_analysis.independent_filter import (
            apply_filter, filter_datasets_to_universe,
        )
        if_active = self._if_check.isChecked()
        if if_active:
            if_threshold = self._if_pct_spin.value() / 100.0
            gene_universe = apply_filter(
                self._datasets,
                threshold=if_threshold,
                stat=META_IF_FILTER_STAT,
                min_studies=self._min_studies_spin.value(),
            )
            datasets_for_pool = filter_datasets_to_universe(
                self._datasets, gene_universe)
            n_pre = sum(len(d['df']) for d in self._datasets) // max(len(self._datasets), 1)
            n_post = len(gene_universe)
            self.log_message.emit(
                f"Independent filter: {n_post:,} genes pass "
                f"(threshold = {if_threshold:.0%}, stat = {META_IF_FILTER_STAT}; "
                f"~{n_pre:,} per-study before filter)"
            )
        else:
            datasets_for_pool = self._datasets

        de_dfs = [d['df'] for d in datasets_for_pool]
        labels = [d['name'] for d in datasets_for_pool]

        self._run_btn.setEnabled(False)
        self._status_label.setText(
            f"Running (DE method: {de_method_key}, "
            f"calibration: {cal_key})...")
        self.log_message.emit(
            f"Starting {self._run_noun().lower()} run: {len(de_dfs)} datasets, "
            f"{', '.join(method_keys)} ({cal_key}, DE={de_method_key})")
        if self.progress_bar:
            self.progress_bar.setRange(0, 0)

        # CC permutation uses fast NB GLM as the DESeq2 equivalent; Welch passes through.
        cc_de_method = {
            'deseq2':    'deseq2',
            'welch_cpm': 'welch_cpm',
            'welch_raw': 'welch_raw',
        }.get(de_method_key, 'deseq2')

        params = {
            'method_keys': method_keys,
            'calibration': cal_key,
            'min_studies': self._min_studies_spin.value(),
            'cc_n_perms': 1000,
            'pseudobulk_paths': pb_paths,
            'de_method':   cc_de_method,
        }

        # LOO worker re-uses these so the user doesn't have to re-pick methods.
        self._last_consensus_params = params
        self._last_consensus_pb_paths = pb_paths
        self._last_consensus_de_dfs = de_dfs
        self._last_consensus_labels = labels

        self._run_btn.setEnabled(False)
        self._status_label.setText("Running...")
        self._last_logged_perm_decile = -1  # so new runs log at each 10%
        self.sidebar_status.emit("Running meta-analysis…", '#2196F3')
        if self.progress_bar:
            self.progress_bar.setRange(0, 0)  # indeterminate spinner

        self._warn_unadjusted(de_dfs)
        # The folder is fixed when the run starts: the selection may
        # change while it runs (ADR-007).
        self._pending_output_selection = self._output_selection
        self._worker = GeneMAWorker(de_dfs, labels, params)
        run_worker(
            self._worker,
            on_finished=self._on_finished,
            on_failed=self._on_failed,
            on_progress=self._on_progress,
        )

    def _warn_unadjusted(self, de_dfs):
        """Name the studies left out of the direction counts."""
        from kosmic.meta_analysis.direction import unadjusted_warning
        warning = unadjusted_warning(de_dfs)
        if warning:
            self.log_message.emit(f"Warning: {warning}")

    def _on_progress(self, message):
        self._status_label.setText(message)
        self.sidebar_status.emit(message, '#2196F3')
        if message.startswith('CC permutation ') and '/' in message:
            try:
                parts = message.split()
                cur, tot = parts[2].split('/')
                cur_i, tot_i = int(cur), int(tot)
                if cur_i > 0 and self.progress_bar:
                    self.progress_bar.setRange(0, 100)
                    self.progress_bar.setValue(int(cur_i / tot_i * 100))
                # Track highest decile seen so parallel batches (~18 perms each)
                # log on every 10% crossed instead of falling between modulo checks.
                decile = int(cur_i / tot_i * 10)
                if cur_i == tot_i or decile > self._last_logged_perm_decile:
                    self.log_message.emit(message)
                    self._last_logged_perm_decile = decile
            except (IndexError, ValueError):
                pass
        else:
            self.log_message.emit(message)

    def _on_failed(self, message):
        self._run_btn.setEnabled(True)
        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)
        self._status_label.setText(f"Failed: {message}")
        self.sidebar_status.emit(
            f"Failed: {message.splitlines()[0]}", '#F44336')
        self.log_message.emit(f"Analysis failed: {message}")
        dialogs.warning(self, "Failed", message)

    def _on_finished(self, data):
        self._run_btn.setEnabled(True)
        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(100)

        message = f"{self._run_noun()} complete"
        self._status_label.setText(message)
        self.log_message.emit(message)
        self._meta_df = data['meta_df']
        self._log_direction_summary()
        self._method_dfs = data['method_dfs']
        self._method_keys = data['method_keys']

        # Annotate with Olink panel membership + (hypothesis mode)
        # the parent pathway each consensus gene came from, before
        # anything downstream consumes the DataFrame.
        self._apply_translational_annotation()
        self._annotate_pathway_membership()

        self._results_by_mode[self._mode] = {
            'meta_df': self._meta_df,
            'method_dfs': self._method_dfs,
            'method_keys': self._method_keys,
        }

        self._populate_volcano()
        self._results_tab.populate(self._meta_df, self._method_keys)
        self._update_gene_completer()
        self._run_output_selection = self._pending_output_selection
        self._auto_save()
        self._record_meta_provenance('meta_gene')

        # Skip mode-hidden tabs so we don't compute layouts the user
        # can't see (see 'set_mode' for hiding rules).
        if self._mode == 'hypothesis':
            self._pathway_forest_tab.draw(
                self._meta_df,
                getattr(self, '_hyp_pathways', {}),
                getattr(self, '_pathway_ma_results', None),
            )
        else:
            self._venn_tab.populate(self._meta_df, self._method_keys)
            self._update_panel_summary_label()
        last = self._last_selected_gene
        if last and self._meta_df is not None:
            if last in self._meta_df['names'].values:
                self._show_gene_forest(last)
                self._results_tab.select_gene(last)
        self.analysis_complete.emit()

    def _log_direction_summary(self):
        """Report how many significant genes have opposite-direction effects."""
        from kosmic.meta_analysis.direction import count_significant_conflicts
        n_sig, n_conflict = count_significant_conflicts(
            self._meta_df, fdr=DEFAULT_FDR)
        if n_conflict is None:
            return
        self.log_message.emit(
            f"{n_sig:,} significant (FDR < {DEFAULT_FDR:g}), {n_conflict:,} "
            f"with opposite-direction effects across studies")

    def _apply_translational_annotation(self):
        """
        Add Olink panel membership boolean columns to self._meta_df.

        Idempotent; safe when panels are missing on disk.
        """
        if self._meta_df is None or 'names' not in self._meta_df.columns:
            return
        if all(f'on_{p}' in self._meta_df.columns for p in DEFAULT_PANELS):
            return
        try:
            self._meta_df = annotate_consensus(
                self._meta_df, gene_col='names', panels=DEFAULT_PANELS)
        except Exception as exc:
            self.log_message.emit(
                f"Translational panel annotation skipped: {exc}")

    def _annotate_pathway_membership(self):
        """Add a 'pathways' column to self._meta_df from self._hyp_pathways.

        Hypothesis-mode only -- builds a gene -> pathway-list reverse map
        from the user's selected pathways and joins it onto the
        consensus DataFrame so the Results tab can show which pathway
        each gene came from. No-op outside hypothesis mode or when no
        pathways are selected.
        """
        if self._meta_df is None or 'names' not in self._meta_df.columns:
            return
        pathways = getattr(self, '_hyp_pathways', None) or {}
        if not pathways:
            return

        gene_to_pathways: dict[str, list[str]] = {}
        for pw_name, genes in pathways.items():
            for gene in genes:
                key = str(gene).upper()
                gene_to_pathways.setdefault(key, []).append(pw_name)

        def _lookup(name):
            hits = gene_to_pathways.get(str(name).upper())
            return '; '.join(hits) if hits else ''

        self._meta_df['pathways'] = self._meta_df['names'].map(_lookup)

    def _update_panel_summary_label(self):
        """Refresh the Olink panel coverage line at the bottom of the Results tab."""
        if self._meta_df is None or len(self._meta_df) == 0:
            self._results_tab.set_panel_summary("")
            return
        try:
            summ = panel_summary(
                self._meta_df, gene_col='names',
                fdr_col='fdr', fdr_threshold=DEFAULT_FDR,
                panels=DEFAULT_PANELS)
        except Exception as exc:
            self._results_tab.set_panel_summary("")
            self.log_message.emit(f"Panel summary skipped: {exc}")
            return

        n_consensus = summ['n_consensus']
        if n_consensus == 0:
            self._results_tab.set_panel_summary(
                "No significant consensus genes (FDR < 0.05).")
            return

        parts = [f"Consensus: {n_consensus} genes (FDR < 0.05)"]
        for key in DEFAULT_PANELS:
            info = summ['panels'].get(key)
            if info is None:
                continue
            parts.append(
                f"{info['label']}: {info['n_on_panel']}/{n_consensus} "
                f"({info['pct']:.1f}%)")
        self._results_tab.set_panel_summary("  |  ".join(parts))

    # --- Volcano ---
    def _populate_volcano(self):
        """Render the consensus volcano via the shared InteractiveVolcano widget."""
        self._volcano_plot.clear_data()
        if self._meta_df is None:
            return

        meta_df = self._meta_df
        names = meta_df['names'].tolist()

        # X-axis: rank methods use mean signed rank; effect-size / p-value methods use logFC.
        keys = self._method_keys or []
        is_rank = (len(keys) == 1
                   and keys[0] in self._RANK_FOREST_KEYS)
        if is_rank and 'mean_signed_rank' in meta_df.columns:
            x = meta_df['mean_signed_rank'].values.astype(float)
            x_label = 'Mean Signed Rank'
        else:
            x = meta_df['logfoldchanges'].values.astype(float)
            x_label = 'Log2 Fold Change'

        # Base method (DL preferred) supplies the y-axis FDR.
        base_key = None
        if self._method_keys:
            for k in self._method_keys:
                if k.startswith('dl'):
                    base_key = k
                    break
            if base_key is None:
                base_key = self._method_keys[0]

        # Count how many statistical families each gene passes at FDR < 0.05.
        # A family agrees if ANY of its selected methods is significant
        # (within-family agreement is redundant -- shared framework).
        from kosmic.gui.meta_analysis.dialogs.ma_settings_dialog import METHOD_FAMILIES
        selected_keys = set(self._method_keys or [])
        family_key_groups = []
        for fam_key, fam_label, methods in METHOD_FAMILIES:
            fam_methods = [k for k, _ in methods if k in selected_keys]
            if fam_methods:
                family_key_groups.append((fam_label, fam_methods))
        n_families = len(family_key_groups)

        n_sig_per_gene = np.zeros(len(meta_df), dtype=int)
        base_fdr = np.ones(len(meta_df))
        for _fam_label, fam_methods in family_key_groups:
            fam_sig = np.zeros(len(meta_df), dtype=bool)
            for key in fam_methods:
                fdr_col = f'fdr_{key}'
                if fdr_col not in meta_df.columns:
                    continue
                fdr_m = meta_df[fdr_col].values.astype(float)
                fam_sig |= (fdr_m < DEFAULT_FDR)
                if key == base_key:
                    base_fdr = fdr_m
            n_sig_per_gene[fam_sig] += 1

        n_methods = n_families  # single method -> levels are 0 or 1

        use_z = (self._yaxis_combo.currentIndex() == 1
                 and 'z_stat' in meta_df.columns)
        if use_z:
            y = np.abs(meta_df['z_stat'].values.astype(float))
        else:
            y = neg_log10(base_fdr)

        # Used by _on_gene_search to find a clicked point.
        self._scatter_data = {
            'x': x, 'y': y, 'names': names,
            'pvals_raw': base_fdr, 'padj': base_fdr,
            'n_sig': n_sig_per_gene,
        }

        # Grey when non-significant; red/blue shaded by family-agreement count.
        def _shade(base_r, base_g, base_b, frac):
            gr = 150
            r = int(gr + (base_r - gr) * frac)
            g = int(gr + (base_g - gr) * frac)
            b = int(gr + (base_b - gr) * frac)
            a = int(80 + 140 * frac)
            return pg.mkBrush(r, g, b, a)

        brushes = []
        for i in range(len(meta_df)):
            level = int(n_sig_per_gene[i])
            if level == 0:
                brushes.append(pg.mkBrush(150, 150, 150, 60))
            else:
                frac = level / max(n_methods, 1)
                if x[i] > 0:
                    brushes.append(_shade(220, 50, 50, frac))
                else:
                    brushes.append(_shade(50, 50, 220, frac))

        y_label = '|Z-statistic|' if use_z else f'-Log10 FDR ({base_key})'
        n_sig = int((n_sig_per_gene >= 1).sum())
        if n_methods == 1:
            method_name = (self._method_keys[0] if self._method_keys
                           else 'method')
            title = f'{method_name}: {n_sig} significant genes'
        else:
            n_all = int((n_sig_per_gene == n_methods).sum())
            title = (f'{n_all} genes in all {n_methods} families '
                     f'({n_sig} in any)')

        # Captured into closure so the hover label survives method switches.
        _xl, _yl = x_label, y_label

        def _fmt(name, gx, gy):
            return f"{name}\n{_xl}: {gx:.3f}\n{_yl}: {gy:.2f}"

        self._volcano_plot.set_scatter(
            x, y, names,
            brushes=brushes, sizes=[6] * len(meta_df),
            groups=n_sig_per_gene,
            x_label=x_label, y_label=y_label, title=title,
            hover_fmt=_fmt,
        )

        # Olink panel overlay rings -- distinct colour/size per panel.
        if (getattr(self, '_olink_cvd_check', None) is not None
                and self._olink_cvd_check.isChecked()
                and 'on_olink_cvd_iii' in meta_df.columns):
            mask = meta_df['on_olink_cvd_iii'].fillna(False).to_numpy(dtype=bool)
            self._volcano_plot.add_overlay_ring(
                mask, color=(255, 200, 0, 220), size=14, pen_width=2)
        if (getattr(self, '_olink_explore_check', None) is not None
                and self._olink_explore_check.isChecked()
                and 'on_olink_explore_3072' in meta_df.columns):
            mask = meta_df['on_olink_explore_3072'].fillna(False).to_numpy(dtype=bool)
            self._volcano_plot.add_overlay_ring(
                mask, color=(0, 200, 255, 160), size=10, pen_width=1)

        if not use_z:
            self._volcano_plot.add_threshold_line(y=-np.log10(DEFAULT_FDR))

        self._volcano_plot.set_symmetric_x()

        # Pre-compute shared forest x-limits from BOTH pooled CIs AND
        # per-study CIs (logFC +/- 1.96*SE). Per-study CIs are typically
        # wider than the pooled CI, so the axis must accommodate them
        # or error bars get clipped.
        self._forest_xlim = None
        extreme = 0.0
        if 'ci_lower' in meta_df.columns and 'ci_upper' in meta_df.columns:
            ci_lo = meta_df['ci_lower'].dropna()
            ci_hi = meta_df['ci_upper'].dropna()
            if len(ci_lo) > 0:
                extreme = max(abs(ci_lo.min()), abs(ci_hi.max()))
        if 'study_effects' in meta_df.columns:
            for effs in meta_df['study_effects']:
                if not isinstance(effs, list):
                    continue
                for eff in effs:
                    lfc = eff.get('logfc', 0.0)
                    se = eff.get('se', 0.0)
                    if np.isfinite(lfc) and np.isfinite(se) and se > 0:
                        lo = abs(lfc - 1.96 * se)
                        hi = abs(lfc + 1.96 * se)
                        extreme = max(extreme, lo, hi)
        if extreme > 0:
            pad = extreme * 0.05
            self._forest_xlim = (-(extreme + pad), extreme + pad)

        self._build_sig_toggles(n_methods, n_sig_per_gene)

    def eventFilter(self, obj, event):
        """Reposition the floating volcano overlays when it resizes."""
        if obj is self._volcano_plot and event.type() == event.Type.Resize:
            btn = self._filter_btn
            if btn.isVisible():
                margin = 6
                btn.move(obj.width() - btn.width() - margin, margin)
            self._place_gene_legend()
        return super().eventFilter(obj, event)

    def _set_gene_legend(self, text):
        """Set the floating readout, or hide it when there is no gene."""
        self._gene_info_label.setText(text or "")
        self._gene_info_label.setVisible(bool(text))
        self._place_gene_legend()

    def _place_gene_legend(self):
        """Sit the legend inside the plot area, clear of the y-axis."""
        lbl = self._gene_info_label
        # Gate on the content, not on isVisible(): a gene picked while
        # another tab is showing reports isVisible() False, and skipping
        # placement there would strand the legend at (0, 0) when the
        # Plots tab came back.
        if not lbl.text():
            return
        margin = 6
        # Anchor inside the plotting area, not the widget: the widget's
        # top strip is the title and its left strip is the y-axis, and
        # sitting at (margin, margin) put the legend across the title.
        left, top = _AXIS_W, 0
        try:
            rect = self._volcano_plot.getPlotItem().vb.geometry()
            left, top = int(rect.left()), int(rect.top())
        except (AttributeError, RuntimeError):
            pass
        # Cap the width so a long readout wraps instead of covering the
        # plot, and leave room for the Filter button on the same line.
        avail = max(120, self._volcano_plot.width() - left
                    - self._filter_btn.width() - margin * 4)
        lbl.setMaximumWidth(min(320, avail))
        lbl.adjustSize()
        lbl.move(left + margin, top + margin)
        lbl.raise_()

    def _build_sig_toggles(self, n_methods, n_sig_per_gene):
        """Populate the volcano filter dropdown with one checkable action per agreement level."""
        self._filter_menu.clear()
        self._sig_toggles = []

        unit = 'family' if n_methods <= 1 else 'families'
        for level in range(n_methods + 1):
            count = int((n_sig_per_gene == level).sum())
            if count == 0:
                continue
            if level == 0:
                label = f"Not significant ({count:,})"
            elif n_methods == 1:
                label = f"Significant ({count:,})"
            elif level == n_methods:
                label = f"All {level} {unit} agree ({count:,})"
            else:
                label = f"{level}/{n_methods} {unit} ({count:,})"

            action = self._filter_menu.addAction(label)
            action.setCheckable(True)
            action.setChecked(True)
            action.toggled.connect(
                lambda checked, lv=level: self._toggle_sig_level(lv, checked))
            self._sig_toggles.append((level, action))

        self._filter_btn.setVisible(bool(self._sig_toggles))
        if self._sig_toggles:
            self._filter_btn.adjustSize()
            margin = 6
            vw = self._volcano_plot.width()
            self._filter_btn.move(
                max(vw - self._filter_btn.width() - margin, margin),
                margin)

    def _toggle_sig_level(self, level, visible):
        """Show or hide dots for a given significance level."""
        self._volcano_plot.set_group_visible(int(level), bool(visible))

    def _update_gene_completer(self):
        """Refresh the gene-search autocomplete list from _meta_df."""
        if self._meta_df is not None and 'names' in self._meta_df.columns:
            genes = sorted(self._meta_df['names'].astype(str).unique())
        else:
            genes = []
        self._gene_completer_model.setStringList(genes)

    def _place_search_strip(self, index=None):
        """Move the gene-search line to whichever tab is being read.

        Plots and Results both want it; every other tab has nothing for
        it to act on, so it goes back under the plots.
        """
        if index is None:
            index = self._tabs.currentIndex()
        if self._tabs.widget(index) is self._results_tab:
            self._results_tab.take_search_widget(self._search_strip)
        else:
            self._plot_controls.insertWidget(0, self._search_strip, 1)
        self._search_strip.show()

    def _on_gene_search(self, text):
        """Highlight matches on the volcano and select one in the table."""
        self._results_tab.select_matching(text)
        if self._highlight_scatter is not None:
            self._volcano_plot.removeItem(self._highlight_scatter)
            self._highlight_scatter = None

        if not text or self._scatter_data is None:
            return

        sd = self._scatter_data
        text_upper = text.upper()
        matches = [i for i, n in enumerate(sd['names'])
                   if text_upper in n.upper()]

        if not matches:
            return

        xi = np.array([sd['x'][i] for i in matches])
        yi = np.array([sd['y'][i] for i in matches])
        ni = [sd['names'][i] for i in matches]

        self._highlight_scatter = pg.ScatterPlotItem(
            x=xi, y=yi, size=12,
            pen=pg.mkPen('white', width=2),
            brush=pg.mkBrush(255, 255, 0, 200),
            symbol='star', data=ni,
            hoverable=True, hoverSize=16,
            tip=lambda x, y, data: f"{data}\nlogFC: {x:.3f}",
        )
        # Direct sigClicked so scene-level clicks can't get blocked by these stars on top.
        self._highlight_scatter.sigClicked.connect(
            self._on_highlight_clicked)
        self._volcano_plot.addItem(self._highlight_scatter)

        exact = [i for i in matches
                 if sd['names'][i].upper() == text_upper]
        if len(exact) == 1:
            gene = sd['names'][exact[0]]
            self._show_gene_forest(gene)
            self._results_tab.select_gene(gene)

    def _on_highlight_clicked(self, scatter, points, ev):
        """Handle click on a highlighted (search-result) star."""
        if not points:
            return
        gene = points[0].data()
        if gene:
            self._show_gene_forest(gene)
            self._results_tab.select_gene(gene)

    def _redraw_volcano(self, _=None):
        if self._meta_df is not None:
            self._populate_volcano()

    def _on_volcano_gene_selected(self, gene: str) -> None:
        """Show the gene-level forest plot and scroll the results table to the gene."""
        self._show_gene_forest(gene)
        self._results_tab.select_gene(gene)

    # --- Forest plot ---
    # Methods whose forest uses the signed-rank renderer.
    _RANK_FOREST_KEYS = frozenset((
        'sumrank', 'sumrank_top50',
    ))

    # Methods whose forest uses the p-value renderer (combined-p, no pooled logFC).
    _PVALUE_FOREST_KEYS = frozenset((
        'fisher', 'rop', 'wop', 'gwop',
    ))

    def _on_gene_search_entered(self):
        """Enter commits the search: show that gene's details in full.

        Typing only highlights and scrolls -- redrawing the forest on
        every keystroke would be wasteful -- so Enter is what fills in
        the readout beside the box.
        """
        text = self._gene_search.text().strip()
        if not text or not self._results_tab.select_matching(text):
            return
        rows = self._results_tab.table.selectionModel().selectedRows()
        if not rows:
            return
        item = self._results_tab.table.item(rows[0].row(), 0)
        if item:
            self._show_gene_forest(item.text())

    def _show_gene_forest(self, gene_name):
        """
        Dispatcher: picks the appropriate forest-plot renderer.

        - Single rank method  -> rank forest (per-study signed ranks).
        - Single p-value method -> p-value forest (per-study -log10 p).
        - Otherwise -> parametric effect-size forest with pooled diamond.

        Also remembers the selected gene for cache restore.
        """
        self._last_selected_gene = gene_name
        self._forest_plot.clear_plot_items()
        self._forest_plot.hide_unavailable_message()
        if self._meta_df is None:
            return
        row = self._meta_df[self._meta_df['names'] == gene_name]
        if row.empty:
            return
        row = row.iloc[0]
        study_effects = row.get('study_effects', [])
        if not study_effects or not isinstance(study_effects, list):
            self._set_gene_legend(
                f"{gene_name}: no per-study data")
            return

        keys = self._method_keys or []
        if len(keys) == 1:
            k = keys[0]
            if k in self._RANK_FOREST_KEYS:
                return self._show_gene_forest_rank(
                    gene_name, row, study_effects)
            if k in self._PVALUE_FOREST_KEYS:
                return self._show_gene_forest_pvalue(
                    gene_name, row, study_effects)
        return self._show_gene_forest_parametric(
            gene_name, row, study_effects)

    def _show_gene_forest_parametric(self, gene_name, row, study_effects):
        """Effect-size forest with per-method p-values in info label."""
        import re

        fg = get_color('fg_primary')
        accent = get_color('accent_primary')
        n = len(study_effects)

        style_pg_plot(self._forest_plot,
                      title=f'{gene_name} -- Forest Plot',
                      left_label='', bottom_label='Log2 Fold Change')

        y_positions = list(range(n))
        x_vals = [eff.get('logfc', 0.0) for eff in study_effects]
        ses = [eff.get('se', float('nan')) for eff in study_effects]

        # Error bars
        for i in range(n):
            se = ses[i]
            if np.isfinite(se) and se > 0:
                lo = x_vals[i] - 1.96 * se
                hi = x_vals[i] + 1.96 * se
                err = pg.PlotDataItem(
                    [lo, hi], [y_positions[i], y_positions[i]],
                    pen=pg.mkPen(fg, width=2.0))
                self._forest_plot.addItem(err)

        weights = np.array([
            1.0 / (se**2) if np.isfinite(se) and se > 0 else 0.0
            for se in ses])
        if weights.max() > 0:
            norm_w = weights / weights.max()
            sizes = 6 + 12 * np.sqrt(norm_w)
        else:
            sizes = np.full(n, 10.0)

        scatter = pg.ScatterPlotItem(
            x=np.array(x_vals), y=np.array(y_positions, dtype=float),
            size=sizes, pen=pg.mkPen(fg), brush=pg.mkBrush(accent),
            symbol='s')
        self._forest_plot.addItem(scatter)

        # Diamond CI is raw 95%; fill = FDR-significant, hollow = CI excludes
        # zero but doesn't survive multiple testing.
        pooled = row.get('logfoldchanges', 0)
        ci_lo = row.get('ci_lower', pooled)
        ci_hi = row.get('ci_upper', pooled)

        # Use consensus FDR if available; else first per-method FDR.
        gene_fdr = row.get('fdr', float('nan'))
        if not np.isfinite(gene_fdr) and self._method_keys:
            for key in self._method_keys:
                f = row.get(f'fdr_{key}', float('nan'))
                if np.isfinite(f):
                    gene_fdr = f
                    break
        fdr_sig = np.isfinite(gene_fdr) and gene_fdr < DEFAULT_FDR

        dy = -1.0
        dh = 0.3
        diamond_x = [ci_lo, pooled, ci_hi, pooled, ci_lo]
        diamond_y = [dy, dy + dh, dy, dy - dh, dy]
        if fdr_sig:
            diamond = pg.PlotDataItem(
                diamond_x, diamond_y,
                pen=pg.mkPen('red', width=2),
                fillLevel=dy, brush=pg.mkBrush(255, 0, 0, 80))
        else:
            diamond = pg.PlotDataItem(
                diamond_x, diamond_y,
                pen=pg.mkPen('#888', width=2, style=Qt.PenStyle.DashLine),
                fillLevel=None)
        self._forest_plot.addItem(diamond)

        zero = pg.InfiniteLine(
            pos=0, angle=90,
            pen=pg.mkPen(fg, width=1, style=Qt.PenStyle.DashLine))
        self._forest_plot.addItem(zero)

        raw_labels = [eff.get('dataset', f'Study_{i+1}')
                      for i, eff in enumerate(study_effects)]
        short_codes = []
        for lbl in raw_labels:
            m = re.match(
                r'(GSE\d+|GSM\d+|SRP\d+|SCP\d+|PRJNA\d+)',
                lbl, re.IGNORECASE)
            short_codes.append(m.group(1) if m else lbl.split('_')[0])

        keys = self._method_keys or []
        pooled_label = keys[0].upper() if len(keys) == 1 else 'Pooled'

        y_ticks = [(float(i), str(i + 1)) for i in range(n)]
        y_ticks.append((dy, pooled_label))
        ax = self._forest_plot.getAxis('left')
        ax.setTicks([y_ticks])
        self._forest_plot.setYRange(dy - 0.8, n - 0.5)

        # Shared x-limits so the axis doesn't jump between genes.
        if self._forest_xlim:
            self._forest_plot.setXRange(*self._forest_xlim, padding=0)

        i2 = row.get('heterogeneity_i2', float('nan'))
        i2_str = f'I\u00b2={i2:.0%}' if np.isfinite(i2) else ''

        if np.isfinite(gene_fdr):
            if fdr_sig:
                fdr_html = (f'<span style="color:#4CAF50; font-weight:bold">'
                            f'p &lt; 0.05 ({gene_fdr:.2e})</span>')
            else:
                ci_excl = (ci_lo > 0 or ci_hi < 0)
                note = (' -- CI excludes zero but does not survive '
                        'multiple testing' if ci_excl else '')
                fdr_html = (f'<span style="color:#F44336">'
                            f'p &gt; 0.05 ({gene_fdr:.2e}){note}</span>')
        else:
            fdr_html = '<span style="color:#999">unavailable</span>'

        # Per-method FDR lines only in consensus mode with >1 method.
        method_lines = []
        keys = self._method_keys or []
        if len(keys) > 1:
            for key in keys:
                fdr_col = f'fdr_{key}'
                f = row.get(fdr_col, float('nan'))
                if np.isfinite(f):
                    colour = '#4CAF50' if f < DEFAULT_FDR else '#F44336'
                    method_lines.append(
                        f'<span style="color:{colour}">{key}: '
                        f'{f:.2e}</span>')

        key_parts = [f'{i+1}={short_codes[i]}' for i in range(n)]

        html = (
            f'<b style="font-size:11pt">{gene_name}</b><br>'
            f'<b>log2FC:</b> {pooled:.3f} '
            f'[{ci_lo:.3f}, {ci_hi:.3f}]<br>'
            f'<b>FDR-adjusted p:</b> {fdr_html}<br>'
        )
        if i2_str:
            html += f'<b>{i2_str}</b> | '
        html += f'<b>Studies:</b> {n}<br>'
        if method_lines:
            html += ' | '.join(method_lines) + '<br>'
        html += (f'<span style="color:#999; font-size:8pt">'
                 f'{", ".join(key_parts)}</span>')

        self._gene_info_label.setTextFormat(Qt.TextFormat.RichText)
        self._set_gene_legend(html)

    def _show_gene_forest_rank(self, gene_name, row, study_effects):
        """SumRank forest: per-study signed ranks on [-1, +1] (no pooled diamond)."""
        import re

        style_pg_plot(self._forest_plot,
                      title=f'{gene_name} -- Per-Study Ranks',
                      left_label='', bottom_label='Signed Rank')

        n = len(study_effects)
        fg = get_color('fg_primary')

        y_positions = list(range(n))
        x_vals = [eff.get('signed_rank', 0.0) for eff in study_effects]

        raw_labels = [eff.get('dataset', f'Study_{i+1}')
                      for i, eff in enumerate(study_effects)]
        short_codes = []
        for lbl in raw_labels:
            m = re.match(
                r'(GSE\d+|GSM\d+|SRP\d+|SCP\d+|PRJNA\d+)',
                lbl, re.IGNORECASE)
            short_codes.append(m.group(1) if m else lbl.split('_')[0])

        colors = []
        for x in x_vals:
            if x > 0:
                colors.append(pg.mkBrush(220, 60, 60, 180))   # up
            elif x < 0:
                colors.append(pg.mkBrush(60, 60, 220, 180))   # down
            else:
                colors.append(pg.mkBrush(150, 150, 150, 150))

        scatter = pg.ScatterPlotItem(
            x=np.array(x_vals), y=np.array(y_positions, dtype=float),
            size=10, pen=pg.mkPen(fg), brush=colors, symbol='s')
        self._forest_plot.addItem(scatter)

        zero_line = pg.InfiniteLine(
            pos=0, angle=90,
            pen=pg.mkPen(fg, width=1, style=Qt.PenStyle.DashLine))
        self._forest_plot.addItem(zero_line)

        y_ticks = [(float(i), str(i + 1)) for i in range(n)]
        ax = self._forest_plot.getAxis('left')
        ax.setTicks([y_ticks])
        self._forest_plot.setYRange(-0.5, n - 0.5)
        self._forest_plot.setXRange(-1.1, 1.1)

        mean_sr = row.get('mean_signed_rank', 0.0)
        pval = row.get('pvals_pooled', float('nan'))
        direction = row.get('direction', '?')
        key_str = '  '.join(
            f'{i+1}={short_codes[i]}' for i in range(n))

        self._set_gene_legend(
            f"{gene_name}  |  mean signed rank = {mean_sr:.3f}  |  "
            f"direction = {direction}  |  "
            f"p = {pval:.2e}  |  k = {n}  |  {key_str}")

    def _show_gene_forest_pvalue(self, gene_name, row, study_effects):
        """
        P-value forest: per-study -log10(p), direction-coloured.

        Used for combined-p methods (Fisher, wOP, gwOP, rOP) where
        per-study p-values are the right axis, not log fold changes.
        """
        import re

        style_pg_plot(self._forest_plot,
                      title=f'{gene_name} -- Per-Study Evidence',
                      left_label='', bottom_label='-log10(p-value)')

        n = len(study_effects)
        fg = get_color('fg_primary')

        y_positions = list(range(n))
        pvals = [eff.get('pval', 1.0) for eff in study_effects]
        x_vals = list(neg_log10(pvals, floor=1e-16))
        logfcs = [eff.get('logfc', 0.0) for eff in study_effects]

        raw_labels = [eff.get('dataset', f'Study_{i+1}')
                      for i, eff in enumerate(study_effects)]
        short_codes = []
        for lbl in raw_labels:
            m = re.match(
                r'(GSE\d+|GSM\d+|SRP\d+|SCP\d+|PRJNA\d+)',
                lbl, re.IGNORECASE)
            short_codes.append(m.group(1) if m else lbl.split('_')[0])

        # Colour by logFC direction
        colors = []
        for lfc in logfcs:
            if lfc > 0:
                colors.append(pg.mkBrush(220, 60, 60, 180))
            elif lfc < 0:
                colors.append(pg.mkBrush(60, 60, 220, 180))
            else:
                colors.append(pg.mkBrush(150, 150, 150, 150))

        # Size by weight (wOP/gwOP) or equal (rOP/Fisher)
        weights = [eff.get('weight', 1.0) for eff in study_effects]
        w_arr = np.array(weights)
        if w_arr.max() > 0 and w_arr.max() != w_arr.min():
            norm_w = w_arr / w_arr.max()
            sizes = 6 + 12 * np.sqrt(norm_w)
        else:
            sizes = np.full(n, 10.0)

        scatter = pg.ScatterPlotItem(
            x=np.array(x_vals), y=np.array(y_positions, dtype=float),
            size=sizes, pen=pg.mkPen(fg), brush=colors, symbol='s')
        self._forest_plot.addItem(scatter)

        ref_line = pg.InfiniteLine(  # -log10(0.05) ~= 1.3
            pos=-np.log10(DEFAULT_FDR), angle=90,
            pen=pg.mkPen(fg, width=1, style=Qt.PenStyle.DashLine))
        self._forest_plot.addItem(ref_line)

        # rOP: highlight the deciding r-th order-statistic study with a gold ring.
        r_used = row.get('r_used', None)
        if r_used is not None:
            try:
                r_used_i = int(r_used)
                sorted_idx = sorted(range(n), key=lambda i: pvals[i])
                if 0 < r_used_i <= len(sorted_idx):
                    deciding_idx = sorted_idx[r_used_i - 1]
                    ring = pg.ScatterPlotItem(
                        x=[x_vals[deciding_idx]],
                        y=[float(deciding_idx)],
                        size=18, pen=pg.mkPen('#FFD700', width=2),
                        brush=pg.mkBrush(0, 0, 0, 0), symbol='o')
                    self._forest_plot.addItem(ring)
            except (TypeError, ValueError):
                pass

        y_ticks = [(float(i), str(i + 1)) for i in range(n)]
        ax = self._forest_plot.getAxis('left')
        ax.setTicks([y_ticks])
        self._forest_plot.setYRange(-0.5, n - 0.5)

        pval = row.get('pvals_pooled', float('nan'))
        direction = row.get('direction', '?')
        r_info = (f"r={row.get('r_used', '?')}"
                  if 'r_used' in row.index else '')
        pooling = row.get('pooling_method',
                          (self._method_keys or [''])[0])
        key_str = '  '.join(
            f'{i+1}={short_codes[i]}' for i in range(n))

        self._set_gene_legend(
            f"{gene_name}  |  {pooling} p = {pval:.2e}  |  "
            f"direction = {direction}  |  {r_info}  |  "
            f"k = {n}  |  {key_str}")

    def _on_results_gene_selected(self, gene: str):
        self._show_gene_forest(gene)
        self._tabs.setCurrentIndex(0)

    # --- Auto-save ---
    def _direction_params(self):
        """Significant and opposite-direction counts for provenance."""
        from kosmic.meta_analysis.direction import (
            count_significant_conflicts, studies_without_adjusted_p)
        if self._meta_df is None:
            return {}
        n_sig, n_conflict = count_significant_conflicts(
            self._meta_df, fdr=DEFAULT_FDR)
        if n_conflict is None:
            return {}
        return {
            'n_significant': n_sig,
            'n_significant_direction_conflict': n_conflict,
            'direction_study_fdr': DEFAULT_FDR,
            'direction_studies_without_adjusted_p': studies_without_adjusted_p(
                getattr(self, '_last_consensus_de_dfs', None) or []),
        }

    def _record_meta_provenance(self, stage):
        """Record the meta-analysis pooling step into a project-level sidecar so
        the combined cross-study methods view can show it. Best-effort."""
        if not self._project_folder:
            return
        try:
            from kosmic.meta_analysis.io import (
                gather_study_provenance, record_meta_stage)
            tokens = {
                label: (rec.get('fingerprint') or {}).get('token')
                for label, rec in gather_study_provenance(self._project_folder)
            }
            if_active = (self._if_check.isChecked()
                         if hasattr(self, '_if_check') else False)
            params = {
                'pooling_methods': list(self._method_keys or []),
                'calibration': ('cc_perm' if self._cal_combo.currentIndex() == 1
                                else 'analytical'),
                'mode': getattr(self, '_mode', None),
                'min_studies': self._min_studies_spin.value(),
                'independent_filter': if_active,
                'n_studies': len(self._run_study_names()),
                'studies': self._run_study_names(),
                'study_tokens': tokens,
            }
            params.update(self._direction_params())
            # Through the shared recorder, so this keeps only the
            # latest run like every other meta stage. Calling
            # record_stage directly here bypassed that and kept
            # appending one entry per Run click.
            record_meta_stage(self._project_folder, stage, params,
                              selection=self._run_output_selection)
        except Exception as e:
            self.log_message.emit(f"Could not record meta provenance: {e}")

    def _auto_save(self):
        """Save the pooled DataFrame + matching settings sidecar JSON.

        Both files share a stem keyed by methods + calibration; re-running
        the same combo overwrites the pair so settings stay in sync with
        the CSV.
        """
        if not self._project_folder or self._meta_df is None:
            return
        try:
            output_dir = meta_output_dir(self._project_folder,
                                         self._run_output_selection)
            output_dir.mkdir(parents=True, exist_ok=True)
            keys_str = '_'.join(self._method_keys) if self._method_keys else 'consensus'
            cal = ('cc' if self._cal_combo.currentIndex() == 1
                   else 'analytical')
            stem = f'consensus_{keys_str}_{cal}'

            csv_path = output_dir / f'{stem}.csv'
            export_df = self._meta_df.drop(
                columns=['study_effects'], errors='ignore')
            export_df.to_csv(csv_path, index=False)

            json_path = output_dir / f'{stem}.json'
            self._write_settings_sidecar(json_path)

            self.log_message.emit(
                f"Results saved to {output_dir.relative_to(self._project_folder).as_posix()}: "
                f"{csv_path.name} + {json_path.name}"
            )
            from kosmic.meta_analysis.io import mixed_selection_warning
            warning = mixed_selection_warning(self._run_output_selection)
            if warning:
                self.log_message.emit(warning)
        except Exception as e:
            self.log_message.emit(f"Auto-save failed: {e}")

    def _write_settings_sidecar(self, path):
        """Write a JSON sidecar capturing the run settings + result summary."""
        from kosmic import META_IF_FILTER_STAT
        from kosmic.meta_analysis.io import write_meta_settings_sidecar

        params = dict(getattr(self, '_last_consensus_params', {}) or {})
        # Drop large transient fields we don't want in the JSON.
        params.pop('pseudobulk_paths', None)

        if_active = self._if_check.isChecked() if hasattr(self, '_if_check') else False
        params['independent_filter'] = {
            'active':    if_active,
            'threshold': (self._if_pct_spin.value() / 100.0
                          if if_active and hasattr(self, '_if_pct_spin') else None),
            'stat':      META_IF_FILTER_STAT if if_active else None,
        }
        params['mode'] = getattr(self, '_mode', None)
        params['n_studies_input'] = len(self._run_study_names())
        params['study_names'] = self._run_study_names()
        if hasattr(self, '_pathway_gene_sets') and self._pathway_gene_sets:
            params['n_pathways_in_geneset'] = len(self._pathway_gene_sets)

        n_sig = (
            int((self._meta_df['fdr'] < DEFAULT_FDR).sum())
            if 'fdr' in self._meta_df.columns else None
        )

        write_meta_settings_sidecar(
            path,
            tool='KOSMIC Gene-Level Meta-Analysis',
            analysis_type='Multi-method consensus meta-analysis',
            parameters=params,
            results_summary={
                'n_genes_tested': len(self._meta_df),
                'n_significant':  n_sig,
                'fdr_threshold':  DEFAULT_FDR,
                'n_significant_direction_conflict': self._direction_params().get(
                    'n_significant_direction_conflict'),
            },
        )

    # --- Enrichment (exploratory mode) ---
    def _run_study_names(self):
        """Studies of the last run, as launched; the selection may since
        have changed (ADR-007)."""
        if self._last_consensus_labels is not None:
            return list(self._last_consensus_labels)
        return [d.get('name', '?') for d in (self._datasets or [])]

    @property
    def run_output_selection(self):
        """Output folder of the last pooling run (ADR-007), or None."""
        return self._run_output_selection

    def cache_loo_result(self, result):
        """Store a LOO result produced by the Validation step.

        The per-mode result cache lives here because this page is what
        restores state when the user switches between discovery and
        hypothesis mode.
        """
        self._loo_result = result
        if self._results_by_mode.get(self._mode) is not None:
            self._results_by_mode[self._mode]['loo_result'] = result

    def get_run_context(self):
        """Everything the Validation step needs from the last pooling run.

        Handed over as a dict rather than letting another page reach
        into this one's attributes: LOO has to refit on the *same*
        studies and settings for the held-out test to mean anything.
        """
        return {
            'meta_df':     self._meta_df,
            'method_keys': self._get_method_keys(),
            'mode':        self._mode,
            'top50':       self._top50_cb.isChecked(),
            'params':      self._last_consensus_params,
            'de_dfs':      self._last_consensus_de_dfs,
            'labels':      self._last_consensus_labels,
            'pb_paths':    self._last_consensus_pb_paths,
        }

    def set_hypothesis_pathways(self, pathways_dict):
        """Receive selected pathways from the PathwayExplorer page."""
        self._hyp_pathways = dict(pathways_dict)
        n_pw = len(pathways_dict)
        n_genes = len({g for genes in pathways_dict.values()
                       for g in genes})
        self.log_message.emit(
            f"Hypothesis: {n_pw} pathways, {n_genes} unique genes")

    def _hyp_selected_genes(self):
        """Union of genes across the selected pathway set."""
        gene_union = set()
        for genes in self._hyp_pathways.values():
            for g in genes:
                gene_union.add(g.upper())
        return gene_union

    def _run_hypothesis(self):
        """Run meta-analysis restricted to genes in the selected pathways."""
        gene_set = self._hyp_selected_genes()
        n_pathways = len(self._hyp_pathways)
        self.log_message.emit(
            f"Hypothesis: {n_pathways} pathways, "
            f"{len(gene_set)} genes")

        if not gene_set:
            dialogs.warning(
                self, "No Genes",
                "No genes found in the selected pathways.\n\n"
                "Go to the 'Select Pathways' step and choose a "
                "gene set first.")
            return

        method_keys = self._get_method_keys()
        self.log_message.emit(
            f"Hypothesis: method_keys = {method_keys}")
        if len(method_keys) < 1:
            self.log_message.emit("Hypothesis: aborting -- no methods selected")
            dialogs.warning(
                self, "No Methods Selected",
                "Select at least one pooling method.")
            return
        de_method_key = self._current_de_method_key()
        n_loaded, n_total_found, fallback = (
            self._reload_datasets_for_de_method(de_method_key))
        self.log_message.emit(
            f"Hypothesis: DE method = {de_method_key}; reloaded "
            f"{n_loaded} studies.")

        if not self._datasets:
            self.log_message.emit("Hypothesis: aborting -- no datasets")
            dialogs.warning(
                self, "No Data",
                f"No DE result files matching DE method "
                f"'{de_method_key}' found. Check that DE has been "
                f"run and saved for each study, or pick a different "
                f"DE method.")
            return

        if fallback:
            self.log_message.emit(
                f"  studies falling back to non-preferred DE method: "
                f"{', '.join(fallback)}")
            n_fb = len(fallback)
            n_match = n_loaded - n_fb
            msg = (
                f"DE method picked: {de_method_key}.\n\n"
                f"Of {n_loaded} loaded studies:\n"
                f"  {n_match} have files matching '{de_method_key}'\n"
                f"  {n_fb} fell back to a different DE method:\n     "
                + "\n     ".join(fallback) + "\n\n"
                "Continuing would mix DE methods across studies. "
                "Continue anyway with mixed methods?")
            if not dialogs.confirm(self, "Mixed DE methods detected", msg):
                self._status_label.setText(
                    f"Aborted: {n_fb} studies missing "
                    f"{de_method_key} DE files.")
                return
        self.log_message.emit(
            f"Hypothesis: {len(self._datasets)} datasets loaded")

        filtered_dfs = []
        for d in self._datasets:
            df = d['df']
            name_col = 'names' if 'names' in df.columns else df.columns[0]
            mask = df[name_col].str.upper().isin(gene_set)
            self.log_message.emit(
                f"  {d.get('name', '?')}: {mask.sum()}/{len(df)} pathway genes matched")
            if mask.sum() > 0:
                fd = d.copy()
                fd['df'] = df[mask].copy()
                filtered_dfs.append(fd)

        if not filtered_dfs:
            self.log_message.emit("Hypothesis: aborting -- no gene overlap")
            dialogs.warning(
                self, "No Overlap",
                "None of the pathway genes were found in the DE results.")
            return

        # Independent filter: drop low-cross-study-pct genes from the
        # pathway-restricted universe before pooling. Recomputes FDR on
        # the survivors (gain power on the surviving pathway genes).
        if self._if_check.isChecked():
            from kosmic import META_IF_FILTER_STAT
            from kosmic.meta_analysis.independent_filter import (
                apply_filter, filter_datasets_to_universe,
            )
            if_threshold = self._if_pct_spin.value() / 100.0
            gene_universe = apply_filter(
                filtered_dfs,
                threshold=if_threshold,
                stat=META_IF_FILTER_STAT,
                min_studies=self._min_studies_spin.value(),
            )
            n_pre = sum(len(d['df']) for d in filtered_dfs) // max(len(filtered_dfs), 1)
            filtered_dfs = filter_datasets_to_universe(
                filtered_dfs, gene_universe)
            self.log_message.emit(
                f"Independent filter: {len(gene_universe):,} genes pass "
                f"(threshold = {if_threshold:.0%}, stat = {META_IF_FILTER_STAT}; "
                f"~{n_pre:,} per-study before filter)"
            )
            if not any(len(d['df']) > 0 for d in filtered_dfs):
                self.log_message.emit(
                    "Hypothesis: aborting -- independent filter removed all genes")
                dialogs.warning(
                    self, "Filter Removed All Genes",
                    "The independent filter removed every gene from the "
                    "pathway universe. Lower the threshold or untick the "
                    "filter.")
                return

        cal_key = ['analytical', 'cc_perm'][self._cal_combo.currentIndex()]
        self.log_message.emit(f"Hypothesis: calibration = {cal_key}")

        # CC perm needs pseudobulk paths (same lookup as exploratory mode).
        pb_paths = []
        if cal_key == 'cc_perm':
            pb_paths = find_pseudobulk_paths(
                self._datasets, self._project_folder,
                log_cb=self.log_message.emit)
            self.log_message.emit(
                f"Hypothesis: found {len(pb_paths)}/{len(self._datasets)} "
                "pseudobulk files")
            if not pb_paths:
                self.log_message.emit(
                    "Hypothesis: aborting -- not all datasets have pseudobulk")
                dialogs.warning(
                    self, "No Pseudobulk",
                    "CC permutation requires pseudobulk CSV files alongside "
                    "the DE results for ALL datasets.")
                self._run_btn.setEnabled(True)
                return

        # Worker computes BH on whatever it receives, so the filtered
        # subset gets its own FDR universe automatically.
        n_pathways = len(self._hyp_pathways)

        self._run_btn.setEnabled(False)
        msg = (f"Hypothesis: pooling {len(gene_set)} genes from "
               f"{n_pathways} pathways using "
               f"{', '.join(method_keys)} ({cal_key})")
        self._status_label.setText(msg)
        self.log_message.emit(msg)
        if self.progress_bar:
            self.progress_bar.setRange(0, 0)  # indeterminate spinner

        de_dfs = [d['df'] for d in filtered_dfs]
        labels = [d['name'] for d in filtered_dfs]

        cc_de_method = {
            'deseq2':    'deseq2',
            'welch_cpm': 'welch_cpm',
            'welch_raw': 'welch_raw',
        }.get(de_method_key, 'deseq2')

        params = {
            'method_keys': method_keys,
            'calibration': cal_key,
            'min_studies': self._min_studies_spin.value(),
            'cc_n_perms': 1000,
            'pseudobulk_paths': pb_paths,
            'de_method':   cc_de_method,
        }

        # LOO worker re-uses these so the user doesn't have to re-pick methods.
        self._last_consensus_params = params
        self._last_consensus_pb_paths = pb_paths
        self._last_consensus_de_dfs = de_dfs
        self._last_consensus_labels = labels

        self._warn_unadjusted(de_dfs)
        # The folder is fixed when the run starts: the selection may
        # change while it runs (ADR-007).
        self._pending_output_selection = self._output_selection
        self._worker = GeneMAWorker(de_dfs, labels, params)
        run_worker(
            self._worker,
            on_finished=self._on_finished,
            on_failed=self._on_failed,
            on_progress=self._on_progress,
        )

