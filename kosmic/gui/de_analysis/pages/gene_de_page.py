# Gene DE page: pseudobulk DE, volcano plot, results table,
# pathway coverage.

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QComboBox, QCheckBox,
    QLineEdit, QCompleter, QListWidget, QListWidgetItem,
    QHeaderView, QApplication,
)
from PyQt6.QtCore import QSettings, Qt, pyqtSignal, QStringListModel, QTimer
from pathlib import Path
import numpy as np

from kosmic.gui.shared.theme import NoScrollDoubleSpinBox, NoScrollComboBox, NoScrollSpinBox
from kosmic.gui.shared.widgets import (
    Column, ResultsTableView, SidebarTabbedPage,
    SecondaryLabel, SectionHeader, PrimaryButton,
    StageAccordion, StageSummaryCard,
)
from kosmic.paths import (
    de_analysis_dir, de_result_path, de_stats_dir, significant_path,
)
from kosmic.de.de_analysis import required_donor_count
from kosmic.numerical import neg_log10
from kosmic.gui.shared.plots import InteractiveHeatmap, InteractiveVolcano
from kosmic.gui.de_analysis.workers import (
    CellTypeBatchWorker, MetabolicDEWorker, MetaExportWorker,
)
from kosmic import (
    DEFAULT_FDR, DEFAULT_LFC_THRESHOLD,
    DE_DETECTION_MIN_PCT, DE_DETECTION_ON, DE_DETECTION_MIN_DONOR_FRAC,
    DE_DESEQ2_INDEPENDENT_FILTER, DE_DESEQ2_COOKS_FILTER,
    DE_FILTER_MIN_COUNT, DE_FILTER_MIN_SAMPLES,
)
from kosmic.gui.shared import dialogs, run_worker


_DE_COLOR_SCHEMES = [
    ((220, 50, 50, 200), (50, 50, 220, 200)),     # red / blue
    ((230, 150, 30, 200), (30, 170, 160, 200)),   # orange / teal
    ((200, 50, 200, 200), (50, 200, 220, 200)),   # magenta / cyan
]


def _compute_de_volcano_data(
    de_results,
    *,
    pathway_genes=None,
    use_raw_pvals=False,
    pval_threshold=DEFAULT_FDR,
    fc_threshold=DEFAULT_LFC_THRESHOLD,
):
    """
    Build the arrays InteractiveVolcano expects from a DE-results frame.

    Returns '(x, y, names, brushes, sizes, y_label)'. Off-pathway
    points (when 'pathway_genes' is given) render small + transparent
    so they don't crowd the layer of interest.
    """
    import pyqtgraph as pg

    x = de_results['logfoldchanges'].values
    if use_raw_pvals and 'pvals' in de_results.columns:
        pval_col = de_results['pvals'].values
        y_label = '-Log10 Raw P-value'
    else:
        pval_col = de_results['pvals_adj'].values
        y_label = '-Log10 Adjusted P-value (FDR)'
    y = neg_log10(pval_col)
    names = de_results['names'].tolist()

    has_pathway_overlay = (pathway_genes is not None
                           and len(pathway_genes) > 0)

    s = QSettings("KOSMIC", "KOSMIC")
    sig_point_size = int(s.value("plot/volcano_point_size", 9))
    color_scheme = int(s.value("plot/volcano_color_scheme", 0))
    up_rgba, down_rgba = _DE_COLOR_SCHEMES[
        min(color_scheme, len(_DE_COLOR_SCHEMES) - 1)]

    # Reuse one brush object per colour category (not one per point) so the
    # scatter can be drawn as a few single-brush layers -- pyqtgraph's fast
    # path -- instead of 26k per-point symbols.
    b_offpath = pg.mkBrush(200, 200, 200, 80)
    b_up = pg.mkBrush(*up_rgba)
    b_down = pg.mkBrush(*down_rgba)
    b_ns = pg.mkBrush(150, 150, 150, 150)
    ns_size = sig_point_size - 2 if has_pathway_overlay else sig_point_size - 1

    brushes = []
    sizes = []
    for i, name in enumerate(names):
        lfc = x[i]
        in_pathway = not has_pathway_overlay or name in pathway_genes
        if not in_pathway:
            brushes.append(b_offpath)
            sizes.append(4)
        elif pval_col[i] < pval_threshold and abs(lfc) > fc_threshold:
            brushes.append(b_up if lfc > 0 else b_down)
            sizes.append(sig_point_size)
        else:
            brushes.append(b_ns)
            sizes.append(ns_size)

    return x, y, names, brushes, sizes, y_label


def _render_de_volcano(widget, de_results, *, pathway_genes=None,
                      use_raw_pvals=False,
                      pval_threshold=DEFAULT_FDR, fc_threshold=DEFAULT_LFC_THRESHOLD,
                      disease_label: str = '',
                      control_label: str = ''):
    """Populate an 'InteractiveVolcano' widget with a DE-results frame."""
    x, y, names, brushes, sizes, y_label = _compute_de_volcano_data(
        de_results,
        pathway_genes=pathway_genes,
        use_raw_pvals=use_raw_pvals,
        pval_threshold=pval_threshold,
        fc_threshold=fc_threshold,
    )
    widget.set_scatter(
        x, y, names, brushes=brushes, sizes=sizes,
        x_label='Log2 Fold Change',
        y_label=y_label,
        title='Volcano Plot',
    )
    widget.add_threshold_line(y=-np.log10(pval_threshold))
    if fc_threshold > 0:
        widget.add_threshold_line(x=fc_threshold)
        widget.add_threshold_line(x=-fc_threshold)
    # Centre the log2FC axis on zero (up and down symmetric) instead of letting
    # autoscale skew it, matching the meta-analysis volcanoes.
    widget.set_symmetric_x()
    # Direction-of-effect arrows pinned to top corners. Forest-plot
    # convention: positive log2FC = "up in disease" on the right.
    if disease_label:
        widget.set_direction_labels(
            up_text=f"Up in {disease_label}",
            down_text=f"Down in {disease_label}",
        )
    else:
        widget.set_direction_labels(None, None)


def _parse_rgba(rgba_str: str) -> tuple:
    """
    Parse an 'rgba(r, g, b, a)' / 'rgb(r, g, b)' theme colour.

    Returns '(r, g, b, a)' ints. Falls back to off-white when the
    string can't be parsed (theme regression safety).
    """
    import re
    m = re.match(r'rgba?\((\d+),\s*(\d+),\s*(\d+),?\s*(\d*)\)', rgba_str or '')
    if m:
        return (
            int(m.group(1)), int(m.group(2)), int(m.group(3)),
            int(m.group(4)) if m.group(4) else 255,
        )
    return (255, 255, 255, 180)


def _draw_grouped_bars_with_dots(
    plot,
    x_idx: int,
    c_vals,
    d_vals,
    *,
    ctrl_color: str,
    dis_color: str,
    fg_color: str,
    dot_rgba: tuple,
    width: float = 0.35,
    seed: int = 0,
):
    """Render a control-vs-disease bar (mean +/- SEM) with per-sample dots."""
    import pyqtgraph as pg

    c_vals = np.asarray(c_vals, dtype=float)
    d_vals = np.asarray(d_vals, dtype=float)

    c_mean = float(c_vals.mean()) if c_vals.size else 0.0
    d_mean = float(d_vals.mean()) if d_vals.size else 0.0
    c_sem = float(c_vals.std(ddof=1) / np.sqrt(len(c_vals))) if len(c_vals) > 1 else 0.0
    d_sem = float(d_vals.std(ddof=1) / np.sqrt(len(d_vals))) if len(d_vals) > 1 else 0.0

    plot.addItem(pg.BarGraphItem(
        x=[x_idx - width / 2], height=[c_mean], width=width,
        brush=pg.mkBrush(ctrl_color + 'cc'),
        pen=pg.mkPen(ctrl_color, width=0.5),
    ))
    plot.addItem(pg.BarGraphItem(
        x=[x_idx + width / 2], height=[d_mean], width=width,
        brush=pg.mkBrush(dis_color + 'cc'),
        pen=pg.mkPen(dis_color, width=0.5),
    ))

    for val, sem, xpos in [
        (c_mean, c_sem, x_idx - width / 2),
        (d_mean, d_sem, x_idx + width / 2),
    ]:
        if sem > 0:
            plot.addItem(pg.ErrorBarItem(
                x=np.array([xpos]), y=np.array([val]),
                top=np.array([sem]), bottom=np.array([sem]),
                beam=width * 0.3,
                pen=pg.mkPen(fg_color, width=1.5),
            ))

    rng = np.random.default_rng(seed)
    for vals, xpos in [
        (c_vals, x_idx - width / 2),
        (d_vals, x_idx + width / 2),
    ]:
        if vals.size > 0:
            jitter = rng.uniform(-width * 0.2, width * 0.2, vals.size)
            plot.addItem(pg.ScatterPlotItem(
                x=xpos + jitter, y=vals,
                size=6, pen=pg.mkPen(fg_color, width=0.5),
                brush=pg.mkBrush(*dot_rgba),
            ))


def attach_gene_completer(line_edit, model, parent):
    """Autocomplete each entry of a comma-separated gene list.

    QCompleter matches against the widget's entire text, so once a comma
    is typed the prefix becomes "POSTN, COL" and nothing ever matches --
    completion silently works for the first gene only. This drives the
    completer manually: the prefix is the text after the last comma, and
    an accepted suggestion replaces just that final entry.

    Returns the completer, which the caller must keep a reference to.
    """
    completer = QCompleter(model, parent)
    completer.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
    completer.setCompletionMode(QCompleter.CompletionMode.PopupCompletion)
    completer.setWidget(line_edit)

    def _on_edited(_text=None):
        prefix = line_edit.text().split(',')[-1].strip()
        if not prefix:
            completer.popup().hide()
            return
        completer.setCompletionPrefix(prefix)
        if completer.completionCount():
            completer.complete()
        else:
            completer.popup().hide()

    def _on_activated(choice):
        parts = line_edit.text().split(',')
        parts[-1] = choice
        line_edit.setText(', '.join(p.strip() for p in parts if p.strip()))
        line_edit.setCursorPosition(len(line_edit.text()))

    # textEdited, not textChanged: setting the text in code (clicking a
    # volcano point, say) must not pop the list open.
    line_edit.textEdited.connect(_on_edited)
    completer.activated.connect(_on_activated)
    return completer


class GeneDEPage(SidebarTabbedPage):
    """Step 2: pseudobulk DE execution, volcano, results table."""

    help_id = "de/gene_de"
    CONTENT_MARGINS = (8, 0, 0, 0)

    analysis_complete = pyqtSignal()  # emitted when DE finishes successfully

    def __init__(self, workspace):
        super().__init__()
        self.ws = workspace
        self.de_worker = None
        self._meta_export_worker = None
        self._last_run_comparison = None  # settings snapshot at last DE run
        self._setup_ui()
        self.ws.mode_changed.connect(self._on_mode_changed)

    def _setup_ui(self):
        left_layout = self.sidebar_layout

        # Sidebar follows the Quality Control pattern: a summary card per
        # concern, expanding in place to its settings. Four settings groups
        # stacked in one scroll made it impossible to see at a glance what a
        # run would do; the collapsed summaries say it in two lines each.
        left_layout.addWidget(SectionHeader("ANALYSIS SETUP"))
        self._accordion = StageAccordion(left_layout)

        model_page = QWidget()
        model_lay = QVBoxLayout(model_page)
        design_page = QWidget()
        design_lay = QVBoxLayout(design_page)
        filters_page = QWidget()
        filters_lay = QVBoxLayout(filters_page)
        fdr_page = QWidget()
        fdr_lay = QVBoxLayout(fdr_page)
        for _lay in (model_lay, design_lay, filters_lay, fdr_lay):
            _lay.setContentsMargins(0, 4, 0, 0)
            _lay.setSpacing(4)


        model_lay.addWidget(QLabel("Replicate unit:"))
        self.unit_combo = NoScrollComboBox()
        self.unit_combo.addItems([
            "Sample / donor (pseudobulk)",
            "Cell (scanpy Wilcoxon)",
        ])
        self.unit_combo.setCurrentIndex(0)  # Default to pseudobulk
        self.unit_combo.setToolTip(
            "Sample / donor (pseudobulk): aggregate cells into one profile per\n"
            "sample and test across samples. The statistically correct unit.\n\n"
            "Cell (scanpy Wilcoxon): rank_genes_groups on individual cells\n"
            "(the marker/dotplot method), disease vs control. Ignores donor\n"
            "grouping, so it is pseudoreplicated and inflates significance --\n"
            "use only when too few samples exist for pseudobulk."
        )
        self.unit_combo.currentIndexChanged.connect(self._on_unit_changed)
        model_lay.addWidget(self.unit_combo)

        self._unit_caveat = SecondaryLabel(
            "Cell mode treats each cell as a replicate (pseudoreplication): "
            "expect far more 'significant' genes than pseudobulk.")
        self._unit_caveat.setWordWrap(True)
        self._unit_caveat.setVisible(False)
        model_lay.addWidget(self._unit_caveat)

        model_lay.addWidget(QLabel("Pseudobulk method:"))
        self.de_method_combo = NoScrollComboBox()
        self.de_method_combo.addItems(["log-CPM (Welch)", "DESeq2 (NB-GLM)"])
        self.de_method_combo.setCurrentIndex(1)  # Default to DESeq2
        self.de_method_combo.setToolTip(
            "Both collapse cells to one value per donor, then test across donors; "
            "they differ in the normalisation + statistical model (each carries "
            "its own normalisation).\n"
            "log-CPM (Welch): pseudobulk sum, CPM-normalised, log2, then Welch's "
            "t-test per gene. Fast and assumption-light -- a simple cross-check.\n"
            "DESeq2 (NB-GLM): negative-binomial GLM on raw count sums with DESeq2 "
            "size factors, Wald test, dispersion + LFC shrinkage. Recommended for "
            "counts.\n"
            "(Ignored when Replicate unit = Cell, which always uses Wilcoxon.)"
        )
        self.de_method_combo.currentIndexChanged.connect(self._on_de_method_changed)
        model_lay.addWidget(self.de_method_combo)

        self.eb_checkbox = QCheckBox("Empirical Bayes moderation")
        self.eb_checkbox.setChecked(False)
        self.eb_checkbox.setEnabled(False)  # Disabled by default (DESeq2 selected)
        self.eb_checkbox.setToolTip(
            "Shrink per-gene variance toward a common prior (limma-style).\n"
            "Reduces false positives from genes with unusually low variance.\n"
            "Only applies to t-test method."
        )
        self.eb_checkbox.toggled.connect(self._on_eb_toggled)
        model_lay.addWidget(self.eb_checkbox)

        filters_lay.addWidget(QLabel("Donor filter: min cells per donor"))
        self.min_cells_spin = NoScrollSpinBox()
        self.min_cells_spin.setRange(1, 1000)
        self.min_cells_spin.setValue(self.ws.min_cells)
        self.min_cells_spin.setToolTip(
            "A donor contributing fewer cells than this is dropped from the\n"
            "comparison -- their pseudobulk profile would be summed from too\n"
            "few cells to mean anything.\n\n"
            "Counted over whatever is loaded. On a whole dataset that is all\n"
            "of a donor's cells; in a per-cell-type run it is only their\n"
            "cells OF THAT TYPE -- which is why this setting decides which\n"
            "cell types are testable at all. A donor with 5,000 cells but 30\n"
            "adipocytes fails it for the adipocyte run.\n\n"
            "10 is the convention (muscat's pbDS default). Above about 20 you\n"
            "start excluding whole cell types rather than thin donors.")
        self.min_cells_spin.valueChanged.connect(self._on_min_cells_changed)
        filters_lay.addWidget(self.min_cells_spin)

        filters_lay.addWidget(QLabel("Donor filter: min transcripts per donor"))
        self.min_counts_spin = NoScrollSpinBox()
        self.min_counts_spin.setRange(0, 10_000_000)
        self.min_counts_spin.setSingleStep(10_000)
        self.min_counts_spin.setGroupSeparatorShown(True)
        self.min_counts_spin.setSpecialValueText("off")
        self.min_counts_spin.setValue(self.ws.min_counts)
        self.min_counts_spin.setToolTip(
            "A donor whose cells sum to fewer transcripts than this is\n"
            "dropped from the comparison. 0 turns the floor off.\n\n"
            "The cell-count filter above does not guard depth: 10 shallow\n"
            "nuclei clear it but can sum to a few thousand transcripts,\n"
            "which is too thin a profile for a count model. Summed over\n"
            "every gene, for the same cells the cell-count filter sees.\n\n"
            "Gao et al. use 50,000 for heart snRNA-seq. Each donor's total\n"
            "is written to the pseudobulk file as 'total_counts', so you\n"
            "can see where a floor would fall before setting one.")
        self.min_counts_spin.valueChanged.connect(self._on_min_counts_changed)
        filters_lay.addWidget(self.min_counts_spin)


        # Built here, shown in the Advanced dialog rather than the sidebar:
        # all three are at standard values and an ordinary run never
        # touches them.
        # Two groups: the gene filter itself, and an optional
        # cell-level floor that is off by default. The dialog
        # heads them separately -- lumping them together was
        # why nobody could tell which was which.
        self._detection_rows = []
        self._gene_filter_rows = []
        self._detection_rows = [("Detected in \u2265 (% of a donor's cells)",
                                 None)]
        self.detection_pct_spin = NoScrollDoubleSpinBox()
        self.detection_pct_spin.setRange(0, 100)
        self.detection_pct_spin.setValue(
            DE_DETECTION_MIN_PCT * 100 if DE_DETECTION_ON else 0.0)
        self.detection_pct_spin.setSingleStep(1)
        self.detection_pct_spin.setSuffix(" %")
        self.detection_pct_spin.setToolTip(
            "Dropout floor: a gene counts as present in a donor when it is\n"
            "non-zero in at least this share of that donor's cells.\n\n"
            "Off by default. Standard pseudobulk DE has no such filter --\n"
            "this is the idea behind Seurat's FindMarkers 'min.pct', which\n"
            "guards a cell-level test. Above about 1% it stops pre-trimming\n"
            "and starts overriding the donor-level gene filter.")
        self._detection_rows[-1] = (self._detection_rows[-1][0],
                                    self.detection_pct_spin)

        self._detection_rows.append(
            ("...in \u2265 (% of a group's donors)", None))
        self.detection_donor_frac_spin = NoScrollDoubleSpinBox()
        self.detection_donor_frac_spin.setRange(0, 100)
        self.detection_donor_frac_spin.setValue(DE_DETECTION_MIN_DONOR_FRAC * 100)
        self.detection_donor_frac_spin.setSingleStep(10)
        self.detection_donor_frac_spin.setSuffix(" %")
        self.detection_donor_frac_spin.setToolTip(
            "How many of a group's donors must each clear the percentage\n"
            "above, counted within disease and within control separately.\n\n"
            "Without this the rate is pooled over a group's cells, so one\n"
            "donor expressing a gene in 30% of their cells outvotes eleven\n"
            "donors at 1% and the gene is tested on the strength of one\n"
            "person. 0 restores that pooled behaviour.")
        self._detection_rows[-1] = (self._detection_rows[-1][0],
                                    self.detection_donor_frac_spin)

        self.detection_study_check = QCheckBox("Filter genes within each study")
        self.detection_study_check.setChecked(True)
        self.detection_study_check.setToolTip(
            "A methodological choice, not a tuning knob -- which is why it\n"
            "is on the page rather than behind Advanced.\n\n"
            "OFF (conventional): the gene filter is applied once to the\n"
            "pooled donors x genes table. This is what a standalone\n"
            "mega-analysis normally does, and it keeps genes that only one\n"
            "cohort measures well.\n\n"
            "ON (default): the filter runs inside each study and a gene\n"
            "must pass in all of them. Use this when the pooled result will\n"
            "be compared against per-study DE -- otherwise the pooled arm\n"
            "tests genes the per-study arm was never offered, and the\n"
            "difference between them is partly a difference in what each\n"
            "was allowed to look at. On the DCM endothelial data that was\n"
            "592 genes.\n\n"
            "Either way it is recorded in the run's provenance.")


        self.filter_min_count_spin = NoScrollDoubleSpinBox()
        self.filter_min_count_spin.setRange(0, 1000)
        self.filter_min_count_spin.setDecimals(0)
        self.filter_min_count_spin.setValue(float(DE_FILTER_MIN_COUNT))
        self.filter_min_count_spin.setSingleStep(1)
        self._gene_filter_rows.append(("Counts in a typical donor", None))
        self.filter_min_count_spin.setToolTip(
            "How much of the gene has to be there. Not a raw count: it sets\n"
            "a CPM cutoff of (this number / median donor library) x 1e6, so\n"
            "the same threshold means the same thing in a deep donor and a\n"
            "shallow one.\n\n"
            "10 against a 3.3M-count library works out to about 3 CPM.\n"
            "0 disables the count cutoff.")


        self.filter_min_samples_spin = NoScrollDoubleSpinBox()
        self.filter_min_samples_spin.setRange(0, 1000)
        self.filter_min_samples_spin.setDecimals(0)
        self.filter_min_samples_spin.setValue(float(DE_FILTER_MIN_SAMPLES))
        self.filter_min_samples_spin.setSingleStep(1)
        # 0 is not a count, it is 'work it out' -- so say so in the box
        # rather than in the caption. Kept current by _refresh_card_summaries.
        self.filter_min_samples_spin.setSpecialValueText("auto")
        self._gene_filter_rows[-1] = (self._gene_filter_rows[-1][0],
                                      self.filter_min_count_spin)
        self._gene_filter_rows.append(
            ("...in at least (donors)", None))
        self.filter_min_samples_spin.setToolTip(
            "How many donors must reach the level above, counted across the\n"
            "whole cohort rather than within a group -- so a gene present in\n"
            "one arm only still passes.\n\n"
            "0 = auto: the smaller group's size, eased for large cohorts\n"
            "by 10 + (n - 10) x 0.7, which is where a fractional requirement\n"
            "like 11.4 of 24 comes from.")

        self._gene_filter_rows[-1] = (self._gene_filter_rows[-1][0],
                                      self.filter_min_samples_spin)

        model_lay.addWidget(QLabel("Count source:"))
        self.count_source_combo = NoScrollComboBox()
        self.count_source_combo.addItems(["Raw counts", "Decontaminated (DecontX)"])
        self.count_source_combo.setCurrentIndex(0)
        self.count_source_combo.setToolTip(
            "Raw counts: use the original counts.\n"
            "Decontaminated: use DecontX-corrected counts "
            "(layers['decontX_counts']).\n"
            "Available only after running the Decontaminate step before Subset."
        )
        model_lay.addWidget(self.count_source_combo)

        filters_lay.addWidget(self.detection_study_check)

        self._advanced_btn = QPushButton("Gene filter settings...")
        self._advanced_btn.setToolTip(
            "The two numbers behind the sentence above, plus an optional\n"
            "cell-level floor that is off by default. All are at standard\n"
            "values; an ordinary run changes none of them.")
        self._advanced_btn.clicked.connect(self._show_advanced_filters)
        filters_lay.addWidget(self._advanced_btn)

        is_deseq2 = self.ws.de_method == 'deseq2'

        # Sample-level covariates for the DESeq2 design. 'study' is the
        # one this exists for: fitting one model over pooled cohorts
        # without it leaves the cohort offset in the residual, inflating
        # dispersion and shrinking every effect toward zero.
        design_lay.addWidget(SecondaryLabel("Adjust for"))
        self.covariate_list = QListWidget()
        self.covariate_list.setMaximumHeight(90)
        self.covariate_list.setEnabled(is_deseq2)
        self.covariate_list.setToolTip(
            "Extra sample-level columns to put in the DESeq2 design,\n"
            "before condition. Tick 'study' when the dataset pools\n"
            "several cohorts." + "\n\n"
            "A column that is constant, unique per sample, or that\n"
            "explains the same split as condition cannot be fitted\n"
            "and is dropped, with a note in the log.")
        design_lay.addWidget(self.covariate_list)

        self.indep_filter_check = QCheckBox("DESeq2 independent filtering")
        self.indep_filter_check.setChecked(DE_DESEQ2_INDEPENDENT_FILTER)
        self.indep_filter_check.setEnabled(is_deseq2)
        self.indep_filter_check.setToolTip(
            "DESeq2 drops low-mean genes from the BH correction to gain power "
            "(Bourgon independent filtering). Filtered genes get a blank FDR and "
            "an 'indep. filter' tag in the results, not a misleading FDR of 1. "
            "Uncheck to correct over every tested gene.")
        fdr_lay.addWidget(self.indep_filter_check)

        self.cooks_filter_check = QCheckBox("Cook's outlier filtering")
        self.cooks_filter_check.setChecked(DE_DESEQ2_COOKS_FILTER)
        self.cooks_filter_check.setEnabled(is_deseq2)
        self.cooks_filter_check.setToolTip(
            "DESeq2 blanks a gene's p-value when one sample is an extreme "
            "Cook's-distance outlier. Tagged \"Cook's\" in the results. Uncheck "
            "to keep those genes (useful with few replicates).")
        fdr_lay.addWidget(self.cooks_filter_check)

        # Any DE-affecting control changing after a run makes the shown results
        # stale -> the Run button flips to "Re-run DE Analysis".
        for _sig in (self.detection_pct_spin.valueChanged,
                     self.detection_donor_frac_spin.valueChanged,
                     self.filter_min_count_spin.valueChanged,
                     self.filter_min_samples_spin.valueChanged,
                     self.count_source_combo.currentIndexChanged,
                     self.detection_study_check.toggled,
                     self.indep_filter_check.toggled,
                     self.cooks_filter_check.toggled):
            _sig.connect(self._refresh_run_dirty)

        # Restate the whole rule in one sentence whenever a piece changes.
        for _sig in (self.min_cells_spin.valueChanged,
                     self.min_counts_spin.valueChanged,
                     self.detection_pct_spin.valueChanged,
                     self.detection_donor_frac_spin.valueChanged,
                     self.filter_min_count_spin.valueChanged,
                     self.filter_min_samples_spin.valueChanged,
                     self.detection_study_check.toggled):
            _sig.connect(self._refresh_filter_summary)

        # Collapsed card summaries track every control they describe.
        for _sig in (self.min_cells_spin.valueChanged,
                     self.min_counts_spin.valueChanged,
                     self.detection_pct_spin.valueChanged,
                     self.detection_donor_frac_spin.valueChanged,
                     self.filter_min_count_spin.valueChanged,
                     self.filter_min_samples_spin.valueChanged,
                     self.detection_study_check.toggled,
                     self.count_source_combo.currentIndexChanged,
                     self.de_method_combo.currentIndexChanged,
                     self.unit_combo.currentIndexChanged,
                     self.eb_checkbox.toggled,
                     self.indep_filter_check.toggled,
                     self.cooks_filter_check.toggled):
            _sig.connect(self._refresh_card_summaries)
        self.covariate_list.itemChanged.connect(self._refresh_card_summaries)

        # The primary action lives at the foot of the Settings group -- no
        # separate box needed.
        self.run_btn = PrimaryButton("Run DE Analysis")
        self.run_btn.setToolTip(
            "Run pseudobulk DE using the selected gene scope.\n"
            "Volcano plot highlights pathway genes; results table shows all tested."
        )
        self.run_btn.clicked.connect(self._run_analysis)
        self.run_btn.setEnabled(False)

        self.batch_btn = QPushButton("Run per cell type...")
        self.batch_btn.setToolTip(
            "Run this same contrast separately for each cell type, with"
            " the settings above, and write one result per type." + "\n\n" +
            "A disease effect in endothelium and one in cardiomyocytes"
            " are different contrasts; pooling them into a single"
            " pseudobulk per donor averages them away.")
        self.batch_btn.clicked.connect(self._run_batch_by_cell_type)
        self.batch_btn.setEnabled(False)

        # Not an output -- a check on the assumption the model rests on.
        # DESeq2's median-of-ratios normalisation assumes the two arms have
        # the same transcriptome size per cell; this says whether they do.
        # It runs automatically after every DE and logs its verdict, so the
        # button only reopens the plot.
        self._transcriptome_btn = QPushButton("Transcriptome size check...")
        self._transcriptome_btn.setToolTip(
            "Does total per-cell expression differ between the two arms?\n"
            "Shows three metrics per donor: total UMI, genes detected,\n"
            "UMI per gene.\n\n"
            "A ratio beyond 0.9-1.1 means CPM and DESeq2 normalisation are\n"
            "standing on a shaky assumption and may mask a real global\n"
            "shift. Runs itself after each DE; this reopens the plot.")
        self._transcriptome_btn.clicked.connect(self._on_transcriptome_btn)
        self._transcriptome_btn.setEnabled(False)

        self.status_label = SecondaryLabel("")
        self.status_label.setWordWrap(True)

        for _lay in (model_lay, design_lay, filters_lay, fdr_lay):
            _lay.addStretch()

        # What the dataset IS comes before every setting, because it
        # decides most of them: whether 'study' belongs in the design,
        # whether detection has to be applied per study, and whether a
        # batch effect is expected or alarming.
        # Two lines, no expander: it states what was loaded and there is
        # nothing behind it to set.
        self._dataset_card = StageSummaryCard("Dataset", icon="dataset")
        self._dataset_card.set_summary(["No data loaded", ""])
        left_layout.addWidget(self._dataset_card)

        self._model_card = StageSummaryCard("Model", icon="bar-chart-2")
        self._model_card.set_settings_widget(model_page)
        left_layout.addWidget(self._model_card)
        self._accordion.add_card('model', self._model_card)

        self._design_card = StageSummaryCard("Design", icon="layers")
        self._design_card.set_settings_widget(design_page)
        left_layout.addWidget(self._design_card)
        self._accordion.add_card('design', self._design_card)

        self._filters_card = StageSummaryCard("Filters", icon="filter-funnel")
        self._filters_card.set_settings_widget(filters_page)
        left_layout.addWidget(self._filters_card)
        self._accordion.add_card('filters', self._filters_card)

        self._fdr_card = StageSummaryCard("FDR handling", icon="circle-dot")
        self._fdr_card.set_settings_widget(fdr_page)
        left_layout.addWidget(self._fdr_card)
        self._accordion.add_card('fdr', self._fdr_card)

        left_layout.addStretch()
        self._accordion.finalize()

        # Actions live below the accordion, always reachable -- collapsing a
        # card must never hide the Run button.
        left_layout.addWidget(self.run_btn)
        left_layout.addWidget(self.batch_btn)
        left_layout.addWidget(self._transcriptome_btn)
        left_layout.addWidget(self.status_label)

        self._refresh_card_summaries()


        # Display filters and gene lookup are built here but live in the
        # control strip under the volcano, not in this sidebar -- they act
        # on what the plot and table show, so they belong beside the plot.
        # Scope filters trigger a BH recompute and are batched behind
        # Apply; find-gene stays live.

        # Detection-based filtering happens pre-DESeq2 (the detection
        # spinboxes in Settings) and again inside DESeq2's own independent
        # filtering, so no post-hoc pct filter lives here.
        self._pathway_only_check = QCheckBox("Pathway genes only")
        self._pathway_only_check.setChecked(True)
        self._pathway_only_check.setToolTip(
            "Display only: hides the grey non-pathway dots on the volcano and "
            "non-pathway rows in the table. Uncheck to show all genes. FDR and "
            "significance are unchanged -- they always reflect the full tested "
            "gene set.")
        self._pathway_only_check.toggled.connect(self._on_display_filter_toggled)

        self._raw_pval_check = QCheckBox("Raw p-values")
        self._raw_pval_check.setToolTip(
            "Show raw (uncorrected) p-values instead of FDR-adjusted. "
            "Applies to the volcano + results table.")
        self._raw_pval_check.toggled.connect(self._on_display_filter_toggled)

        self._apply_filters_btn = QPushButton("Filters up to date")
        self._apply_filters_btn.setEnabled(False)
        self._apply_filters_btn.setToolTip(
            "Apply pending filter changes to the volcano and results table.")
        self._apply_filters_btn.clicked.connect(self._apply_filters)
        self._filters_dirty = False

        self._filter_count_label = SecondaryLabel("No DE results yet.")

        self._volcano_gene_search = QLineEdit()
        self._volcano_gene_search.setPlaceholderText("Gene name(s), comma-separated")
        self._volcano_gene_search.setMaximumWidth(240)
        self._volcano_gene_search.setToolTip(
            "Highlight matching genes on the volcano with a magenta "
            "star and filter the results table to just those genes. "
            "Comma-separated for multiple (e.g. VEGFA, KDR, FLT1). "
            "Applied live -- no Apply click needed.")
        self._volcano_gene_search_model = QStringListModel(self)
        self._volcano_gene_completer = attach_gene_completer(
            self._volcano_gene_search, self._volcano_gene_search_model, self)
        self._volcano_gene_search.textChanged.connect(self._on_gene_lookup_changed)

        left_layout.addStretch()

        results_tabs = self.tabs

        # --- Volcano tab ---
        volcano_widget = QWidget()
        volcano_layout = QVBoxLayout(volcano_widget)
        self.interactive_volcano = InteractiveVolcano()
        self.interactive_volcano.gene_selected.connect(self._on_gene_selected)
        volcano_layout.addWidget(self.interactive_volcano, 1)

        # The plot's own usage hint moves onto the plot, so the strip below
        # it has room for the controls that act on what it shows.
        self.interactive_volcano.setToolTip(
            "Hover a point for details; click to select that gene in the "
            "results table.")

        volcano_controls = QHBoxLayout()
        volcano_controls.setSpacing(8)

        # After a per-cell-type batch there are N result sets and only one
        # volcano. Without this the run wrote its CSVs and the plot went on
        # showing whatever preceded it.
        self._batch_view_label = QLabel("Cell type:")
        self._batch_view_combo = NoScrollComboBox()
        self._batch_view_combo.setMinimumWidth(180)
        self._batch_view_combo.setToolTip(
            "Which cell type's results to show. Set by the last per-cell-type "
            "run; every cell type was also written to its own CSV.")
        self._batch_view_combo.currentIndexChanged.connect(
            self._on_batch_view_changed)
        self._batch_view_label.hide()
        self._batch_view_combo.hide()
        volcano_controls.addWidget(self._batch_view_label)
        volcano_controls.addWidget(self._batch_view_combo)
        volcano_controls.addSpacing(12)

        volcano_controls.addWidget(QLabel("Find gene:"))
        volcano_controls.addWidget(self._volcano_gene_search)
        volcano_controls.addSpacing(12)
        volcano_controls.addWidget(self._pathway_only_check)
        volcano_controls.addWidget(self._raw_pval_check)
        volcano_controls.addWidget(self._apply_filters_btn)
        volcano_controls.addStretch()
        volcano_controls.addWidget(self._filter_count_label)
        volcano_controls.addSpacing(12)

        # Explicit reset button -- pyqtgraph's "View All" right-click is
        # not discoverable.
        self._reset_volcano_btn = QPushButton("Reset view")
        self._reset_volcano_btn.setToolTip(
            "Auto-range the volcano back to fit all points "
            "(undoes any zoom/pan).")
        self._reset_volcano_btn.clicked.connect(
            lambda: self.interactive_volcano.autoRange())
        volcano_controls.addWidget(self._reset_volcano_btn)

        volcano_layout.addLayout(volcano_controls)

        results_tabs.addTab(volcano_widget, "Volcano Plot")

        # --- Results table tab ---
        results_widget = QWidget()
        results_layout_inner = QVBoxLayout(results_widget)

        filter_layout = QHBoxLayout()
        filter_layout.addWidget(QLabel("Filter:"))

        self.pval_filter = NoScrollDoubleSpinBox()
        self.pval_filter.setRange(0, 1)
        self.pval_filter.setValue(DEFAULT_FDR)
        self.pval_filter.setSingleStep(0.01)
        self.pval_filter.setPrefix("FDR < ")
        self.pval_filter.setToolTip(
            "Threshold applied to FDR or raw p-values depending on the "
            "'Raw p-values' toggle next to the volcano. Sig flag flips "
            "to use whichever is active.")
        self.pval_filter.valueChanged.connect(self._mark_filters_dirty)
        filter_layout.addWidget(self.pval_filter)

        self.fc_filter = NoScrollDoubleSpinBox()
        self.fc_filter.setRange(0, 10)
        self.fc_filter.setValue(DEFAULT_LFC_THRESHOLD)
        self.fc_filter.setSingleStep(0.1)
        self.fc_filter.setPrefix("|log2FC| > ")
        self.fc_filter.setToolTip(
            "Keep genes whose ESTIMATED |log2FC| exceeds this."
            + "\n\n" +
            "A display cut on the point estimate, not a test: the FDR"
            " beside it answers 'did this gene change?', not"
            " 'did it change by more than this much?'."
            + "\n\n" +
            "Read the 95% CI column before describing a gene as"
            " changing by more than N-fold -- an estimate of 2.9 with"
            " an interval of [1.4, 4.3] supports a weaker claim than"
            " the 2.9 suggests.")
        self.fc_filter.valueChanged.connect(self._mark_filters_dirty)
        filter_layout.addWidget(self.fc_filter)

        self.sig_count_label = QLabel("")
        filter_layout.addWidget(self.sig_count_label)
        filter_layout.addStretch()
        results_layout_inner.addLayout(filter_layout)

        self.results_table = ResultsTableView()
        self.results_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.results_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self.results_table.setAlternatingRowColors(True)
        self.results_table.setSortingEnabled(True)
        results_layout_inner.addWidget(self.results_table)

        results_tabs.addTab(results_widget, "DE Results")

        # --- Heatmap tab ---
        heatmap_widget = QWidget()
        heatmap_layout = QVBoxLayout(heatmap_widget)

        heatmap_controls = QHBoxLayout()
        heatmap_controls.addWidget(QLabel("Rows:"))
        # Was "Pathway:", which left the tab dead in discovery mode --
        # there are no pathway sets, so the combo was empty and the
        # heatmap never drew. The rows are genes; where they come from
        # is the choice.
        self.heatmap_pathway_combo = QComboBox()
        # No stretch: its contents are short and fixed, whereas the gene
        # list beside it is the thing that grows. Giving the combo the
        # stretch left the gene box at its size hint -- a stub you could
        # not read three symbols in.
        self.heatmap_pathway_combo.setMinimumWidth(180)
        self.heatmap_pathway_combo.setMaximumWidth(260)
        self.heatmap_pathway_combo.setToolTip(
            "Which genes to draw as rows. Top-N comes from the current "
            "results, ordered by FDR; pathway entries appear when gene "
            "sets are loaded.")
        self.heatmap_pathway_combo.currentTextChanged.connect(
            self._generate_heatmap)
        heatmap_controls.addWidget(self.heatmap_pathway_combo)

        heatmap_controls.addWidget(QLabel("or genes:"))
        self._heatmap_gene_input = QLineEdit()
        self._heatmap_gene_input.setPlaceholderText(
            "POSTN, COL1A1, SFRP4 ...")
        self._heatmap_gene_input.setMinimumWidth(260)
        self._heatmap_gene_input.setToolTip(
            "Draw these genes instead. Comma-separated; overrides the "
            "selection on the left while it has text in it.")
        self._heatmap_gene_model = QStringListModel(self)
        self._heatmap_gene_completer = attach_gene_completer(
            self._heatmap_gene_input, self._heatmap_gene_model, self)
        self._heatmap_gene_input.editingFinished.connect(
            lambda: self._generate_heatmap(
                self.heatmap_pathway_combo.currentText()))
        heatmap_controls.addWidget(self._heatmap_gene_input, 1)

        heatmap_layout.addLayout(heatmap_controls)

        self.interactive_heatmap = InteractiveHeatmap()
        # The colormap hint used to sit on the controls row, where its
        # ~200px of text squeezed the gene list down to its minimum.
        # It belongs on the plot it describes.
        self.interactive_heatmap.setToolTip(
            "Click a gene to select it. Colormap and clipping are under "
            "View -> Plot Settings.")
        self.interactive_heatmap.gene_selected.connect(self._on_gene_selected)
        heatmap_layout.addWidget(self.interactive_heatmap, 1)
        results_tabs.addTab(heatmap_widget, "Gene Heatmap")

        # --- Gene Expression tab ---
        ge_widget = QWidget()
        ge_layout = QVBoxLayout(ge_widget)

        ge_controls = QHBoxLayout()
        ge_controls.addWidget(QLabel("Gene:"))
        self._ge_gene_input = QLineEdit()
        self._ge_gene_input.setPlaceholderText("e.g. GAPDH, HK1, VEGFA")
        self._ge_gene_input.setToolTip(
            "Enter one or more gene symbols, comma-separated.\n"
            "Starts completing from the genes in the current results, so\n"
            "you can find one without leaving for the volcano.")
        # Same autocomplete as the volcano's Find gene: typing a gene was
        # otherwise a blind guess against 12,000 names, and the only way
        # to fill this box was to click a point on another tab.
        self._ge_gene_model = QStringListModel(self)
        self._ge_gene_completer = attach_gene_completer(
            self._ge_gene_input, self._ge_gene_model, self)
        self._ge_gene_input.returnPressed.connect(self._run_gene_expression)
        ge_controls.addWidget(self._ge_gene_input, 1)
        self._ge_run_btn = QPushButton("Plot")
        self._ge_run_btn.clicked.connect(self._run_gene_expression)
        self._ge_run_btn.setEnabled(False)
        ge_controls.addWidget(self._ge_run_btn)
        ge_layout.addLayout(ge_controls)

        from kosmic.gui.shared.plots import InteractivePlot
        self._ge_mean_plot = InteractivePlot(
            title='Expression per Donor',
            left_label='Expression',
            unavailable_message='Click a gene in the volcano to see its expression.',
        )
        self._ge_mean_plot.setMinimumHeight(200)
        ge_layout.addWidget(self._ge_mean_plot)

        self._ge_pct_plot = InteractivePlot(
            title='% Cells Expressing Gene',
            left_label='% Expressing',
            unavailable_message='Click a gene in the volcano to see expression coverage.',
        )
        self._ge_pct_plot.setMinimumHeight(200)
        ge_layout.addWidget(self._ge_pct_plot)

        results_tabs.addTab(ge_widget, "Gene Expression")

        # --- Pathway Gene Expression tab ---
        pw_gene_widget = QWidget()
        pw_gene_layout = QVBoxLayout(pw_gene_widget)
        pw_gene_controls = QHBoxLayout()
        pw_gene_controls.addWidget(QLabel("Pathway:"))
        self._pw_gene_combo = NoScrollComboBox()
        self._pw_gene_combo.currentTextChanged.connect(self._on_pw_gene_combo_changed)
        pw_gene_controls.addWidget(self._pw_gene_combo, 1)
        pw_gene_controls.addStretch()
        pw_gene_layout.addLayout(pw_gene_controls)

        self._pw_gene_bar_plot = InteractivePlot(
            title='Per-Gene Expression',
            left_label='Expression',
            unavailable_message='Pick a pathway above to see per-gene expression.',
        )
        pw_gene_layout.addWidget(self._pw_gene_bar_plot, 1)
        self._pw_gene_tab_idx = results_tabs.addTab(
            pw_gene_widget, "Pathway Genes")

        # --- Top DE Genes tab ---
        from kosmic.gui.shared.plots import InteractiveTopDEPlot
        topde_tab = QWidget()
        topde_layout = QVBoxLayout(topde_tab)

        self._topde_plot = InteractiveTopDEPlot()
        self._topde_plot.gene_selected.connect(self._on_gene_selected)
        self._topde_plot.level_changed.connect(self._on_topde_level_changed)
        topde_layout.addWidget(self._topde_plot, 1)

        topde_controls = QHBoxLayout()
        topde_controls.addWidget(QLabel("Top N:"))
        self._topde_n_spin = NoScrollSpinBox()
        self._topde_n_spin.setRange(5, 50)
        self._topde_n_spin.setValue(15)
        self._topde_n_spin.valueChanged.connect(self._on_topde_n_changed)
        topde_controls.addWidget(self._topde_n_spin)

        self._topde_back_btn = QPushButton("Back to Pathways")
        self._topde_back_btn.clicked.connect(lambda: self._topde_plot.go_back())
        self._topde_back_btn.setVisible(False)
        topde_controls.addWidget(self._topde_back_btn)

        reset_topde_btn = QPushButton("Reset View")
        reset_topde_btn.clicked.connect(self._topde_plot.reset_view)
        topde_controls.addWidget(reset_topde_btn)
        topde_controls.addStretch()
        info_topde = SecondaryLabel("Click pathway to drill into genes, click gene to select")
        topde_controls.addWidget(info_topde)
        topde_layout.addLayout(topde_controls)
        self._topde_tab_idx = results_tabs.addTab(topde_tab, "Top DE Genes")

        self._results_tabs = results_tabs
        # Draw the expression-matrix-backed tabs only when opened, so their
        # cost is never paid on page entry.
        results_tabs.currentChanged.connect(self._on_results_tab_changed)

    def _on_results_tab_changed(self, index):
        """Lazily draw the Heatmap / Pathway Genes tabs on first view (they
        compute the per-donor expression matrix, which shouldn't run on entry)."""
        title = self._results_tabs.tabText(index)
        if title == "Gene Heatmap" and self.heatmap_pathway_combo.count():
            self._generate_heatmap(self.heatmap_pathway_combo.currentText())
        elif title == "Pathway Genes" and self._pw_gene_combo.count():
            self._on_pw_gene_combo_changed(self._pw_gene_combo.currentText())

    # --- Apply-button gating for filter changes ---

    def _mark_filters_dirty(self, *_args):
        """Flag pending filter changes and light up the Apply button."""
        self._filters_dirty = True
        if hasattr(self, '_apply_filters_btn') and self._apply_filters_btn is not None:
            self._apply_filters_btn.setEnabled(True)
            self._apply_filters_btn.setText("Apply filters ●")
            # Bold + accent so the pending state is visually obvious.
            self._apply_filters_btn.setProperty("role", "primary")
            style = self._apply_filters_btn.style()
            if style is not None:
                style.unpolish(self._apply_filters_btn)
                style.polish(self._apply_filters_btn)

    def _apply_filters(self):
        """Commit pending filter changes to the volcano + results table.

        The significant-gene CSV is defined by the FDR / log2FC
        thresholds, so it is rewritten here too -- a stale file on disk
        claiming different thresholds is worse than none.
        """
        self._filter_results()
        self._on_pval_toggle()
        self._auto_export_significant()
        # _on_pval_toggle re-rendered the scatter, so re-apply the highlight.
        if hasattr(self, '_volcano_gene_search'):
            self._on_volcano_gene_search(self._volcano_gene_search.text())
        self._filters_dirty = False
        if hasattr(self, '_apply_filters_btn') and self._apply_filters_btn is not None:
            self._apply_filters_btn.setEnabled(False)
            self._apply_filters_btn.setText("Filters up to date")
            self._apply_filters_btn.setProperty("role", "")
            style = self._apply_filters_btn.style()
            if style is not None:
                style.unpolish(self._apply_filters_btn)
                style.polish(self._apply_filters_btn)

    def _on_gene_lookup_changed(self, _text=None):
        """
        Live gene-name lookup: refilter the table + restamp highlight stars
        (no scatter re-render, no FDR recompute).
        """
        self._filter_results()
        if hasattr(self, '_volcano_gene_search'):
            self._on_volcano_gene_search(self._volcano_gene_search.text())

    # --- Pathway-UI visibility ---

    def _refresh_count_source(self):
        """Enable the DecontX count source only when the layer is present."""
        adata = self.ws.current_adata
        has_decontx = (adata is not None
                       and 'decontX_counts' in getattr(adata, 'layers', {}))
        model = self.count_source_combo.model()
        model.item(1).setEnabled(has_decontx)
        if not has_decontx and self.count_source_combo.currentIndex() == 1:
            self.count_source_combo.setCurrentIndex(0)

    def _selected_counts_layer(self):
        """Return the counts layer name to feed DE, or None for raw counts."""
        if self.count_source_combo.currentIndex() == 1:
            return 'decontX_counts'
        return None

    def _refresh_pathway_ui_visibility(self):
        """Hide the pathway-driven views when there are no pathways.

        The pathway-only display filter usefully narrows a genome-wide
        discovery result; in hypothesis mode the results are already
        restricted to the committed list, making the toggle a no-op.

        Two whole tabs are pathway-driven and are hidden outright without
        gene sets. 'Pathway Genes' is self-explanatory; 'Top DE Genes' is
        less obviously so -- it is a pathway drill-down (click a pathway
        to open its genes, with a 'Back to Pathways' control), not a
        top-N gene list, so in discovery mode it has nothing to draw.
        """
        has_pathways = bool(getattr(self.ws, 'pathway_gene_sets', None))
        discovery = getattr(self.ws, 'analysis_mode', None) != 'scoring'
        if hasattr(self, '_pathway_only_check'):
            self._pathway_only_check.setVisible(has_pathways and discovery)

        for attr in ('_pw_gene_tab_idx', '_topde_tab_idx'):
            idx = getattr(self, attr, None)
            if idx is not None and self.tabs is not None:
                self.tabs.setTabVisible(idx, has_pathways)

    def _on_mode_changed(self, mode):
        """Reset the pathway-only volcano filter to the mode's natural default.

        Discovery shows the whole transcriptome, hypothesis restricts to the
        committed panel. Without this, switching hypothesis -> discovery would
        leave the filter on and keep the volcano restricted to the old panel
        genes even after a genome-wide rerun.
        """
        self._refresh_pathway_ui_visibility()
        if not hasattr(self, '_pathway_only_check'):
            return
        want = (mode == 'scoring')
        if self._pathway_only_check.isChecked() != want:
            self._pathway_only_check.blockSignals(True)
            self._pathway_only_check.setChecked(want)
            self._pathway_only_check.blockSignals(False)
            if self.ws.de_results is not None:
                self._on_pval_toggle()

    # --- Activation ---

    def on_activated(self):
        """Called when page becomes visible. Auto-loads previous results if available."""
        has_data = self.ws.current_adata is not None
        self.run_btn.setEnabled(has_data)
        self.batch_btn.setEnabled(has_data)
        self._ge_run_btn.setEnabled(has_data)
        self._refresh_study_control()
        self._refresh_dataset_card()
        self._refresh_filter_summary()
        self._refresh_card_summaries()
        self._refresh_pathway_ui_visibility()
        self._populate_covariates()

        self.unit_combo.blockSignals(True)
        self.unit_combo.setCurrentIndex(1 if self.ws.unit == 'cell' else 0)
        self.unit_combo.blockSignals(False)

        self.de_method_combo.blockSignals(True)
        self.de_method_combo.setCurrentIndex(1 if self.ws.de_method == 'deseq2' else 0)
        self.de_method_combo.blockSignals(False)
        self.eb_checkbox.setEnabled(self.ws.de_method != 'deseq2')
        is_deseq2 = self.ws.de_method == 'deseq2'
        self.indep_filter_check.setEnabled(is_deseq2)
        self.cooks_filter_check.setEnabled(is_deseq2)
        self._apply_unit_ui_state()

        self.min_cells_spin.blockSignals(True)
        self.min_cells_spin.setValue(self.ws.min_cells)
        self.min_cells_spin.blockSignals(False)
        self.min_counts_spin.blockSignals(True)
        self.min_counts_spin.setValue(self.ws.min_counts)
        self.min_counts_spin.blockSignals(False)

        self._refresh_count_source()

        if self.ws.pathway_gene_sets and not self._pw_gene_combo.count():
            self._populate_pw_gene_combo()

        # Populating the 26k-gene table + volcano is fast but not free; defer
        # it one event-loop tick so the tab switch paints immediately. Surface
        # the load clearly (sidebar label + status bar) while it happens.
        if has_data and self.results_table.rowCount() == 0 and self.ws.de_results is not None:
            self.status_label.setText("Loading gene DE results...")
            self.ws.status_message.emit("Loading gene DE results...")
        QTimer.singleShot(0, self._populate_on_activate)
        # Catch setup-page edits (groups, columns) made since the last run.
        self._refresh_run_dirty()
        # Catch upstream data changes (re-QC, roles, decontX) since this DE ran.
        self._check_upstream_stale()

    def _populate_on_activate(self):
        """Deferred results populate, so switching to the page is instant."""
        try:
            if self.ws.de_results is not None:
                if self.results_table.rowCount() == 0:
                    self._populate_results()
                if not self.heatmap_pathway_combo.count() and self.ws.pathway_gene_sets:
                    self._populate_heatmap_combo()
                if not self._pw_gene_combo.count() and self.ws.pathway_gene_sets:
                    self._populate_pw_gene_combo()
                self._update_topde()
                self.status_label.setText("")
                return
            self._try_load_previous_results()
        finally:
            self.ws.status_message.emit("Ready")

    def _current_comparison(self):
        """Return a dict of the DE settings that require a re-run if changed."""
        full_genome = True  # per-cell normalisation needs the full matrix
        return {
            'control': self.ws.control_label,
            'disease': self.ws.disease_label,
            'condition_col': self.ws.condition_col,
            'sample_col': self.ws.sample_col,
            'min_cells': self.ws.min_cells,
            'min_counts': self.ws.min_counts,
            'de_method': self.ws.de_method,
            'unit': self.ws.unit,
            'moderate': self.eb_checkbox.isChecked(),
            'detection_min_pct': round(self.detection_pct_spin.value(), 4),
            'detection_min_donor_frac': round(
                self.detection_donor_frac_spin.value(), 4),
            'gene_filter_per_study': bool(self._per_study_filter_col()),
            'counts_layer': self._selected_counts_layer() or 'raw',
            'indep_filter': self.indep_filter_check.isChecked(),
            'cooks_filter': self.cooks_filter_check.isChecked(),
            'full_genome': full_genome,
            'geneset_label': self.ws.geneset_label if not full_genome else '',
        }

    def _refresh_run_dirty(self, *_args):
        """Flag when a DE-affecting setting changed since the last run, so the
        results on screen are stale and the Run button prompts a re-run."""
        if not hasattr(self, 'run_btn'):
            return
        last = getattr(self, '_last_run_comparison', None)
        if self.ws.de_results is None or last is None:
            self.run_btn.setText("Run DE Analysis")
            return
        if self._current_comparison() != last:
            self.run_btn.setText("Re-run DE Analysis")
            self.status_label.setText(
                "Settings changed since the last run. Re-run to update the "
                "results shown.")
        else:
            self.run_btn.setText("Run DE Analysis")
            if self.status_label.text().startswith("Settings changed"):
                self.status_label.setText("")

    def _study_dir(self):
        """Directory holding the current study's h5ad + provenance sidecar."""
        h5ad = getattr(self.ws, 'h5ad_path', None)
        return Path(h5ad).parent if h5ad else None

    def _record_de_provenance(self, stage):
        """Append a DE stage to the study's provenance, tagged with the upstream
        data token it was computed against. Best-effort; never blocks."""
        study_dir = self._study_dir()
        if study_dir is None or self.ws.current_adata is None:
            return
        try:
            from kosmic import provenance
            upstream = provenance.compute_fingerprint(
                self.ws.current_adata,
                h5ad_path=getattr(self.ws, 'h5ad_path', None))['token']
            params = {
                'method': self.ws.de_method,
                'unit': self.ws.unit,
                'min_cells': self.ws.min_cells,
                'min_counts': self.ws.min_counts,
                'detection_min_pct': round(self.detection_pct_spin.value() / 100.0, 4),
                'detection_min_donor_frac': round(
                    self.detection_donor_frac_spin.value() / 100.0, 4),
                'gene_filter_per_study': bool(self._per_study_filter_col()),
                'filter_min_count': self.filter_min_count_spin.value(),
                'filter_min_samples': int(self.filter_min_samples_spin.value()),
                'count_source': self._selected_counts_layer() or 'raw',
                'independent_filter': self.indep_filter_check.isChecked(),
                'cooks_filter': self.cooks_filter_check.isChecked(),
                'fdr_threshold': self.pval_filter.value(),
                'log2fc_threshold': self.fc_filter.value(),
                'mode': getattr(self.ws, 'analysis_mode', None),
            }
            provenance.record_stage(
                study_dir, Path(self.ws.h5ad_path).stem, stage, params,
                consumed_token=upstream)  # DE does not change the data
        except Exception as e:
            self.ws.log_message.emit(f"Could not record DE provenance: {e}")

    def _check_upstream_stale(self):
        """Flag when the data changed since the loaded DE ran (re-QC, roles,
        decontX upstream) -- separate from the in-session settings dirty check."""
        study_dir = self._study_dir()
        if (study_dir is None or self.ws.de_results is None
                or self.ws.current_adata is None):
            return
        try:
            from kosmic import provenance
            de_stage = provenance.last_stage(study_dir, 'gene_de')
            consumed = de_stage.get('consumed_token') if de_stage else None
            if not consumed:
                return
            live = provenance.compute_fingerprint(
                self.ws.current_adata,
                h5ad_path=getattr(self.ws, 'h5ad_path', None))['token']
            if consumed != live:
                self.run_btn.setText("Re-run DE Analysis")
                self.status_label.setText(
                    "Upstream data changed since this DE ran. Re-run to update.")
        except Exception:
            pass

    def _try_load_previous_results(self):
        """Load previous DE results -- prefer adata.uns, fall back to CSV on disk."""
        import pandas as pd

        method = self.ws.de_method
        moderate = self.eb_checkbox.isChecked()
        if method == 'ttest' and moderate:
            label = 'welch_cpm_eb'
        elif method == 'ttest':
            label = 'welch_cpm'
        else:
            label = 'deseq2'

        adata = self.ws.current_adata
        if adata is not None and 'de_results' in adata.uns:
            uns_results = adata.uns['de_results']
            if isinstance(uns_results, dict) and label in uns_results:
                de_df = uns_results[label]
                if isinstance(de_df, pd.DataFrame) and not de_df.empty and 'names' in de_df.columns:
                    self._apply_loaded_results(de_df, f"adata.uns['{label}']")
                    return

        if not self.ws.project_dir:
            return

        stats_dir = de_stats_dir(self.ws.project_dir)
        if not stats_dir.is_dir():
            return

        gse_id = self.ws.gse_accession or (
            Path(self.ws.h5ad_path).stem if self.ws.h5ad_path else None
        )
        if not gse_id:
            return

        path = de_result_path(self.ws.project_dir, gse_id, label)
        if not path.is_file():
            return
        try:
            de_results = pd.read_csv(path)
            if de_results.empty or 'names' not in de_results.columns:
                return
            self._apply_loaded_results(de_results, path.name)
        except Exception as e:
            self.ws.log_message.emit(f"Could not load previous DE results: {e}")

    def _apply_loaded_results(self, de_results, source_name):
        """Apply loaded DE results to the workspace and refresh the UI."""
        import pandas as pd

        self.ws.de_results = de_results
        # Disk results carry no preserved genome-wide ranking; clear any stale
        # one so fgsea falls back to these loaded (genome-wide) results.
        self.ws.gene_de_genome_wide = de_results.attrs.get('genome_wide_ranking')

        if 'pvals_adj' in de_results.columns and 'abs_logfoldchange' in de_results.columns:
            sig_mask = (de_results['pvals_adj'] < DEFAULT_FDR) & (de_results['abs_logfoldchange'] > DEFAULT_LFC_THRESHOLD)
            self.ws.significant_genes = de_results[sig_mask]
        else:
            self.ws.significant_genes = pd.DataFrame()

        if self.ws.pathway_gene_sets and self.ws.current_adata is not None:
            from kosmic.de.de_analysis import prepare_gene_coverage
            from kosmic.scrna.inspect.detection import detect_species
            species = detect_species(list(self.ws.current_adata.var_names))
            var_names = set(self.ws.current_adata.raw.var_names) if self.ws.current_adata.raw is not None else set(self.ws.current_adata.var_names)
            _, self.ws.pathway_coverage = prepare_gene_coverage(
                self.ws.pathway_gene_sets, var_names, species
            )

        if self.ws.sample_df is None and self.ws.current_adata is not None:
            adata = self.ws.current_adata
            sample_col = self.ws.sample_col
            condition_col = self.ws.condition_col
            if sample_col and condition_col and sample_col in adata.obs.columns:
                obs = adata.obs
                sample_info = obs.groupby(sample_col).agg(
                    condition=(condition_col, 'first'),
                    n_cells=(condition_col, 'size'),
                ).reset_index()
                sample_info.columns = ['sample', 'condition', 'n_cells']
                self.ws.sample_df = sample_info

        self._populate_results()

        # Default the pathway-only filter to the mode (off for discovery so the
        # genome-wide volcano shows, on for hypothesis) before rendering.
        if hasattr(self, '_pathway_only_check'):
            self._pathway_only_check.blockSignals(True)
            self._pathway_only_check.setChecked(self.ws.analysis_mode == 'scoring')
            self._pathway_only_check.blockSignals(False)
        self._refresh_pathway_ui_visibility()

        # Render the volcano honouring the pathway-only / pct checkbox states.
        self._on_pval_toggle()

        self._populate_heatmap_combo()
        self._populate_pw_gene_combo()
        self._update_topde()

        if self.ws.analysis_mode == 'scoring':
            self.ws.mark_step_complete(3)  # scoring step 3 = Gene DE
        elif self.ws.analysis_mode == 'discovery':
            self.ws.mark_step_complete(2)  # discovery step 2 = Gene DE

        # Loaded results are treated as matching the current settings, so the
        # Run button is not flagged stale until something changes.
        self._last_run_comparison = self._current_comparison()
        self._refresh_run_dirty()

        n_sig = len(self.ws.significant_genes) if self.ws.significant_genes is not None else 0
        self.ws.log_message.emit(
            f"Loaded previous DE results: {source_name} "
            f"({len(de_results)} genes, {n_sig} significant)"
        )
        self.ws.status_message.emit(
            f"Loaded previous results -- {len(de_results)} genes ({n_sig} significant)"
        )
        self.status_label.setText(f"Loaded: {source_name}")

    # --- Settings handlers ---

    def _on_de_method_changed(self, index):
        self.ws.de_method = 'deseq2' if index == 1 else 'ttest'
        is_ttest = index == 0
        # In cell mode the method picker is inert (always Wilcoxon).
        self.eb_checkbox.setEnabled(is_ttest and self.ws.unit != 'cell')
        if not is_ttest:
            self.eb_checkbox.setChecked(False)
        # DESeq2 filter toggles apply only to the DESeq2 engine.
        self.indep_filter_check.setEnabled(not is_ttest)
        self.cooks_filter_check.setEnabled(not is_ttest)
        self._reload_for_current_method()

    def _on_unit_changed(self, index):
        self.ws.unit = 'cell' if index == 1 else 'sample'
        self._apply_unit_ui_state()
        self._reload_for_current_method()

    def _apply_unit_ui_state(self):
        """Grey out the pseudobulk-only controls when cell mode is active."""
        is_cell = self.ws.unit == 'cell'
        self._unit_caveat.setVisible(is_cell)
        self.de_method_combo.setEnabled(not is_cell)
        self.eb_checkbox.setEnabled(
            (not is_cell) and self.ws.de_method == 'ttest')

    def _on_eb_toggled(self, checked):
        # EB only applies to t-test (toggles welch_cpm vs welch_cpm_eb storage key).
        if self.ws.de_method == 'ttest':
            self._reload_for_current_method()

    def _reload_for_current_method(self):
        """Clear current results and try to load stored results for the active method."""
        self.ws.de_results = None
        self.ws.significant_genes = None
        self._try_load_previous_results()
        self._refresh_run_dirty()

    def _on_min_cells_changed(self, value):
        self.ws.min_cells = value
        self._refresh_run_dirty()

    def _on_min_counts_changed(self, value):
        self.ws.min_counts = value
        self._refresh_run_dirty()

    # --- Run DE ---

    #: Columns never worth offering as a covariate -- either the thing
    #: being tested, the sample identifier itself, or per-cell metrics.
    _COVARIATE_SKIP = {
        'condition', '_role', 'sample', 'cell_type', 'cell_type_atlas',
        'cell_type_fine', 'cell_type_auto', 'leiden', 'leiden_atlas',
    }

    def _populate_covariates(self):
        """Offer sample-level obs columns as DESeq2 design terms.

        A covariate has to be constant within a sample to mean anything
        at the pseudobulk level, and has to have at least two levels but
        fewer than one per sample to be estimable. Filtering here keeps
        the list short and honest rather than offering 40 columns of
        which two are usable.
        """
        adata = self.ws.current_adata
        if adata is None or not hasattr(self, 'covariate_list'):
            return
        previously = set(self._selected_covariates())
        self.covariate_list.clear()
        if 'sample' not in adata.obs.columns:
            return

        n_samples = adata.obs['sample'].nunique()
        grouped = adata.obs.groupby('sample', observed=True)
        candidates = []
        for col in adata.obs.columns:
            if col in self._COVARIATE_SKIP or col.startswith('_'):
                continue
            series = adata.obs[col]
            if not (str(series.dtype) in ('category', 'object')
                    or hasattr(series, 'cat')):
                continue
            n_levels = series.nunique(dropna=True)
            if n_levels < 2 or n_levels >= n_samples:
                continue
            # Must not vary within a sample.
            if int(grouped[col].nunique().max()) != 1:
                continue
            candidates.append((col, n_levels))

        # Drop anything the design could not fit anyway. 'group', 'region'
        # and 'state' are often the condition under another name; offering
        # them means a user ticks a term that is silently dropped at run
        # time, and reads a design that was never fitted.
        candidates = self._estimable_only(adata, candidates)
        self._n_estimable_covariates = len(candidates)

        for col, n_levels in candidates:
            item = QListWidgetItem(f"{col}  ({n_levels} levels)")
            item.setData(Qt.ItemDataRole.UserRole, col)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            # 'study' is why this control exists, so tick it by default
            # when the dataset actually has more than one.
            default_on = (col == 'study') or col in previously
            item.setCheckState(Qt.CheckState.Checked if default_on
                               else Qt.CheckState.Unchecked)
            self.covariate_list.addItem(item)

        if not self.covariate_list.count():
            # An empty box reads as "broken"; say why it is empty.
            placeholder = QListWidgetItem(
                "Nothing to adjust for in this dataset")
            placeholder.setFlags(Qt.ItemFlag.NoItemFlags)
            self.covariate_list.addItem(placeholder)

        # Nothing estimable means nothing to open the card for.
        if hasattr(self, '_design_card'):
            self._design_card.set_expandable(
                self._n_estimable_covariates > 0)

    def _detection_study_col(self):
        """The column to stratify detection by, or None.

        Only meaningful on a combined object. A single-study dataset has
        nothing to stratify, so the control is disabled rather than
        silently doing nothing.
        """
        adata = self.ws.current_adata
        if adata is None or 'study' not in adata.obs.columns:
            return None
        if adata.obs['study'].astype(str).nunique() < 2:
            return None
        return 'study'

    def _per_study_filter_col(self):
        """The study column to filter within, honouring the user's choice.

        None when the dataset is a single study, or when the choice is
        off -- in which case the gene filter runs once on the pooled
        table, which is what a standalone mega-analysis conventionally
        does.
        """
        if not self.detection_study_check.isChecked():
            return None
        return self._detection_study_col()

    def _refresh_study_control(self):
        """Enable the per-study option only when there is more than one."""
        if not hasattr(self, 'detection_study_check'):
            return
        col = self._detection_study_col()
        self.detection_study_check.setEnabled(col is not None)
        self.detection_study_check.setVisible(col is not None)
        if col is not None:
            n = self.ws.current_adata.obs[col].astype(str).nunique()
            self.detection_study_check.setText(
                f"Filter genes within each study  ({n} studies)")

    def _auto_min_samples(self):
        """The donor count 'auto' resolves to, or None if not yet knowable."""
        adata = self.ws.current_adata
        if adata is None or '_role' not in adata.obs.columns:
            return None, None
        if 'sample' not in adata.obs.columns:
            return None, None
        per_donor = adata.obs.groupby('sample', observed=True)['_role'].first()
        per_donor = per_donor.astype(str).str.lower()
        n_d = int((per_donor == 'disease').sum())
        n_c = int((per_donor == 'control').sum())
        if not n_d or not n_c:
            return None, None
        return required_donor_count(
            n_d, n_c,
            min_samples=int(self.filter_min_samples_spin.value())), n_d + n_c

    def _cohort_shape(self):
        """Cached cohort description for the loaded dataset."""
        from kosmic.scrna.inspect.batch import cohort_shape

        adata = self.ws.current_adata
        if adata is None:
            return None
        token = (id(adata), adata.n_obs, self.ws.sample_col)
        if getattr(self, '_cohort_token', None) != token:
            self._cohort_token = token
            self._cohort_cache = cohort_shape(
                adata, sample_col=self.ws.sample_col)
        return self._cohort_cache

    def _refresh_dataset_card(self):
        """Say plainly whether this is one study or several pooled."""
        from kosmic.scrna.inspect.batch import cohort_headline

        if not hasattr(self, '_dataset_card'):
            return
        shape = self._cohort_shape()
        if shape is None:
            self._dataset_card.set_summary(
                ["No data loaded", "Load a dataset from the Project page"])
            return
        self._dataset_card.set_summary(list(cohort_headline(shape)))

    def _refresh_card_summaries(self, _=None):
        """Two lines per card saying what a run would actually do.

        The point of the collapsed card is that you should not have to
        open it to know the answer, so these mirror the controls inside
        rather than restating their labels.
        """
        if not hasattr(self, '_model_card'):
            return
        self._refresh_dataset_card()

        unit = ("per donor (pseudobulk)" if self.ws.unit == 'sample'
                else "per cell (pseudoreplicated)")
        engine = ("Wilcoxon" if self.ws.unit == 'cell'
                  else self.de_method_combo.currentText())
        counts = self.count_source_combo.currentText()
        model_lines = [engine, f"{unit} \u00b7 {counts}"]
        if self.ws.unit == 'sample' and self.eb_checkbox.isChecked():
            model_lines[1] += " \u00b7 EB"
        self._model_card.set_summary(model_lines)

        covariates = self._selected_covariates()
        design = "~" + " + ".join([*covariates, "condition"])
        shape = self._cohort_shape()
        if shape and shape['is_atlas'] and 'study' not in covariates:
            # The one combination worth shouting about: several cohorts in
            # one model with no term for which cohort a donor came from.
            second = f"{shape['n_studies']} studies NOT adjusted for"
        elif shape and shape['is_atlas']:
            second = f"adjusting for {shape['n_studies']} studies"
        elif not shape:
            second = "no data loaded"
        elif not getattr(self, '_n_estimable_covariates', 0):
            # Say what is true rather than blaming the study count: a
            # single study can still have a usable covariate (sex, run),
            # and an atlas can have none beyond 'study'.
            second = "nothing estimable to adjust for"
        elif covariates:
            second = f"adjusting for {', '.join(covariates)}"
        else:
            second = (f"{self._n_estimable_covariates} available, "
                      f"none used")
        self._design_card.set_summary([design, second])

        auto_n, total = self._auto_min_samples()
        min_samples = int(self.filter_min_samples_spin.value())
        if min_samples > 0:
            donors = f"\u2265{min_samples} donors"
        elif auto_n is not None:
            donors = f"\u2265{auto_n:g} of {total} donors"
        else:
            donors = "the smaller group"
        # Two lines, one per filter, each naming its own subject.
        # Quoting the count threshold here put a "10 counts" directly
        # under a "10 cells" -- two unrelated tens, read as one setting.
        per_study = (self._per_study_requirements()
                     if self._per_study_filter_col() else [])
        if per_study:
            spread = " + ".join(f"{req:g}/{n}" for _name, req, n in per_study)
            second = f"Genes: carried by {spread} donors, per study"
        else:
            second = f"Genes: carried by {donors}"
        # Spell out what 'auto' resolves to, inside the box itself.
        auto_txt = "auto" if auto_n is None else f"auto ({auto_n:g})"
        if self.filter_min_samples_spin.specialValueText() != auto_txt:
            self.filter_min_samples_spin.setSpecialValueText(auto_txt)
        if self.detection_pct_spin.value() > 0:
            second += f", {self.detection_pct_spin.value():g}%+ detection"
        donors = f"Donors: \u2265{self.min_cells_spin.value()} cells each"
        if self.min_counts_spin.value():
            donors += f", \u2265{self.min_counts_spin.value():,} transcripts"
        self._filters_card.set_summary([donors, second])

        if self.ws.de_method == 'deseq2' and self.ws.unit == 'sample':
            indep = ("independent filtering on"
                     if self.indep_filter_check.isChecked()
                     else "independent filtering off")
            cooks = ("Cook's outliers blanked"
                     if self.cooks_filter_check.isChecked()
                     else "Cook's filtering off")
            self._fdr_card.set_summary([indep, cooks])
        else:
            self._fdr_card.set_summary(
                ["DESeq2 only", "not used by this engine"])

    def _show_advanced_filters(self):
        """The gene filter's numbers, plus the optional cell-level floor.

        One QFormLayout for both groups so there is a single label
        column -- two separate forms sized their columns independently
        and the rows did not line up. Headers use 'subsection_title'
        (13px bold, primary) rather than the sidebar's dim 11px
        'section_header', which looked washed out beside the captions.
        """
        from PyQt6.QtWidgets import (
            QDialog, QDialogButtonBox, QFormLayout, QVBoxLayout,
        )

        def _heading(text):
            label = QLabel(text)
            label.setProperty("role", "subsection_title")
            return label

        dlg = QDialog(self)
        dlg.setWindowTitle("Gene filter settings")
        dlg.setFixedWidth(400)
        layout = QVBoxLayout(dlg)

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight
                               | Qt.AlignmentFlag.AlignVCenter)
        form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        form.setHorizontalSpacing(10)
        form.setVerticalSpacing(8)

        def _rows(heading, rows):
            form.addRow(_heading(heading))
            for caption, widget in rows:
                if widget is None:
                    continue
                widget.setMaximumWidth(110)
                form.addRow(caption + ":", widget)

        _rows("The gene filter", self._gene_filter_rows)

        live = SecondaryLabel("")
        live.setWordWrap(True)
        form.addRow("", live)

        _rows("Optional: cell-level floor  (not standard)",
              self._detection_rows)
        layout.addLayout(form)

        def _live(_=None):
            live.setText(self._filter_sentence())

        spins = (self.filter_min_count_spin, self.filter_min_samples_spin,
                 self.detection_pct_spin, self.detection_donor_frac_spin)
        for spin in spins:
            spin.valueChanged.connect(_live)
        _live()

        layout.addSpacing(8)
        reset = QPushButton("Reset to standard values")
        reset.clicked.connect(self._reset_advanced_filters)
        reset.clicked.connect(_live)
        layout.addWidget(reset)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(dlg.reject)
        buttons.accepted.connect(dlg.accept)
        layout.addWidget(buttons)
        dlg.exec()

        for spin in spins:
            spin.valueChanged.disconnect(_live)
            spin.setParent(None)
        self._refresh_card_summaries()

    def _reset_advanced_filters(self):
        """Back to the documented defaults."""
        self.detection_pct_spin.setValue(
            DE_DETECTION_MIN_PCT * 100 if DE_DETECTION_ON else 0.0)
        self.detection_donor_frac_spin.setValue(
            DE_DETECTION_MIN_DONOR_FRAC * 100)
        self.filter_min_count_spin.setValue(float(DE_FILTER_MIN_COUNT))
        self.filter_min_samples_spin.setValue(float(DE_FILTER_MIN_SAMPLES))

    def _per_study_requirements(self):
        """[(study, required, n_donors)] when the filter runs per study.

        The pooled figure is wrong once stratification is on: a gene has
        to clear each study's own requirement, computed from that study's
        own arms, not from the combined cohort.
        """
        shape = self._cohort_shape()
        if not shape or not shape.get('is_atlas'):
            return []
        rows = []
        for study in shape['studies']:
            required = required_donor_count(
                study['n_disease'], study['n_control'],
                min_samples=int(self.filter_min_samples_spin.value()))
            if required is None:
                continue
            rows.append((study['study'], required, study['n_donors']))
        return rows

    def _filter_sentence(self):
        """The gene rule as a short phrase, for the dialog's live readout.

        Deliberately terse. Every spin box carries a tooltip with the
        mechanics; a settings dialog that explains itself in paragraphs
        is unreadable at the moment you actually want to change a value.
        """
        min_count = self.filter_min_count_spin.value()
        min_samples = int(self.filter_min_samples_spin.value())

        per_study = (self._per_study_requirements()
                     if self._per_study_filter_col() else [])
        if per_study:
            spread = " and ".join(f"{req:g} of {n} ({name})"
                                  for name, req, n in per_study)
            text = f"\u2265{min_count:g} counts in {spread} donors"
        else:
            auto_n, total = self._auto_min_samples()
            if min_samples > 0:
                text = f"\u2265{min_count:g} counts in \u2265{min_samples} donors"
            elif auto_n is not None:
                text = (f"\u2265{min_count:g} counts in \u2265{auto_n:g} of "
                        f"the {total} donors (the smaller arm)")
            else:
                text = f"\u2265{min_count:g} counts in the smaller arm's donors"

        pct = self.detection_pct_spin.value()
        if pct > 0:
            text += f", detected in \u2265{pct:g}% of cells"
        return text + "."

    def _refresh_filter_summary(self, _=None):
        """Kept for the signal wiring; the card no longer shows a sentence."""
        self._refresh_card_summaries()

    def _estimable_only(self, adata, candidates):
        """Keep only covariates the DESeq2 design could actually fit.

        Delegates to the pipeline's own rule so the list offered here and
        the list accepted at run time cannot drift apart.
        """
        from kosmic.de.de_analysis import usable_covariates

        if not candidates or '_role' not in adata.obs.columns:
            return candidates
        cols = [c for c, _n in candidates]
        per_donor = adata.obs.groupby('sample', observed=True).first()
        roles = per_donor['_role'].astype(str).str.lower().values
        keep = set(usable_covariates(per_donor, cols, roles))
        return [(c, n) for c, n in candidates if c in keep]

    def _selected_covariates(self):
        if not hasattr(self, 'covariate_list'):
            return []
        return [self.covariate_list.item(i).data(Qt.ItemDataRole.UserRole)
                for i in range(self.covariate_list.count())
                if self.covariate_list.item(i).checkState()
                == Qt.CheckState.Checked]

    def _run_analysis(self):
        if self.ws.current_adata is None:
            dialogs.warning(self, "Error", "No data loaded.")
            return
        if not self._confirm_mixed_cell_types():
            return

        sample_col = self.ws.sample_col
        condition_col = self.ws.condition_col
        control_label = self.ws.control_label
        disease_label = self.ws.disease_label

        if not all([sample_col, condition_col, control_label, disease_label]):
            dialogs.warning(self, "Error", "Please configure all analysis parameters in Setup.")
            return

        if control_label == disease_label:
            dialogs.warning(self, "Error", "Control and disease groups must be different.")
            return

        min_cells = self.ws.min_cells
        de_method = self.ws.de_method

        full_genome = True  # per-cell normalisation needs the full matrix

        fdr_genes = self._committed_fdr_genes()

        self.de_worker = MetabolicDEWorker(
            self.ws.current_adata, sample_col, condition_col,
            control_label, disease_label,
            self.ws.pathway_gene_sets, min_cells,
            min_counts=self.ws.min_counts,
            de_method=de_method, full_genome=full_genome,
            moderate=self.eb_checkbox.isChecked(),
            unit=self.ws.unit,
            counts_layer=self._selected_counts_layer(),
            detection_min_pct=self.detection_pct_spin.value() / 100.0,
            detection_min_donor_frac=(
                self.detection_donor_frac_spin.value() / 100.0),
            detection_study_col=self._per_study_filter_col(),
            fdr_genes=fdr_genes,
            deseq2_independent_filter=self.indep_filter_check.isChecked(),
            deseq2_cooks_filter=self.cooks_filter_check.isChecked(),
            filter_min_count=self.filter_min_count_spin.value(),
            filter_min_samples=int(self.filter_min_samples_spin.value()),
            covariates=self._selected_covariates(),
        )
        run_worker(
            self.de_worker,
            on_finished=self._on_analysis_finished,
            on_failed=self._on_analysis_failed,
            on_progress=self._on_progress,
            on_progress_pct=self._on_progress_pct,
        )
        self.run_btn.setEnabled(False)
        engine = ("Wilcoxon" if self.ws.unit == 'cell'
                  else self.de_method_combo.currentText())
        self.ws.status_message.emit(f"Running {engine} DE...")
        if self.ws.progress_bar:
            self.ws.progress_bar.setValue(0)

    # --- Run DE once per cell type ---

    #: obs columns most likely to hold cell-type labels, best first. Same
    #: ranking the Subset step uses, so the two agree on what to offer.
    _CELL_TYPE_KEYWORDS = [
        'cell_type', 'celltype', 'cell type',
        'cell_annotation', 'annotation', 'cluster', 'leiden', 'louvain',
    ]

    def _cell_type_columns(self, adata):
        """Categorical obs columns, cell-type-looking ones first."""
        ranked, other = [], []
        for col in adata.obs.columns:
            if col.startswith('_'):
                continue
            series = adata.obs[col]
            is_cat = (str(series.dtype) in ('category', 'object')
                      or hasattr(series, 'cat'))
            if any(kw in col.lower() for kw in self._CELL_TYPE_KEYWORDS):
                ranked.append(col)
            elif is_cat:
                other.append(col)

        def _rank(col):
            low = col.lower()
            for i, kw in enumerate(self._CELL_TYPE_KEYWORDS):
                if low == kw:
                    return i
            for i, kw in enumerate(self._CELL_TYPE_KEYWORDS):
                if kw in low:
                    return i + len(self._CELL_TYPE_KEYWORDS)
            return len(self._CELL_TYPE_KEYWORDS) * 2

        ranked.sort(key=_rank)
        return ranked, other

    def _cell_type_multiplicity(self):
        """(n_cell_types, column) for the loaded object, or (None, None).

        The smallest count across the candidate columns, because a
        subset dataset often keeps a stale finer-grained column: an
        endothelial subset can carry 'cell_type' with 6 leftover labels
        and 'cell_type_atlas' with 1. One column saying "this is a
        single cell type" is enough.
        """
        adata = self.ws.current_adata
        if adata is None:
            return None, None
        ranked, _other = self._cell_type_columns(adata)
        best = (None, None)
        for col in ranked:
            n = adata.obs[col].astype(str).nunique()
            if best[0] is None or n < best[0]:
                best = (n, col)
        return best

    def _confirm_mixed_cell_types(self):
        """Warn before pooling several cell types into one contrast.

        Not a hard block: whole-tissue pseudobulk is a real analysis, and
        the check cannot know it was not intended. But it is far more
        often a mistake -- the donor profile becomes a mix, and a shift
        in cell composition then reads as a change in expression.
        """
        n_types, col = self._cell_type_multiplicity()
        if not n_types or n_types < 2:
            return True

        from PyQt6.QtWidgets import QMessageBox

        box = QMessageBox(self)
        box.setIcon(QMessageBox.Icon.Warning)
        box.setWindowTitle("Several cell types loaded")
        box.setText(
            f"This dataset holds {n_types} cell types (column '{col}').")
        box.setInformativeText(
            "Running one analysis over all of them sums every cell of a "
            "donor into a single profile, so cardiomyocytes and endothelium "
            "are averaged together. A shift in cell composition between "
            "your groups then looks like a change in expression.\n\n"
            "Run per cell type unless you specifically want a whole-tissue "
            "result.")
        per_type = box.addButton("Run per cell type...",
                                 QMessageBox.ButtonRole.AcceptRole)
        anyway = box.addButton("Run over all cell types",
                               QMessageBox.ButtonRole.DestructiveRole)
        box.addButton(QMessageBox.StandardButton.Cancel)
        box.setDefaultButton(per_type)
        box.exec()

        if box.clickedButton() is per_type:
            self._run_batch_by_cell_type()
            return False
        return box.clickedButton() is anyway

    def _run_batch_by_cell_type(self):
        """Pick a cell-type column and a set of types, then run them all."""
        from PyQt6.QtWidgets import (
            QDialog, QDialogButtonBox, QVBoxLayout, QHBoxLayout,
        )

        from kosmic.de.batch import plan_cell_types

        adata = self.ws.current_adata
        if adata is None:
            dialogs.warning(self, "Error", "No data loaded.")
            return
        if not all([self.ws.sample_col, self.ws.condition_col]):
            dialogs.warning(
                self, "Error",
                "Configure the sample and condition columns in Setup first.")
            return
        if '_role' not in adata.obs.columns:
            dialogs.warning(
                self, "Roles not set",
                "Set Disease / Control / Exclude in the Inspect tab "
                "(Sample Setup) and save, then come back.")
            return

        n_types, col = self._cell_type_multiplicity()
        if n_types == 1:
            dialogs.info(
                self, "Only one cell type",
                f"Every cell here has the same '{col}' label, so a "
                f"per-cell-type run would produce exactly one result -- "
                f"the same contrast Run DE Analysis gives you.\n\n"
                f"Use Run DE Analysis instead, or load a dataset that "
                f"still has several cell types.")
            return

        ranked, other = self._cell_type_columns(adata)
        if not ranked and not other:
            dialogs.warning(
                self, "No cell-type column",
                "No categorical obs column to split on. Annotate the "
                "dataset first, or propagate labels from an atlas.")
            return

        dlg = QDialog(self)
        dlg.setWindowTitle("Run DE per cell type")
        dlg.setMinimumWidth(560)
        layout = QVBoxLayout(dlg)

        layout.addWidget(SecondaryLabel(
            "Each cell type is analysed on its own, using the settings on "
            "this page. Results are written per type and appear in the "
            "Meta-Analysis workspace as separate entries, so pooling "
            "happens across studies within a cell type."))

        col_row = QHBoxLayout()
        col_row.addWidget(QLabel("Cell type column:"))
        col_combo = NoScrollComboBox()
        for col in ranked:
            col_combo.addItem(f"★ {col}", col)
        for col in other:
            col_combo.addItem(col, col)
        col_row.addWidget(col_combo, 1)
        layout.addLayout(col_row)

        type_list = QListWidget()
        layout.addWidget(type_list, 1)

        summary = SecondaryLabel("")
        summary.setWordWrap(True)
        layout.addWidget(summary)

        btn_row = QHBoxLayout()
        all_btn = QPushButton("Select all eligible")
        none_btn = QPushButton("Clear")
        btn_row.addWidget(all_btn)
        btn_row.addWidget(none_btn)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Run")
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        layout.addWidget(buttons)

        min_cells = self.min_cells_spin.value()
        min_counts = self.min_counts_spin.value()
        counts_layer = self._selected_counts_layer()

        def _checked_types():
            return [type_list.item(i).data(Qt.ItemDataRole.UserRole)
                    for i in range(type_list.count())
                    if type_list.item(i).flags() != Qt.ItemFlag.NoItemFlags
                    and type_list.item(i).checkState()
                    == Qt.CheckState.Checked]

        def _update_summary():
            n = len(_checked_types())
            covariates = self._selected_covariates()
            design = "~" + " + ".join([*covariates, "condition"])
            summary.setText(
                f"{n} cell type(s) selected · design {design} · "
                f"method {self.de_method_combo.currentText()}")

        def _repopulate():
            type_list.clear()
            col = col_combo.currentData()
            try:
                plans = plan_cell_types(
                    adata, col, self.ws.sample_col, min_cells=min_cells,
                    min_counts=min_counts, counts_layer=counts_layer)
            except ValueError as exc:
                summary.setText(str(exc))
                return
            for plan in plans:
                text = (f"{plan.cell_type}  --  {plan.n_cells:,} cells, "
                        f"{plan.n_disease_samples} disease / "
                        f"{plan.n_control_samples} control donors")
                if not plan.eligible:
                    text += f"  [{plan.reason}]"
                item = QListWidgetItem(text)
                item.setData(Qt.ItemDataRole.UserRole, plan.cell_type)
                if plan.eligible:
                    item.setFlags(item.flags()
                                  | Qt.ItemFlag.ItemIsUserCheckable)
                    item.setCheckState(Qt.CheckState.Checked)
                else:
                    # Ineligible types are shown, not hidden: "why is
                    # pericyte missing" is a question worth answering in
                    # place rather than by silence.
                    item.setFlags(Qt.ItemFlag.NoItemFlags)
                type_list.addItem(item)
            _update_summary()

        def _set_all(checked):
            state = (Qt.CheckState.Checked if checked
                     else Qt.CheckState.Unchecked)
            for i in range(type_list.count()):
                item = type_list.item(i)
                if item.flags() != Qt.ItemFlag.NoItemFlags:
                    item.setCheckState(state)
            _update_summary()

        col_combo.currentIndexChanged.connect(_repopulate)
        type_list.itemChanged.connect(lambda _i: _update_summary())
        all_btn.clicked.connect(lambda: _set_all(True))
        none_btn.clicked.connect(lambda: _set_all(False))
        _repopulate()

        if dlg.exec() != QDialog.DialogCode.Accepted:
            return

        cell_types = _checked_types()
        if not cell_types:
            dialogs.warning(self, "Nothing selected",
                            "Tick at least one cell type.")
            return

        self._start_batch(col_combo.currentData(), cell_types)

    def _committed_fdr_genes(self):
        """The gene list BH should be corrected over, or None.

        Hypothesis mode narrows the FDR denominator to the committed
        sets; discovery corrects genome-wide. The DESeq2 fit stays
        genome-wide either way -- it needs the whole transcriptome for
        size factors and the dispersion trend -- so this only changes
        the denominator, and the results subset that follows it.

        Shared by the single run and the per-cell-type batch: the batch
        used to run discovery whatever mode was selected, because it
        never asked for this.
        """
        if (getattr(self.ws, 'analysis_mode', None) != 'scoring'
                or not self.ws.pathway_gene_sets):
            return None
        genes = set()
        for members in self.ws.pathway_gene_sets.values():
            genes.update(members)
        return genes or None

    def _start_batch(self, cell_type_col, cell_types):
        gse_id = self.ws.gse_accession or (
            Path(self.ws.h5ad_path).stem if self.ws.h5ad_path else "unknown")
        self._batch_worker = CellTypeBatchWorker(
            self.ws.current_adata, cell_type_col, self.ws.sample_col,
            self.ws.condition_col, self.ws.project_dir, gse_id,
            cell_types=cell_types,
            de_method=self.ws.de_method,
            moderate=self.eb_checkbox.isChecked(),
            covariates=self._selected_covariates(),
            # Same source the dialog's eligibility counts came from, so
            # what it promised is what actually runs.
            min_cells=self.min_cells_spin.value(),
            min_counts=self.min_counts_spin.value(),
            pathway_gene_sets=self.ws.pathway_gene_sets,
            fdr_genes=self._committed_fdr_genes(),
            counts_layer=self._selected_counts_layer(),
            detection_min_pct=self.detection_pct_spin.value() / 100.0,
            detection_min_donor_frac=(
                self.detection_donor_frac_spin.value() / 100.0),
            detection_study_col=self._per_study_filter_col(),
            deseq2_independent_filter=self.indep_filter_check.isChecked(),
            deseq2_cooks_filter=self.cooks_filter_check.isChecked(),
            filter_min_count=self.filter_min_count_spin.value(),
            filter_min_samples=int(self.filter_min_samples_spin.value()),
        )
        run_worker(
            self._batch_worker,
            on_finished=self._on_batch_finished,
            on_failed=self._on_batch_failed,
            on_progress=self._on_progress,
            on_progress_pct=self._on_progress_pct,
        )
        self.run_btn.setEnabled(False)
        self.batch_btn.setEnabled(False)
        self.ws.status_message.emit(
            f"Running DE across {len(cell_types)} cell types...")
        if self.ws.progress_bar:
            self.ws.progress_bar.setValue(0)

    def _on_batch_finished(self, payload):
        runs, _summary, message = payload
        self.run_btn.setEnabled(True)
        self.batch_btn.setEnabled(True)
        self.ws.log_message.emit(message)
        self.status_label.setText(message)
        self.ws.status_message.emit("Ready")
        if self.ws.progress_bar:
            self.ws.progress_bar.setValue(0)

        lines = []
        for run in runs:
            if run.status == 'ok':
                lines.append(
                    f"{run.cell_type}: {run.n_significant:,} significant of "
                    f"{run.n_genes_tested:,} tested ({run.n_samples} donors)")
            else:
                lines.append(f"{run.cell_type}: {run.status} -- {run.message}")
        # Each cell type is filtered on its own donors, so their gene
        # universes differ. Say so rather than let a later comparison
        # conflate "not significant" with "never tested".
        tail = ("\n\nEach cell type was filtered against its own donors, so "
                "their gene universes differ. Restrict to the shared genes "
                "before comparing counts between them.")

        self._batch_runs = [r for r in runs if r.status == 'ok'
                            and r.de_results is not None]
        self._populate_batch_view()
        if self._batch_runs:
            tail += ("\n\nUse the 'Cell type' picker under the volcano to "
                     "look at each one.")
        dialogs.info(self, "Batch DE complete",
                     message + "\n\n" + "\n".join(lines) + tail)

    def _populate_batch_view(self):
        """Offer the finished cell types, and show the first."""
        runs = getattr(self, '_batch_runs', [])
        combo = self._batch_view_combo
        combo.blockSignals(True)
        combo.clear()
        for run in runs:
            combo.addItem(
                f"{run.cell_type}  ({run.n_significant:,} sig)", run.cell_type)
        combo.blockSignals(False)

        show = bool(runs)
        self._batch_view_label.setVisible(show)
        combo.setVisible(show)
        if show:
            combo.setCurrentIndex(0)
            self._on_batch_view_changed(0)

    def _on_batch_view_changed(self, index):
        """Render the chosen cell type's results in the volcano + table.

        Only the DE-results views are switched. The Gene Expression and
        heatmap tabs read the loaded AnnData, which still holds every
        cell type -- so they are not showing this contrast's cells. The
        status line says so rather than letting the two quietly disagree.
        """
        runs = getattr(self, '_batch_runs', [])
        if not runs or not (0 <= index < len(runs)):
            return
        from kosmic.de.de_analysis import significant_subset

        run = runs[index]
        self.ws.de_results = run.de_results
        self.ws.significant_genes = significant_subset(
            run.de_results, fdr=self.pval_filter.value(),
            lfc=self.fc_filter.value())
        self.ws.gene_de_genome_wide = None

        self._populate_results()
        self._on_pval_toggle()
        self.status_label.setText(
            f"Showing {run.cell_type}: {run.n_significant:,} significant of "
            f"{run.n_genes_tested:,} tested, {run.n_samples} donors. "
            f"Per-gene expression tabs still show the whole loaded dataset.")

    def _on_batch_failed(self, message):
        self.run_btn.setEnabled(True)
        self.batch_btn.setEnabled(True)
        self.ws.log_message.emit(message)
        self.status_label.setText("Batch DE failed -- see log.")
        self.ws.status_message.emit("Ready")
        if self.ws.progress_bar:
            self.ws.progress_bar.setValue(0)

    def _on_progress(self, message):
        self.status_label.setText(message)
        self.ws.log_message.emit(message)

    def _on_progress_pct(self, pct):
        if self.ws.progress_bar:
            self.ws.progress_bar.setValue(pct)

    def _on_analysis_finished(self, payload):
        de_results, significant_genes, pathway_coverage, sample_df, message = payload
        self.run_btn.setEnabled(True)
        self.ws.log_message.emit(message)
        if self.ws.progress_bar:
            self.ws.progress_bar.setRange(0, 100)
            self.ws.progress_bar.setValue(0)

        # A whole-dataset run replaces whatever the batch picker showed.
        self._batch_runs = []
        if hasattr(self, '_batch_view_combo'):
            self._populate_batch_view()

        self.ws.de_results = de_results
        self.ws.significant_genes = significant_genes
        self.ws.pathway_coverage = pathway_coverage
        self.ws.sample_df = sample_df
        # Snapshot the settings this run used, so later edits flag a re-run.
        self._last_run_comparison = self._current_comparison()
        self._refresh_run_dirty()
        # Record this DE run + the upstream data token it ran against.
        self._record_de_provenance('gene_de')
        # In hypothesis mode de_results is restricted to the committed genes;
        # keep the genome-wide ranking for consumers whose null is the whole
        # transcriptome (e.g. fgsea). None in discovery mode (de_results is
        # already genome-wide).
        self.ws.gene_de_genome_wide = de_results.attrs.get('genome_wide_ranking')

        self._populate_results()

        self._auto_export_significant()
        self._auto_export_pathway_stats()
        self._auto_export_gene_level_stats()

        # Render through the shared toggle path so the initial volcano honours
        # the 'Pathway genes only' / 'Raw p-values' checkbox states.
        self._on_pval_toggle()
        self._populate_heatmap_combo()

        n_sig = len(significant_genes) if significant_genes is not None else 0
        self.ws.status_message.emit(
            f"Gene DE complete -- {len(de_results)} genes ({n_sig} significant)"
        )

        self._run_transcriptome_size_test()
        self.analysis_complete.emit()

    def _on_analysis_failed(self, message: str):
        self.run_btn.setEnabled(True)
        self.ws.log_message.emit(message)
        if self.ws.progress_bar:
            self.ws.progress_bar.setRange(0, 100)
            self.ws.progress_bar.setValue(0)
        self.ws.status_message.emit("Gene DE failed")
        dialogs.warning(self, "Analysis Failed", message)

    # --- Transcriptome size test ---

    def _run_transcriptome_size_test(self):
        """Check if total per-cell expression differs between conditions."""
        if not hasattr(self.ws, 'current_adata') or self.ws.current_adata is None:
            return
        if not self.ws.sample_col or not self.ws.condition_col:
            return
        if not self.ws.disease_label or not self.ws.control_label:
            return

        try:
            from kosmic.scrna.inspect.transcriptome_size import test_transcriptome_size_adata
            result = test_transcriptome_size_adata(
                self.ws.current_adata,
                sample_col=self.ws.sample_col,
                condition_col=self.ws.condition_col,
                disease_label=self.ws.disease_label,
                control_label=self.ws.control_label,
            )
            if not result:
                return

            umi = result['metrics']['mean_total_umi']
            ratio_pct = abs(umi['ratio'] - 1) * 100
            direction = 'smaller' if umi['ratio'] < 1 else 'larger'
            sig = ' *' if umi['pval'] < DEFAULT_FDR else ''

            self.ws.log_message.emit(
                f"Transcriptome size: disease ({result['n_disease']} donors) "
                f"is {ratio_pct:.1f}% {direction} than control "
                f"({result['n_control']} donors) -- "
                f"mean UMI/cell: {umi['mean_disease']:.0f} vs "
                f"{umi['mean_control']:.0f} "
                f"(ratio={umi['ratio']:.3f}, p={umi['pval']:.4f}{sig})")

            if umi['ratio'] < 0.9 or umi['ratio'] > 1.1:
                self.ws.log_message.emit(
                    "WARNING: Per-cell transcriptome size differs >10% between "
                    "conditions. CPM/DESeq2 normalization assumes equal "
                    "transcriptome size and may mask real global changes.")

            self._transcriptome_size_result = result
            self._transcriptome_btn.setEnabled(True)
        except Exception as e:
            self.ws.log_message.emit(f"Transcriptome size test error: {e}")

    def _on_transcriptome_btn(self):
        """Show the transcriptome size plot from the stored result."""
        if not hasattr(self, '_transcriptome_size_result') or self._transcriptome_size_result is None:
            dialogs.info(self, "No Data", "Run DE analysis first.")
            return
        self._show_transcriptome_size_plot(self._transcriptome_size_result)

    def _show_transcriptome_size_plot(self, result):
        """Show a three-panel plot of transcriptome metrics per donor."""
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from PyQt6.QtWidgets import QDialog, QVBoxLayout
        import pandas as pd
        from kosmic.gui.shared.theme import style_mpl_plot

        per_sample = pd.DataFrame(result['per_sample'])
        disease_label = self.ws.disease_label
        control_label = self.ws.control_label

        metrics = result['metrics']

        metric_keys = [
            ('mean_total_umi', 'Transcripts per cell\n(UMI counts)'),
            ('mean_genes_detected', 'Genes detected\nper cell'),
            ('mean_umi_per_gene', 'UMI per gene\nper cell'),
        ]
        metric_keys = [(k, label) for k, label in metric_keys if k in per_sample.columns]

        n_panels = len(metric_keys)
        fig = Figure(figsize=(4.5 * n_panels, 4), dpi=100)

        rng = np.random.default_rng(42)

        for panel_i, (col, ylabel) in enumerate(metric_keys):
            ax = fig.add_subplot(1, n_panels, panel_i + 1)

            c_vals = per_sample.loc[
                per_sample['condition'] == control_label, col].values
            d_vals = per_sample.loc[
                per_sample['condition'] == disease_label, col].values

            if len(c_vals) == 0 or len(d_vals) == 0:
                continue

            bp = ax.boxplot(
                [c_vals, d_vals],
                labels=[f'{control_label}\n(n={len(c_vals)})',
                        f'{disease_label}\n(n={len(d_vals)})'],
                widths=0.5,
                patch_artist=True,
            )
            bp['boxes'][0].set_facecolor('#4CAF50')
            bp['boxes'][0].set_alpha(0.4)
            bp['boxes'][1].set_facecolor('#F44336')
            bp['boxes'][1].set_alpha(0.4)

            for i, (vals, color) in enumerate(
                    [(c_vals, '#4CAF50'), (d_vals, '#F44336')]):
                jitter = rng.uniform(-0.1, 0.1, len(vals))
                ax.scatter(np.full(len(vals), i + 1) + jitter, vals,
                           c=color, alpha=0.7, s=30, zorder=3,
                           edgecolors='white', linewidths=0.5)

            m = metrics.get(col, {})
            ratio = m.get('ratio', d_vals.mean() / max(c_vals.mean(), 1e-10))
            pval = m.get('pval', 1.0)
            ratio_pct = abs(ratio - 1) * 100
            direction = 'smaller' if ratio < 1 else 'larger'
            sig_str = f"p={pval:.4f}"
            if pval < DEFAULT_FDR:
                sig_str += ' *'

            style_mpl_plot(fig, ax,
                           title=f'{disease_label} is {ratio_pct:.1f}% '
                                 f'{direction}\n({sig_str})',
                           left_label=ylabel)

        fig.suptitle('Transcriptome Size Test', fontsize=12, y=1.02)
        fig.tight_layout()

        if self.ws.project_dir:
            out_dir = de_analysis_dir(self.ws.project_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
            for fmt in ('png', 'pdf'):
                fig.savefig(out_dir / f'transcriptome_size_test.{fmt}',
                            dpi=300, bbox_inches='tight')
            self.ws.log_message.emit(
                f"Saved: {out_dir / 'transcriptome_size_test.png'}")

        dlg = QDialog(self)
        dlg.setWindowTitle('Transcriptome Size Test')
        dlg.setMinimumSize(400 * n_panels, 450)
        layout = QVBoxLayout(dlg)
        canvas = FigureCanvas(fig)
        layout.addWidget(canvas)
        dlg.exec()

    # --- Meta-analysis export (background) ---

    def _export_for_meta(self):
        """Export the chosen method's DE result + pseudobulk for meta-analysis.

        The other DE methods are not recomputed and the h5ad is not rewritten;
        the meta-analysis reads the per-study CSV directly.
        """
        gse_id = self.ws.gse_accession or (
            Path(self.ws.h5ad_path).stem if self.ws.h5ad_path else "unknown")
        self._meta_export_worker = MetaExportWorker(
            self.ws.current_adata, self.ws.sample_col, self.ws.condition_col,
            self.ws.min_cells, self.ws.de_method, self.eb_checkbox.isChecked(),
            self.ws.de_results, self.ws.project_dir, gse_id,
            min_counts=self.ws.min_counts,
        )
        if self.ws.progress_bar:
            self.ws.progress_bar.setValue(0)
        self.status_label.setText("Exporting DE for meta-analysis...")
        self.ws.status_message.emit("Exporting DE for meta-analysis...")
        run_worker(
            self._meta_export_worker,
            on_finished=self._on_meta_export_finished,
            on_failed=self._on_meta_export_failed,
            on_progress=self._on_meta_export_progress,
            on_progress_pct=self._on_meta_export_pct,
        )

    def _on_meta_export_progress(self, msg):
        self.ws.log_message.emit(f"[Meta-export] {msg}")
        self.status_label.setText(msg)

    def _on_meta_export_pct(self, pct):
        if self.ws.progress_bar:
            self.ws.progress_bar.setValue(pct)

    def _on_meta_export_finished(self, message):
        self.ws.log_message.emit(f"[Meta-export] {message}")
        self.status_label.setText(message)
        self.ws.status_message.emit("Ready")
        if self.ws.progress_bar:
            self.ws.progress_bar.setValue(0)

    def _on_meta_export_failed(self, message: str):
        self.ws.log_message.emit(f"[Meta-export] {message}")
        self.status_label.setText(message)
        self.ws.status_message.emit("Ready")
        if self.ws.progress_bar:
            self.ws.progress_bar.setValue(0)

    # --- Results display ---

    def _populate_results(self):
        self._filter_results()

    def _filter_results(self, _=None):
        if self.ws.de_results is None:
            return

        pval_threshold = self.pval_filter.value()
        fc_threshold = self.fc_filter.value()
        pathway_only = self._pathway_only_check.isChecked()

        df = self.ws.de_results.copy()

        use_raw_pvals = (getattr(self, '_raw_pval_check', None) is not None
                         and self._raw_pval_check.isChecked())

        gene_search_text = ''
        if getattr(self, '_volcano_gene_search', None) is not None:
            gene_search_text = self._volcano_gene_search.text().strip()

        # 'Pathway genes only' is a pure display filter: it hides non-pathway
        # rows but never touches FDR. Significance stays as the DE engine called
        # it over the full tested gene set, whatever subset is on screen.
        pathway_subset_active = False
        if pathway_only and self.ws.pathway_coverage:
            pathway_genes = set()
            for info in self.ws.pathway_coverage.values():
                pathway_genes.update(info['genes'])
            if pathway_genes:
                df = df[df['names'].isin(pathway_genes)].copy()
                pathway_subset_active = True

        # Captured before the name filter narrows further.
        fdr_scope = len(df)

        if gene_search_text:
            needles = [g.strip().upper() for g in gene_search_text.replace(';', ',').split(',') if g.strip()]
            if needles and 'names' in df.columns:
                names_upper = df['names'].astype(str).str.upper()
                mask = np.zeros(len(df), dtype=bool)
                for needle in needles:
                    mask |= names_upper.str.contains(needle, regex=False, na=False)
                df = df[mask].copy()

        sig_pcol = 'pvals' if use_raw_pvals else 'pvals_adj'
        df['_significant'] = (
            (df[sig_pcol].astype(float) < pval_threshold)
            & (df['logfoldchanges'].astype(float).abs() > fc_threshold))

        schema = [
            Column('Gene',        'names',          's'),
            Column('Pathway',     'pathways',       's'),
            Column('log2FC',      'logfoldchanges', '.3f'),
        ]
        # The 95% interval on the fold change, beside the fold change.
        # The FDR answers "did this change?"; it says nothing about how
        # big the change is, and a point estimate on its own invites the
        # reader to treat it as exact. 'log2FC 2.88 [1.42, 4.33]' makes
        # the uncertainty impossible to miss and needs no new statistics
        # -- it is lfc +/- 1.96 x se, already on every row.
        if 'ci_lower' in df.columns and 'ci_upper' in df.columns:
            df = df.copy()
            lo = df['ci_lower'].astype(float)
            hi = df['ci_upper'].astype(float)
            df['_ci'] = [
                '' if not (np.isfinite(a) and np.isfinite(b))
                else f"[{a:.2f}, {b:.2f}]"
                for a, b in zip(lo, hi)]
            schema.append(Column('95% CI', '_ci', 's'))
        schema.extend([
            Column('p-value',     'pvals',          '.2e'),
            Column('FDR',         'pvals_adj',      '.2e'),
        ])
        if 'pct_disease' in df.columns and 'pct_control' in df.columns:
            df = df.copy()
            df['pct_disease_pct'] = df['pct_disease'].astype(float) * 100
            df['pct_control_pct'] = df['pct_control'].astype(float) * 100
            schema.extend([
                Column('% Disease', 'pct_disease_pct', '.1f'),
                Column('% Control', 'pct_control_pct', '.1f'),
            ])
        # Surface DESeq2's independent / Cook's filtering: filtered genes show
        # a blank FDR, so a tag tells them apart from tested-but-non-significant.
        if 'filter_status' in df.columns and (df['filter_status'] != 'tested').any():
            df = df.copy()
            _flabel = {'tested': '', 'independent_filter': 'indep. filter',
                       'cooks_outlier': "Cook's"}
            df['_filter_label'] = df['filter_status'].map(_flabel).fillna('')
            schema.append(Column('DESeq2 filter', '_filter_label', 's'))
        # 'Passes' rather than 'Significant': the column is the AND of an
        # FDR cut and a fold-change cut, and only the first is a test.
        schema.append(Column('Passes', '_significant', 'bool'))
        self.results_table.set_schema(schema)
        self.results_table.set_data(df)

        sig_count = int(df['_significant'].sum())
        total = len(df)
        if use_raw_pvals:
            scope_note = "  (raw p-values, no MTC)"
        else:
            scope_note = ""
        self.sig_count_label.setText(
            f"{sig_count} significant / {total} genes{scope_note}")

        self._update_count_readout(pval_threshold, fc_threshold, use_raw_pvals)

        active = []
        if pathway_subset_active:
            active.append("pathway")
        if gene_search_text:
            active.append(f"name='{gene_search_text}'")
        if use_raw_pvals:
            active.append("raw-p")
        filter_str = ", ".join(active) if active else "none"
        if self.ws.unit == 'cell':
            engine = "Wilcoxon"
        elif self.ws.de_method == 'deseq2':
            engine = "DESeq2"
        else:
            engine = "Welch's t-test"
        full_total = len(self.ws.de_results) if self.ws.de_results is not None else 0
        self.ws.log_message.emit(
            f"DE filter: {full_total:,} in results by {engine} → "
            f"{fdr_scope:,} shown → "
            f"{sig_count:,} significant  filters=[{filter_str}]"
        )

    def _update_count_readout(self, pval_threshold, fc_threshold, use_raw_pvals):
        """Mode-aware gene funnel shown in the Display-filters panel.

        Hypothesis mode: Universe -> Pathway genes -> Tested -> Significant
        (the pathway is the focus). Discovery mode: Universe -> Tested ->
        Significant genome-wide, with Pathway genes as a trailing context
        number. 'Tested' counts genes that got a valid result (a real p-value,
        not NaN'd out by DESeq2's independent / Cook's filtering) -- fewer than
        the universe merely run through DESeq2 for normalisation.
        """
        if self._filter_count_label is None or self.ws.de_results is None:
            return
        de = self.ws.de_results
        attrs = de.attrs
        universe = attrs.get('n_genes_input')

        pathway_present = None
        cov = getattr(self.ws, 'pathway_coverage', None)
        if cov:
            pw = set()
            for c in cov.values():
                pw.update(c.get('genes', []))
            pathway_present = len(pw)

        # Valid-tested: a real result, not NaN'd by DESeq2's filtering.
        sig_pcol = 'pvals' if use_raw_pvals else 'pvals_adj'
        if 'filter_status' in de.columns:
            valid = (de['filter_status'] == 'tested').to_numpy()
        elif sig_pcol in de.columns:
            valid = de[sig_pcol].notna().to_numpy()
        else:
            valid = np.ones(len(de), dtype=bool)
        n_tested = int(valid.sum())

        if sig_pcol in de.columns:
            p = de[sig_pcol].to_numpy(dtype=float)
            lfc = de['logfoldchanges'].to_numpy(dtype=float)
            n_sig = int((valid & (p < pval_threshold)
                         & (np.abs(lfc) > fc_threshold)).sum())
        else:
            n_sig = 0

        hypothesis = getattr(self.ws, 'analysis_mode', None) == 'scoring'
        lines = []
        if universe:
            lines.append(f"Universe genes: {universe:,}")
        if hypothesis and pathway_present is not None:
            lines.append(f"Pathway genes: {pathway_present:,}")
        lines.append(f"Tested: {n_tested:,}")
        lines.append(f"Significant: {n_sig:,}")
        if not hypothesis and pathway_present is not None:
            lines.append(f"Pathway genes: {pathway_present:,}")

        fcounts = attrs.get('deseq2_filter_counts', {})
        if fcounts.get('independent_filter', 0):
            lines.append(f"Indep-filtered: {fcounts['independent_filter']:,}")
        if fcounts.get('cooks_outlier', 0):
            lines.append(f"Cook's: {fcounts['cooks_outlier']:,}")

        self._filter_count_label.setText("  ·  ".join(lines))
        scope = "pathway genes" if hypothesis else "all genes"
        self._filter_count_label.setToolTip(
            "Universe genes: the whole transcriptome run through DESeq2 for "
            "normalisation.\n"
            "Pathway genes: genes in the selected pathway set(s) present in the "
            "data.\n"
            f"Tested: {scope} that got a valid test -- passed the detection "
            "pre-filter and returned a real p-value (not dropped by DESeq2's "
            "independent / Cook's filtering).\n"
            "Significant: tested genes with FDR < threshold and |log2FC| > "
            "threshold.\n"
            "Indep-filtered / Cook's: genes DESeq2 excluded from the FDR "
            "correction.")

    def _on_gene_selected(self, gene_name):
        self.ws.log_message.emit(f"Selected gene: {gene_name}")
        self.results_table.select_first_by_text(0, gene_name)
        if hasattr(self, '_ge_gene_input'):
            self._ge_gene_input.setText(gene_name)

    def _on_display_filter_toggled(self, _checked=None):
        """A Display-filters checkbox changed: refilter the table + volcano."""
        self._filter_results()
        self._on_pval_toggle()

    def _on_pval_toggle(self, _=None):
        """Refresh the volcano with the current checkbox states."""
        if self.ws.de_results is None:
            return
        use_raw = self._raw_pval_check.isChecked()
        pathway_only = self._pathway_only_check.isChecked()

        pathway_genes = set()
        if self.ws.pathway_coverage:
            for info in self.ws.pathway_coverage.values():
                pathway_genes.update(info['genes'])

        # 'Pathway genes only' both filters to the panel and drives the pathway
        # overlay (panel genes coloured by direction, others greyed). With it
        # off, show every gene coloured by significance the standard way -- no
        # overlay -- so discovery mode is a normal genome-wide volcano.
        de_data = self.ws.de_results
        overlay = None
        if pathway_only and pathway_genes:
            de_data = de_data[de_data['names'].isin(pathway_genes)].copy()
            overlay = pathway_genes

        _render_de_volcano(
            self.interactive_volcano, de_data,
            pathway_genes=overlay,
            use_raw_pvals=use_raw,
            disease_label=self.ws.disease_label,
            control_label=self.ws.control_label,
        )
        self._refresh_volcano_gene_completer(de_data)
        # Re-apply the active highlight after re-render.
        self._on_volcano_gene_search(self._volcano_gene_search.text())

    def _refresh_volcano_gene_completer(self, de_results):
        """Sync every gene autocomplete to the current results.

        Three inputs take gene names -- the volcano's Find gene, the
        heatmap's row list and the Gene Expression box -- and they all
        draw from the same set, so they are refreshed together rather
        than drifting apart.
        """
        models = [self._volcano_gene_search_model]
        for attr in ('_heatmap_gene_model', '_ge_gene_model'):
            model = getattr(self, attr, None)
            if model is not None:
                models.append(model)

        if de_results is None or 'names' not in de_results.columns:
            for model in models:
                model.setStringList([])
            return
        genes = sorted(de_results['names'].astype(str).unique())
        for model in models:
            model.setStringList(genes)

    def _on_volcano_gene_search(self, text):
        """
        Highlight matching genes on the volcano with a magenta star.

        Accepts a comma/semicolon-separated panel (e.g. "VEGFA, KDR, FLT1");
        each token is a case-insensitive substring against visible genes.
        """
        text = (text or '').strip()
        if not text:
            self.interactive_volcano.clear_highlight()
            return
        needles = [g.strip() for g in text.replace(';', ',').split(',') if g.strip()]
        if not needles:
            self.interactive_volcano.clear_highlight()
            return

        self.interactive_volcano.highlight_genes(
            needles,
            symbol='star',
            color=(255, 255, 255),                # white outline
            fill_color=(255, 64, 192, 255),       # hot magenta fill
            size=18,
            pen_width=2.5,
            with_labels=True,
        )

    # --- Plot generation ---

    def _populate_heatmap_combo(self):
        """Populate the heatmap row-source combo.

        Top-N entries come first and are what discovery mode uses; the
        pathway entries only appear when gene sets are loaded. The
        heatmap itself is drawn lazily when the tab is opened (it
        computes the expression matrix), so this never blocks page entry.
        """
        self.heatmap_pathway_combo.blockSignals(True)
        self.heatmap_pathway_combo.clear()
        if self.ws.de_results is not None:
            for n in (20, 30, 50, 100):
                self.heatmap_pathway_combo.addItem(f"Top {n} by FDR")
        if self.ws.pathway_gene_sets:
            for pw_name in sorted(self.ws.pathway_gene_sets.keys()):
                self.heatmap_pathway_combo.addItem(pw_name)
        # Hypothesis runs are about the committed sets, so start there.
        if (self.ws.pathway_gene_sets
                and getattr(self.ws, 'analysis_mode', None) == 'scoring'):
            self.heatmap_pathway_combo.setCurrentIndex(
                self.heatmap_pathway_combo.count()
                - len(self.ws.pathway_gene_sets))
        self.heatmap_pathway_combo.blockSignals(False)

    def _expr_normalization(self):
        """Viz normalisation matching the selected DE method: VST for DESeq2,
        log-CPM for the Welch engine, so plots show the quantity the test used."""
        return 'vst' if self.ws.de_method == 'deseq2' else 'cpm'

    def _bars_normalization(self):
        """Absolute-expression transform for the bar charts, matching the DE
        method: log2 size-factor-normalised counts for DESeq2, log2 CPM for
        Welch. Both are depth-corrected and put zero counts at zero (unlike
        VST), so undetected genes read as empty bars."""
        return 'deseq2' if self.ws.de_method == 'deseq2' else 'cpm'

    def _bars_label(self):
        return ('log2 norm. counts' if self._bars_normalization() == 'deseq2'
                else 'log2 CPM')

    def _expr_matrix(self, norm):
        """Cached per-donor expression matrix under the given normalisation.

        The heatmap requests the DE-method-matched transform (VST for DESeq2,
        log-CPM for Welch); the expression bars always request 'cpm' (log2 CPM)
        because that is an interpretable absolute unit with zeros at zero, where
        VST's arbitrary scale + non-zero floor would mislead. Computed once per
        (dataset, layer, sample/condition col, normalisation) and cached.
        """
        adata = self.ws.current_adata
        if adata is None or not self.ws.sample_col or not self.ws.condition_col:
            return None
        counts_layer = self._selected_counts_layer()
        key = (id(adata), counts_layer, self.ws.sample_col,
               self.ws.condition_col, norm)
        cache = getattr(self, '_expr_cache', None)
        if cache is None:
            cache = self._expr_cache = {}
        if key in cache:
            return cache[key]

        from kosmic.de.de_analysis import pseudobulk_expression_matrix
        nice = 'VST' if norm == 'vst' else 'log-CPM'
        self.status_label.setText(f"Computing {nice} expression for plots...")
        self.ws.status_message.emit(f"Computing {nice} expression for plots...")
        if self.ws.progress_bar:
            self.ws.progress_bar.setRange(0, 0)  # indeterminate spinner
        # Flush the busy state to screen before the blocking compute so the
        # user sees "Computing..." + spinner rather than an apparent freeze.
        QApplication.processEvents()
        try:
            df = pseudobulk_expression_matrix(
                adata, self.ws.sample_col, self.ws.condition_col,
                normalization=norm, counts_layer=counts_layer,
                min_cells=self.ws.min_cells, min_counts=self.ws.min_counts)
        finally:
            self.status_label.setText("")
            self.ws.status_message.emit("Ready")
            if self.ws.progress_bar:
                self.ws.progress_bar.setRange(0, 100)
                self.ws.progress_bar.setValue(0)
        cache[key] = df
        return df

    def _heatmap_row_genes(self, selection):
        """Genes to draw, from the typed list, a pathway, or the top hits."""
        typed = (self._heatmap_gene_input.text().strip()
                 if hasattr(self, '_heatmap_gene_input') else '')
        if typed:
            wanted = [g.strip().upper() for g in typed.split(',') if g.strip()]
            de = self.ws.de_results
            if de is None or 'names' not in de.columns:
                return []
            upper = {str(g).upper(): str(g) for g in de['names']}
            return [upper[g] for g in wanted if g in upper]

        if selection and selection.startswith('Top '):
            de = self.ws.de_results
            if de is None or 'pvals_adj' not in de.columns:
                return []
            try:
                n = int(selection.split()[1])
            except (IndexError, ValueError):
                n = 30
            ranked = de.dropna(subset=['pvals_adj']).nsmallest(n, 'pvals_adj')
            return [str(g) for g in ranked['names']]

        return list(self.ws.pathway_gene_sets.get(selection, [])
                    if self.ws.pathway_gene_sets else [])

    def _generate_heatmap(self, pathway_name):
        """Draw the gene heatmap for whichever row source is selected."""
        pathway_genes = self._heatmap_row_genes(pathway_name)
        if not pathway_genes:
            return

        try:
            expr_df = self._expr_matrix(self._expr_normalization())
            if expr_df is None or expr_df.empty:
                return

            # Restrict to genes that were actually in the DE (passed the
            # detection filter). Undetectable genes are dropped so their
            # z-scored rows -- pure standardised noise -- don't appear.
            de = self.ws.de_results
            tested = (set(de['names'].astype(str))
                      if de is not None and 'names' in de.columns else None)
            available_genes = [
                g for g in pathway_genes
                if g in expr_df.columns and (tested is None or g in tested)]
            if len(available_genes) < 2:
                return

            # Significance markers (up / down in disease) for the retained rows.
            gene_markers = {}
            if (de is not None
                    and {'names', 'pvals_adj', 'logfoldchanges'} <= set(de.columns)):
                pth = self.pval_filter.value()
                fth = self.fc_filter.value()
                sub = de[de['names'].isin(available_genes)]
                p = sub['pvals_adj'].astype(float)
                lfc = sub['logfoldchanges'].astype(float)
                sig = sub[(p < pth) & (lfc.abs() > fth)]
                for name, lv in zip(sig['names'].astype(str),
                                    sig['logfoldchanges'].astype(float)):
                    gene_markers[name] = 'up' if lv > 0 else 'down'

            sample_gene_df = expr_df[['sample', 'condition'] + available_genes].copy()

            from kosmic.gui.shared.theme import get_plot_settings
            ps = get_plot_settings()
            self.interactive_heatmap.set_data(
                sample_gene_df, available_genes, pathway_name,
                self.ws.control_label, self.ws.disease_label,
                cmap=ps['heatmap_cmap'],
                vmin_pct=ps['heatmap_vmin'],
                vmax_pct=ps['heatmap_vmax'],
                gene_markers=gene_markers,
            )
        except Exception as e:
            self.ws.log_message.emit(f"Heatmap error for {pathway_name}: {e}")

    # --- Export methods ---

    # --- Auto-export helpers (pathway stats, gene stats, settings) ---

    def _get_output_dir(self):
        if self.ws.project_dir:
            return de_analysis_dir(self.ws.project_dir)
        return Path.cwd() / "de_analysis_results"

    def _auto_export_significant(self):
        """Write the significant-gene table beside the full one.

        Every run already writes '{accession}_DE_{method}.csv'; this is
        the same table cut to the FDR and log2FC thresholds currently set
        on the results tab, which is what the counts in the app report.
        It used to be a button, which meant the file existed only if
        someone remembered to press it.
        """
        from kosmic.de.de_analysis import significant_subset

        if self.ws.de_results is None or self.ws.de_results.empty:
            return
        if not self.ws.project_dir:
            return

        gse_id = self.ws.gse_accession or (
            Path(self.ws.h5ad_path).stem if self.ws.h5ad_path else "unknown")
        try:
            sig = significant_subset(
                self.ws.de_results,
                fdr=self.pval_filter.value(),
                lfc=self.fc_filter.value())
            out = significant_path(self.ws.project_dir, gse_id)
            out.parent.mkdir(parents=True, exist_ok=True)
            sig.to_csv(out, index=False)
            self.ws.log_message.emit(
                f"Exported {len(sig):,} genes with FDR<"
                f"{self.pval_filter.value():g} and estimated |log2FC|>"
                f"{self.fc_filter.value():g}: {out.name}")
        except Exception as e:
            self.ws.log_message.emit(
                f"Warning: could not export significant genes: {e}")

    def _auto_export_pathway_stats(self):
        """Export pathway-level statistics for downstream meta-analysis."""
        from kosmic.de.sample_stats_export import compute_pathway_stats

        if self.ws.current_adata is None or self.ws.sample_df is None:
            return
        if not all([self.ws.sample_col, self.ws.condition_col,
                    self.ws.control_label, self.ws.disease_label]):
            return

        gse_id = self.ws.gse_accession or (
            Path(self.ws.h5ad_path).stem if self.ws.h5ad_path else "unknown")

        try:
            norm_df, raw_df = compute_pathway_stats(
                self.ws.current_adata, self.ws.sample_df,
                self.ws.sample_col, self.ws.condition_col,
                self.ws.control_label, self.ws.disease_label,
                self.ws.pathway_gene_sets, gse_id,
            )
            if norm_df.empty:
                return

            stats_dir = de_stats_dir(self.ws.project_dir)
            stats_dir.mkdir(parents=True, exist_ok=True)
            norm_path = stats_dir / f"{gse_id}_PathwayStats_Norm.csv"
            raw_path = stats_dir / f"{gse_id}_PathwayStats_Raw.csv"
            norm_df.to_csv(norm_path, index=False)
            raw_df.to_csv(raw_path, index=False)
            self.ws.log_message.emit(
                f"Exported pathway stats: {norm_path.name}, {raw_path.name}")
        except Exception as e:
            self.ws.log_message.emit(
                f"Warning: Could not export pathway stats: {str(e)}")

    def _auto_export_gene_level_stats(self):
        """Export gene-level statistics for downstream meta-analysis."""
        from kosmic.de.sample_stats_export import compute_gene_level_stats

        if self.ws.current_adata is None or self.ws.sample_df is None:
            return
        if not all([self.ws.sample_col, self.ws.condition_col,
                    self.ws.control_label, self.ws.disease_label]):
            return

        gse_id = self.ws.gse_accession or (
            Path(self.ws.h5ad_path).stem if self.ws.h5ad_path else "unknown")

        try:
            norm_df, raw_df = compute_gene_level_stats(
                self.ws.current_adata, self.ws.sample_df,
                self.ws.sample_col, self.ws.condition_col,
                self.ws.control_label, self.ws.disease_label,
                self.ws.pathway_gene_sets, gse_id,
            )
            if norm_df.empty:
                self.ws.log_message.emit(
                    "Warning: No pathway genes found in data")
                return

            self.ws.log_message.emit(
                f"Exporting gene-level stats for {len(norm_df)} genes...")
            stats_dir = de_stats_dir(self.ws.project_dir)
            stats_dir.mkdir(parents=True, exist_ok=True)
            norm_path = stats_dir / f"{gse_id}_GeneStats_Norm.csv"
            raw_path = stats_dir / f"{gse_id}_GeneStats_Raw.csv"
            norm_df.to_csv(norm_path, index=False)
            raw_df.to_csv(raw_path, index=False)
            self.ws.log_message.emit(
                f"Exported gene stats: {norm_path.name}, {raw_path.name}")
        except Exception as e:
            self.ws.log_message.emit(
                f"Warning: Could not export gene-level stats: {str(e)}")
            import traceback
            traceback.print_exc()

    # --- Gene Expression tab ---

    def _draw_expression_bars(self, plot, expr_df, genes):
        """Draw per-donor log-CPM expression bars (control vs disease) for
        'genes'. Gene names are drawn as angled TextItems below the axis
        (pyqtgraph cannot rotate tick labels). Returns the genes drawn."""
        import pyqtgraph as pg
        from PyQt6.QtGui import QFont
        from kosmic.gui.shared.theme import get_color, get_font_sizes
        genes_used = [g for g in genes if g in expr_df.columns]
        if not genes_used:
            return []
        ctrl_color = get_color('plot_control')
        dis_color = get_color('plot_disease')
        fg = get_color('fg_primary')
        dot_rgba = _parse_rgba(get_color('plot_dot'))
        cond = expr_df['condition'].astype(str).values
        control_mask = cond == str(self.ws.control_label)
        disease_mask = cond == str(self.ws.disease_label)

        plot.clear_plot_items()
        plot.hide_unavailable_message()
        ymax = 0.0
        for i, gene in enumerate(genes_used):
            vals = expr_df[gene].to_numpy(dtype=float)
            _draw_grouped_bars_with_dots(
                plot, i, vals[control_mask], vals[disease_mask],
                ctrl_color=ctrl_color, dis_color=dis_color,
                fg_color=fg, dot_rgba=dot_rgba, seed=42 + i,
            )
            if vals.size:
                ymax = max(ymax, float(np.nanmax(vals)))
        ymax = ymax if ymax > 0 else 1.0

        # Angled gene labels below the axis (pyqtgraph can't rotate tick text).
        plot.getAxis('bottom').setTicks([[]])
        lbl_pt = max(6, int(get_font_sizes()['tick']) - (2 if len(genes_used) > 12 else 0))
        y_lbl = -0.04 * ymax
        for i, gene in enumerate(genes_used):
            t = pg.TextItem(gene, color=fg, anchor=(1, 0.5), angle=30)
            t.setFont(QFont('Arial', lbl_pt))
            t.setPos(i + 0.2, y_lbl)
            plot.addItem(t)

        plot.setLabel('left', self._bars_label())
        plot.setXRange(-0.6, len(genes_used) - 0.4, padding=0)
        plot.setYRange(-0.34 * ymax, ymax * 1.08, padding=0)
        return genes_used

    def _run_gene_expression(self):
        """
        Plot per-condition expression for the gene(s) in the input field.

        Two panels: per-donor normalised expression (VST or log-CPM, matching
        the DE method) and per-donor % expressing (from the cell-level matrix).
        """
        from kosmic.gui.shared.theme import get_color
        from scipy import sparse

        gene_text = self._ge_gene_input.text().strip()
        if not gene_text or self.ws.current_adata is None:
            return
        if not self.ws.sample_col or not self.ws.condition_col:
            return

        adata = self.ws.current_adata
        genes = [g.strip() for g in gene_text.replace(';', ',').split(',') if g.strip()]
        missing = [g for g in genes if g not in adata.var_names]
        genes = [g for g in genes if g in adata.var_names]
        if missing:
            self.ws.log_message.emit(f"{len(missing)} gene(s) not found: {', '.join(missing[:5])}")
        if not genes:
            return

        expr_df = self._expr_matrix(self._bars_normalization())
        if expr_df is None or expr_df.empty:
            return
        genes_used = self._draw_expression_bars(self._ge_mean_plot, expr_df, genes)
        if not genes_used:
            return

        # % expressing panel: fraction of cells with a non-zero count per donor
        # (detection, normalisation-free), over the same donors as the bars.
        ctrl_color = get_color('plot_control')
        dis_color = get_color('plot_disease')
        fg = get_color('fg_primary')
        dot_color = _parse_rgba(get_color('plot_dot'))
        cond = expr_df['condition'].astype(str).values
        control_mask = cond == str(self.ws.control_label)
        disease_mask = cond == str(self.ws.disease_label)
        sample_ids = expr_df['sample'].astype(str).values
        obs_sample = adata.obs[self.ws.sample_col].astype(str).values

        self._ge_pct_plot.clear_plot_items()
        self._ge_pct_plot.hide_unavailable_message()
        for i, gene in enumerate(genes_used):
            var_idx = list(adata.var_names).index(gene)
            col = adata.X[:, var_idx]
            col = col.toarray().flatten() if sparse.issparse(col) else np.asarray(col).flatten()
            c_pcts, d_pcts = [], []
            for s_idx, sid in enumerate(sample_ids):
                cell_mask = obs_sample == sid
                if cell_mask.sum() == 0:
                    continue
                pct = (col[cell_mask] > 0).mean() * 100
                if control_mask[s_idx]:
                    c_pcts.append(pct)
                elif disease_mask[s_idx]:
                    d_pcts.append(pct)
            _draw_grouped_bars_with_dots(
                self._ge_pct_plot, i,
                np.asarray(c_pcts), np.asarray(d_pcts),
                ctrl_color=ctrl_color, dis_color=dis_color,
                fg_color=fg, dot_rgba=dot_color, seed=42 + i,
            )
        self._ge_pct_plot.getAxis('bottom').setTicks(
            [[(i, g) for i, g in enumerate(genes_used)]])
        self._results_tabs.setCurrentWidget(self._results_tabs.widget(
            self._results_tabs.indexOf(self._ge_mean_plot.parent().parent())
        ))
        self.ws.log_message.emit(f"Gene expression plotted: {', '.join(genes_used)}")

    # --- Pathway Genes tab ---

    def _on_pw_gene_combo_changed(self, pathway_name):
        """Plot per-gene expression bars for the selected pathway."""
        from kosmic.gui.shared.theme import get_color

        if not pathway_name or self.ws.current_adata is None:
            return
        genes = list(self.ws.pathway_gene_sets.get(pathway_name, []))
        if not genes:
            return

        adata = self.ws.current_adata
        if not self.ws.sample_col or not self.ws.condition_col:
            return

        genes = [g for g in genes if g in adata.var_names]
        if not genes:
            self.ws.log_message.emit(f"No genes from '{pathway_name}' found in dataset")
            return

        expr_df = self._expr_matrix(self._bars_normalization())
        if expr_df is None or expr_df.empty:
            return
        genes_used = self._draw_expression_bars(self._pw_gene_bar_plot, expr_df, genes)
        if not genes_used:
            return

        self._pw_gene_bar_plot.setTitle(
            f"{pathway_name} -- {len(genes_used)} genes",
            color=get_color('fg_primary'))

    def _populate_pw_gene_combo(self):
        """Fill the pathway combo from workspace gene sets. The per-gene bars
        are drawn lazily when the Pathway Genes tab is opened (they compute the
        expression matrix), so this never blocks page entry."""
        self._pw_gene_combo.blockSignals(True)
        self._pw_gene_combo.clear()
        if self.ws.pathway_gene_sets:
            for pw_name in self.ws.pathway_gene_sets.keys():
                self._pw_gene_combo.addItem(pw_name)
        self._pw_gene_combo.blockSignals(False)

    # --- Top DE Genes tab ---

    def _on_topde_n_changed(self, value):
        if self.ws.de_results is not None:
            self._topde_plot.top_n = value
            if self._topde_plot.current_level == "gene" and self._topde_plot.current_pathway:
                self._topde_plot.show_gene_level(self._topde_plot.current_pathway)
            elif self._topde_plot.current_level == "pathway":
                self._topde_plot.show_pathway_level()

    def _on_topde_level_changed(self, level):
        self._topde_back_btn.setVisible(level == "gene")

    def _update_topde(self):
        """Update the Top DE Genes plot from current DE results."""
        if self.ws.de_results is not None and self.ws.pathway_gene_sets:
            self._topde_plot.set_data(
                self.ws.de_results, self.ws.pathway_gene_sets,
                self._topde_n_spin.value()
            )

    def refresh_theme(self):
        """Re-apply theme to interactive plots."""
        from kosmic.gui.shared.theme import style_pg_plot
        if hasattr(self, 'interactive_volcano'):
            self.interactive_volcano.refresh_theme()
        if hasattr(self, 'interactive_heatmap'):
            self.interactive_heatmap.refresh_theme()
        if hasattr(self, '_ge_mean_plot'):
            style_pg_plot(self._ge_mean_plot, title='Mean Expression per Condition',
                          left_label='Mean Expression')
            style_pg_plot(self._ge_pct_plot, title='% Cells Expressing Gene',
                          left_label='% Expressing')
        if hasattr(self, '_pw_gene_bar_plot'):
            style_pg_plot(self._pw_gene_bar_plot, title='Per-Gene Expression',
                          left_label='Mean Expression')
        if hasattr(self, '_topde_plot'):
            self._topde_plot.refresh_theme()
