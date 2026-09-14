# Marker Check tab
# Sense-check an annotation: plot known markers for a few of your cell
# types against those same types. Rows are markers grouped by the type
# they mark, columns are the types you picked, so a correct annotation
# reads as a diagonal.
#
# The picker lists *your* cell types, not the marker database's. Those
# vocabularies rarely agree -- CellTypist says "Ventricular
# Cardiomyocyte" where a curated set says "Cardiomyocyte" -- and picking
# on the database's axis means picking on one nobody thinks in.
#
# Defaults to four types, because the readable version of this plot is
# small. This replaced a gene-family browser that drew the union of a
# whole collection (818 genes for GPCRs) against every group at once.
# Those families are still reachable: a bundled .gmt loads through
# "From file...", because a GMT is already {name: [genes]}, the same
# shape a marker panel has.

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import (
    QComboBox, QFileDialog, QHBoxLayout, QLabel, QListWidget,
    QListWidgetItem, QSpinBox,
)

from kosmic.gui.shared import dialogs, run_worker
from kosmic.gui.shared.plots import GeneGroupDotPlot
from kosmic.gui.shared.widgets import (
    BaseWorker, CaptionLabel, Column, ResultsTable, SecondaryButton,
    SecondaryLabel, SectionHeader, SidebarPage,
)
from kosmic.scrna.annotate.score import DEFAULT_MARKERS
from kosmic.scrna.inspect.gene_group import (
    compare_annotations, marker_panel, match_marker_types,
)

_MIN_GROUPS = 2
_MAX_GROUPS = 60

_SOURCE_BUILTIN = "Built-in curated"
_SOURCE_PANGLAO = "PanglaoDB"
_SOURCE_CELLMARKER = "CellMarker2"
_SOURCE_FILE = "From file..."

#: How many cell types to tick on first arrival. Four types at up to six
#: markers each is ~24 rows, which fits on screen and can be read.
_DEFAULT_PICKED = 4


class _SummaryWorker(BaseWorker):
    """Compute a GeneGroupSummary off the UI thread."""

    def __init__(self, adata, genes, group_col, parent=None):
        super().__init__(parent)
        self._adata = adata
        self._genes = genes
        self._group_col = group_col

    def _run(self):
        from kosmic.scrna.inspect.gene_group import summarise_gene_group
        return summarise_gene_group(
            self._adata, self._genes, group_col=self._group_col)


class GeneGroupTab(SidebarPage):
    """Check an annotation against known cell-type markers."""

    help_id = "scrna/marker_check"

    log_message = pyqtSignal(str)

    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self._workers = set()
        self._gen = 0
        self._loading = False
        self._last_adata_version = -1
        self._markers = dict(DEFAULT_MARKERS)
        self._owner = {}
        self._matched = {}
        self._concordance = {}
        self._picked_groups = []
        self._setup_ui()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------
    def _setup_ui(self):
        side = self.sidebar_layout

        side.addWidget(SectionHeader("MARKERS"))
        self._source_combo = QComboBox()
        self._source_combo.addItems(
            [_SOURCE_BUILTIN, _SOURCE_PANGLAO, _SOURCE_CELLMARKER,
             _SOURCE_FILE])
        self._source_combo.setToolTip(
            "Where the marker genes come from. 'From file...' reads any "
            "GMT, so your own panel works here too.")
        self._source_combo.currentIndexChanged.connect(self._on_source_changed)
        side.addWidget(self._source_combo)

        self._tissue_combo = QComboBox()
        self._tissue_combo.setVisible(False)
        self._tissue_combo.currentIndexChanged.connect(self._on_tissue_changed)
        side.addWidget(self._tissue_combo)

        row = QHBoxLayout()
        row.addWidget(SecondaryLabel("Markers per type"))
        self._per_type_spin = QSpinBox()
        self._per_type_spin.setRange(1, 20)
        self._per_type_spin.setValue(6)
        self._per_type_spin.valueChanged.connect(lambda _v: self._recompute())
        row.addWidget(self._per_type_spin)
        row.addStretch()
        side.addLayout(row)

        side.addSpacing(8)
        side.addWidget(SectionHeader("CELL TYPES"))
        side.addWidget(CaptionLabel(
            "Tick a few. The plot shows only these and the markers for "
            "them -- the whole annotation at once is unreadable."))

        self._type_list = QListWidget()
        self._type_list.itemChanged.connect(self._on_type_toggled)
        side.addWidget(self._type_list, 1)

        btns = QHBoxLayout()
        none_btn = SecondaryButton("None")
        none_btn.clicked.connect(lambda: self._set_all_types(False))
        btns.addWidget(none_btn)
        all_btn = SecondaryButton("All")
        all_btn.clicked.connect(lambda: self._set_all_types(True))
        btns.addWidget(all_btn)
        btns.addStretch()
        side.addLayout(btns)

        side.addSpacing(8)
        row = QHBoxLayout()
        row.addWidget(SecondaryLabel("Group by"))
        self._group_combo = QComboBox()
        self._group_combo.currentIndexChanged.connect(
            self._on_group_col_changed)
        row.addWidget(self._group_combo, 1)
        side.addLayout(row)

        self._concordance_label = SecondaryLabel("")
        self._concordance_label.setWordWrap(True)
        self._concordance_label.hide()
        side.addWidget(self._concordance_label)

        self._status_label = SecondaryLabel("")
        self._status_label.setWordWrap(True)
        side.addWidget(self._status_label)

        # ---- Content ----
        layout = self.content_layout

        header = QLabel("Marker Check")
        header.setFont(QFont("Segoe UI", 14, QFont.Weight.Bold))
        layout.addWidget(header)

        self._coverage_label = SecondaryLabel("")
        self._coverage_label.setWordWrap(True)
        layout.addWidget(self._coverage_label)

        self._dotplot = GeneGroupDotPlot()
        layout.addWidget(self._dotplot, 3)

        layout.addWidget(CaptionLabel(
            "One row per marker, in the same order as the plot above. "
            "'Peak group' should match 'Marks' -- where it does not, that "
            "type is not expressing what it should."))

        self._table = ResultsTable()
        self._table.set_schema([
            Column("Gene", "gene", "s"),
            Column("Marks", "marks", "s",
                   tooltip="The cell type this gene is a marker for."),
            Column("Peak group", "peak_group", "s",
                   tooltip="Which group expresses it most. For a good "
                           "annotation this matches 'Marks'."),
            Column("Mean expr", "mean_expr", ".2f"),
            Column("% expressing", "pct_expressing", ".0%"),
        ])
        self._table.setAlternatingRowColors(False)
        # Sorting off: the point of this table is that its rows line up
        # with the plot's, and sorting silently breaks that.
        self._table.setSortingEnabled(False)
        self._table.verticalHeader().setVisible(False)
        layout.addWidget(self._table, 2)

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------
    def on_tab_activated(self):
        adata = self.main_window.current_adata
        if adata is None:
            self._status_label.setText(
                "No data loaded. Load, cluster and annotate first.")
            return
        version = getattr(self.main_window, '_adata_version', 0)
        if version != self._last_adata_version:
            self._last_adata_version = version
            self._populate_group_columns(adata)
            self._populate_types(adata)
            self._recompute()
        self._refresh_concordance(adata)

    def set_adata(self, adata):
        self._last_adata_version = -1

    def reset_state(self):
        self._gen += 1
        for w in list(self._workers):
            if w.isRunning():
                w.quit()
                w.wait(2000)
        self._workers.clear()
        self._last_adata_version = -1
        self._table.setRowCount(0)
        self._type_list.clear()
        self._coverage_label.setText("")
        self._status_label.setText("")

    # ------------------------------------------------------------------
    # Cell types
    # ------------------------------------------------------------------
    def _group_values(self, adata):
        col = self._group_combo.currentText()
        if not col or col not in adata.obs.columns:
            return []
        counts = adata.obs[col].value_counts(dropna=True)
        return [str(v) for v, n in counts.items() if n > 0]

    def _populate_types(self, adata):
        """List the cell types in the data, matched to marker sets."""
        values = self._group_values(adata)
        self._matched = match_marker_types(values, self._markers)

        self._loading = True
        self._type_list.clear()
        picked = 0
        for name in values:
            marker_set = self._matched.get(name)
            label = (f"{name}  ({marker_set})" if marker_set
                     else f"{name}  -- no markers")
            item = QListWidgetItem(label)
            item.setData(Qt.ItemDataRole.UserRole, name)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            # Only types with markers are worth ticking, and only a few
            # of those: this plot is readable at about four.
            take = bool(marker_set) and picked < _DEFAULT_PICKED
            item.setCheckState(Qt.CheckState.Checked if take
                               else Qt.CheckState.Unchecked)
            picked += take
            if not marker_set:
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEnabled)
            self._type_list.addItem(item)
        self._loading = False

        n_unmatched = len(values) - len(self._matched)
        if n_unmatched:
            self._status_label.setText(
                f"{n_unmatched} of {len(values)} cell type(s) have no "
                f"markers in this source and cannot be checked. Try "
                f"another source, or supply your own via 'From file...'.")
        else:
            self._status_label.setText("")

    def _selected_types(self):
        return [self._type_list.item(i).data(Qt.ItemDataRole.UserRole)
                for i in range(self._type_list.count())
                if self._type_list.item(i).checkState()
                == Qt.CheckState.Checked]

    def _set_all_types(self, checked: bool):
        self._loading = True
        for i in range(self._type_list.count()):
            item = self._type_list.item(i)
            if item.flags() & Qt.ItemFlag.ItemIsEnabled:
                item.setCheckState(Qt.CheckState.Checked if checked
                                   else Qt.CheckState.Unchecked)
        self._loading = False
        self._recompute()

    def _on_type_toggled(self, _item):
        if not self._loading:
            self._recompute()

    def _on_group_col_changed(self):
        if self._loading:
            return
        adata = self.main_window.current_adata
        if adata is not None:
            self._populate_types(adata)
            self._recompute()

    def _populate_group_columns(self, adata):
        self._loading = True
        self._group_combo.clear()
        for col in adata.obs.columns:
            series = adata.obs[col]
            is_cat = (str(series.dtype) in ('category', 'object')
                      or hasattr(series, 'cat'))
            if not is_cat:
                continue
            if _MIN_GROUPS <= series.nunique(dropna=True) <= _MAX_GROUPS:
                self._group_combo.addItem(col)
        for pref in ('cell_type', 'cell_type_atlas', 'leiden'):
            idx = self._group_combo.findText(pref)
            if idx >= 0:
                self._group_combo.setCurrentIndex(idx)
                break
        self._loading = False

    # ------------------------------------------------------------------
    # Independent vs atlas labels
    # ------------------------------------------------------------------
    def _refresh_concordance(self, adata):
        self._concordance = compare_annotations(adata)
        rate = self._concordance.get('agreement')
        if rate is None:
            self._concordance_label.hide()
            return
        n_unlab = self._concordance['n_unlabelled']
        note = (f", {n_unlab:,} without an atlas label" if n_unlab else "")
        self._concordance_label.setText(
            f"Independent vs atlas labels: {rate * 100:.1f}% agree "
            f"({self._concordance['n_compared']:,} cells compared{note}). "
            f"Switch 'Group by' between cell_type and cell_type_atlas to "
            f"check each.")
        self._concordance_label.show()

    # ------------------------------------------------------------------
    # Marker sources
    # ------------------------------------------------------------------
    def _on_source_changed(self):
        source = self._source_combo.currentText()
        self._tissue_combo.setVisible(
            source in (_SOURCE_PANGLAO, _SOURCE_CELLMARKER))
        if source == _SOURCE_BUILTIN:
            self._set_markers(dict(DEFAULT_MARKERS))
        elif source == _SOURCE_FILE:
            self._load_from_file()
        else:
            self._populate_tissues(source)

    def _populate_tissues(self, source):
        try:
            if source == _SOURCE_PANGLAO:
                from kosmic.reference.markers.panglaodb import get_organs
                names = get_organs()
            else:
                from kosmic.reference.markers.cellmarker2 import get_tissues
                names = get_tissues()
        except Exception as exc:
            self._status_label.setText(f"Could not read {source}: {exc}")
            return
        self._loading = True
        self._tissue_combo.clear()
        for name in sorted(names):
            self._tissue_combo.addItem(name)
        idx = self._tissue_combo.findText("Heart")
        if idx >= 0:
            self._tissue_combo.setCurrentIndex(idx)
        self._loading = False
        self._on_tissue_changed()

    def _on_tissue_changed(self):
        if self._loading:
            return
        source = self._source_combo.currentText()
        tissue = self._tissue_combo.currentText()
        if not tissue:
            return
        self._status_label.setText(f"Loading {source} markers ({tissue})...")
        try:
            if source == _SOURCE_PANGLAO:
                from kosmic.reference.markers.panglaodb import (
                    get_cell_types, get_markers,
                )
                types = get_cell_types(organ=tissue)
                markers = get_markers(cell_types=types, organ=tissue)
            else:
                from kosmic.reference.markers.cellmarker2 import (
                    get_cell_types, get_markers,
                )
                types = get_cell_types(tissue=tissue)
                markers = get_markers(cell_types=types, tissue=tissue)
        except Exception as exc:
            self._status_label.setText(f"Could not read {source}: {exc}")
            return
        self._set_markers({k: v for k, v in markers.items() if v})

    def _load_from_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open marker set", "",
            "Gene sets (*.gmt *.tsv *.csv *.txt);;All files (*)")
        if not path:
            self._source_combo.setCurrentText(_SOURCE_BUILTIN)
            return
        try:
            from kosmic.reference.pathways.formats import load_gmt
            markers = load_gmt(path)
        except Exception as exc:
            dialogs.warning(self, "Could not read file", str(exc))
            self._source_combo.setCurrentText(_SOURCE_BUILTIN)
            return
        if not markers:
            dialogs.warning(self, "Nothing loaded",
                            "No gene sets were found in that file.")
            self._source_combo.setCurrentText(_SOURCE_BUILTIN)
            return
        self._set_markers(markers)

    def _set_markers(self, markers: dict):
        self._markers = markers
        adata = self.main_window.current_adata
        if adata is not None:
            self._populate_types(adata)
        self._recompute()

    # ------------------------------------------------------------------
    # Compute
    # ------------------------------------------------------------------
    def _recompute(self):
        if self._loading:
            return
        adata = self.main_window.current_adata
        if adata is None:
            return

        picked = self._selected_types()
        if not picked:
            self._table.setRowCount(0)
            self._dotplot.show_unavailable_message(
                "Tick one or more cell types on the left.")
            self._coverage_label.setText("")
            return

        # Markers are keyed by the database's names; the user picked ours.
        by_marker_name, seen = {}, set()
        for observed in picked:
            marker_set = self._matched.get(observed)
            if marker_set and marker_set not in seen:
                seen.add(marker_set)
                by_marker_name[marker_set] = self._markers[marker_set]

        genes, owner = marker_panel(
            by_marker_name, list(by_marker_name),
            max_per_type=self._per_type_spin.value())
        self._owner = owner
        if not genes:
            self._table.setRowCount(0)
            self._dotplot.show_unavailable_message(
                "No markers for the selected types.")
            return

        self._picked_groups = picked
        self._status_label.setText("Summarising markers...")
        self._gen += 1
        gen = self._gen
        worker = _SummaryWorker(adata, genes, self._group_combo.currentText())
        self._workers.add(worker)
        worker.finished.connect(lambda w=worker: self._workers.discard(w))

        def on_done(summary, g=gen):
            if g == self._gen:
                self._on_summary(summary)

        def on_fail(err, g=gen):
            if g == self._gen:
                self._on_failed(err)

        run_worker(worker, on_finished=on_done, on_failed=on_fail)

    def _on_failed(self, err):
        self._status_label.setText(f"Error: {err}")
        self.log_message.emit(f"Marker summary failed: {err}")

    def _on_summary(self, summary):
        n_present = len(summary.present)
        n_missing = len(summary.missing)
        picked = self._picked_groups

        missing_note = ""
        if summary.missing:
            preview = ", ".join(list(summary.missing)[:6])
            if n_missing > 6:
                preview += f" (+{n_missing - 6} more)"
            missing_note = f"  Not detected: {preview}."
        self._coverage_label.setText(
            f"{len(picked)} cell type(s), {n_present} of "
            f"{n_present + n_missing} markers detected." + missing_note)

        mean = summary.mean_matrix
        pct = summary.pct_matrix
        if summary.groups and not mean.empty:
            # Only the picked groups as columns: the whole annotation at
            # once is what made this unreadable.
            cols = [c for c in mean.columns if str(c) in picked]
            if cols:
                mean, pct = mean[cols], pct[cols]
            self._dotplot.set_data(mean, pct)
        else:
            self._dotplot.show_unavailable_message(
                "Annotate cell types to see the dot plot.")

        # Table rows in plot order, so the two line up.
        df = summary.per_gene.reset_index()
        gene_col = 'gene' if 'gene' in df.columns else df.columns[0]
        if 'peak_group' not in df.columns:
            df['peak_group'] = ''
        df['marks'] = df[gene_col].map(self._owner).fillna('')
        order = {g: i for i, g in enumerate(mean.index)}
        df = (df[df[gene_col].isin(order)]
              .assign(_o=lambda d: d[gene_col].map(order))
              .sort_values('_o')
              .drop(columns='_o'))
        self._table.set_data(df)

        self._status_label.setText("")
        self.main_window.mark_step_complete(5)
        self.log_message.emit(
            f"Marker check: {n_present} markers across {len(picked)} types")

    def refresh_theme(self):
        self._dotplot.refresh_theme()
