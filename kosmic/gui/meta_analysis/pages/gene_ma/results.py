# Results tab: per-gene consensus table, sig-count label, panel summary.
# A checkbox limits the table to genes whose studies disagree in
# direction ('direction_conflict', see kosmic.meta_analysis.direction),
# and a row of level checkboxes (0/k ... k/k, conflict: the studies that
# call a gene significant in the same direction, see
# kosmic.meta_analysis.pooling.vote_count) limits it to chosen levels.
#
# Parent drives via 'populate'; 'gene_selected' fires on row click.
from __future__ import annotations

from typing import Optional

import pandas as pd
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QAbstractItemView, QCheckBox, QHBoxLayout, QLabel, QVBoxLayout, QWidget,
)

from kosmic import DEFAULT_FDR
from kosmic.gui.shared.widgets import Column, ResultsTable, SecondaryLabel
from kosmic.meta_analysis.direction import count_significant_conflicts


_PANEL_HEADERS = {
    'on_olink_cvd_iii':      'Olink CVD III',
    'on_olink_explore_3072': 'Olink Explore',
}


class ResultsTab(QWidget):
    """Consensus results table + sig-count + panel-coverage summary."""

    gene_selected = pyqtSignal(str)

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._df: Optional[pd.DataFrame] = None
        self._build_ui()

    # -- UI construction -----------------------------------------------

    def _build_ui(self) -> None:
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)

        _selectable = (Qt.TextInteractionFlag.TextSelectableByMouse
                       | Qt.TextInteractionFlag.TextSelectableByKeyboard)

        self.sig_label = SecondaryLabel("")
        self.sig_label.setTextInteractionFlags(_selectable)
        lay.addWidget(self.sig_label)

        self.conflict_only_cb = QCheckBox(
            "Show only genes with opposite-direction effects")
        self.conflict_only_cb.setToolTip(
            "Genes that are significant (study FDR < 0.05) with a positive "
            "log2FC in at least one study and a negative log2FC in "
            "another. P-value methods such as Fisher's can call these "
            "genes significant because they do not use the sign.")
        self.conflict_only_cb.toggled.connect(self._apply_filter)
        self.conflict_only_cb.setEnabled(False)
        lay.addWidget(self.conflict_only_cb)

        # One checkbox per level, rebuilt by populate (k varies by run).
        self._level_row = QHBoxLayout()
        self._level_row.setContentsMargins(0, 0, 0, 0)
        self._level_row.setSpacing(10)
        lay.addLayout(self._level_row)
        self.level_checks: dict = {}
        self._levels: Optional[pd.Series] = None

        self.table = ResultsTable()
        self.table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setSortingEnabled(True)
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self._on_selection_changed)
        lay.addWidget(self.table)

        # The page's gene-search line is parked here while this tab is
        # showing, so it sits under the table exactly as it sits under
        # the volcano on the Plots tab.
        self._search_slot = QHBoxLayout()
        self._search_slot.setContentsMargins(0, 0, 0, 0)
        lay.addLayout(self._search_slot)

        # Olink translational-potential summary -- clinical-measurability
        # flag for the consensus gene list. NOT a validation: panel
        # membership says nothing about whether the protein is changed.
        self.panel_summary_label = QLabel(
            "Run a consensus analysis to see Olink panel coverage.")
        self.panel_summary_label.setWordWrap(True)
        self.panel_summary_label.setProperty("role", "padded_status")
        self.panel_summary_label.setToolTip(
            "Which of your consensus genes can be measured today on a "
            "routine commercial Olink proteomics assay. Clinical-"
            "measurability flag, NOT a validation: panel membership "
            "says nothing about whether the protein is actually changed.")
        self.panel_summary_label.setTextInteractionFlags(_selectable)
        lay.addWidget(self.panel_summary_label)

    # -- public API ----------------------------------------------------

    def populate(
        self,
        meta_df: Optional[pd.DataFrame],
        method_keys: Optional[list],
        n_studies: Optional[int] = None,
    ) -> None:
        """Build the schema + load rows from the consensus result.

        'n_studies' (k) labels the agreement levels; without it the
        largest per-gene study count is used.
        """
        if meta_df is None:
            self._df = None
            self.table.setRowCount(0)
            return

        df = meta_df

        cols = [Column('Gene', 'names', 's')]

        # Hypothesis-mode parent annotates genes with their containing
        # pathways; surface that here so users can see which pathway a
        # consensus gene came from without leaving the Results tab.
        if 'pathways' in df.columns:
            cols.append(Column('Pathway', 'pathways', 's'))

        cols.append(Column('log2FC', 'logfoldchanges', '.3f'))

        # Vote counting has no confidence interval: run on its own, its CI
        # and I² columns come from the inverse-variance fit kept for the
        # plots, so they are not shown as if they were its own.
        votes = 'votes' in df.columns
        vote_only = list(method_keys or []) == ['vote_count']
        if 'ci_lower' in df.columns and not vote_only:
            cols += [Column('CI lo', 'ci_lower', '.3f'),
                     Column('CI hi', 'ci_upper', '.3f')]
        if 'heterogeneity_i2' in df.columns and not vote_only:
            cols.append(Column('I²', 'heterogeneity_i2', '.1%'))

        # Per-method FDR (precomputed in worker).
        if method_keys:
            for key in method_keys:
                fcol = f'fdr_{key}'
                if fcol in df.columns:
                    cols.append(Column(f'FDR({key})', fcol, '.2e'))

        cols.append(Column('Consensus FDR', 'fdr', '.2e'))
        cols.append(Column('k', 'n_studies', 'd'))
        if votes:                       # vote counting was one of the methods
            cols.append(Column('Votes', 'votes', 'd',
                               tooltip='Studies with study FDR < 0.05 in the '
                                       'more frequent direction'))

        has_direction = 'direction_conflict' in df.columns
        if has_direction:
            cols += [
                Column('Up', 'n_up', 'd',
                       tooltip='Studies with study FDR < 0.05 and log2FC > 0'),
                Column('Down', 'n_down', 'd',
                       tooltip='Studies with study FDR < 0.05 and log2FC < 0'),
                Column('Conflict', 'direction_conflict', 'bool',
                       tooltip='Significant up in one study and down in another'),
            ]

        # Panel-membership cols (added by translational annotation
        # before populate is called) -- present only if shipped.
        for panel_col, header in _PANEL_HEADERS.items():
            if panel_col in df.columns:
                cols.append(Column(header, panel_col, 'bool'))

        self.table.set_schema(cols)
        self._df = df
        self.conflict_only_cb.setEnabled(has_direction)
        if not has_direction:
            self.conflict_only_cb.setChecked(False)
        self._build_levels(df, n_studies)
        self._apply_filter()

        sig_count, n_conflict = count_significant_conflicts(
            df, fdr=DEFAULT_FDR)
        text = f"{sig_count} consensus genes / {len(df)} total (FDR < 0.05)"
        if n_conflict:
            text += (f"; {n_conflict} with opposite-direction effects "
                     f"across studies")
        self.sig_label.setText(text)

    def _build_levels(self, df: pd.DataFrame, n_studies: Optional[int]) -> None:
        """One checkbox per agreement level, labelled with its gene count."""
        from kosmic.meta_analysis.pooling.vote_count import vote_levels, vote_level

        while self._level_row.count():
            item = self._level_row.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self.level_checks = {}
        self._levels = None
        if not {'n_up', 'n_down'} <= set(df.columns) or df.empty:
            return
        k = int(n_studies or df['n_studies'].max())
        self._levels = pd.Series(
            vote_level(df['n_up'], df['n_down'], k).to_numpy(), index=df.index)
        counts = vote_levels(df, k)['total']
        label = QLabel("Studies agreeing:")
        label.setToolTip("Studies that call the gene significant (study FDR "
                         "< 0.05) in the same direction; 'conflict' means "
                         "significant up in one study and down in another.")
        self._level_row.addWidget(label)
        for level, n in counts.items():
            cb = QCheckBox(f"{level} ({n:,})")
            cb.setChecked(True)
            cb.toggled.connect(self._apply_filter)
            self.level_checks[level] = cb
            self._level_row.addWidget(cb)
        self._level_row.addStretch(1)

    def _apply_filter(self) -> None:
        """Load the table, limited to direction conflicts and to the
        ticked agreement levels."""
        if self._df is None:
            return
        df = self._df
        if (self.conflict_only_cb.isChecked()
                and 'direction_conflict' in df.columns):
            df = df[df['direction_conflict'].fillna(False).astype(bool)]
        if self._levels is not None and self.level_checks:
            shown = [lv for lv, cb in self.level_checks.items() if cb.isChecked()]
            df = df[self._levels.reindex(df.index).isin(shown)]
        self.table.set_data(df)

    def clear(self) -> None:
        """
        Reset the table and sig-count label; keep the panel summary
        as a placeholder (the parent updates it via
        'set_panel_summary').
        """
        self._df = None
        self.table.setRowCount(0)
        self.table.setColumnCount(0)
        self.sig_label.setText("")
        self.panel_summary_label.setText(
            "Run a consensus analysis to see Olink panel coverage.")

    def set_panel_summary(self, text: str) -> None:
        self.panel_summary_label.setText(text)

    def take_search_widget(self, widget: QWidget) -> None:
        """Park the page's shared search line under this tab's table."""
        self._search_slot.addWidget(widget)

    def select_gene(self, gene_name: str) -> bool:
        """
        Find a gene row by name and select+scroll to it without
        re-emitting 'gene_selected'. Returns 'True' if found.
        """
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 0)
            if item and item.text() == gene_name:
                self.table.blockSignals(True)
                self.table.selectRow(row)
                self.table.scrollToItem(item)
                self.table.blockSignals(False)
                return True
        return False

    def select_matching(self, text: str) -> bool:
        """
        Select the best row for a partial gene name, as typed into the
        page's search box. Exact name wins, then a prefix, then any
        substring -- so typing "COL1" lands on COL1A1 rather than on
        whichever gene happens to sit highest in the table.
        """
        needle = (text or "").strip().upper()
        if not needle:
            return False
        names = []
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 0)
            names.append(item.text() if item else "")
        upper = [n.upper() for n in names]
        for test in (lambda n: n == needle,
                     lambda n: n.startswith(needle),
                     lambda n: needle in n):
            for row, name in enumerate(upper):
                if name and test(name):
                    return self.select_gene(names[row])
        return False

    # -- internal: emit selection signal -------------------------------

    def _on_selection_changed(self) -> None:
        rows = self.table.selectionModel().selectedRows()
        if not rows:
            return
        item = self.table.item(rows[0].row(), 0)
        if item:
            self.gene_selected.emit(item.text())


__all__ = ["ResultsTab"]
