# Cross-Dataset Reproducibility tab.
#
# Three sub-tabs: Pairwise Matrices / Per-Study Orphans / Per-Gene.
# Owns its own state cache; parent calls 'populate' on worker finish.
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QAbstractItemView, QLabel, QSplitter, QTableWidgetItem, QTabWidget,
    QVBoxLayout, QWidget,
)

from kosmic.gui.shared import borderless
from kosmic.gui.shared.widgets import (
    Column, CopyableTableWidget, ResultsTable, SectionHeader, SecondaryLabel,
)


_DEFAULT_STATUS = (
    "Cross-dataset reproducibility BEFORE any meta-analysis pooling. "
    "Three complementary views in the sub-tabs below. All use the "
    "testable denominator convention (genes not measured in a "
    "comparison study don't count as failures).")


class ReproducibilityTab(QWidget):
    """Pre-meta-analysis cross-dataset reproducibility view."""

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        # Last-rendered state (read by callers if drilling in further).
        self.matrices = None
        self.repro_df: Optional[pd.DataFrame] = None
        self.orphan_df: Optional[pd.DataFrame] = None

        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)

        # Headline / status text is selectable so users can copy
        # the reproducibility numbers straight into a paper draft.
        _selectable = (Qt.TextInteractionFlag.TextSelectableByMouse
                       | Qt.TextInteractionFlag.TextSelectableByKeyboard)

        self.status_label = SecondaryLabel(_DEFAULT_STATUS)
        self.status_label.setWordWrap(True)
        self.status_label.setTextInteractionFlags(_selectable)
        lay.addWidget(self.status_label)

        self.headline_label = QLabel("")
        self.headline_label.setWordWrap(True)
        self.headline_label.setProperty("role", "headline")
        self.headline_label.setTextInteractionFlags(_selectable)
        lay.addWidget(self.headline_label)

        self._inner_tabs = QTabWidget()
        self._inner_tabs.addTab(
            self._build_matrix_subtab(), "Pairwise Matrices")
        self._inner_tabs.addTab(
            self._build_orphan_subtab(), "Per-Study Orphans")
        self._inner_tabs.addTab(
            self._build_per_gene_subtab(), "Per-Gene Reproducibility")
        lay.addWidget(self._inner_tabs, stretch=1)

    # -- sub-tab builders ----------------------------------------------

    def _build_matrix_subtab(self) -> QWidget:
        # Strict + direction-only side by side. Loose is omitted here
        # (redundant between the other two); it still shows in the
        # orphan + per-gene sub-tabs.
        widget = QWidget()
        m_lay = QVBoxLayout(widget)
        m_lay.setContentsMargins(2, 2, 2, 2)
        intro = SecondaryLabel(
            "K x K pairwise agreement. Row i = significant genes from "
            "study i; column j = study they're compared against. Each "
            "cell = % of row i's DEGs that are also called in column j "
            "with the same direction, restricted to genes measured in j "
            "(testable denominator). Rows and columns do not have to "
            "match -- agreement(i -> j) can differ from agreement(j -> i).")
        intro.setWordWrap(True)
        m_lay.addWidget(intro)

        m_split = QSplitter(Qt.Orientation.Vertical)

        strict_wrap = QWidget()
        strict_lay = borderless(QVBoxLayout, strict_wrap)
        strict_lay.addWidget(SectionHeader(
            "Strict mode -- 'significant' = FDR < 0.05 AND same direction"))
        self.matrix_table_strict = CopyableTableWidget()
        self.matrix_table_strict.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectItems)
        self.matrix_table_strict.setSortingEnabled(False)
        self.matrix_table_strict.setAlternatingRowColors(True)
        strict_lay.addWidget(self.matrix_table_strict)
        m_split.addWidget(strict_wrap)

        dir_wrap = QWidget()
        dir_lay = borderless(QVBoxLayout, dir_wrap)
        dir_lay.addWidget(SectionHeader(
            "Direction-only mode -- 'significant' = any non-zero "
            "log2FC in the same direction (no p-value test)"))
        self.matrix_table_dir = CopyableTableWidget()
        self.matrix_table_dir.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectItems)
        self.matrix_table_dir.setSortingEnabled(False)
        self.matrix_table_dir.setAlternatingRowColors(True)
        dir_lay.addWidget(self.matrix_table_dir)
        m_split.addWidget(dir_wrap)

        m_split.setSizes([275, 275])
        m_lay.addWidget(m_split)
        return widget

    def _build_orphan_subtab(self) -> QWidget:
        widget = QWidget()
        lay = QVBoxLayout(widget)
        lay.setContentsMargins(2, 2, 2, 2)
        intro = SecondaryLabel(
            "One row per study. For each study's own DEG list, how "
            "reproducible are those calls?\n"
            "  * Orphans = called here, called nowhere else with the "
            "same direction.\n"
            "  * >= Half = also called in at least ceil(K/2) total "
            "studies (consensus core).\n"
            "Three modes shown side-by-side: (S) Strict / (L) Loose / "
            "(D) Direction-only. Higher orphan % = lower reproducibility "
            "for that study's DEGs.")
        intro.setWordWrap(True)
        lay.addWidget(intro)
        self.orphan_table = ResultsTable()
        self.orphan_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self.orphan_table.setSortingEnabled(True)
        self.orphan_table.setAlternatingRowColors(True)
        lay.addWidget(self.orphan_table)
        return widget

    def _build_per_gene_subtab(self) -> QWidget:
        widget = QWidget()
        lay = QVBoxLayout(widget)
        lay.setContentsMargins(2, 2, 2, 2)
        intro = SecondaryLabel(
            "One row per gene observed in any study. Shows how many "
            "studies call the gene significant AND agree on the majority "
            "direction, in each of the three modes: (S) Strict / (L) "
            "Loose / (D) Direction-only. A gene called up in 4 studies "
            "and down in 1 reports n_studies = 4, not 5 -- reproduced "
            "genes have both significance and direction concordance "
            "baked in. Sort by any 'n studies' column to find the most "
            "widely-called genes.")
        intro.setWordWrap(True)
        lay.addWidget(intro)
        self.repro_table = ResultsTable()
        self.repro_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self.repro_table.setSortingEnabled(True)
        self.repro_table.setAlternatingRowColors(True)
        lay.addWidget(self.repro_table)
        return widget

    # -- public API ----------------------------------------------------

    def set_status(self, text: str) -> None:
        self.status_label.setText(text)

    def set_headline(self, text: str) -> None:
        self.headline_label.setText(text)

    def clear(self) -> None:
        """Reset everything to the empty state (called on dataset reload)."""
        self.matrices = None
        self.repro_df = None
        self.orphan_df = None
        for table in (self.matrix_table_strict, self.matrix_table_dir,
                      self.orphan_table, self.repro_table):
            table.setRowCount(0)
            table.setColumnCount(0)
        self.headline_label.setText("")
        self.status_label.setText(_DEFAULT_STATUS)

    def populate(self, result: dict) -> None:
        """Render the whole tab from a 'ReproducibilityWorker' result."""
        matrices = result['matrices']
        pw_summaries = result['pw_summaries']
        repro_df = result['repro_df']
        orphan_df = result['orphan_df']
        summary = result['summary']

        self.matrices = matrices
        self.repro_df = repro_df
        self.orphan_df = orphan_df

        K = summary['n_studies']
        half = summary['half_threshold']

        m_strict = summary['modes'].get('strict', {})
        m_loose = summary['modes'].get('loose', {})
        m_dir = summary['modes'].get('direction', {})

        pw_s = pw_summaries.get('strict', {})
        pw_d = pw_summaries.get('direction_only', {})

        def _pw_line(label, s):
            if not s or not np.isfinite(s.get('mean', float('nan'))):
                return f"  {label}: n/a"
            return (f"  {label}: mean {s['mean']*100:.1f}% "
                    f"(range {s['min']*100:.1f}-{s['max']*100:.1f}%)")

        def _rep_line(label, m):
            n = m.get('n_genes_in_geq_half', 0)
            tot = m.get('n_unique_DEGs', 0)
            pct = m.get('pct_genes_in_geq_half', 0.0)
            return f"  {label}: {n:,} / {tot:,} ({pct:.1f}%)"

        self.set_headline(
            f"K = {K} studies.\n"
            f"Pairwise agreement (mean off-diagonal):\n"
            f"{_pw_line('Strict', pw_s)}\n"
            f"{_pw_line('Direction-only', pw_d)}\n"
            f"Genes reproduced in >= {half}/{K} studies "
            f"(concordant direction):\n"
            f"{_rep_line('Strict', m_strict)}\n"
            f"{_rep_line('Loose', m_loose)}\n"
            f"{_rep_line('Direction-only', m_dir)}")

        self._populate_matrix(matrices['strict'], self.matrix_table_strict)
        self._populate_matrix(
            matrices['direction_only'], self.matrix_table_dir)
        self._populate_orphan(orphan_df)
        self._populate_per_gene(repro_df)

    # -- internal populate helpers -------------------------------------

    @staticmethod
    def _populate_matrix(matrix: pd.DataFrame, table: CopyableTableWidget) -> None:
        """
        Fill 'table' with the K x K agreement 'matrix' as
        row-wise percentages.
        """
        K = matrix.shape[0]
        arr = matrix.to_numpy(dtype=float)
        table.setColumnCount(K)
        table.setRowCount(K)
        table.setHorizontalHeaderLabels(list(matrix.columns))
        table.setVerticalHeaderLabels(list(matrix.index))
        for i in range(K):
            for j in range(K):
                v = arr[i, j]
                if not np.isfinite(v):
                    item = QTableWidgetItem('n/a')
                    item.setData(Qt.ItemDataRole.UserRole, float('nan'))
                else:
                    item = QTableWidgetItem(f'{v:.1%}')
                    item.setData(Qt.ItemDataRole.UserRole, float(v))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                table.setItem(i, j, item)
        table.resizeColumnsToContents()

    def _populate_orphan(self, orphan_df: pd.DataFrame) -> None:
        self.orphan_table.set_schema([
            Column('Study',       'study',                   's',
                   "Study label."),
            Column('DEGs (S)',    'n_DEGs_strict',           'd',
                   "Significant DEGs in this study under STRICT mode "
                   "(pvals_adj < 0.05 + concordant with majority direction)."),
            Column('Orphans (S)', 'n_orphans_strict',        'd',
                   "Strict DEGs from this study called in NO other study "
                   "(same direction). Paper-style 'failed to reproduce' count."),
            Column('Orphan % (S)', 'orphan_pct_strict',      '.1f',
                   "100 * Orphans / DEGs for strict. Higher = lower "
                   "reproducibility for this study's DEGs."),
            Column('>=Half (S)',  'n_in_geq_half_strict',    'd',
                   "Strict DEGs from this study also called in at least "
                   "ceil(K/2) total studies -- the consensus core slice."),
            Column('DEGs (L)',    'n_DEGs_loose',            'd',
                   "Same as DEGs (S) but using LOOSE mode (raw p < 0.05)."),
            Column('Orphan % (L)', 'orphan_pct_loose',       '.1f',
                   "Orphan rate under loose mode."),
            Column('DEGs (D)',    'n_DEGs_direction',        'd',
                   "DEGs under DIRECTION-ONLY mode (any non-zero log2FC). "
                   "Tests direction agreement only, regardless of p-value."),
            Column('Orphan % (D)', 'orphan_pct_direction',   '.1f',
                   "Orphan rate under direction-only mode. The honest "
                   "agreement-ceiling number when individual studies are "
                   "underpowered to clear FDR alone."),
        ])
        self.orphan_table.set_data(orphan_df)

    def _populate_per_gene(self, repro_df: pd.DataFrame) -> None:
        self.repro_table.set_schema([
            Column('Gene',           'gene',                          's',
                   "Gene symbol."),
            Column('n studies (S)',  'n_studies_strict',              'd',
                   "Number of studies that called this gene strict-significant "
                   "AND agree with the majority direction. "
                   "Out of K total studies."),
            Column('n studies (L)',  'n_studies_loose',               'd',
                   "Same as n studies (S) but using loose mode."),
            Column('n studies (D)',  'n_studies_direction',           'd',
                   "Same as n studies (S) but using direction-only mode."),
            Column('Direction',      'dominant_direction_strict',     'd',
                   "+1 = up across the majority of studies that called it; "
                   "-1 = down; 0 = tied."),
            Column('Concordant (S)', 'direction_concordant_strict',   'bool',
                   "Did EVERY study calling this gene under strict mode agree "
                   "on direction? False means some called it up and some down."),
            Column('Mean log2FC',    'mean_log2FC',                   '.3f',
                   "Mean log2FC across all studies that measured this gene "
                   "(includes non-significant calls)."),
            Column('n measured',     'n_measured',                    'd',
                   "How many studies measured this gene (denominator for the "
                   "n studies counts is K, not this; this just records gene "
                   "presence)."),
        ])
        self.repro_table.set_data(repro_df)


__all__ = ["ReproducibilityTab"]
