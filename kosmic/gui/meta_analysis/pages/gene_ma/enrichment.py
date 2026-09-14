# Enrichment tab: status line + Results sub-tab + Dot Plot sub-tab.
#
# The tab is a pure renderer; the parent 'GeneMAPage' runs the
# worker and calls 'EnrichmentTab.set_results' on completion.
from __future__ import annotations

import contextlib
from typing import Optional

import numpy as np
import pandas as pd
from PyQt6.QtWidgets import (
    QAbstractItemView, QTabWidget, QVBoxLayout, QWidget,
)

from kosmic import DEFAULT_FDR
from kosmic.gui.shared import borderless
from kosmic.gui.shared.widgets import (
    Column, ResultsTable, SecondaryLabel,
)
from kosmic.numerical import neg_log10


_TABLE_SCHEMA = [
    Column('Term', 'Term', 's'),
    Column('P-value', 'P_value', '.2e'),
    Column('FDR', 'FDR', '.2e'),
    Column('Fold', 'Fold_Enrichment', '.1f'),
    Column('Genes', 'Gene_Count', 'd'),
    Column('Size', 'Term_Size', 'd'),
]


class EnrichmentTab(QWidget):
    """Tabbed enrichment results: 'Results' table + matplotlib dot plot."""

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)

        # Matplotlib imports deferred until needed -- keeps import-time
        # lightweight for tests that don't render plots.
        import matplotlib
        matplotlib.use('QtAgg')
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import (
            FigureCanvasQTAgg as FigureCanvas)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)

        self.status_label = SecondaryLabel("")
        outer.addWidget(self.status_label)

        self._subtabs = QTabWidget()

        # Results table tab
        table_tab = QWidget()
        tt_lay = borderless(QVBoxLayout, table_tab)
        self.table = ResultsTable()
        self.table.set_schema(_TABLE_SCHEMA)
        self.table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSortingEnabled(True)
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setStretchLastSection(True)
        tt_lay.addWidget(self.table)
        self._subtabs.addTab(table_tab, "Results")

        # Dot plot tab
        dotplot_tab = QWidget()
        dp_lay = borderless(QVBoxLayout, dotplot_tab)
        self.fig = Figure(figsize=(7, 5), tight_layout=True)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        self._cbar = None
        dp_lay.addWidget(self.canvas)
        self._subtabs.addTab(dotplot_tab, "Dot Plot")

        outer.addWidget(self._subtabs)

    # -- public API ----------------------------------------------------

    def set_status(self, text: str) -> None:
        self.status_label.setText(text)

    def set_results(self, result_df: Optional[pd.DataFrame]) -> None:
        """
        Populate the table + dot plot from an enrichment result.

        Empty / None results render an empty table and dot plot with
        a "no terms" message; callers should also update the status
        label to explain why.
        """
        if result_df is None or len(result_df) == 0:
            self.table.set_data(pd.DataFrame(columns=[c.key for c in _TABLE_SCHEMA]))
            self._draw_dotplot(None)
            return

        self.table.set_data(result_df)
        self._draw_dotplot(result_df)

    def clear(self) -> None:
        """Reset everything to an empty state (called on consensus rerun)."""
        self.table.setRowCount(0)
        self.table.setColumnCount(0)
        self.set_status("")
        self._draw_dotplot(None)

    # -- dot plot ------------------------------------------------------

    def _draw_dotplot(self, result_df: Optional[pd.DataFrame]) -> None:
        """
        Top-20 significant terms by FDR.

        X = fold enrichment, Y = term name (categorical),
        size = gene count, colour = -log10(FDR).
        """
        ax = self.ax
        ax.clear()
        if self._cbar is not None:
            with contextlib.suppress(Exception):
                self._cbar.remove()
            self._cbar = None

        if result_df is None or len(result_df) == 0:
            self.canvas.draw_idle()
            return

        df = result_df.copy()
        if 'FDR' in df.columns:
            df = df[df['FDR'] < DEFAULT_FDR].sort_values('FDR').head(20)
        else:
            df = df.head(20)
        df = df.iloc[::-1].reset_index(drop=True)

        n = len(df)
        if n == 0:
            ax.text(0.5, 0.5, 'No significant terms (FDR < 0.05)',
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=10, color='#888')
            self.canvas.draw_idle()
            return

        fold = df['Fold_Enrichment'].values.astype(float)
        fdr = df['FDR'].values.astype(float)
        gene_count = df['Gene_Count'].values.astype(float)
        terms = df['Term'].tolist()

        neg_log_fdr = neg_log10(fdr, floor=1e-16)
        sizes = 30 + 250 * (gene_count / max(gene_count.max(), 1))

        sc = ax.scatter(fold, range(n), s=sizes, c=neg_log_fdr,
                        cmap='Reds', edgecolors='#333',
                        linewidths=0.5, zorder=3)

        ax.set_yticks(range(n))
        ax.set_yticklabels(terms, fontsize=7.5)
        ax.set_xlabel('Fold Enrichment', fontsize=9)
        ax.set_title(f'Top {n} significant terms (FDR < 0.05)',
                     fontsize=10)
        ax.grid(axis='x', alpha=0.3, zorder=0)
        ax.set_axisbelow(True)
        ax.margins(y=0.02)

        self._cbar = self.fig.colorbar(
            sc, ax=ax, label='-log10(FDR)', shrink=0.7, pad=0.02)

        # Gene-count size legend
        gc_vals = sorted(set([int(gene_count.min()),
                              int(np.median(gene_count)),
                              int(gene_count.max())]))
        gc_vals = [v for v in gc_vals if v > 0]
        handles = [ax.scatter([], [], s=30 + 250 * v / max(gene_count.max(), 1),
                              c='grey', edgecolors='#333', linewidths=0.5)
                   for v in gc_vals]
        if handles:
            ax.legend(handles, [str(v) for v in gc_vals],
                      title='Genes', loc='lower right',
                      fontsize=7, title_fontsize=8, framealpha=0.9)

        self.fig.tight_layout()
        self.canvas.draw_idle()


__all__ = ["EnrichmentTab"]
