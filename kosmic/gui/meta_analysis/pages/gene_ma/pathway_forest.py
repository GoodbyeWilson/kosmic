# Pathway Forest tab (hypothesis mode): one panel per sub-pathway,
# filtered by a "Significant / All / pick one" combo. Pure renderer --
# parent passes inputs to 'PathwayForestTab.draw'.
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import pyqtgraph as pg
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QComboBox, QFrame, QHBoxLayout, QLabel, QScrollArea,
    QVBoxLayout, QWidget,
)

from kosmic.gui.shared import borderless
from kosmic.gui.shared.plots import InteractivePlot
from kosmic.gui.shared.theme import get_color, style_pg_plot


class PathwayForestTab(QWidget):
    """Per-pathway forest plots, filtered by a top-of-tab combo."""

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._panels: dict = {}            # pw_name -> InteractivePlot
        self._sig_pathways: set = set()
        self._build_ui()

    # -- UI construction -----------------------------------------------

    def _build_ui(self) -> None:
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)

        top = QHBoxLayout()
        top.addWidget(QLabel("Pathway:"))
        self.combo = QComboBox()
        self.combo.addItem("All")
        self.combo.currentIndexChanged.connect(self._on_combo_changed)
        top.addWidget(self.combo, 1)
        top.addStretch()
        outer.addLayout(top)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        self._scroll_content = QWidget()
        self._scroll_layout = borderless(QVBoxLayout, self._scroll_content)
        self._scroll_layout.setSpacing(8)
        self._scroll_layout.addStretch()
        scroll.setWidget(self._scroll_content)
        outer.addWidget(scroll)

    # -- public API ----------------------------------------------------

    def clear(self) -> None:
        """Drop all panels + reset the combo."""
        layout = self._scroll_layout
        while layout.count() > 1:  # keep the trailing stretch
            item = layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        self._panels.clear()
        self._sig_pathways = set()
        self.combo.blockSignals(True)
        self.combo.clear()
        self.combo.addItem("All")
        self.combo.blockSignals(False)

    def draw(
        self,
        meta_df: Optional[pd.DataFrame],
        pathways: Optional[dict],
        pathway_ma_results: Optional[pd.DataFrame],
    ) -> None:
        """
        Render one panel per pathway in 'pathways'.

        Parameters
        ----------
        meta_df
            Consensus meta-analysis result -- must contain 'logfoldchanges',
            'ci_lower', 'ci_upper', 'fdr', 'names' columns.
        pathways
            '{pathway_name: [gene_symbols]}' -- the hypothesis-mode
            sub-pathways selected by the user.
        pathway_ma_results
            Optional pathway-level meta-analysis result; used only to
            mark which pathways are significant in the combo filter.
        """
        self.clear()

        if meta_df is None or len(meta_df) == 0:
            return
        df = meta_df.copy()
        for col in ('logfoldchanges', 'ci_lower', 'ci_upper', 'fdr', 'names'):
            if col not in df.columns:
                return

        if not pathways:
            return

        fg = get_color('fg_primary')
        gene_lookup = df.set_index('names')

        # Shared x-axis range across all panels.
        all_ci = df[['ci_lower', 'ci_upper']].values.flatten()
        all_ci = all_ci[np.isfinite(all_ci)]
        if len(all_ci) == 0:
            return
        extreme = max(abs(all_ci.min()), abs(all_ci.max()), 0.1)
        pad = extreme * 0.05
        x_range = (-(extreme + pad), extreme + pad)

        # Significant pathways (from pathway-level MA, if any).
        sig_pathways: set = set()
        if pathway_ma_results is not None and len(pathway_ma_results) > 0:
            pw_ma = pathway_ma_results
            name_col = 'names' if 'names' in pw_ma.columns else pw_ma.columns[0]
            fdr_col = next((c for c in ('fdr', 'FDR', 'pvals_adj')
                            if c in pw_ma.columns), None)
            if fdr_col:
                sig_mask = pw_ma[fdr_col].astype(float) < 0.05
                sig_pathways = set(pw_ma.loc[sig_mask, name_col].values)

        # Build pathway groups: only pathways with any of their genes
        # present in the meta-analysis result get a panel.
        groups = []
        for pw_name in sorted(pathways.keys()):
            pw_genes = []
            for g in pathways[pw_name]:
                key = g.upper() if g.upper() in gene_lookup.index else g
                if key in gene_lookup.index:
                    r = gene_lookup.loc[key]
                    if isinstance(r, pd.DataFrame):
                        r = r.iloc[0]
                    pw_genes.append((g, r))
            if pw_genes:
                pw_genes.sort(key=lambda x: x[1]['logfoldchanges'])
                groups.append((pw_name, pw_genes))

        # One panel per group.
        for pw_name, gene_rows in groups:
            n = len(gene_rows)
            height = max(25 * n + 60, 120)

            panel = InteractivePlot(
                title=pw_name,
                bottom_label='Log2 Fold Change',
                left_label='Gene',
                unavailable_message='',
            )
            panel.LEFT_AXIS_WIDTH = 100
            panel.getAxis('left').setWidth(100)
            panel.setMouseEnabled(x=True, y=False)
            panel.setFixedHeight(height)
            panel.hide_unavailable_message()

            n_sig = 0
            for i, (gene_name, r) in enumerate(gene_rows):
                lfc = float(r['logfoldchanges'])
                ci_lo = float(r['ci_lower'])
                ci_hi = float(r['ci_upper'])
                fdr_val = float(r['fdr'])

                if np.isfinite(ci_lo) and np.isfinite(ci_hi):
                    bar = pg.PlotDataItem(
                        [ci_lo, ci_hi], [float(i), float(i)],
                        pen=pg.mkPen(fg, width=1.5))
                    panel.addItem(bar)

                if np.isfinite(lfc):
                    sig = fdr_val < 0.05
                    if sig:
                        n_sig += 1
                        colour = (220, 50, 50) if lfc > 0 else (50, 50, 220)
                        brush = pg.mkBrush(*colour, 200)
                        pen_dot = pg.mkPen(*colour, width=1)
                    else:
                        brush = pg.mkBrush(0, 0, 0, 0)
                        pen_dot = pg.mkPen('#888', width=1)
                    dot = pg.ScatterPlotItem(
                        x=[lfc], y=[float(i)], size=8,
                        pen=pen_dot, brush=brush, symbol='o',
                        hoverable=True,
                        tip=lambda x, y, data, t=gene_name, f=fdr_val, lf=lfc: (
                            "{}\nlogFC: {:+.3f}\nFDR: {:.2e}".format(t, lf, f)),
                    )
                    panel.addItem(dot)

            zero = pg.InfiniteLine(
                pos=0, angle=90,
                pen=pg.mkPen(fg, width=1, style=Qt.PenStyle.DashLine))
            panel.addItem(zero)

            ticks = [(float(i), gene_rows[i][0]) for i in range(n)]
            ax = panel.getAxis('left')
            ax.setTicks([ticks])
            ax.setStyle(tickFont=pg.QtGui.QFont('Segoe UI', 7))
            panel.setYRange(-0.5, n - 0.5)
            panel.setXRange(*x_range)

            title = "{} ({} genes, {} sig)".format(
                pw_name.replace('_', ' '), n, n_sig)
            style_pg_plot(panel, title=title,
                          left_label='',
                          bottom_label='Pooled Log2 Fold Change')

            # Insert before the trailing stretch.
            self._scroll_layout.insertWidget(
                self._scroll_layout.count() - 1, panel)
            self._panels[pw_name] = panel

        # Combo: default to "Significant" if pathway-MA had any hits.
        self.combo.blockSignals(True)
        self.combo.clear()
        if sig_pathways:
            n_sig_groups = sum(1 for pw, _ in groups if pw in sig_pathways)
            self.combo.addItem(
                f"Significant ({n_sig_groups})", userData='__significant__')
        self.combo.addItem("All", userData='__all__')
        for pw_name, _ in groups:
            suffix = ('' if pw_name in sig_pathways or not sig_pathways
                      else ' (NS)')
            self.combo.addItem(
                pw_name.replace('_', ' ') + suffix, userData=pw_name)
        self.combo.blockSignals(False)
        self._sig_pathways = sig_pathways

        self._on_combo_changed()

    # -- internal: combo-driven panel filter ---------------------------

    def _on_combo_changed(self) -> None:
        key = self.combo.currentData()
        sig = self._sig_pathways
        for pw_name, panel in self._panels.items():
            if key == '__all__' or key is None:
                panel.setVisible(True)
            elif key == '__significant__':
                panel.setVisible(pw_name in sig)
            else:
                panel.setVisible(pw_name == key)


__all__ = ["PathwayForestTab"]
