# Interactive two-level Top-DE-Genes plot.
#
# Pathway-level summary first; clicking a pathway drills down to the per-gene
# view inside that pathway. PyQtGraph-based.
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtGui import QColor, QFont

from kosmic.gui.shared.plots.interactive_plot import (
    InteractivePlot, make_tooltip_label,
)
from kosmic.gui.shared.theme import get_color, get_font_sizes, style_pg_plot
from kosmic import DEFAULT_FDR, DEFAULT_LFC_THRESHOLD
from kosmic.numerical import neg_log10 as _neg_log10

pg.setConfigOptions(antialias=True)


def _semi_transparent(theme_key: str, alpha: int = 200) -> QColor:
    """Theme color with the given alpha, ready for 'pg.mkBrush'."""
    c = QColor(get_color(theme_key))
    c.setAlpha(alpha)
    return c


class InteractiveTopDEPlot(InteractivePlot):
    """Two-level Top-DE-Genes plot: pathway summary -> per-pathway drill-down."""

    TITLE = "Top DE Genes by Pathway"
    UNAVAILABLE_MESSAGE = "Run differential expression to see top DE genes."

    pathway_selected = pyqtSignal(str)
    gene_selected = pyqtSignal(str)
    level_changed = pyqtSignal(str)  # "pathway" or "gene"

    def __init__(self):
        super().__init__()
        self.setMouseEnabled(x=True, y=True)

        self.current_level = "pathway"
        self.current_pathway = None
        self.de_results = None
        self.pathway_gene_sets = None
        self.top_n = 15
        self.bar_data = None

        self.tooltip_label = make_tooltip_label()
        self.addItem(self.tooltip_label)

        self.scene().sigMouseMoved.connect(self.on_mouse_move)
        self.scene().sigMouseClicked.connect(self.on_mouse_click)

    def set_data(self, de_results, pathway_gene_sets, top_n=15):
        """Store data and show pathway-level view."""
        self.de_results = de_results
        self.pathway_gene_sets = pathway_gene_sets
        self.top_n = top_n
        self.show_pathway_level()

    def _get_sig_genes(self):
        df = self.de_results
        return df[(df['pvals_adj'] < DEFAULT_FDR)
                  & (df['abs_logfoldchange'] > DEFAULT_LFC_THRESHOLD)].copy()

    def show_pathway_level(self):
        """Draw split horizontal bars showing sig up/down counts per pathway."""
        self.clear_plot_items()
        self.hide_unavailable_message()
        self.current_level = "pathway"
        self.current_pathway = None
        self.level_changed.emit("pathway")

        self.tooltip_label = make_tooltip_label()
        self.addItem(self.tooltip_label)

        if self.de_results is None or self.pathway_gene_sets is None:
            return

        sig = self._get_sig_genes()
        if len(sig) == 0:
            self.setTitle('No significant DE genes',
                          color=get_color('fg_primary'),
                          size=f'{get_font_sizes()["title"]:.0f}pt')
            return

        pathway_stats = []
        for pathway_name, genes in self.pathway_gene_sets.items():
            genes_upper = {g.upper() for g in genes}
            pw_genes = sig[sig['names'].str.upper().isin(genes_upper)]
            sig_up = int((pw_genes['logfoldchanges'] > 0).sum())
            sig_down = int((pw_genes['logfoldchanges'] < 0).sum())
            if sig_up + sig_down > 0:
                pathway_stats.append({
                    'pathway': pathway_name,
                    'display': pathway_name.replace('_', ' '),
                    'sig_up': sig_up,
                    'sig_down': sig_down,
                    'total_sig': sig_up + sig_down,
                })

        if not pathway_stats:
            self.setTitle('No pathways with significant genes',
                          color=get_color('fg_primary'),
                          size=f'{get_font_sizes()["title"]:.0f}pt')
            return

        import pandas as pd
        pw_df = (pd.DataFrame(pathway_stats)
                   .sort_values('total_sig', ascending=False)
                   .head(self.top_n)
                   .sort_values('total_sig', ascending=True)
                   .reset_index(drop=True))

        y_positions = np.arange(len(pw_df))
        self.bar_data = {
            'y': y_positions,
            'names': pw_df['pathway'].tolist(),
            'display': pw_df['display'].tolist(),
            'sig_up': pw_df['sig_up'].values,
            'sig_down': pw_df['sig_down'].values,
        }

        up_bar = pg.BarGraphItem(
            x0=np.zeros(len(pw_df)), x1=pw_df['sig_up'].values.astype(float),
            y=y_positions, height=0.6,
            brush=_semi_transparent('plot_upregulated'),
            pen=pg.mkPen('k', width=1),
        )
        self.addItem(up_bar)

        down_bar = pg.BarGraphItem(
            x0=-pw_df['sig_down'].values.astype(float),
            x1=np.zeros(len(pw_df)),
            y=y_positions, height=0.6,
            brush=_semi_transparent('plot_downregulated'),
            pen=pg.mkPen('k', width=1),
        )
        self.addItem(down_bar)

        y_axis = self.getAxis('left')
        ticks = [(int(i), name) for i, name in enumerate(pw_df['display'].values)]
        y_axis.setTicks([ticks])
        y_axis.setWidth(None)

        zero_line = pg.InfiniteLine(
            pos=0, angle=90, pen=pg.mkPen('k', width=0.8))
        self.addItem(zero_line)

        self.setYRange(-0.5, len(pw_df) - 0.5)
        self.autoRange()

        style_pg_plot(
            self, title='Top DE Genes by Pathway (click to view genes)',
            bottom_label='Upregulated →          ← Downregulated')

    def show_gene_level(self, pathway_name):
        """Draw horizontal bars for top N genes in the selected pathway."""
        self.clear_plot_items()
        self.hide_unavailable_message()
        self.current_level = "gene"
        self.current_pathway = pathway_name
        self.level_changed.emit("gene")

        self.tooltip_label = make_tooltip_label()
        self.addItem(self.tooltip_label)

        if self.de_results is None or self.pathway_gene_sets is None:
            return

        genes = self.pathway_gene_sets.get(pathway_name, [])
        genes_upper = {g.upper() for g in genes}
        sig = self._get_sig_genes()
        pw_sig = sig[sig['names'].str.upper().isin(genes_upper)].copy()

        if len(pw_sig) == 0:
            self.setTitle(
                f'{pathway_name.replace("_", " ")} — No significant genes',
                color=get_color('fg_primary'),
                size=f'{get_font_sizes()["title"]:.0f}pt')
            return

        pw_sig = (pw_sig.sort_values('pvals_adj', ascending=True)
                       .head(self.top_n)
                       .sort_values('pvals_adj', ascending=False)
                       .reset_index(drop=True))

        n = len(pw_sig)
        y_positions = np.arange(n)
        neg_log10 = _neg_log10(pw_sig['pvals_adj'])
        lfc_vals = pw_sig['logfoldchanges'].values

        self.bar_data = {
            'y': y_positions,
            'names': pw_sig['names'].tolist(),
            'neg_log10': neg_log10,
            'lfc': lfc_vals,
            'pvals_adj': pw_sig['pvals_adj'].values,
        }

        for i in range(n):
            up = lfc_vals[i] > 0
            bar = pg.BarGraphItem(
                x0=[0], x1=[neg_log10[i]], y=[y_positions[i]], height=0.6,
                brush=_semi_transparent(
                    'plot_upregulated' if up else 'plot_downregulated'),
                pen=pg.mkPen('k', width=1),
            )
            self.addItem(bar)

            direction = '↑' if up else '↓'
            lfc_color = get_color(
                'plot_upregulated' if up else 'plot_downregulated')
            lfc_label = pg.TextItem(
                text=f'{direction} log2FC: {lfc_vals[i]:+.2f}',
                color=lfc_color, anchor=(0, 0.5))
            lfc_label.setFont(QFont(
                'Arial', int(get_font_sizes()['axis_label']),
                QFont.Weight.Bold))
            lfc_label.setPos(neg_log10[i] + 0.1, y_positions[i])
            self.addItem(lfc_label)

        y_axis = self.getAxis('left')
        ticks = [(int(i), pw_sig['names'].iloc[i]) for i in range(n)]
        y_axis.setTicks([ticks])
        y_axis.setWidth(None)

        self.setYRange(-0.5, n - 0.5)
        self.autoRange()

        display_name = pathway_name.replace('_', ' ')
        style_pg_plot(
            self,
            title=f'{display_name} — Top DE Genes (click gene to select)',
            bottom_label='-log10(adjusted p-value)')

    def go_back(self):
        self.show_pathway_level()

    def on_mouse_move(self, pos):
        if self.bar_data is None:
            return

        mouse_point = self.plotItem.vb.mapSceneToView(pos)
        x, y = mouse_point.x(), mouse_point.y()
        y_idx = int(round(y))

        if 0 <= y_idx < len(self.bar_data['names']):
            if self.current_level == "pathway":
                name = self.bar_data['display'][y_idx]
                up = self.bar_data['sig_up'][y_idx]
                down = self.bar_data['sig_down'][y_idx]
                self.tooltip_label.setText(
                    f'{name}\nUp: {up}  Down: {down}  Total: {up + down}')
            else:
                name = self.bar_data['names'][y_idx]
                lfc = self.bar_data['lfc'][y_idx]
                padj = self.bar_data['pvals_adj'][y_idx]
                self.tooltip_label.setText(
                    f'{name}\nlog2FC: {lfc:.3f}\nadj p: {padj:.2e}')
            self.tooltip_label.setPos(x, y)
            self.tooltip_label.show()
        else:
            self.tooltip_label.hide()

    def on_mouse_click(self, event):
        if self.bar_data is None:
            return

        pos = event.scenePos()
        mouse_point = self.plotItem.vb.mapSceneToView(pos)
        y_idx = int(round(mouse_point.y()))

        if 0 <= y_idx < len(self.bar_data['names']):
            if self.current_level == "pathway":
                pathway_name = self.bar_data['names'][y_idx]
                self.pathway_selected.emit(pathway_name)
                self.show_gene_level(pathway_name)
            else:
                gene_name = self.bar_data['names'][y_idx]
                self.gene_selected.emit(gene_name)

    def reset_view(self):
        self.autoRange()

    def refresh_theme(self):
        style_pg_plot(self)
