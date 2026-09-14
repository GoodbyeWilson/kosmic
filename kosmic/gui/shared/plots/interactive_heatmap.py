# Interactive z-score gene heatmap (PyQtGraph ImageItem).
#
# Features: hover shows gene/sample/value, zoom/pan, condition colour bars,
# percentile-clipped colormap, manual colorbar with min/max labels, themed
# fonts.
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtGui import QFont

from kosmic.gui.shared.plots.interactive_plot import (
    InteractivePlot, make_tooltip_label,
)
from kosmic.gui.shared.theme import get_color, get_font_sizes, style_pg_plot

pg.setConfigOptions(antialias=True)


class InteractiveHeatmap(InteractivePlot):
    """Interactive z-score heatmap using 'pyqtgraph.ImageItem'."""

    TITLE = "Gene Heatmap"
    UNAVAILABLE_MESSAGE = "Click a pathway to see its z-score heatmap."

    gene_selected = pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.setMouseEnabled(x=True, y=True)

        self.tooltip_label = make_tooltip_label()
        self.addItem(self.tooltip_label)

        self._heatmap_data = None
        self._gene_names: list[str] = []
        self._sample_names: list[str] = []
        self._conditions: list[str] = []
        self._img_item = None
        self._cbar_item = None
        self._condition_bars: list = []

        self.scene().sigMouseMoved.connect(self._on_mouse_move)
        self.scene().sigMouseClicked.connect(self._on_mouse_click)

    def set_data(self, sample_gene_df, genes, pathway_name,
                 control_label, disease_label,
                 cmap='viridis', vmin_pct=5, vmax_pct=95,
                 gene_markers=None):
        """Render a z-score heatmap from sample-level gene expression.

        Parameters
        ----------
        sample_gene_df : pandas.DataFrame
            Columns: 'sample', 'condition', + gene columns.
        genes : list of str
            Gene names to display (must be columns in sample_gene_df).
        pathway_name : str
            For the title.
        control_label, disease_label : str
            Condition labels (control plotted first).
        cmap : str
            Matplotlib colormap name.
        vmin_pct, vmax_pct : float
            Percentile clipping (0-100).
        gene_markers : dict, optional
            '{gene: suffix}' appended to the row label only (e.g. ' up'/' down'
            for significant genes). The clicked/emitted gene name stays clean.
        """
        gene_markers = gene_markers or {}
        import pandas as pd

        # Remove every item except the empty-state overlay so we can
        # decide below whether to keep it visible (no genes) or hide it
        # (real heatmap rendered).
        self.clear_plot_items()

        self.tooltip_label = make_tooltip_label()
        # Draw above the ImageItem / colorbar / bars added below, otherwise the
        # tooltip renders behind the heatmap and only shows where it overhangs.
        self.tooltip_label.setZValue(1000)
        self.addItem(self.tooltip_label)

        available = [g for g in genes if g in sample_gene_df.columns]
        if len(available) < 2:
            self.show_unavailable_message(
                "Need at least 2 genes from this pathway to render a heatmap.")
            return

        # Real data: drop the empty-state overlay before rendering.
        self.hide_unavailable_message()

        df = sample_gene_df.copy()

        if not np.any(df[available].values < 0):
            for gene in available:
                m = df[gene].mean()
                s = df[gene].std()
                df[gene] = (df[gene] - m) / s if s > 0 else 0

        control_mask = df['condition'] == control_label
        disease_mask = df['condition'] == disease_label
        other_mask = ~control_mask & ~disease_mask
        df_sorted = pd.concat(
            [df[control_mask], df[disease_mask], df[other_mask]])

        self._gene_names = available
        self._sample_names = df_sorted['sample'].tolist()
        self._conditions = df_sorted['condition'].tolist()

        matrix = df_sorted[available].values
        self._heatmap_data = matrix

        vmin = np.nanpercentile(matrix, vmin_pct)
        vmax = np.nanpercentile(matrix, vmax_pct)
        if vmin == vmax:
            vmax = vmin + 1

        import matplotlib.pyplot as plt
        mpl_cmap = plt.get_cmap(cmap)
        lut = (mpl_cmap(np.linspace(0, 1, 256)) * 255).astype(np.uint8)

        self._img_item = pg.ImageItem()
        self._img_item.setImage(matrix, levels=(vmin, vmax))
        self._img_item.setLookupTable(lut)
        n_samples = len(self._sample_names)
        n_genes = len(available)
        self._img_item.setRect(0, 0, n_samples, n_genes)
        self.addItem(self._img_item)

        # Manual colorbar.
        cbar_data = np.linspace(vmin, vmax, 256).reshape(1, 256)
        self._cbar_img = pg.ImageItem()
        self._cbar_img.setImage(cbar_data, levels=(vmin, vmax))
        self._cbar_img.setLookupTable(lut)
        cbar_x = n_samples + 0.5
        cbar_w = max(1, n_samples * 0.03)
        self._cbar_img.setRect(cbar_x, 0, cbar_w, n_genes)
        self.addItem(self._cbar_img)

        cbar_label_font = QFont(
            'Arial', int(get_font_sizes()['colorbar_tick']))
        cbar_label_x = cbar_x + cbar_w + 0.1

        vmin_lbl = pg.TextItem(text=f'{vmin:.1f}',
                               color=get_color('fg_secondary'),
                               anchor=(0, 0.5))
        vmin_lbl.setFont(cbar_label_font)
        vmin_lbl.setPos(cbar_label_x, 0)
        self.addItem(vmin_lbl)

        vmax_lbl = pg.TextItem(text=f'{vmax:.1f}',
                               color=get_color('fg_secondary'),
                               anchor=(0, 0.5))
        vmax_lbl.setFont(cbar_label_font)
        vmax_lbl.setPos(cbar_label_x, n_genes)
        self.addItem(vmax_lbl)

        zscore_lbl = pg.TextItem(text='Z-score',
                                 color=get_color('fg_secondary'),
                                 anchor=(0.5, 0.5), angle=-90)
        zscore_lbl.setFont(QFont('Arial', int(get_font_sizes()['axis_label'])))
        zscore_lbl.setPos(cbar_label_x + 2.5, n_genes / 2)
        self.addItem(zscore_lbl)

        n_ctrl = int(control_mask.sum())
        n_dis = int(disease_mask.sum())
        cond_colors = {
            control_label: get_color('plot_control'),
            disease_label: get_color('plot_disease'),
        }

        bar_y = len(available)
        bar_h = 0.4
        if n_ctrl > 0:
            ctrl_bar = pg.BarGraphItem(
                x0=[0], x1=[n_ctrl], y=[bar_y + bar_h / 2], height=bar_h,
                brush=pg.mkBrush(cond_colors[control_label]),
                pen=pg.mkPen(None))
            self.addItem(ctrl_bar)
        if n_dis > 0:
            dis_bar = pg.BarGraphItem(
                x0=[n_ctrl], x1=[n_ctrl + n_dis],
                y=[bar_y + bar_h / 2], height=bar_h,
                brush=pg.mkBrush(cond_colors[disease_label]),
                pen=pg.mkPen(None))
            self.addItem(dis_bar)

        fs = get_font_sizes()
        legend_font = QFont('Arial', int(fs['annotation']))
        ctrl_lbl = pg.TextItem(text=control_label,
                               color=cond_colors[control_label],
                               anchor=(0.5, 1))
        ctrl_lbl.setFont(legend_font)
        ctrl_lbl.setPos(n_ctrl / 2, bar_y + bar_h + 0.15)
        self.addItem(ctrl_lbl)

        dis_lbl = pg.TextItem(text=disease_label,
                              color=cond_colors[disease_label],
                              anchor=(0.5, 1))
        dis_lbl.setFont(legend_font)
        dis_lbl.setPos(n_ctrl + n_dis / 2, bar_y + bar_h + 0.15)
        self.addItem(dis_lbl)

        display_name = pathway_name.replace('_', ' ')
        title_size = f'{fs["title"]:.0f}pt'
        self.setTitle(f'{display_name} Gene Expression  '
                      f'({control_label} → {disease_label})',
                      color=get_color('fg_primary'), size=title_size)

        y_axis = self.getAxis('left')
        # Clean gene names as tick labels (no markers, so nothing inflates the
        # label metrics and forces pyqtgraph to cull rows).
        gene_ticks = [(i + 0.5, g) for i, g in enumerate(available)]
        y_axis.setTicks([gene_ticks])
        y_axis.setPen(pg.mkPen(None))
        # Shrink the tick font as the gene count grows so pyqtgraph does not
        # cull overlapping labels -- otherwise many rows render unlabelled.
        tick_pt = int(get_font_sizes()['tick'])
        if n_genes > 40:
            tick_pt = max(5, tick_pt - 3)
        elif n_genes > 25:
            tick_pt = max(6, tick_pt - 2)
        y_axis.setStyle(tickLength=-1, tickTextOffset=4,
                        tickFont=QFont('Arial', tick_pt))
        y_axis.label.hide()

        # Significance markers as separate items in the gutter left of the
        # heatmap, so gene-name labels stay clean (and uncalled). 'up'/'down'
        # = significant up/down in disease.
        for i, g in enumerate(available):
            direction = gene_markers.get(g)
            if not direction:
                continue
            sym = '▲' if direction == 'up' else '▼'
            colour = (get_color('plot_disease') if direction == 'up'
                      else get_color('plot_control'))
            mark = pg.TextItem(sym, color=colour, anchor=(0.5, 0.5))
            mark.setFont(QFont('Arial', max(6, tick_pt)))
            mark.setPos(-0.35, i + 0.5)
            self.addItem(mark)

        x_axis = self.getAxis('bottom')
        sample_ticks = [(i + 0.5, str(i + 1)) for i in range(n_samples)]
        x_axis.setTicks([sample_ticks])
        x_axis.setPen(pg.mkPen(None))
        x_axis.setStyle(tickLength=-1, tickTextOffset=4)
        x_axis.label.hide()

        self.showAxis('top', False)
        self.showAxis('right', False)

        vb = self.plotItem.getViewBox()
        vb.setDefaultPadding(0)

        self.setYRange(-0.5, bar_y + bar_h + 0.8, padding=0)
        # Left margin (-0.7) leaves room for the significance markers.
        self.setXRange(-0.7, n_samples + cbar_w + 4, padding=0)
        self.plotItem.layout.setContentsMargins(0, 0, 0, 5)

    def _on_mouse_move(self, pos):
        if self._heatmap_data is None:
            return

        mouse_point = self.plotItem.vb.mapSceneToView(pos)
        x, y = mouse_point.x(), mouse_point.y()

        col = int(x)
        row = int(y)

        if (0 <= row < len(self._gene_names)
                and 0 <= col < len(self._sample_names)):
            gene = self._gene_names[row]
            sample = self._sample_names[col]
            cond = self._conditions[col]
            val = self._heatmap_data[col, row]
            self.tooltip_label.setText(
                f'{gene}\n{sample} ({cond})\nz-score: {val:.2f}')
            self.tooltip_label.setPos(x, y)
            self.tooltip_label.show()
        else:
            self.tooltip_label.hide()

    def _on_mouse_click(self, event):
        if self._heatmap_data is None:
            return

        pos = event.scenePos()
        mouse_point = self.plotItem.vb.mapSceneToView(pos)
        row = int(mouse_point.y())

        if 0 <= row < len(self._gene_names):
            self.gene_selected.emit(self._gene_names[row])

    def reset_view(self):
        self.autoRange()

    def refresh_theme(self):
        style_pg_plot(self)
