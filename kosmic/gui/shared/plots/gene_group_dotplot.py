# Gene-group dot plot (PyQtGraph ScatterPlotItem).
#
# Rows are genes, columns are groups (clusters / cell types). Dot size
# encodes the fraction of cells expressing the gene; dot colour encodes
# mean expression. Mirrors InteractiveHeatmap's theming, colorbar and
# hover conventions.
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt6.QtGui import QFont

from kosmic.gui.shared.plots.interactive_plot import (
    InteractivePlot, make_tooltip_label,
)
from kosmic.gui.shared.theme import get_color, get_font_sizes, style_pg_plot

pg.setConfigOptions(antialias=True)

_MIN_DOT = 4.0   # px radius at 0% expressing
_MAX_DOT = 22.0  # px radius at 100% expressing


class GeneGroupDotPlot(InteractivePlot):
    """Genes x groups dot plot: size = percent expressing, colour = mean."""

    TITLE = "Gene Group Expression"
    UNAVAILABLE_MESSAGE = "Choose a gene group to see its expression across groups."

    def __init__(self):
        super().__init__()
        self.setMouseEnabled(x=True, y=True)
        self.tooltip_label = make_tooltip_label()
        self.addItem(self.tooltip_label)
        self._genes: list[str] = []
        self._groups: list[str] = []
        self._mean = None   # genes x groups
        self._pct = None    # genes x groups
        self._scatter = None
        self.scene().sigMouseMoved.connect(self._on_mouse_move)

    def set_data(self, mean_matrix, pct_matrix, cmap='viridis'):
        """Render from genes x groups mean / percent DataFrames.

        Parameters
        ----------
        mean_matrix, pct_matrix : pandas.DataFrame
            Same index (genes) and columns (groups). 'pct_matrix' is a
            fraction in 0-1.
        cmap : str
            Matplotlib colormap name for mean expression.
        """
        self.clear_plot_items()
        self.tooltip_label = make_tooltip_label()
        self.addItem(self.tooltip_label)

        genes = list(mean_matrix.index)
        groups = list(mean_matrix.columns)
        if len(genes) == 0 or len(groups) == 0:
            self.show_unavailable_message(
                "No detected genes in this group to plot.")
            return
        self.hide_unavailable_message()

        self._genes = genes
        self._groups = groups
        self._mean = mean_matrix.values.astype(float)
        self._pct = pct_matrix.reindex(index=genes, columns=groups).values.astype(float)

        n_genes = len(genes)
        n_groups = len(groups)

        vmin = float(np.nanmin(self._mean))
        vmax = float(np.nanmax(self._mean))
        if vmin == vmax:
            vmax = vmin + 1.0

        import matplotlib.pyplot as plt
        mpl_cmap = plt.get_cmap(cmap)

        def colour_for(val):
            t = 0.0 if vmax == vmin else (val - vmin) / (vmax - vmin)
            t = float(np.clip(t, 0, 1))
            r, g, b, _ = mpl_cmap(t)
            return pg.mkBrush(int(r * 255), int(g * 255), int(b * 255))

        spots = []
        for gi, _gene in enumerate(genes):
            y = n_genes - 1 - gi  # first gene at top
            for gj in range(n_groups):
                m = self._mean[gi, gj]
                p = self._pct[gi, gj]
                if not np.isfinite(m):
                    continue
                size = _MIN_DOT + (0.0 if not np.isfinite(p) else p) * (_MAX_DOT - _MIN_DOT)
                spots.append({
                    'pos': (gj + 0.5, y + 0.5),
                    'size': size,
                    'brush': colour_for(m),
                    'pen': pg.mkPen(get_color('fg_secondary'), width=0.5),
                })

        self._scatter = pg.ScatterPlotItem(spots=spots, pxMode=True)
        self.addItem(self._scatter)

        self._draw_colorbar(vmin, vmax, mpl_cmap, n_groups, n_genes)
        self._draw_size_legend(n_groups, n_genes)
        self._set_axes(genes, groups, n_groups, n_genes)

    def _draw_colorbar(self, vmin, vmax, mpl_cmap, n_groups, n_genes):
        lut = (mpl_cmap(np.linspace(0, 1, 256)) * 255).astype(np.uint8)
        cbar = pg.ImageItem()
        cbar.setImage(np.linspace(vmin, vmax, 256).reshape(1, 256), levels=(vmin, vmax))
        cbar.setLookupTable(lut)
        cbar_x = n_groups + 0.6
        cbar_w = max(0.4, n_groups * 0.04)
        cbar.setRect(cbar_x, 0, cbar_w, n_genes)
        self.addItem(cbar)

        fs = get_font_sizes()
        font = QFont('Arial', int(fs['colorbar_tick']))
        lx = cbar_x + cbar_w + 0.1
        for val, y in ((vmin, 0), (vmax, n_genes)):
            lbl = pg.TextItem(text=f'{val:.1f}', color=get_color('fg_secondary'),
                              anchor=(0, 0.5))
            lbl.setFont(font)
            lbl.setPos(lx, y)
            self.addItem(lbl)
        title = pg.TextItem(text='Mean expr', color=get_color('fg_secondary'),
                            anchor=(0.5, 0.5), angle=-90)
        title.setFont(QFont('Arial', int(fs['axis_label'])))
        title.setPos(lx + 2.4, n_genes / 2)
        self.addItem(title)
        self._cbar_w = cbar_w

    def _draw_size_legend(self, n_groups, n_genes):
        fs = get_font_sizes()
        font = QFont('Arial', int(fs['annotation']))
        x = n_groups + 0.9
        spots, labels = [], []
        for frac in (0.25, 0.5, 1.0):
            y = -1.0 - (1.0 - frac) * 1.2
            spots.append({
                'pos': (x, y),
                'size': _MIN_DOT + frac * (_MAX_DOT - _MIN_DOT),
                'brush': pg.mkBrush(get_color('fg_secondary')),
                'pen': pg.mkPen(None),
            })
            labels.append((x + 0.7, y, f'{int(frac * 100)}%'))
        self.addItem(pg.ScatterPlotItem(spots=spots, pxMode=True))
        for lx, ly, text in labels:
            lbl = pg.TextItem(text=text, color=get_color('fg_secondary'),
                              anchor=(0, 0.5))
            lbl.setFont(font)
            lbl.setPos(lx, ly)
            self.addItem(lbl)
        cap = pg.TextItem(text='% expressing', color=get_color('fg_secondary'),
                          anchor=(0, 0.5))
        cap.setFont(font)
        cap.setPos(x - 0.5, -0.2)
        self.addItem(cap)

    def _set_axes(self, genes, groups, n_groups, n_genes):
        y_axis = self.getAxis('left')
        y_axis.setTicks([[(n_genes - 1 - i + 0.5, g) for i, g in enumerate(genes)]])
        y_axis.setPen(pg.mkPen(None))
        y_axis.setStyle(tickLength=-1, tickTextOffset=4)
        y_axis.label.hide()

        x_axis = self.getAxis('bottom')
        x_axis.setTicks([[(j + 0.5, str(g)) for j, g in enumerate(groups)]])
        x_axis.setPen(pg.mkPen(None))
        x_axis.setStyle(tickLength=-1, tickTextOffset=4)
        x_axis.label.hide()

        self.showAxis('top', False)
        self.showAxis('right', False)
        vb = self.plotItem.getViewBox()
        vb.setDefaultPadding(0.02)
        self.setYRange(-2.6, n_genes + 0.3, padding=0)
        self.setXRange(0, n_groups + getattr(self, '_cbar_w', 1) + 4, padding=0)

    def _on_mouse_move(self, pos):
        if self._mean is None:
            return
        pt = self.plotItem.vb.mapSceneToView(pos)
        col = int(pt.x())
        row_from_top = len(self._genes) - 1 - int(pt.y())
        if 0 <= row_from_top < len(self._genes) and 0 <= col < len(self._groups):
            gene = self._genes[row_from_top]
            grp = self._groups[col]
            m = self._mean[row_from_top, col]
            p = self._pct[row_from_top, col]
            self.tooltip_label.setText(
                f'{gene} @ {grp}\nmean: {m:.2f}\n% expr: {p * 100:.0f}%')
            self.tooltip_label.setPos(pt.x(), pt.y())
            self.tooltip_label.show()
        else:
            self.tooltip_label.hide()

    def reset_view(self):
        self.autoRange()

    def refresh_theme(self):
        style_pg_plot(self)
