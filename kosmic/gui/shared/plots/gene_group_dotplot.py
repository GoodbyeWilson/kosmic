# Marker dot plot widget, styled after the Figures workspace's marker dot
# plot (kosmic.visualisation.scrna.markers.create_marker_dotplot).
#
# Columns are marker genes, consecutive for each cell type they mark
# (the table beneath names which); rows are the data's cell types, placed by 'order_rows' under
# the marker set they express most, so a correct annotation reads as a
# diagonal. Dot area is proportional to the fraction of cells expressing
# the gene; colour is mean expression on the figure's grey-to-red scale,
# capped at the 99.5th percentile.
#
# Drawn with QPainter rather than pyqtgraph: the layout needs rotated
# gene and panel labels and a legend at a fixed pixel size, none of which
# a pyqtgraph view box provides. Text follows the application theme.
from __future__ import annotations

import numpy as np
from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import (
    QBrush, QColor, QFont, QFontMetrics, QLinearGradient, QPainter, QPen,
)
from PyQt6.QtWidgets import QSizePolicy, QToolTip, QWidget

from kosmic.gui.shared.theme import get_color, get_font_sizes

_MAX_CELL = 44      # px: the largest row height / column width
_DOT_FILL = 0.9     # diameter of a 100% dot as a fraction of the cell
_MIN_DOT = 2.0      # px: a gene detected in no cells is still drawn
_LEGEND_W = 150     # px
_PAD = 12           # px: outer margin and gaps
_SIZE_KEY = (20, 40, 60, 80, 100)


def _cmap():
    """The figure's grey-red-dark red scale, its grey taken from the theme:
    the figure's light grey is for a white page and glares on a dark one."""
    from matplotlib.colors import LinearSegmentedColormap

    from kosmic.visualisation.style import UP_COLOR
    return LinearSegmentedColormap.from_list(
        'marker_check', [get_color('fg_tertiary'), UP_COLOR, '#7A1F1F'])


class GeneGroupDotPlot(QWidget):
    """Marker genes (columns) against cell types (rows)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setMinimumSize(300, 200)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self._message = "No markers to show."
        self._genes: list[str] = []
        self._rows: list[str] = []
        self._owner: dict = {}
        self._mean = None   # rows x genes
        self._pct = None    # rows x genes
        self._vmax = 1.0
        self._clipped = False
        self._lut = []
        self._geom = None   # (x0, y0, cell) of the dot grid, set by paint

    # -- data ------------------------------------------------------------

    def set_data(self, mean_matrix, pct_matrix, owner):
        """Draw genes x groups mean and fraction-expressing tables.

        Parameters
        ----------
        mean_matrix, pct_matrix : pandas.DataFrame
            genes x groups, as from ``summarise_gene_group``; the fraction
            is in 0-1. Genes are drawn in ``mean_matrix``'s order.
        owner : dict
            ``{gene: marker set}``, used to order the rows.
        """
        from kosmic.visualisation.scrna.markers import order_rows

        genes = list(mean_matrix.index)
        if not genes or mean_matrix.shape[1] == 0:
            self.show_unavailable_message("No detected markers to plot.")
            return
        self._genes = genes
        self._owner = owner
        cols = list(order_rows(mean_matrix, owner))
        self._rows = [str(c) for c in cols]
        self._mean = mean_matrix.loc[genes, cols].to_numpy(float).T
        self._pct = pct_matrix.loc[genes, cols].to_numpy(float).T
        self._vmax = float(np.nanpercentile(self._mean, 99.5)) or 1.0
        self._clipped = float(np.nanmax(self._mean)) > self._vmax
        self._lut = [QColor(*(int(c * 255) for c in _cmap()(t)[:3]))
                     for t in np.linspace(0, 1, 256)]
        self.update()

    def show_unavailable_message(self, message=None):
        if message is not None:
            self._message = message
        self._mean = None
        self._genes, self._rows = [], []
        self.update()

    def _colour(self, value):
        t = 0.0 if not np.isfinite(value) else min(max(value / self._vmax, 0.0), 1.0)
        return self._lut[int(round(t * 255))]

    # -- painting ----------------------------------------------------------

    def _fonts(self):
        fs = get_font_sizes()
        tick = QFont()
        tick.setPointSizeF(fs['tick'])
        note = QFont()
        note.setPointSizeF(fs['annotation'])
        return tick, note

    def paintEvent(self, _event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        fg = QColor(get_color('fg_primary'))
        if self._mean is None:
            p.setPen(QColor(get_color('fg_tertiary')))
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, self._message)
            self._geom = None
            return

        tick, note = self._fonts()
        fm_t, fm_n = QFontMetrics(tick), QFontMetrics(note)
        n_r, n_g = len(self._rows), len(self._genes)

        left = max(fm_t.horizontalAdvance(r) for r in self._rows) + _PAD + 8
        bottom = max(fm_t.horizontalAdvance(g) for g in self._genes) + _PAD + 8
        top = _PAD + fm_n.height()
        avail_h = self.height() - top - bottom
        # The legend's size key is a row of dots drawn at the plot's dot
        # sizes, so the legend widens with the cells; one pass settles it.
        legend_w = _LEGEND_W
        for _ in range(2):
            avail_w = self.width() - left - legend_w - 2 * _PAD
            cell = max(4.0, min(_MAX_CELL, avail_w / n_g, avail_h / n_r))
            legend_w = max(_LEGEND_W, len(_SIZE_KEY) * (cell * _DOT_FILL + 6) + _PAD)
        grid_w, grid_h = cell * n_g, cell * n_r
        # Centre the plot and legend together.
        x0 = left + max(0.0, (avail_w - grid_w) / 2)
        y0 = top + max(0.0, (avail_h - grid_h) / 2)
        self._geom = (x0, y0, cell)

        # Axes lines on the left and bottom, as in the figure.
        p.setPen(QPen(QColor(get_color('fg_secondary')), 1))
        p.drawLine(QPointF(x0, y0), QPointF(x0, y0 + grid_h))
        p.drawLine(QPointF(x0, y0 + grid_h), QPointF(x0 + grid_w, y0 + grid_h))

        # Row labels.
        p.setFont(tick)
        p.setPen(fg)
        for i, name in enumerate(self._rows):
            r = QRectF(0, y0 + i * cell, x0 - 8, cell)
            p.drawText(r, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, name)

        # Gene labels, rotated, below the grid.
        for j, gene in enumerate(self._genes):
            p.save()
            p.translate(x0 + (j + 0.5) * cell, y0 + grid_h + 6)
            p.rotate(-90)
            p.drawText(QRectF(-fm_t.horizontalAdvance(gene), -cell / 2,
                              fm_t.horizontalAdvance(gene), cell),
                       Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                       gene)
            p.restore()

        # Dots.
        p.setPen(Qt.PenStyle.NoPen)
        max_d = cell * _DOT_FILL
        for i in range(n_r):
            for j in range(n_g):
                frac = self._pct[i, j]
                d = max(_MIN_DOT, max_d * np.sqrt(frac if np.isfinite(frac) else 0.0))
                p.setBrush(QBrush(self._colour(self._mean[i, j])))
                p.drawEllipse(QPointF(x0 + (j + 0.5) * cell, y0 + (i + 0.5) * cell),
                              d / 2, d / 2)

        self._paint_legend(p, x0 + grid_w + 2 * _PAD, y0, legend_w - _PAD,
                           max_d, tick, note, fg)

    def _paint_legend(self, p, lx, ly, w, max_d, tick, note, fg):
        fm_t, fm_n = QFontMetrics(tick), QFontMetrics(note)

        # Size key: a row of dots, labelled below.
        p.setFont(note)
        p.setPen(fg)
        p.drawText(QRectF(lx, ly, w, fm_n.height()), Qt.AlignmentFlag.AlignHCenter,
                   'Cells expressing (%)')
        y = ly + fm_n.height() + 6 + max_d / 2
        step = w / len(_SIZE_KEY)
        p.setFont(tick)
        for k, pct in enumerate(_SIZE_KEY):
            cx = lx + (k + 0.5) * step
            d = max(_MIN_DOT, max_d * np.sqrt(pct / 100))
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QBrush(QColor('#7A7A7A')))
            p.drawEllipse(QPointF(cx, y), d / 2, d / 2)
            p.setPen(fg)
            p.drawText(QRectF(cx - step / 2, y + max_d / 2 + 2, step, fm_t.height()),
                       Qt.AlignmentFlag.AlignHCenter, str(pct))

        # Colour bar: horizontal, titled above, labelled below.
        y += max_d / 2 + fm_t.height() + 18
        p.setFont(note)
        p.drawText(QRectF(lx, y, w, fm_n.height()), Qt.AlignmentFlag.AlignHCenter,
                   'Mean expression')
        y += fm_n.height() + 4
        bar = QRectF(lx + 4, y, w - 8, 10)
        grad = QLinearGradient(bar.left(), 0, bar.right(), 0)
        for t in np.linspace(0, 1, 9):
            grad.setColorAt(float(t), self._lut[int(round(t * 255))])
        p.fillRect(bar, QBrush(grad))
        p.setPen(QPen(fg, 0.8))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(bar)
        p.setFont(tick)
        top_label = f'{self._vmax:.1f}' + ('+' if self._clipped else '')
        p.drawText(QRectF(bar.left() - 20, bar.bottom() + 2, 40, fm_t.height()),
                   Qt.AlignmentFlag.AlignHCenter, '0')
        p.drawText(QRectF(bar.right() - 30, bar.bottom() + 2, 60, fm_t.height()),
                   Qt.AlignmentFlag.AlignHCenter, top_label)

    # -- hover -------------------------------------------------------------

    def mouseMoveEvent(self, event):
        if self._geom is None:
            return
        x0, y0, cell = self._geom
        pos = event.position()
        j = int((pos.x() - x0) // cell)
        i = int((pos.y() - y0) // cell)
        if 0 <= i < len(self._rows) and 0 <= j < len(self._genes) \
                and pos.x() >= x0 and pos.y() >= y0:
            QToolTip.showText(
                event.globalPosition().toPoint(),
                f"{self._genes[j]} in {self._rows[i]}\n"
                f"mean: {self._mean[i, j]:.2f}\n"
                f"expressing: {self._pct[i, j] * 100:.0f}%",
                self)
        else:
            QToolTip.hideText()
