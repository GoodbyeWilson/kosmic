# Area-proportional Venn diagram widget. Qt-native, no matplotlib.
#
# Accepts data via 'set_data(sig_sets: dict[str, set])' where keys
# are method names and values are sets of gene names.

from __future__ import annotations

from typing import Dict, Set

import numpy as np
from PyQt6.QtWidgets import QWidget, QSizePolicy
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import (
    QPainter, QPen, QColor, QBrush, QFont, QFontMetrics,
    QPainterPath,
)


from kosmic.gui.shared.theme import get_color, get_font_sizes
# Venn geometry -- the area-proportional layout, shared with the exported
# Methods Comparison Venn.
from kosmic.visualisation.venn_layout import (  # noqa: F401
    _distance_for_overlap, _layout_venn2, _layout_venn3, _layout_venn4, _lens_area,
)

# Set colours (semi-transparent fill, opaque border)
_SET_COLORS = [
    QColor(230, 80, 80, 70),    # red
    QColor(80, 130, 230, 70),   # blue
    QColor(80, 200, 120, 70),   # green
    QColor(220, 180, 50, 70),   # gold
]

_SET_BORDER_COLORS = [
    QColor(200, 50, 50, 200),
    QColor(50, 100, 200, 200),
    QColor(50, 170, 90, 200),
    QColor(190, 150, 30, 200),
]


# Path construction helpers

def _make_ellipse_path(cx, cy, rx, ry, rot_deg):
    """Return a QPainterPath for a (possibly rotated) ellipse."""
    path = QPainterPath()
    # Build axis-aligned ellipse centred at origin
    path.addEllipse(QPointF(0, 0), rx, ry)
    # Rotate then translate
    from PyQt6.QtGui import QTransform
    xf = QTransform()
    xf.translate(cx, cy)
    xf.rotate(rot_deg)
    return xf.map(path)


def _region_centroid(path: QPainterPath) -> QPointF:
    """
    Approximate the visual centroid of a QPainterPath region.

    Samples a grid of points inside the path's bounding rect and
    returns the mean of those that are contained.  Falls back to
    bounding-rect centre if the region is too thin to hit.
    """
    br = path.boundingRect()
    if br.width() < 1 or br.height() < 1:
        return br.center()

    n_samples = 20  # per axis
    xs = np.linspace(br.left(), br.right(), n_samples)
    ys = np.linspace(br.top(), br.bottom(), n_samples)

    px, py, count = 0.0, 0.0, 0
    for x in xs:
        for y in ys:
            if path.contains(QPointF(x, y)):
                px += x
                py += y
                count += 1

    if count == 0:
        return br.center()
    return QPointF(px / count, py / count)


# VennWidget

class VennWidget(QWidget):
    """
    QPainter-based Venn diagram for 2, 3, or 4 sets.

    Uses QPainterPath boolean operations to compute the exact geometric
    region for each membership pattern, then labels at the region centroid.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._sig_sets: Dict[str, Set[str]] = {}
        self._keys: list[str] = []
        self._cached_geom = None   # (geom, label_pos) -- reused on repaint
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMinimumSize(300, 250)

    def set_data(self, sig_sets: Dict[str, Set[str]]):
        self._sig_sets = sig_sets
        self._keys = list(sig_sets.keys())
        self._cached_geom = self._compute_layout()
        self.update()

    def _compute_layout(self):
        """
        Compute (geom, label_pos) for the current data.

        Runs the area-proportional solver once; cached until set_data()
        is called again.
        """
        n = min(len(self._keys), 4)
        keys = self._keys[:n]
        ss = self._sig_sets

        if n == 1:
            # Single circle centred
            r = 0.30
            geom = [(0.5, 0.5, r, r, 0)]
            labels = [(0.5, 0.5 - r - 0.06)]
            return geom, labels
        elif n == 2:
            sizes = [len(ss[keys[0]]), len(ss[keys[1]])]
            overlap = len(ss[keys[0]] & ss[keys[1]])
            circles, labels = _layout_venn2(sizes, overlap)
            geom = [(cx, cy, r, r, 0) for cx, cy, r in circles]
            return geom, labels
        elif n == 3:
            sizes = [len(ss[k]) for k in keys]
            overlaps_2 = [
                len(ss[keys[0]] & ss[keys[1]]),
                len(ss[keys[0]] & ss[keys[2]]),
                len(ss[keys[1]] & ss[keys[2]]),
            ]
            circles, labels = _layout_venn3(sizes, overlaps_2)
            geom = [(cx, cy, r, r, 0) for cx, cy, r in circles]
            return geom, labels
        elif n >= 4:
            ss_4 = {keys[i]: ss[keys[i]] for i in range(4)}
            return _layout_venn4(ss_4)
        return None

    def paintEvent(self, event):
        if not self._keys or not self._sig_sets or not self._cached_geom:
            return

        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        n = min(len(self._keys), 4)
        keys = self._keys[:n]
        ss = self._sig_sets

        w, h = self.width(), self.height()
        # Square drawing box centred in widget
        s = min(w, h) * 0.92
        ox = (w - s) / 2
        oy = (h - s) / 2

        geom, label_pos = self._cached_geom

        # --- Build QPainterPaths ---
        set_paths = []
        for i, (cx_f, cy_f, rx_f, ry_f, rot) in enumerate(geom):
            path = _make_ellipse_path(
                ox + cx_f * s, oy + cy_f * s, rx_f * s, ry_f * s, rot)
            set_paths.append(path)

        # --- Draw filled ellipses ---
        for i, path in enumerate(set_paths):
            painter.setPen(QPen(_SET_BORDER_COLORS[i], 2.0))
            painter.setBrush(QBrush(_SET_COLORS[i]))
            painter.drawPath(path)

        # --- Compute exclusive regions via path boolean ops ---
        fg = QColor(get_color('fg_primary'))
        fs = get_font_sizes()
        count_font = QFont()
        count_font.setPointSize(
            int(max(8, fs.get('annotation', 10) - (1 if n == 4 else 0))))
        count_font.setBold(True)
        painter.setFont(count_font)
        painter.setPen(fg)

        for mask in range(1, 1 << n):
            combo = tuple(keys[i] for i in range(n) if mask & (1 << i))
            excluded = tuple(keys[i] for i in range(n)
                             if not (mask & (1 << i)))
            gene_set = set.intersection(*(ss[k] for k in combo))
            for k in excluded:
                gene_set = gene_set - ss[k]
            count = len(gene_set)
            if count == 0:
                continue

            in_indices = [i for i in range(n) if mask & (1 << i)]
            out_indices = [i for i in range(n) if not (mask & (1 << i))]

            region = set_paths[in_indices[0]]
            for idx in in_indices[1:]:
                region = region.intersected(set_paths[idx])
            for idx in out_indices:
                region = region.subtracted(set_paths[idx])

            centroid = _region_centroid(region)
            self._draw_count(painter, centroid.x(), centroid.y(), count)

        # --- Set name labels ---
        label_font = QFont()
        label_font.setPointSize(
            int(max(8, fs.get('axis_label', 11) - (1 if n == 4 else 0))))
        label_font.setBold(True)
        painter.setFont(label_font)
        for i, key in enumerate(keys):
            lx, ly = label_pos[i]
            self._draw_label(painter, ox + lx * s, oy + ly * s, key,
                             _SET_BORDER_COLORS[i])

        painter.end()

    # ------------------------------------------------------------------
    # Drawing helpers
    # ------------------------------------------------------------------
    def _draw_count(self, painter: QPainter, x: float, y: float, count: int):
        """Draw a count number centred at (x, y)."""
        text = str(count)
        fm = QFontMetrics(painter.font())
        tw = fm.horizontalAdvance(text)
        th = fm.height()
        painter.drawText(int(x - tw / 2), int(y + th / 4), text)

    def _draw_label(self, painter: QPainter, x: float, y: float,
                    text: str, color: QColor):
        """Draw a set label at (x, y) in the given colour."""
        painter.setPen(color)
        fm = QFontMetrics(painter.font())
        tw = fm.horizontalAdvance(text)
        th = fm.height()
        painter.drawText(int(x - tw / 2), int(y + th / 4), text)
        # Restore fg colour
        painter.setPen(QColor(get_color('fg_primary')))


