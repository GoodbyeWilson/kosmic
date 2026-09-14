# Area-proportional Venn diagram widget. Qt-native, no matplotlib.
#
# Accepts data via 'set_data(sig_sets: dict[str, set])' where keys
# are method names and values are sets of gene names.

from __future__ import annotations

import math
from typing import Dict, Set

import numpy as np
from PyQt6.QtWidgets import QWidget, QSizePolicy
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import (
    QPainter, QPen, QColor, QBrush, QFont, QFontMetrics,
    QPainterPath,
)


from kosmic.gui.shared.theme import get_color, get_font_sizes


# Helpers


# Venn geometry -- area-proportional layout

# Edwards rotation pattern for 4-set Venn (two pairs of rotated ellipses)
_VENN4_ROTATIONS_RAD = [
    math.radians(-45), math.radians(-45),
    math.radians(45), math.radians(45),
]
_VENN4_ASPECT = 1.5   # ry / rx


def _layout_venn4(sig_sets_4):
    """
    Compute area-proportional layout for 4 ellipses.

    Optimises all 20 parameters (cx, cy, rx, ry, rotation per ellipse)
    so that the geometric area of each of the 15 exclusive regions is
    proportional to its gene count.

    Returns (geom, label_pos) in [0..1] normalised coordinates.
    """
    from scipy.optimize import minimize as sp_minimize

    keys = list(sig_sets_4.keys())[:4]
    sets = [sig_sets_4[k] for k in keys]

    # --- Target area fractions for all 15 exclusive regions ---
    region_counts = {}
    all_genes = set.union(*sets)
    total = max(len(all_genes), 1)
    for gene in all_genes:
        mask = 0
        for i in range(4):
            if gene in sets[i]:
                mask |= (1 << i)
        region_counts[mask] = region_counts.get(mask, 0) + 1

    target_fracs = np.zeros(16)
    for mask in range(1, 16):
        target_fracs[mask] = region_counts.get(mask, 0) / total

    # --- Sampling grid ---
    gn = 50
    gx = np.linspace(0.04, 0.96, gn)
    gy = np.linspace(0.04, 0.96, gn)
    GX, GY = np.meshgrid(gx, gy)
    px = GX.ravel()
    py = GY.ravel()
    n_pts = len(px)

    # Initial radii proportional to sqrt(|set|)
    raw = np.array([math.sqrt(max(len(s), 1)) for s in sets])
    init_scale = 0.18 / raw.max()
    init_rx = raw * init_scale
    init_ry = init_rx * _VENN4_ASPECT

    def _objective(params):
        # 20 params: 4 x (cx, cy, rx, ry, rot_rad)
        p = params.reshape(4, 5)

        point_masks = np.zeros(n_pts, dtype=np.int32)
        for i in range(4):
            cx, cy, rx, ry, rot = p[i]
            rx = max(abs(rx), 0.03)
            ry = max(abs(ry), 0.03)
            cr = math.cos(-rot)
            sr = math.sin(-rot)
            dx = px - cx
            dy = py - cy
            lx = dx * cr - dy * sr
            ly = dx * sr + dy * cr
            inside = (lx / rx) ** 2 + (ly / ry) ** 2 <= 1.0
            point_masks[inside] |= (1 << i)

        n_in_any = max((point_masks > 0).sum(), 1)

        err = 0.0
        for mask in range(1, 16):
            actual = (point_masks == mask).sum() / n_in_any
            target = target_fracs[mask]
            # Chi-squared: naturally weights large regions more
            if target > 0.001:
                err += (actual - target) ** 2 / target
            elif actual > 0.001:
                err += actual * 10  # penalise false regions

        # Keep centres in bounds
        for i in range(4):
            for d in range(2):
                v = p[i, d]
                if v < 0.08:
                    err += (0.08 - v) ** 2 * 10
                elif v > 0.92:
                    err += (v - 0.92) ** 2 * 10
            # Keep radii reasonable
            for d in range(2, 4):
                v = abs(p[i, d])
                if v < 0.03:
                    err += (0.03 - v) ** 2 * 10
                elif v > 0.45:
                    err += (v - 0.45) ** 2 * 10
        return err

    # --- Constructive layout: lens-area distances + MDS ---
    # Treat each ellipse as a circle with equivalent area: r = sqrt(rx*ry)
    eq_r = np.sqrt(init_rx * init_ry)

    # For each pair, compute the centre distance that gives the right
    # overlap area using the circle lens-area formula
    pair_dists = np.zeros((4, 4))
    for i in range(4):
        for j in range(i + 1, 4):
            overlap = len(sets[i] & sets[j])
            size_i = max(len(sets[i]), 1)
            size_j = max(len(sets[j]), 1)
            # Target overlap area as fraction of total ellipse area
            # Scale: overlap / |union| * total_diagram_area
            area_i = math.pi * eq_r[i] ** 2
            area_j = math.pi * eq_r[j] ** 2
            # Target lens area proportional to overlap size
            union = size_i + size_j - overlap
            target_area = overlap / max(union, 1) * min(area_i, area_j)
            d = _distance_for_overlap(eq_r[i], eq_r[j], target_area)
            pair_dists[i, j] = pair_dists[j, i] = d

    # Classical MDS on the distance matrix
    D2 = pair_dists ** 2
    H = np.eye(4) - np.ones((4, 4)) / 4
    B = -0.5 * H @ D2 @ H
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1][:2]
    coords = eigvecs[:, idx] * np.sqrt(np.maximum(eigvals[idx], 0))

    # Normalise positions so ellipses fit in [0.08, 0.92] box
    min_x = min(coords[:, 0] - np.maximum(init_rx, init_ry))
    max_x = max(coords[:, 0] + np.maximum(init_rx, init_ry))
    min_y = min(coords[:, 1] - np.maximum(init_rx, init_ry))
    max_y = max(coords[:, 1] + np.maximum(init_rx, init_ry))
    span = max(max_x - min_x, max_y - min_y, 1e-9)

    margin = 0.10
    usable = 1.0 - 2 * margin
    s = usable / span
    for i in range(4):
        coords[i, 0] = margin + (coords[i, 0] - min_x) * s + \
            (usable - (max_x - min_x) * s) / 2
        coords[i, 1] = margin + (coords[i, 1] - min_y) * s + \
            (usable - (max_y - min_y) * s) / 2
        init_rx[i] *= s
        init_ry[i] *= s

    # Initial rotations: spread to minimise nesting
    init_rots = np.linspace(-0.7, 0.7, 4)

    # Build initial x0
    x0 = np.zeros(20)
    for i in range(4):
        x0[i*5+0] = coords[i, 0]
        x0[i*5+1] = coords[i, 1]
        x0[i*5+2] = init_rx[i]
        x0[i*5+3] = init_ry[i]
        x0[i*5+4] = init_rots[i]

    # Multi-start: lens-area guess + random perturbations
    best_res = None
    best_cost = float('inf')
    rng = np.random.RandomState(42)

    for trial in range(3):
        x_try = x0.copy()
        if trial > 0:
            x_try[0::5] += rng.uniform(-0.08, 0.08, 4)
            x_try[1::5] += rng.uniform(-0.08, 0.08, 4)
            x_try[2::5] *= rng.uniform(0.8, 1.3, 4)
            x_try[3::5] *= rng.uniform(0.8, 1.3, 4)
            x_try[4::5] += rng.uniform(-0.3, 0.3, 4)

        res = sp_minimize(_objective, x_try, method='Powell',
                          options={'maxiter': 500, 'ftol': 1e-8})
        if res.fun < best_cost:
            best_cost = res.fun
            best_res = res

    opt = best_res.x.reshape(4, 5)

    # --- Build geom + labels ---
    geom = []
    for i in range(4):
        rx = max(abs(opt[i, 2]), 0.03)
        ry = max(abs(opt[i, 3]), 0.03)
        geom.append((
            float(opt[i, 0]), float(opt[i, 1]),
            float(rx), float(ry),
            math.degrees(opt[i, 4]),
        ))

    cx_mean = np.mean([g[0] for g in geom])
    cy_mean = np.mean([g[1] for g in geom])
    labels = []
    for i in range(4):
        cx, cy, rx, ry, _ = geom[i]
        dx = cx - cx_mean
        dy = cy - cy_mean
        dist = math.sqrt(dx ** 2 + dy ** 2) or 1e-9
        extent = max(rx, ry)
        lx = cx + dx / dist * (extent + 0.06)
        ly = cy + dy / dist * (extent + 0.06)
        lx = max(0.02, min(0.98, lx))
        ly = max(0.02, min(0.98, ly))
        labels.append((float(lx), float(ly)))

    return geom, labels

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


def _lens_area(r1: float, r2: float, d: float) -> float:
    """
    Area of intersection (lens) of two circles with radii r1, r2
    whose centres are distance d apart.
    """
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2
    # Standard lens area formula
    part1 = r1 ** 2 * math.acos((d ** 2 + r1 ** 2 - r2 ** 2) / (2 * d * r1))
    part2 = r2 ** 2 * math.acos((d ** 2 + r2 ** 2 - r1 ** 2) / (2 * d * r2))
    part3 = 0.5 * math.sqrt(
        (-d + r1 + r2) * (d + r1 - r2) *
        (d - r1 + r2) * (d + r1 + r2))
    return part1 + part2 - part3


def _distance_for_overlap(r1: float, r2: float,
                          target_area: float) -> float:
    """
    Binary search for the centre distance that gives the target
    lens area between two circles of radii r1, r2.
    """
    if target_area <= 0:
        return r1 + r2 + 1.0  # no overlap
    max_area = math.pi * min(r1, r2) ** 2
    if target_area >= max_area:
        return abs(r1 - r2)   # full containment

    lo, hi = abs(r1 - r2), r1 + r2
    for _ in range(60):
        mid = (lo + hi) / 2
        if _lens_area(r1, r2, mid) > target_area:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _layout_venn2(sizes, overlap_01):
    """
    Compute (cx, cy, r) for 2 area-proportional circles.

    Returns list of (cx, cy, r) in a normalised [0..1] coordinate system
    and label positions outside each circle.
    """
    # Radii proportional to sqrt(size) so area ~ size
    r0 = math.sqrt(sizes[0])
    r1 = math.sqrt(sizes[1])

    # Normalise so the larger radius is 0.30 of the box
    max_r = max(r0, r1, 1e-9)
    scale = 0.30 / max_r
    r0 *= scale
    r1 *= scale

    # Target overlap area: since circle areas are proportional to set
    # sizes (r = sqrt(size) * scale → area = pi*r² ∝ size), the
    # overlap area should be proportional to the overlap count using
    # the same scale factor.
    # area_per_count = pi * scale² (from r = sqrt(size) * scale)
    area_per_count = math.pi * scale ** 2
    target_area = overlap_01 * area_per_count
    d = _distance_for_overlap(r0, r1, target_area)

    # Centre the pair horizontally at y=0.5
    cx0 = 0.5 - d / 2
    cx1 = 0.5 + d / 2
    cy = 0.50

    circles = [(cx0, cy, r0), (cx1, cy, r1)]
    labels = [(cx0 - r0 * 0.5, cy - r0 - 0.06),
              (cx1 + r1 * 0.5, cy - r1 - 0.06)]
    return circles, labels


def _layout_venn3(sizes, overlaps_2):
    """
    Compute (cx, cy, r) for 3 area-proportional circles.

    sizes: [|A|, |B|, |C|]
    overlaps_2: [|A&B|, |A&C|, |B&C|]

    Places circles in a triangular arrangement with pairwise distances
    determined by desired overlap areas.
    """
    r = [math.sqrt(s) for s in sizes]
    max_r = max(r[0], r[1], r[2], 1e-9)
    scale = 0.27 / max_r
    r = [ri * scale for ri in r]

    # Pairwise distances: target overlap area proportional to overlap
    # count using the same area-per-count scale as the circle areas.
    area_per_count = math.pi * scale ** 2
    pairs = [(0, 1), (0, 2), (1, 2)]
    dists = []
    for idx, (i, j) in enumerate(pairs):
        target_area = overlaps_2[idx] * area_per_count
        d = _distance_for_overlap(r[i], r[j], target_area)
        dists.append(d)

    d01, d02, d12 = dists

    # Place circle 0 at origin, circle 1 to the right
    cx0, cy0 = 0.0, 0.0
    cx1, cy1 = d01, 0.0

    # Circle 2 by triangulation from circles 0 and 1
    if d01 > 1e-9:
        cos_a = (d02 ** 2 + d01 ** 2 - d12 ** 2) / (2 * d02 * d01)
        cos_a = max(-1.0, min(1.0, cos_a))
        sin_a = math.sqrt(1 - cos_a ** 2)
        cx2 = d02 * cos_a
        cy2 = d02 * sin_a
    else:
        cx2, cy2 = 0.0, d02

    # Centre and normalise into [0..1] box
    all_cx = [cx0, cx1, cx2]
    all_cy = [cy0, cy1, cy2]
    min_x = min(c - ri for c, ri in zip(all_cx, r))
    max_x = max(c + ri for c, ri in zip(all_cx, r))
    min_y = min(c - ri for c, ri in zip(all_cy, r))
    max_y = max(c + ri for c, ri in zip(all_cy, r))

    span = max(max_x - min_x, max_y - min_y, 1e-9)
    margin = 0.10
    usable = 1.0 - 2 * margin
    s = usable / span

    circles = []
    for i in range(3):
        nx = margin + (all_cx[i] - min_x) * s + (usable - (max_x - min_x) * s) / 2
        ny = margin + (all_cy[i] - min_y) * s + (usable - (max_y - min_y) * s) / 2
        circles.append((nx, ny, r[i] * s))

    # Labels: outside each circle, away from the centroid
    centroid_x = sum(c[0] for c in circles) / 3
    centroid_y = sum(c[1] for c in circles) / 3
    labels = []
    for cx, cy, ri in circles:
        dx = cx - centroid_x
        dy = cy - centroid_y
        dist = math.sqrt(dx ** 2 + dy ** 2) or 1e-9
        lx = cx + dx / dist * (ri + 0.06)
        ly = cy + dy / dist * (ri + 0.06)
        labels.append((lx, ly))

    return circles, labels


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


