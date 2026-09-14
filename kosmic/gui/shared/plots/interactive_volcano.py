# Shared interactive volcano-style scatter widget.
#
# Used by DE ('gene_de_page') and Meta ('gene_ma_page',
# 'pathway_ma_page'). Each page computes its own significance
# colouring and axis choice; this widget owns the scatter layers, hover
# tooltip, click-to-select signal, threshold lines, and optional overlay
# rings (e.g. Olink panel membership).
from __future__ import annotations

from typing import Callable, Optional, Sequence

import numpy as np
import pyqtgraph as pg
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QFont

from kosmic.gui.shared.plots.interactive_plot import InteractivePlot
from kosmic.gui.shared.theme import get_color, style_pg_plot, get_font_sizes

# Above this many distinct brush+size styles, a per-point-brush scatter is
# drawn as one item rather than split into per-style layers (per-gene
# shading would otherwise create one layer per point).
_MAX_SCATTER_LAYERS = 64


class InteractiveVolcano(InteractivePlot):
    """pyqtgraph volcano scatter with hover + click + overlays.

    The widget is deliberately low-level: callers pre-compute per-gene
    colours (usually derived from FDR + direction) and pass them in.
    The widget handles rendering, hover tooltips, click-to-emit, and
    optional overlay rings.

    Signals
    -------
    gene_selected : str
        Emitted when the user clicks within 'hover_threshold' of a
        point. Payload is the gene name.
    """

    TITLE = "Volcano"
    LEFT_LABEL = "-Log10 FDR"
    BOTTOM_LABEL = "Log2 Fold Change"
    UNAVAILABLE_MESSAGE = "Run differential expression to see the volcano."

    gene_selected = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseEnabled(x=True, y=True)

        # Data cache used by hover/click
        self._xy: Optional[tuple[np.ndarray, np.ndarray]] = None
        self._names: list = []
        self._hover_threshold: float = 0.15
        self._hover_fmt: Optional[Callable[[str, float, float], str]] = None

        # Item bookkeeping so clear_data can tear down cleanly
        self._scatter_items: list = []
        self._overlay_items: list = []
        self._threshold_items: list = []
        # group_key -> ScatterPlotItem; populated when 'groups=' is passed
        # to set_scatter, used by set_group_visible() to toggle layers.
        self._group_items: dict = {}
        # Optional highlight scatter (e.g. gene search) and its labels
        self._highlight_item = None
        self._highlight_label_items: list = []

        # Optional direction-of-effect labels (sticky to top corners)
        self._direction_label_items: list = []  # [(side, TextItem), ...]
        self._direction_label_signal_connected = False

        # Tooltip text item (always present; toggled via show/hide)
        self._tooltip = pg.TextItem(
            text='', color=get_color('fg_primary'), anchor=(0, 1))
        self._tooltip.setFont(
            QFont('Arial', int(get_font_sizes()['tick'])))
        self.addItem(self._tooltip)
        self._tooltip.hide()

        self.scene().sigMouseMoved.connect(self._on_mouse_move)
        self.scene().sigMouseClicked.connect(self._on_mouse_click)
        # Title + axis labels + empty-state overlay applied by InteractivePlot.

    # Data API

    def clear_data(self) -> None:
        """Remove scatter, overlays, threshold lines; keep widget chrome."""
        for item in (self._scatter_items
                     + self._overlay_items
                     + self._threshold_items):
            self.removeItem(item)
        self._scatter_items.clear()
        self._overlay_items.clear()
        self._threshold_items.clear()
        self._group_items.clear()
        if self._highlight_item is not None:
            self.removeItem(self._highlight_item)
            self._highlight_item = None
        self._xy = None
        self._names = []
        self._tooltip.hide()
        self.show_unavailable_message()

    def _add_scatter_layer(self, x, y, size, brush) -> None:
        """Add one hoverable single-brush ScatterPlotItem layer."""
        item = pg.ScatterPlotItem(
            x=x, y=y, size=size, brush=brush,
            pen=pg.mkPen(None), hoverable=True, hoverSize=10,
            hoverPen=pg.mkPen('k', width=2),
        )
        self.addItem(item)
        self._scatter_items.append(item)

    def set_scatter(
        self,
        x: np.ndarray,
        y: np.ndarray,
        names: Sequence[str],
        *,
        brushes: Optional[Sequence] = None,
        sizes: Optional[Sequence] = None,
        groups: Optional[np.ndarray] = None,
        x_label: str = 'Log2 Fold Change',
        y_label: str = '-Log10 FDR',
        title: Optional[str] = None,
        hover_fmt: Optional[Callable[[str, float, float], str]] = None,
    ) -> None:
        """Render the main scatter.

        Parameters
        ----------
        x, y : ndarray
            Same length; one entry per gene.
        names : sequence of str
        brushes : sequence of pg brushes, optional
            Per-gene fill colour. Required unless 'groups' is given.
        sizes : sequence of int, optional
            Per-gene point size. Defaults to 6.
        groups : ndarray of int, optional
            If given, points with the same group are drawn as one
            ScatterPlotItem so the caller can toggle whole layers later
            via 'set_group_visible'.
        x_label, y_label, title : axis + title strings
        hover_fmt : callable(name, x, y) -> str
            Custom tooltip formatter. Defaults to a 3-line label.
        """
        self.clear_data()
        if len(x) == 0:
            style_pg_plot(self, title=title or 'Volcano',
                          left_label=y_label, bottom_label=x_label)
            return

        # Real data: drop the empty-state overlay before rendering.
        self.hide_unavailable_message()

        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        self._xy = (x, y)
        self._names = list(names)

        self._hover_fmt = hover_fmt
        # Data-range-scaled hover threshold
        x_range = float(np.ptp(x)) if len(x) > 1 else 1.0
        y_range = float(np.ptp(y)) if len(y) > 1 else 1.0
        self._hover_threshold = max(x_range, y_range) * 0.02

        default_size = 6
        if sizes is None:
            sizes_arr = np.full(len(x), default_size, dtype=int)
        else:
            sizes_arr = np.asarray(sizes)

        if groups is None:
            if brushes is None:
                self._add_scatter_layer(
                    x, y, sizes_arr, pg.mkBrush(150, 150, 150, 150))
            else:
                # Batch points sharing a brush object + size into one
                # single-brush ScatterPlotItem each, so pyqtgraph uses its
                # fast path instead of drawing one symbol per point. Hover
                # and highlight are position-based ('_xy'), so the split is
                # transparent. Smaller points draw first, keeping larger
                # significant points on top of the grey cloud. Callers that
                # style every point differently (e.g. per-gene shading) fall
                # back to a single per-point-brush item.
                brush_list = list(brushes)
                buckets: dict = {}
                for i in range(len(x)):
                    buckets.setdefault(
                        (id(brush_list[i]), int(sizes_arr[i])), []).append(i)
                if len(buckets) > _MAX_SCATTER_LAYERS:
                    self._add_scatter_layer(x, y, sizes_arr, brush_list)
                else:
                    for (_bid, size), idxs in sorted(buckets.items(),
                                                     key=lambda kv: kv[0][1]):
                        idx = np.asarray(idxs)
                        self._add_scatter_layer(
                            x[idx], y[idx], size, brush_list[idx[0]])
        else:
            groups_arr = np.asarray(groups)
            for g in np.unique(groups_arr):
                mask = groups_arr == g
                if not mask.any():
                    continue
                layer_brushes = (
                    [brushes[i] for i in np.where(mask)[0]]
                    if brushes is not None else
                    pg.mkBrush(150, 150, 150, 150))
                item = pg.ScatterPlotItem(
                    x=x[mask], y=y[mask], size=sizes_arr[mask],
                    brush=layer_brushes,
                    pen=pg.mkPen(None), hoverable=True, hoverSize=10,
                    hoverPen=pg.mkPen('k', width=2),
                    data=[self._names[i] for i in np.where(mask)[0]],
                )
                self.addItem(item)
                self._scatter_items.append(item)
                # Hashable group key for set_group_visible().
                try:
                    key = int(g)
                except (TypeError, ValueError):
                    key = g
                self._group_items[key] = item

        # Tooltip must be on top of the scatter items
        self.removeItem(self._tooltip)
        self.addItem(self._tooltip)

        style_pg_plot(
            self, title=title or 'Volcano',
            left_label=y_label, bottom_label=x_label)

    def add_overlay_ring(
        self,
        mask: np.ndarray,
        *,
        color: tuple,
        size: int = 12,
        pen_width: float = 1.5,
    ) -> None:
        """Draw hollow rings at the (x, y) of every gene where 'mask' is True.

        Used for translational-panel overlays (e.g. Olink CVD III, Explore
        3072). Safe to call multiple times to add distinct panels.
        """
        if self._xy is None:
            return
        x, y = self._xy
        mask = np.asarray(mask, dtype=bool)
        if not mask.any() or len(mask) != len(x):
            return
        item = pg.ScatterPlotItem(
            x=x[mask], y=y[mask], size=size,
            brush=pg.mkBrush(None),
            pen=pg.mkPen(color, width=pen_width),
            hoverable=False,
        )
        self.addItem(item)
        self._overlay_items.append(item)

    def add_threshold_line(
        self,
        *,
        y: Optional[float] = None,
        x: Optional[float] = None,
        color=None,
        style: str = 'dash',
        width: int = 1,
    ) -> None:
        """Add a horizontal ('y=') or vertical ('x=') threshold line."""
        if y is None and x is None:
            return
        style_map = {
            'dash': Qt.PenStyle.DashLine,
            'dot': Qt.PenStyle.DotLine,
            'solid': Qt.PenStyle.SolidLine,
        }
        pen = pg.mkPen(
            color or get_color('fg_secondary'),
            style=style_map.get(style, Qt.PenStyle.DashLine),
            width=width)
        if y is not None:
            line = pg.InfiniteLine(pos=y, angle=0, pen=pen)
        else:
            line = pg.InfiniteLine(pos=x, angle=90, pen=pen)
        self.addItem(line)
        self._threshold_items.append(line)

    def set_group_visible(self, group_key, visible: bool) -> None:
        """Show or hide a scatter layer by its group key.

        Only meaningful when 'set_scatter(..., groups=...)' was used.
        Silently no-ops for unknown keys.
        """
        item = self._group_items.get(group_key)
        if item is not None:
            item.setVisible(visible)

    def highlight_genes(
        self,
        names: Sequence[str],
        *,
        color=None,
        fill_color=None,
        size: int = 14,
        pen_width: float = 2.0,
        symbol: str = 'o',
        with_labels: bool = False,
    ) -> int:
        """Draw a single highlight scatter at every matching gene.

        Case-insensitive substring match against the name list passed to
        'set_scatter'. Replaces any previous highlight. Returns the
        number of matched genes (0 when names is empty or no match).

        Parameters
        ----------
        names : sequence of str
            Substrings to match against the data point names.
        color : pyqtgraph color, optional
            Outline colour. Defaults to the theme's foreground primary.
        fill_color : pyqtgraph color, optional
            Fill colour. 'None' (default) gives a hollow ring; pass a
            colour (e.g. yellow) to get a filled marker.
        size : int
            Marker size in pixels.
        pen_width : float
            Outline thickness.
        symbol : str
            Pyqtgraph symbol code: 'o' (default circle), 'star',
            't' (triangle), 's' (square), '+', etc.
        """
        # Remove previous
        if self._highlight_item is not None:
            self.removeItem(self._highlight_item)
            self._highlight_item = None
        for label in self._highlight_label_items:
            self.removeItem(label)
        self._highlight_label_items = []

        if not names or self._xy is None:
            return 0
        needles = [n.upper() for n in names if n]
        if not needles:
            return 0

        matches = [
            i for i, n in enumerate(self._names)
            if any(needle in str(n).upper() for needle in needles)
        ]
        if not matches:
            return 0
        x, y = self._xy
        brush = pg.mkBrush(fill_color) if fill_color is not None else pg.mkBrush(None)
        item = pg.ScatterPlotItem(
            x=x[matches], y=y[matches], size=size,
            symbol=symbol,
            brush=brush,
            pen=pg.mkPen(color or get_color('fg_primary'), width=pen_width),
            hoverable=True,
        )
        # Lift the highlight above the base scatter dots so click /
        # hover events reach the star, not the underlying point. Z=10
        # is well above pyqtgraph's default (~0).
        item.setZValue(10)
        self.addItem(item, ignoreBounds=True)
        self._highlight_item = item

        if with_labels:
            # Match the fill colour when available so star + label read
            # as one unit. Bold + slightly larger so they stand out from
            # axis-tick text.
            label_color = fill_color if fill_color is not None else (
                color or get_color('fg_primary'))
            label_size = int(get_font_sizes().get('tick', 10)) + 2
            for idx in matches:
                label = pg.TextItem(
                    text=str(self._names[idx]),
                    color=label_color,
                    anchor=(0, 1),  # bottom-left of text sits at the point
                )
                font = QFont('Arial', label_size)
                font.setBold(True)
                label.setFont(font)
                label.setPos(x[idx], y[idx])
                # Lift labels above the highlight stars (Z=10) and far
                # above the base scatter (Z=0); ignoreBounds keeps the
                # text out of viewBox auto-range computation.
                label.setZValue(11)
                self.addItem(label, ignoreBounds=True)
                self._highlight_label_items.append(label)

        return len(matches)

    def clear_highlight(self) -> None:
        """Remove the current highlight overlay, if any."""
        self.highlight_genes([])

    def set_direction_labels(
        self,
        up_text: Optional[str] = None,
        down_text: Optional[str] = None,
    ) -> None:
        """Pin "up in X" / "down in X" markers to the top corners.

        Mirrors the convention forest plots use to disambiguate left
        vs right of zero. Pass 'None' for both to clear existing
        labels.

        The labels are sticky: a sigRangeChanged hook keeps them in
        the corners across zoom/pan so they never drift off-screen.
        """
        # Tear down any previous direction labels first
        for _side, item in self._direction_label_items:
            self.removeItem(item)
        self._direction_label_items = []

        if not up_text and not down_text:
            return

        fg = get_color('fg_primary')
        size = int(get_font_sizes().get('axis', 11))

        if down_text:
            left = pg.TextItem(
                text=f"← {down_text}",
                color=fg,
                anchor=(0, 0),  # top-left of text aligns to point
            )
            font = QFont('Arial', size)
            font.setBold(True)
            left.setFont(font)
            # ignoreBounds keeps the TextItem out of viewBox bounds
            # computation -- otherwise auto-range expands to fit it,
            # which moves the corner, which moves the label, which
            # expands the bounds again. Infinite zoom-out spiral.
            self.addItem(left, ignoreBounds=True)
            self._direction_label_items.append(('left', left))

        if up_text:
            right = pg.TextItem(
                text=f"{up_text} →",
                color=fg,
                anchor=(1, 0),  # top-right of text aligns to point
            )
            font = QFont('Arial', size)
            font.setBold(True)
            right.setFont(font)
            self.addItem(right, ignoreBounds=True)
            self._direction_label_items.append(('right', right))

        self._reposition_direction_labels()
        if not self._direction_label_signal_connected:
            self.getViewBox().sigRangeChanged.connect(
                self._reposition_direction_labels)
            self._direction_label_signal_connected = True

    def _reposition_direction_labels(self, *_args) -> None:
        """Stick the direction labels to the top-left / top-right corners."""
        if not self._direction_label_items:
            return
        vb = self.getViewBox()
        (xmin, xmax), (ymin, ymax) = vb.viewRange()
        pad_x = (xmax - xmin) * 0.01
        pad_y = (ymax - ymin) * 0.02
        for side, item in self._direction_label_items:
            if side == 'left':
                item.setPos(xmin + pad_x, ymax - pad_y)
            else:
                item.setPos(xmax - pad_x, ymax - pad_y)

    def set_symmetric_x(self, pad_frac: float = 0.05) -> None:
        """Symmetric x-axis range around 0, padded from data extremes."""
        if self._xy is None:
            return
        x, _ = self._xy
        if len(x) == 0:
            return
        x_abs_max = max(abs(float(x.min())), abs(float(x.max())))
        pad = max(x_abs_max * pad_frac, 0.1)
        self.setXRange(-(x_abs_max + pad), x_abs_max + pad)

    # Mouse handlers

    def _nearest_index(self, scene_pos) -> tuple[int, float]:
        """Return (nearest_index, distance) in data coords.

        Points DESeq2 couldn't test (Cook's-outlier / independent-filtered
        genes) carry a NaN p-value, so their plotted y is NaN too. Plain
        'np.argmin' on a distance array containing any NaN always returns
        that NaN's index -- regardless of where the click actually was --
        so every click/hover would silently snap to the same untested gene.
        'nanargmin' skips them; only bail out if every point is NaN.
        """
        if self._xy is None:
            return -1, float('inf')
        mouse_point = self.plotItem.vb.mapSceneToView(scene_pos)
        mx, my = mouse_point.x(), mouse_point.y()
        x, y = self._xy
        dx = x - mx
        dy = y - my
        distances = np.sqrt(dx * dx + dy * dy)
        if not np.any(np.isfinite(distances)):
            return -1, float('inf')
        idx = int(np.nanargmin(distances))
        return idx, float(distances[idx])

    def _on_mouse_move(self, pos) -> None:
        if self._xy is None:
            return
        idx, dist = self._nearest_index(pos)
        if idx < 0 or dist >= self._hover_threshold:
            self._tooltip.hide()
            return
        x, y = self._xy
        name = self._names[idx]
        gx, gy = float(x[idx]), float(y[idx])
        if self._hover_fmt is not None:
            text = self._hover_fmt(name, gx, gy)
        else:
            text = f"{name}\n{gx:.3f}\n{gy:.2f}"
        self._tooltip.setText(text)
        self._tooltip.setPos(gx, gy)
        self._tooltip.show()

    def _on_mouse_click(self, event) -> None:
        if self._xy is None:
            return
        idx, dist = self._nearest_index(event.scenePos())
        if idx < 0 or dist >= self._hover_threshold:
            return
        self.gene_selected.emit(self._names[idx])

    def refresh_theme(self) -> None:
        """Re-apply themed styling (call after a theme toggle)."""
        style_pg_plot(
            self,
            title=getattr(self.plotItem.titleLabel, 'text', 'Volcano'),
            left_label=self.plotItem.getAxis('left').labelText or '-Log10 FDR',
            bottom_label=self.plotItem.getAxis('bottom').labelText or 'Log2 Fold Change')
