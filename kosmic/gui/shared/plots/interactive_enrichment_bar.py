# Interactive enrichment bar plot (pyqtgraph).
#
# Two modes sharing one widget on the DE Enrichment tab:
#   - GSEA: diverging NES bars per gene set, top down + top up, coloured by
#     direction with opacity encoding FDR significance.
#   - ORA: -log10(P) bars per term, coloured by majority gene direction, with
#     a significance line at -log10(alpha).
# Matplotlib stays the Figure Export renderer; this is the live in-app plot.
from __future__ import annotations

import numpy as np
import pandas as pd
import pyqtgraph as pg
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor

from kosmic.gui.shared.plots.interactive_plot import (
    InteractivePlot, make_tooltip_label,
)
from kosmic.gui.shared.theme import get_color, style_pg_plot
from kosmic.numerical import neg_log10 as _neg_log10
from kosmic.de.fgsea import select_top_nes

pg.setConfigOptions(antialias=True)

_SIG_ALPHA = 235
_NS_ALPHA = 80
_MAX_LABEL = 52


def _brush(theme_key: str, alpha: int) -> QColor:
    """Theme colour with the given opacity, ready for 'pg.mkBrush'."""
    c = QColor(get_color(theme_key))
    c.setAlpha(alpha)
    return c


def _trim(name: str) -> str:
    s = str(name).replace('_', ' ')
    return s if len(s) <= _MAX_LABEL else s[:_MAX_LABEL - 1] + '…'


class InteractiveEnrichmentBar(InteractivePlot):
    """Native NES / -log10(P) bar plot for the Enrichment tab."""

    TITLE = ""
    UNAVAILABLE_MESSAGE = "Run enrichment to see the bar plot."
    LEFT_AXIS_WIDTH = 0  # let long pathway labels drive the axis width

    def __init__(self):
        super().__init__()
        self.setMouseEnabled(x=True, y=True)
        self.bar_data = None
        self.tooltip_label = make_tooltip_label()
        self.addItem(self.tooltip_label)
        self.scene().sigMouseMoved.connect(self._on_mouse_move)

    # -- shared draw helpers -------------------------------------------

    def _begin(self):
        self.clear_plot_items()
        self.hide_unavailable_message()
        self.tooltip_label = make_tooltip_label()
        self.addItem(self.tooltip_label)

    def _add_group(self, mask, x0, x1, y, theme_key, alpha, outline):
        """Draw one homogeneous bar group (same colour + opacity)."""
        if not np.any(mask):
            return
        pen = pg.mkPen(get_color('fg_primary'), width=1) if outline \
            else pg.mkPen(None)
        self.addItem(pg.BarGraphItem(
            x0=x0[mask], x1=x1[mask], y=y[mask], height=0.68,
            brush=_brush(theme_key, alpha), pen=pen))

    def _set_y_labels(self, names):
        axis = self.getAxis('left')
        axis.setTicks([[(i, _trim(n)) for i, n in enumerate(names)]])
        axis.setWidth(None)

    # -- GSEA (NES) ----------------------------------------------------

    def set_gsea(self, res, alpha):
        """Diverging NES bars, adaptively selected: all sets for a small library,
        else the top gene sets by NES in each direction."""
        self._begin()
        sel = select_top_nes(res)
        if sel is None or sel.empty:
            self.show_unavailable_message("No gene sets to plot.")
            return
        sel = sel.copy()
        sel['pvals_adj'] = pd.to_numeric(sel['pvals_adj'], errors='coerce')

        n = len(sel)
        y = np.arange(n)
        nes = sel['nes'].to_numpy(dtype=float)
        padj = sel['pvals_adj'].to_numpy(dtype=float)
        sig = np.isfinite(padj) & (padj < alpha)
        is_up = nes > 0
        x0 = np.zeros(n)

        self._add_group(is_up & sig, x0, nes, y, 'plot_upregulated', _SIG_ALPHA, True)
        self._add_group(is_up & ~sig, x0, nes, y, 'plot_upregulated', _NS_ALPHA, False)
        self._add_group(~is_up & sig, x0, nes, y, 'plot_downregulated', _SIG_ALPHA, True)
        self._add_group(~is_up & ~sig, x0, nes, y, 'plot_downregulated', _NS_ALPHA, False)

        le = sel.get('leading_edge', pd.Series([''] * n)).astype(str)
        le_n = le.str.split(';').apply(
            lambda gs: len([g for g in gs if g.strip()])).to_numpy()

        self.addItem(pg.InfiniteLine(
            pos=0, angle=90, pen=pg.mkPen(get_color('fg_primary'), width=0.8)))
        self._set_y_labels(sel['names'].tolist())
        self.setYRange(-0.6, n - 0.4)
        self.autoRange()

        self.bar_data = {
            'mode': 'gsea', 'y': y, 'names': sel['names'].tolist(),
            'nes': nes, 'padj': padj, 'le_n': le_n}
        style_pg_plot(
            self, title=f"Preranked GSEA — solid bars: FDR < {alpha:g}",
            bottom_label="Normalised enrichment score (NES)")

    # -- ORA (-log10 P) ------------------------------------------------

    def set_ora(self, df, alpha, top_n=30):
        """-log10(P) bars per term, coloured by majority gene direction."""
        self._begin()
        if df is None or df.empty:
            self.show_unavailable_message("No terms to plot.")
            return

        sig = df[df['FDR'] < alpha].copy()
        sel = (sig if not sig.empty else df).head(top_n).copy()
        sel = sel.iloc[::-1].reset_index(drop=True)  # largest at top of axis

        n = len(sel)
        y = np.arange(n)
        val = np.asarray(_neg_log10(sel['P_value']), dtype=float)
        x0 = np.zeros(n)
        direction = sel.get('Direction', pd.Series([''] * n)).astype(str).tolist()

        up_mask = np.array([d.endswith(' up') for d in direction])
        down_mask = np.array([d.endswith(' down') for d in direction])
        mixed_mask = ~(up_mask | down_mask)

        self._add_group(up_mask, x0, val, y, 'plot_upregulated', _SIG_ALPHA, True)
        self._add_group(down_mask, x0, val, y, 'plot_downregulated', _SIG_ALPHA, True)
        self._add_group(mixed_mask, x0, val, y, 'plot_insignificant', _SIG_ALPHA, False)

        thresh = float(_neg_log10([alpha])[0])
        self.addItem(pg.InfiniteLine(
            pos=thresh, angle=90,
            pen=pg.mkPen(get_color('plot_insignificant'),
                         width=1, style=Qt.PenStyle.DashLine)))

        self._set_y_labels(sel['Term'].tolist())
        self.setYRange(-0.6, n - 0.4)
        self.autoRange()

        self.bar_data = {
            'mode': 'ora', 'y': y, 'names': sel['Term'].tolist(),
            'val': val, 'fdr': sel['FDR'].to_numpy(dtype=float),
            'direction': direction}
        style_pg_plot(
            self, title="Enriched terms — colour = majority gene direction",
            bottom_label="-log10(P-value)")

    # -- hover ---------------------------------------------------------

    def _on_mouse_move(self, pos):
        if self.bar_data is None:
            return
        pt = self.plotItem.vb.mapSceneToView(pos)
        i = int(round(pt.y()))
        names = self.bar_data['names']
        if not (0 <= i < len(names)):
            self.tooltip_label.hide()
            return
        d = self.bar_data
        if d['mode'] == 'gsea':
            text = (f"{names[i]}\nNES: {d['nes'][i]:+.2f}\n"
                    f"FDR: {d['padj'][i]:.2e}\nLeading edge: {int(d['le_n'][i])}")
        else:
            text = (f"{names[i]}\n-log10(P): {d['val'][i]:.2f}\n"
                    f"FDR: {d['fdr'][i]:.2e}\n{d['direction'][i]}")
        self.tooltip_label.setText(text)
        self.tooltip_label.setPos(pt.x(), pt.y())
        self.tooltip_label.show()

    def refresh_theme(self):
        style_pg_plot(self)
