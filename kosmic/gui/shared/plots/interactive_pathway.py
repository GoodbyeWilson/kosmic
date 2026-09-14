# Interactive horizontal bar plot of pathway effect sizes, on PyQtGraph.
#
# Used by the pathway DE and pathway meta-analysis pages. Builds on
# 'InteractivePlot' for hover tooltips, click-to-select and export.
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt6.QtCore import pyqtSignal

from kosmic.gui.shared.plots.interactive_plot import (
    InteractivePlot, make_tooltip_label,
)
from kosmic.gui.shared.theme import get_color, style_pg_plot

pg.setConfigOptions(antialias=True)


class InteractivePathwayPlot(InteractivePlot):
    """Horizontal bar plot of pathway effect sizes with hover tooltips."""

    TITLE = "Pathway Activity"
    BOTTOM_LABEL = "Effect Size"
    UNAVAILABLE_MESSAGE = "Run pathway scoring to see pathway activity."

    pathway_selected = pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.setMouseEnabled(x=True, y=True)

        self.pathway_names = []
        self.bar_items = []

        self.tooltip_label = make_tooltip_label()
        self.addItem(self.tooltip_label)

        self.scene().sigMouseMoved.connect(self.on_mouse_move)

    def set_data(self, stats_df):
        """Set pathway data from stats DataFrame.

        Accepts either column-name convention: (Pathway, Effect_Size)
        or (names, logfoldchanges).
        """
        self.clear_plot_items()

        if stats_df is None or len(stats_df) == 0:
            self.show_unavailable_message()
            return

        if 'names' in stats_df.columns and 'Pathway' not in stats_df.columns:
            stats_df = stats_df.rename(columns={'names': 'Pathway'})
        if ('logfoldchanges' in stats_df.columns
                and 'Effect_Size' not in stats_df.columns):
            stats_df = stats_df.rename(columns={'logfoldchanges': 'Effect_Size'})

        if ('Pathway' not in stats_df.columns
                or 'Effect_Size' not in stats_df.columns):
            self.show_unavailable_message()
            return

        # Real data: drop the empty-state overlay before rendering.
        self.hide_unavailable_message()

        self.tooltip_label = make_tooltip_label()
        self.addItem(self.tooltip_label)

        stats_df = (stats_df.sort_values('Effect_Size', ascending=True)
                            .reset_index(drop=True))

        self.pathway_names = stats_df['Pathway'].tolist()
        y_positions = np.arange(len(stats_df))
        effect_sizes = stats_df['Effect_Size'].values

        def _numcol(col):
            if col not in stats_df.columns:
                return None
            try:
                return np.asarray(stats_df[col].to_numpy(), dtype=float)
            except (TypeError, ValueError):
                return None

        # Uncertainty + significance (present for the deseq2 / cpm results).
        se_vals = _numcol('se')
        p_vals = _numcol('pvals_adj')
        if p_vals is None:
            p_vals = _numcol('pvals')

        self.bar_data = {
            'y': y_positions, 'x': effect_sizes, 'names': self.pathway_names,
        }

        fg = get_color('fg_primary')
        for i, (y, es, name) in enumerate(
                zip(y_positions, effect_sizes, self.pathway_names)):
            # Colour by significance when available, else by effect magnitude.
            if p_vals is not None and np.isfinite(p_vals[i]):
                if p_vals[i] < 0.05:
                    color = (50, 50, 220, 220) if es < 0 else (220, 50, 50, 220)
                else:
                    color = (150, 150, 150, 150)          # not significant
            elif es > 0.3:
                color = (220, 50, 50, 200)
            elif es > 0.1:
                color = (255, 127, 0, 200)
            elif es > -0.1:
                color = (50, 180, 50, 200)
            else:
                color = (50, 50, 220, 200)

            bar = pg.BarGraphItem(
                x0=[0], x1=[es], y=[y], height=0.6,
                brush=color, pen=pg.mkPen('k', width=1))
            self.addItem(bar)
            self.bar_items.append(bar)

            # 95% CI error bar (horizontal) when a per-pathway SE exists.
            if se_vals is not None and np.isfinite(se_vals[i]) and se_vals[i] > 0:
                ci = 1.96 * float(se_vals[i])
                self.addItem(pg.ErrorBarItem(
                    x=np.array([float(es)]), y=np.array([float(y)]),
                    left=np.array([ci]), right=np.array([ci]),
                    beam=0.2, pen=pg.mkPen(fg, width=1.2)))

            label = pg.TextItem(
                text=name.replace('_', ' '),
                color=get_color('fg_primary'),
                anchor=(1, 0.5) if es >= 0 else (0, 0.5))
            label.setPos(0, y)
            self.addItem(label)

        zero_line = pg.InfiniteLine(
            pos=0, angle=90, pen=pg.mkPen('k', width=1))
        self.addItem(zero_line)

        self.setYRange(-0.5, len(stats_df) - 0.5)

        style_pg_plot(self, title='Pathway Activity', bottom_label='Effect Size')

    def on_mouse_move(self, pos):
        """Hover tooltip showing pathway name + effect size."""
        if not hasattr(self, 'bar_data') or self.bar_data is None:
            return

        mouse_point = self.plotItem.vb.mapSceneToView(pos)
        x, y = mouse_point.x(), mouse_point.y()

        y_idx = int(round(y))
        if 0 <= y_idx < len(self.bar_data['names']):
            name = self.bar_data['names'][y_idx]
            es = self.bar_data['x'][y_idx]
            self.tooltip_label.setText(
                f'{name.replace("_", " ")}\nEffect Size: {es:.3f}')
            self.tooltip_label.setPos(x, y)
            self.tooltip_label.show()
        else:
            self.tooltip_label.hide()

    def reset_view(self):
        self.autoRange()

    def refresh_theme(self):
        style_pg_plot(self, title='Pathway Activity',
                      left_label='Score', bottom_label='Pathway')
