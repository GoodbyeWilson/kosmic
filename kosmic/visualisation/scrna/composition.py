"""
Cell-Type Composition
=====================
The percentage of cells of each type in each group (condition or sample),
as 100% stacked bars.
"""

import numpy as np
import pandas as pd

from kosmic.visualisation.labels import draw_bracket


def composition_table(cell_types, groups, group_order=None):
    """Percentage of each cell type within each group.

    Returns a DataFrame with one row per group (in *group_order* when
    given) and one column per cell type, rows summing to 100. Columns run
    from the most to the least abundant type overall (the mean of the
    groups' percentages, so each group counts equally).
    """
    counts = pd.crosstab(pd.Series(np.asarray(groups), name='group'),
                         pd.Series(np.asarray(cell_types), name='cell_type'))
    if group_order is not None:
        counts = counts.reindex([g for g in group_order if g in counts.index])
    pct = counts.div(counts.sum(axis=1), axis=0) * 100
    return pct[pct.mean().sort_values(ascending=False).index]


def create_composition_plot(pct, colours, *, blocks=None, title='Cell-type composition',
                            figsize=None, font_sizes=None):
    """100% stacked bars of *pct* (from ``composition_table``), one bar per
    row, the most abundant type at the bottom.

    *colours* maps each cell type to its colour. *blocks* is an optional
    list of (label, n_bars) bracketing consecutive bars, such as samples
    grouped by condition. The legend lists the types in stack order, top
    first.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    n = len(pct)
    many = n > 12
    fig, ax = plt.subplots(figsize=figsize or ((min(3 + 0.18 * n, 16), 6) if many
                                               else (4.0 + 0.9 * n, 6)))
    x = np.arange(n)
    bottom = np.zeros(n)
    for ct in pct.columns:
        v = pct[ct].to_numpy(dtype=float)
        ax.bar(x, v, bottom=bottom, width=0.85 if many else 0.75, color=colours[ct],
               edgecolor='white', linewidth=0 if many else 0.5, label=ct)
        bottom += v

    ax.set_xticks(x)
    ax.set_xticklabels(pct.index.astype(str), rotation=90 if many else 0,
                       fontsize=font_sizes['tick'] * (0.75 if many else 1))
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylim(0, 100)
    ax.set_ylabel('Cells (%)', fontsize=font_sizes['axis_label'])
    ax.tick_params(axis='y', labelsize=font_sizes['tick'])
    ax.tick_params(axis='x', length=0)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)

    pad = 6
    if blocks:
        start = 0
        for label, size in blocks:
            draw_bracket(ax, start, start + size - 1, label,
                         fontsize=font_sizes['axis_label'])
            start += size
        pad = font_sizes['axis_label'] * 2 + 8
    ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold', pad=pad)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1.02, 0.5),
              frameon=False, fontsize=font_sizes['legend'], handlelength=1.0,
              handleheight=1.0, labelspacing=0.5)
    fig.tight_layout()
    return fig
