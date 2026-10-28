"""
Study Agreement (Vote Levels)
=============================
Genes by the number of studies that call them significant in the same
direction (0/k ... k/k), split into up- and down-regulated, with genes
called up in one study and down in another as a separate 'conflict' bar.
Optionally overlays how many genes at each level another method calls
significant, to show what that method adds to, or confirms in, the
studies' own calls.
"""

import numpy as np

from kosmic.visualisation.style import DOWN_COLOR, FG, NS_COLOR, UP_COLOR


def create_vote_levels_plot(levels, title='Study agreement', significant_label=None,
                            figsize=None, font_sizes=None):
    """Bar chart of genes per agreement level.

    Parameters
    ----------
    levels : DataFrame
        From 'kosmic.meta_analysis.pooling.vote_count.vote_levels': one row
        per level with 'up', 'down', 'total' and optionally 'significant'.
    title : str
    significant_label : str, optional
        Legend label for the 'significant' overlay (e.g. 'REML FDR < 0.05');
        the overlay is drawn only when given and the column is present.
    figsize : tuple, optional
    font_sizes : dict, optional

    Returns
    -------
    matplotlib Figure, or None for an empty table.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    if levels is None or levels.empty:
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    labels = [str(lv) for lv in levels.index]
    x = np.arange(len(labels))
    conflict = np.array([lv == 'conflict' for lv in labels])
    up = levels['up'].to_numpy(float)
    down = levels['down'].to_numpy(float)
    total = levels['total'].to_numpy(float)
    # The 0/k bar has no direction; draw its total in grey like conflicts.
    neutral = conflict | ((up + down) == 0)

    fig, ax = plt.subplots(figsize=figsize or (1.3 + 0.9 * len(labels), 4.8))
    width = 0.7
    ax.bar(x[~neutral], down[~neutral], width, color=DOWN_COLOR, edgecolor=FG,
           linewidth=0.6, label='Down')
    ax.bar(x[~neutral], up[~neutral], width, bottom=down[~neutral], color=UP_COLOR,
           edgecolor=FG, linewidth=0.6, label='Up')
    ax.bar(x[neutral], total[neutral], width, color=NS_COLOR, edgecolor=FG,
           linewidth=0.6)

    handles = [Patch(facecolor=UP_COLOR, edgecolor=FG, label='Up'),
               Patch(facecolor=DOWN_COLOR, edgecolor=FG, label='Down'),
               Patch(facecolor=NS_COLOR, edgecolor=FG, label='None / conflict')]
    show_sig = significant_label and 'significant' in levels.columns
    if show_sig:
        sig = levels['significant'].to_numpy(float)
        ax.bar(x, sig, width * 0.45, facecolor='none', edgecolor='black',
               hatch='///', linewidth=0.9)
        handles.append(Patch(facecolor='none', edgecolor='black', hatch='///',
                             label=significant_label))

    top = total.max() if len(total) else 0
    for xi, t in zip(x, total):
        ax.text(xi, t + top * 0.01, f'{int(t):,}', ha='center', va='bottom',
                fontsize=font_sizes['tick'])
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=font_sizes['tick'])
    ax.set_ylim(0, top * 1.12 if top else 1)
    ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold')
    ax.set_xlabel('Studies significant in the same direction (FDR < 0.05)',
                  fontsize=font_sizes['axis_label'], fontweight='bold')
    ax.set_ylabel('Genes', fontsize=font_sizes['axis_label'], fontweight='bold')
    ax.tick_params(direction='out', labelsize=font_sizes['tick'])
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.legend(handles=handles, frameon=False, fontsize=font_sizes['tick'])
    fig.tight_layout()
    return fig
