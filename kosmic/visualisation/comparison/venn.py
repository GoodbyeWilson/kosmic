"""
Methods Venn diagram
====================
The genes two or three methods call, as a Venn diagram of equal circles:
the number in each region is the genes called by exactly those methods.
The circles are not drawn to scale, so the overlap is read from the
numbers, not the areas. A note under the diagram gives the genes compared
(those every method of the comparison tested) and how many of them none
of the methods shown calls.
"""

from kosmic.visualisation.style import CATEGORY_COLORS, FG, FG_SECONDARY

# Circle centres and label positions, in units of the radius (1).
_LAYOUT = {
    2: {'centres': [(-0.55, 0.0), (0.55, 0.0)],
        'labels': [(-1.05, 1.12, 'center'), (1.05, 1.12, 'center')],
        'regions': {(1, 0): (-0.95, 0.0), (0, 1): (0.95, 0.0), (1, 1): (0.0, 0.0)},
        'xlim': (-1.75, 1.75), 'ylim': (-1.2, 1.4)},
    3: {'centres': [(-0.55, 0.32), (0.55, 0.32), (0.0, -0.62)],
        'labels': [(-1.15, 1.45, 'center'), (1.15, 1.45, 'center'), (0.0, -1.83, 'center')],
        'regions': {(1, 0, 0): (-0.95, 0.55), (0, 1, 0): (0.95, 0.55),
                    (0, 0, 1): (0.0, -1.1), (1, 1, 0): (0.0, 0.72),
                    (1, 0, 1): (-0.55, -0.35), (0, 1, 1): (0.55, -0.35),
                    (1, 1, 1): (0.0, 0.05)},
        'xlim': (-1.75, 1.75), 'ylim': (-2.05, 1.7)},
}


def create_venn_plot(sets, labels=None, n_compared=None, title='Genes called',
                     figsize=None, font_sizes=None):
    """Venn diagram of two or three sets of genes.

    Parameters
    ----------
    sets : dict
        ``{method: set of genes}``, two or three methods in drawing order.
    labels : dict, optional
        ``{method: label}`` shown beside each circle; the key by default.
    n_compared : int, optional
        The genes compared; with it, the note gives the number none of the
        methods shown calls.
    title : str
    figsize : tuple, optional
    font_sizes : dict, optional

    Returns
    -------
    matplotlib Figure, or None unless there are two or three sets.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    from kosmic.comparison.agreement import region_counts

    if sets is None or len(sets) not in _LAYOUT:
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()
    methods = list(sets)
    labels = {m: (labels or {}).get(m, m) for m in methods}
    layout = _LAYOUT[len(methods)]

    fig, ax = plt.subplots(figsize=figsize or (5.0, 4.6 if len(methods) == 3 else 3.6))
    for i, m in enumerate(methods):
        colour = CATEGORY_COLORS[i % len(CATEGORY_COLORS)]
        ax.add_patch(Circle(layout['centres'][i], 1.0, facecolor=colour, alpha=0.45,
                            edgecolor='none'))
        ax.add_patch(Circle(layout['centres'][i], 1.0, facecolor='none', edgecolor=FG,
                            linewidth=0.8))
        x, y, ha = layout['labels'][i]
        ax.text(x, y, f"{labels[m]}\n{len(sets[m]):,} genes", ha=ha,
                va='bottom' if y > 0 else 'top', color=FG,
                fontsize=font_sizes['axis_label'], fontweight='bold')
    for flags, n in region_counts(sets).items():
        x, y = layout['regions'][flags]
        ax.text(x, y, f"{n:,}", ha='center', va='center', color=FG,
                fontsize=font_sizes['annotation'])

    if n_compared is not None:
        called = len(set().union(*sets.values()))
        note = (f"{n_compared:,} genes compared; {max(n_compared - called, 0):,} called by "
                "none of these methods.")
        fig.text(0.01, 0.01, note, ha='left', va='bottom', color=FG_SECONDARY,
                 fontsize=font_sizes['annotation'])
    ax.set_xlim(*layout['xlim'])
    ax.set_ylim(*layout['ylim'])
    ax.set_aspect('equal')
    ax.axis('off')
    fig.suptitle(title, x=0.01, ha='left', color=FG, fontsize=font_sizes['title'],
                 fontweight='bold')
    fig.tight_layout(rect=(0, 0.05 if n_compared is not None else 0, 1, 1))
    return fig


__all__ = ['create_venn_plot']
