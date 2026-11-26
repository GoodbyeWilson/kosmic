"""
Methods Venn diagram
====================
The genes two, three or four methods call, as an area-proportional Venn
diagram: each shape's area is proportional to the number of genes the
method calls, and each overlap to the genes the methods share, so a method
that calls few genes is drawn small, and one whose calls are all called by
another lies inside it. The number in each region is the genes called by
exactly those methods. Two or three methods are drawn as circles; four as
ellipses fitted so that every region approaches its share of the genes.
The layout is the one the Meta-Analysis workspace's Venn uses
(``kosmic.visualisation.venn_layout``).

A note under the diagram gives the call rule, the genes compared (those
every method of the comparison tested) and how many of them none of the
methods shown calls.
"""

from kosmic.visualisation.style import CATEGORY_COLORS, FG, FG_SECONDARY

_GRID = 400   # sampling points per axis for the region label positions


def _inside(geom, xs, ys):
    """Boolean masks, one per shape, of the points (xs, ys) inside it."""
    import math

    import numpy as np
    masks = []
    for cx, cy, rx, ry, rot in geom:
        a = math.radians(rot)
        dx, dy = xs - cx, ys - cy
        lx = dx * math.cos(-a) - dy * math.sin(-a)
        ly = dx * math.sin(-a) + dy * math.cos(-a)
        masks.append((lx / max(rx, 1e-9)) ** 2 + (ly / max(ry, 1e-9)) ** 2 <= 1.0)
    return np.array(masks)


def region_positions(geom, n_sets):
    """``{membership flags: (x, y)}``: where to write each region's count, in
    layout coordinates: the point of the region furthest from its edges,
    found on a grid. (A region's mean point can fall outside it: the
    centre of a crescent lies in the shape it wraps.) A region too small
    to hit is left out."""
    import numpy as np
    from scipy.ndimage import distance_transform_edt
    g = np.linspace(0, 1, _GRID)
    gx, gy = np.meshgrid(g, g)
    masks = _inside(geom, gx.ravel(), gy.ravel()).reshape(len(geom), _GRID, _GRID)
    out = {}
    for code in range(1, 1 << n_sets):
        flags = tuple((code >> i) & 1 for i in range(n_sets))
        hit = np.ones((_GRID, _GRID), bool)
        for i, f in enumerate(flags):
            hit &= masks[i] if f else ~masks[i]
        if hit.any():
            r, c = np.unravel_index(np.argmax(distance_transform_edt(hit)), hit.shape)
            out[flags] = (float(gx[r, c]), float(gy[r, c]))
    return out


def create_venn_plot(sets, labels=None, n_compared=None, title='Genes called',
                     rule_text=None, figsize=None, font_sizes=None):
    """Area-proportional Venn diagram of two, three or four sets of genes.

    Parameters
    ----------
    sets : dict
        ``{method: set of genes}``, two to four methods in drawing order.
    labels : dict, optional
        ``{method: label}`` shown in the key; the method by default.
    n_compared : int, optional
        The genes compared; with it, the note gives the number none of the
        methods shown calls.
    title : str
    rule_text : str, optional
        How a gene is called ('adjusted p < 0.05 ...'), added to the note.
    figsize : tuple, optional
    font_sizes : dict, optional

    Returns
    -------
    matplotlib Figure, or None unless there are two to four sets.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse, Patch

    from kosmic.comparison.agreement import region_counts
    from kosmic.visualisation.venn_layout import venn_layout

    if sets is None or not 2 <= len(sets) <= 4:
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()
    methods = list(sets)
    labels = {m: (labels or {}).get(m, m) for m in methods}
    geom, _labels = venn_layout({m: set(sets[m]) for m in methods})

    fig, ax = plt.subplots(figsize=figsize or (7.5, 5.0))
    # The layout's y runs down the page; the axes' runs up.
    ax.set_xlim(0, 1)
    ax.set_ylim(1, 0)
    for i, (cx, cy, rx, ry, rot) in enumerate(geom):
        colour = CATEGORY_COLORS[i % len(CATEGORY_COLORS)]
        ax.add_patch(Ellipse((cx, cy), 2 * rx, 2 * ry, angle=rot, facecolor=colour,
                             edgecolor='none', alpha=0.45))
        ax.add_patch(Ellipse((cx, cy), 2 * rx, 2 * ry, angle=rot, facecolor='none',
                             edgecolor=FG, linewidth=0.8))
    counts = region_counts(sets)
    for flags, (x, y) in region_positions(geom, len(methods)).items():
        n = counts.get(flags, 0)
        if n:
            ax.text(x, y, f"{n:,}", ha='center', va='center', color=FG,
                    fontsize=font_sizes['annotation'])
    # A key rather than labels by the shapes: a method whose calls lie
    # inside another's is drawn inside it, where a label would cover it.
    handles = [Patch(facecolor=CATEGORY_COLORS[i % len(CATEGORY_COLORS)], alpha=0.6,
                     edgecolor=FG, linewidth=0.8,
                     label=f"{labels[m]} ({len(sets[m]):,} genes)")
               for i, m in enumerate(methods)]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False,
              fontsize=font_sizes['legend'], labelcolor=FG)

    notes = []
    if rule_text:
        notes.append(f"Called: {rule_text}.")
    if n_compared is not None:
        called = len(set().union(*sets.values()))
        notes.append(f"{n_compared:,} genes compared; {max(n_compared - called, 0):,} "
                     "called by none of these methods.")
    if notes:
        fig.text(0.01, 0.01, '\n'.join(notes), ha='left', va='bottom', color=FG_SECONDARY,
                 fontsize=font_sizes['annotation'])
    ax.set_aspect('equal')
    ax.axis('off')
    fig.suptitle(title, x=0.01, ha='left', color=FG, fontsize=font_sizes['title'],
                 fontweight='bold')
    fig.tight_layout(rect=(0, 0.035 * len(notes) + 0.015 if notes else 0, 1, 1))
    return fig


__all__ = ['create_venn_plot', 'region_positions']
