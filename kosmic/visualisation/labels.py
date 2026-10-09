"""
Label Placement
===============
Places text labels at points on an axes without labels overlapping, for
figures that label a few points (volcano genes, embedding clusters).
"""

from kosmic.visualisation.style import FG

# Candidate positions, nearest first: distance from the point in points,
# and the eight directions around it (right before left, up before down).
_RADII = (4, 10, 18, 28, 40, 55)
_DIRECTIONS = ((1, 1), (1, -1), (-1, 1), (-1, -1), (1, 0), (-1, 0), (0, 1), (0, -1))


def place_labels(ax, labels, fontsize, *, centred=False, fontweight='normal',
                 avoid_points=True, within='axes'):
    """Label points without overlap, in order of priority.

    *labels* is a list of (text, x, y) in data coordinates, most important
    first. Each label takes the nearest candidate position that lies inside
    the axes and clears every label already placed (and, with
    *avoid_points*, every labelled point). With *centred*, the first
    candidate is the label centred on its point. A label moved away from
    its point gets a thin leader line; a label with no free position is
    left out. *within* is 'axes' (labels stay inside the plot area) or
    'figure' (labels may extend past the axes but stay on the page). Call
    after the figure layout is final, because placement works in page
    coordinates. Returns the number placed.
    """
    renderer = ax.figure.canvas.get_renderer()
    axes_box = (ax.get_window_extent(renderer) if within == 'axes'
                else ax.figure.bbox)
    taken = []
    if avoid_points:
        for _, x, y in labels:
            px, py = ax.transData.transform((x, y))
            taken.append((px - 3, py - 3, px + 3, py + 3))

    def free(b):
        inside = (b.x0 >= axes_box.x0 and b.x1 <= axes_box.x1
                  and b.y0 >= axes_box.y0 and b.y1 <= axes_box.y1)
        return inside and not any(
            b.x0 < t[2] and b.x1 > t[0] and b.y0 < t[3] and b.y1 > t[1]
            for t in taken)

    candidates = [(0, 0, 0)] if centred else []
    candidates += [(r, dx, dy) for r in _RADII for dx, dy in _DIRECTIONS]
    # A leader line only once the label has clearly left its point.
    leader_from = _RADII[1] if centred else _RADII[0]

    placed = 0
    for text, x, y in labels:
        for r, dx, dy in candidates:
            ha = 'left' if dx > 0 else 'right' if dx < 0 else 'center'
            va = 'bottom' if dy > 0 else 'top' if dy < 0 else 'center'
            leader = (dict(arrowstyle='-', color=FG, lw=0.5, shrinkA=0, shrinkB=1)
                      if r > leader_from else None)
            ann = ax.annotate(text, (x, y), xytext=(dx * r, dy * r),
                              textcoords='offset points', ha=ha, va=va,
                              fontsize=fontsize, fontweight=fontweight,
                              color=FG, arrowprops=leader)
            b = ann.get_window_extent(renderer).padded(1)
            if free(b):
                taken.append((b.x0, b.y0, b.x1, b.y1))
                placed += 1
                break
            ann.remove()
    return placed


def draw_bracket(ax, x0, x1, text, fontsize, *, y=1.01, height=0.012, rotation=0):
    """Bracket over the data x range *x0*..*x1* (inclusive positions) just
    above the axes, with *text* over its centre. *y* and *height* are in
    axes fractions; the bracket is drawn outside the axes, so the figure's
    tight layout keeps room for it."""
    from matplotlib.transforms import blended_transform_factory

    tr = blended_transform_factory(ax.transData, ax.transAxes)
    a, b = x0 - 0.35, x1 + 0.35
    ax.plot([a, a, b, b], [y, y + height, y + height, y], color=FG, lw=0.8,
            transform=tr, clip_on=False)
    ax.text((x0 + x1) / 2, y + height * 1.6, text, transform=tr, fontsize=fontsize,
            ha='center', va='bottom', rotation=rotation, color=FG, clip_on=False)
