# Figure export settings and the single export entry point (ADR-009).
#
# This module sits in kosmic/visualisation/ rather than kosmic/gui/ so that
# analysis code can use it without importing Qt. DPI is therefore passed in
# by the caller instead of being read from the GUI's settings.

from __future__ import annotations

from pathlib import Path

# Type 42 (TrueType) keeps PDF and PostScript text as text; matplotlib's
# default Type 3 fonts are imported into vector graphics editors as outlines.
# svg.fonttype 'none' writes SVG text as <text> elements instead of paths.
EXPORT_RCPARAMS = {
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
}

#: Scatter collections with at least this many points are rasterised on
#: export. Smaller ones stay vector so each point can be selected.
RASTERIZE_THRESHOLD = 20_000

FORMATS = ('png', 'pdf', 'svg')

# House palette for exported figures, fixed so a figure looks the same
# whichever app theme it was exported under: up / down in the disease
# group, not significant, and text and lines.
UP_COLOR = '#E06666'
DOWN_COLOR = '#6D9EEB'
NS_COLOR = '#B0B0B0'
FG = '#333333'
# Condition groups: the disease group takes the up colour, control the down.
DISEASE_COLOR = UP_COLOR
CONTROL_COLOR = DOWN_COLOR
# Categories that carry no direction (e.g. enriched terms): soft colours
# with no red, blue or grey, so they never read as up, down or ns.
CATEGORY_COLORS = ['#8DD3C7', '#BEBADA', '#FDB462', '#B3DE69',
                   '#FCCDE5', '#FFED6F', '#CCEBC5', '#BC80BD']


def diverging_cmap():
    """Blue (down) - white - red (up) colour map for values centred on zero,
    such as log2 fold changes and z-scores."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list('kosmic_diverging',
                                             [DOWN_COLOR, '#F7F7F7', UP_COLOR])


def sequential_cmap():
    """Light grey - red - dark red colour map for values that start at zero,
    such as expression and module scores. The low end stays visible on a
    white background."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list('kosmic_sequential',
                                             ['#D9D9D9', UP_COLOR, '#7A1F1F'])


def apply_export_rcparams() -> None:
    """Set the export rcParams for the whole process."""
    import matplotlib

    matplotlib.rcParams.update(EXPORT_RCPARAMS)


def rasterize_dense_collections(fig, threshold: int = RASTERIZE_THRESHOLD) -> int:
    """Rasterise scatter collections with at least *threshold* points.

    Text, axes and smaller collections stay vector. Returns the number of
    collections changed.
    """
    n = 0
    for ax in fig.get_axes():
        for coll in ax.collections:
            if len(coll.get_offsets()) >= threshold and not coll.get_rasterized():
                coll.set_rasterized(True)
                n += 1
    return n


def save_figure(fig, path, dpi: int, close: bool = True) -> Path:
    """Save *fig* to *path* in the format given by its suffix.

    The background is always white, so a figure drawn under the dark GUI
    theme does not export with a dark background. *dpi* sets the PNG
    resolution and the resolution of rasterised parts of a PDF or SVG.
    Returns the path written.
    """
    import matplotlib.pyplot as plt

    path = Path(path)
    fmt = path.suffix.lower().lstrip('.')
    if fmt not in FORMATS:
        raise ValueError(f"Unsupported format {path.suffix!r}; use one of {FORMATS}")

    apply_export_rcparams()
    rasterize_dense_collections(fig)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.savefig(path, format=fmt, dpi=dpi, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
    finally:
        if close:
            plt.close(fig)
    return path
