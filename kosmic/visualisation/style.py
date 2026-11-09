# Figure export settings and the single export entry point (ADR-009).
#
# This module sits in kosmic/visualisation/ rather than kosmic/gui/ so that
# analysis code can use it without importing Qt. DPI is therefore passed in
# by the caller instead of being read from the GUI's settings.

from __future__ import annotations

from contextlib import contextmanager, nullcontext
from pathlib import Path

# Type 42 (TrueType) keeps PDF and PostScript text as text; matplotlib's
# default Type 3 fonts are imported into vector graphics editors as outlines.
# svg.fonttype 'none' writes SVG text as <text> elements instead of paths.
# Figures use Arial; matplotlib's own default is DejaVu Sans. Helvetica and
# Liberation Sans (metric-compatible with Arial) stand in where Arial is not
# installed.
EXPORT_RCPARAMS = {
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
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
# Secondary text (footnotes, keys), and cells or marks with no value.
FG_SECONDARY = '#666666'
MISSING_COLOR = '#E6E6E6'
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


#: Axis lines (spines) and tick marks in every figure: black at 1 pt, applied
#: by 'style_axes' to the preview and the export, so each plot function does
#: not have to set them and none can drift from the others.
AXIS_COLOR = '#000000'
AXIS_LINEWIDTH = 1.0


def style_axes(fig) -> None:
    """Draw every axis line and tick mark of *fig* in AXIS_COLOR at
    AXIS_LINEWIDTH. Tick labels and other text keep their colour."""
    for ax in fig.get_axes():
        for spine in ax.spines.values():
            spine.set_edgecolor(AXIS_COLOR)
            spine.set_linewidth(AXIS_LINEWIDTH)
        ax.tick_params(which='both', color=AXIS_COLOR, width=AXIS_LINEWIDTH)


def _display_points(artist, renderer):
    """Points (n, 2) in page (display) coordinates that bound *artist*, or
    None when they cannot be found. Scatter points count by their centres:
    a marker at the edge of the axes is drawn whole either way."""
    import numpy as np
    from matplotlib.collections import Collection, PathCollection
    from matplotlib.lines import Line2D

    empty = np.empty((0, 2))
    try:
        if isinstance(artist, Line2D):
            xy = np.asarray(artist.get_xydata(), dtype=float)
            return artist.get_transform().transform(xy) if len(xy) else empty
        if isinstance(artist, PathCollection):      # scatter: one marker per offset
            offsets = np.asarray(artist.get_offsets(), dtype=float)
            return artist.get_offset_transform().transform(offsets) if len(offsets) else empty
        if isinstance(artist, Collection):          # areas, line segments: their vertices
            verts = [p.vertices for p in artist.get_paths() if len(p.vertices)]
            return artist.get_transform().transform(np.vstack(verts)) if verts else empty
        box = artist.get_window_extent(renderer)
        return np.array([[box.x0, box.y0], [box.x1, box.y1]])
    except Exception:                     # an artist that cannot say: keep its clip
        return None


def unclip_contained(fig, tolerance: float = 1.0) -> int:
    """Switch off clipping for artists that lie entirely inside their axes.

    A clipped artist exports with a clip path; in PDF and SVG a vector
    editor turns that into a clipping mask (one per point of a scatter --
    thousands in a volcano -- and one per bar, line or label). An artist
    wholly inside the axes is not cut by the clip, so removing it changes
    nothing on the page. The test is made in page coordinates, so it holds
    for log axes, categorical axes and lines spanning the axes (axhline).
    Artists that reach past the axes by more than *tolerance* pixels, and
    rasterised ones, keep their clip. Returns the number changed.
    """
    import numpy as np

    renderer = fig.canvas.get_renderer() if hasattr(fig.canvas, 'get_renderer') \
        else fig._get_renderer()
    n = 0
    for ax in fig.get_axes():
        box = ax.get_window_extent(renderer)
        artists = [*ax.collections, *ax.lines, *ax.patches, *ax.images, *ax.texts]
        for artist in artists:
            if (not artist.get_clip_on() or artist.get_rasterized()
                    or not artist.get_visible() or artist is ax.patch):
                continue
            pts = _display_points(artist, renderer)
            if pts is None:
                continue
            pts = pts[np.isfinite(pts).all(axis=1)]
            if not len(pts) or (
                    pts[:, 0].min() >= box.x0 - tolerance and pts[:, 0].max() <= box.x1 + tolerance
                    and pts[:, 1].min() >= box.y0 - tolerance
                    and pts[:, 1].max() <= box.y1 + tolerance):
                artist.set_clip_on(False)
                n += 1
    return n


@contextmanager
def _markers_as_paths():
    """Write scatter points and line markers as plain paths in PDF and SVG.

    Both backends draw a repeated marker once and place it by reference (a
    PDF form object, an SVG <use>), which a vector editor imports as one
    group with its own clipping mask per point. The base renderer draws
    each as an ordinary path instead: same picture, slightly larger file.
    """
    from matplotlib.backend_bases import RendererBase
    from matplotlib.backends.backend_pdf import RendererPdf
    from matplotlib.backends.backend_svg import RendererSVG

    saved = []
    for cls in (RendererPdf, RendererSVG):
        for name in ('draw_path_collection', 'draw_markers'):
            saved.append((cls, name, cls.__dict__.get(name)))
            setattr(cls, name, getattr(RendererBase, name))
    try:
        yield
    finally:
        for cls, name, original in saved:
            if original is None:
                delattr(cls, name)
            else:
                setattr(cls, name, original)


#: OpenType features switched off in exported text. Matplotlib 3.11 shapes
#: text with HarfBuzz, which joins 'fi' and 'fl' into one ligature glyph in
#: fonts that have one (DejaVu Sans, Liberation Sans). A PDF writes that
#: glyph as a separate piece of text, so the word opens in a vector editor
#: split around it ('signi fi cant').
NO_LIGATURES = ('-liga', '-clig', '-dlig')


def disable_ligatures(fig) -> None:
    """Switch off ligatures in every text object of *fig*."""
    from matplotlib.text import Text

    for text in fig.findobj(Text):
        own = [f for f in (text.get_fontfeatures() or ()) if f not in NO_LIGATURES]
        text.set_fontfeatures([*own, *NO_LIGATURES])


def save_figure(fig, path, dpi: int, close: bool = True) -> Path:
    """Save *fig* to *path* in the format given by its suffix.

    The background is always white, so a figure drawn under the dark GUI
    theme does not export with a dark background. Axis lines and tick marks
    are drawn black at 1 pt ('style_axes'). In PDF and SVG, artists wholly
    inside their axes are not clipped and scatter points are plain paths,
    so a vector editor opens them without a clipping mask per point. *dpi*
    sets the PNG resolution and the resolution of rasterised parts of a
    PDF or SVG. Returns the path written.
    """
    import matplotlib.pyplot as plt

    path = Path(path)
    fmt = path.suffix.lower().lstrip('.')
    if fmt not in FORMATS:
        raise ValueError(f"Unsupported format {path.suffix!r}; use one of {FORMATS}")

    apply_export_rcparams()
    rasterize_dense_collections(fig)
    disable_ligatures(fig)
    style_axes(fig)
    vector = fmt in ('pdf', 'svg')
    if vector:
        unclip_contained(fig)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with _markers_as_paths() if vector else nullcontext():
            fig.savefig(path, format=fmt, dpi=dpi, bbox_inches='tight',
                        facecolor='white', edgecolor='none')
    finally:
        if close:
            plt.close(fig)
    return path
