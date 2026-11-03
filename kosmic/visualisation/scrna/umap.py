"""
UMAP/t-SNE Publication Plots
==============================
Publication-quality embedding scatter plots.
"""

import numpy as np


# Distinct, print-friendly default (Tableau-style, saturated hues first).
_DEFAULT_PALETTE = [
    '#4E79A7', '#F28E2B', '#59A14F', '#E15759', '#B07AA1',
    '#76B7B2', '#EDC948', '#FF9DA7', '#9C755F', '#B6992D',
    '#499894', '#D37295', '#8CD17D', '#A0CBE8', '#FFBE7D',
    '#86BCB6', '#F1CE63', '#D4A6C8', '#79706E', '#BAB0AC',
    '#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E',
    '#E6AB02', '#A6761D', '#666666', '#386CB0', '#F0027F',
]
# Selectable categorical schemes: friendly name -> matplotlib qualitative cmap.
_CMAP_SCHEMES = {
    'Tableau 20': 'tab20', 'Bold': 'Set1', 'Dark': 'Dark2',
    'Paired': 'Paired', 'Pastel': 'Set3', 'Classic 10': 'tab10',
}
# Order shown in the UI.
PALETTE_NAMES = ['Default'] + list(_CMAP_SCHEMES.keys())


def _soft_palette(n, scheme='Default'):
    """Return n distinct colours from the named scheme (cycled if n exceeds it)."""
    if scheme in _CMAP_SCHEMES:
        import matplotlib.pyplot as plt
        import matplotlib.colors as mc
        base = [mc.to_hex(c) for c in plt.get_cmap(_CMAP_SCHEMES[scheme]).colors]
    else:
        base = _DEFAULT_PALETTE
    return [base[i % len(base)] for i in range(n)]


def category_colours(labels, scheme='Default', categories=None):
    """``{category: colour}`` as the embedding plots assign them (sorted
    categories, in palette order), so another figure of the same column
    uses the same colours.

    *categories* are further values to assign colours over, such as the
    values of the same column in the project's other studies. A category
    then has the same colour in every study, whichever of them it holds.
    """
    values = set(np.asarray(labels, dtype=str))
    if categories is not None:
        values |= {str(c) for c in categories}
    values = sorted(values)
    return dict(zip(values, _soft_palette(len(values), scheme)))


def _corner_axes(ax, embedding_name, fontsize):
    """Plain left and bottom axis lines with the axis names at the corner,
    no ticks."""
    from kosmic.visualisation.style import FG
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color(FG)
    ax.set_xlabel(f'{embedding_name}1', loc='left', fontsize=fontsize)
    ax.set_ylabel(f'{embedding_name}2', loc='bottom', fontsize=fontsize)


def _label_centres(ax, coords, labels, values, fontsize):
    """Write each value's name in bold at the median position of its cells,
    largest groups first, without overlap. Call once the layout is final."""
    from kosmic.visualisation.labels import place_labels
    sizes = {v: int((labels == v).sum()) for v in values}
    centres = [(v, *np.median(coords[labels == v], axis=0))
               for v in sorted(values, key=sizes.get, reverse=True) if sizes[v]]
    # Labels may run past the axes, so a long name stays on its cluster.
    place_labels(ax, centres, fontsize, centred=True, fontweight='bold',
                 avoid_points=False, within='figure')



def create_embedding_plot(coords, labels, title='', embedding_name='UMAP',
                          figsize=(12, 10), point_size=3, alpha=0.7,
                          dpi=150, rasterized=True, font_sizes=None,
                          scheme='Default', labels_on_plot=False, colours=None):
    """Create a single embedding scatter plot colored by categorical labels.

    Parameters
    ----------
    coords : numpy.ndarray
        (n_cells, 2) embedding coordinates.
    labels : array-like
        Categorical labels for each cell.
    title : str
        Plot title.
    embedding_name : str
        'UMAP' or 't-SNE' for axis labels.
    figsize : tuple
        Figure size in inches.
    point_size : float
        Scatter point size.
    alpha : float
        Point transparency.
    dpi : int
        Resolution for rasterized elements.
    rasterized : bool
        Rasterize scatter for publication PDFs.
    labels_on_plot : bool
        Write each category's name in bold at the median position of its
        cells, without overlap, instead of drawing a legend; the axes then
        show plain left and bottom axis lines.
    colours : dict, optional
        ``{category: colour}`` (see 'category_colours'); by default the
        categories present are coloured in sorted order.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    fig, ax = plt.subplots(figsize=figsize)

    labels = np.asarray(labels, dtype=str)
    unique_vals = sorted(set(labels))
    n_unique = len(unique_vals)
    if colours is None:
        colours = category_colours(labels, scheme)
    palette = [colours[v] for v in unique_vals]

    # Honour the requested point size exactly (may be sub-1 for dense
    # embeddings); large cell counts need genuinely tiny dots.
    s = float(point_size)

    for j, val in enumerate(unique_vals):
        mask = labels == val
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            c=[palette[j]], label=f"{val} ({mask.sum():,})",
            s=s, alpha=alpha, rasterized=rasterized, edgecolors='none',
            linewidths=0,
        )

    if title:
        ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold')
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])

    if labels_on_plot:
        # Category names on the plot instead of a legend.
        _corner_axes(ax, embedding_name, font_sizes['legend'])
        plt.tight_layout()
        _label_centres(ax, coords, labels, unique_vals, font_sizes['annotation'])
        return fig

    ax.set_xlabel(f'{embedding_name} 1', fontsize=font_sizes['axis_label'])
    ax.set_ylabel(f'{embedding_name} 2', fontsize=font_sizes['axis_label'])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Adaptive legend — scale from the configured legend size. Legend markers
    # are enlarged independently of the (possibly tiny) plotted point size so
    # they stay visible.
    base_legend = font_sizes['legend']
    ncol = 1 if n_unique <= 20 else 2 if n_unique <= 40 else 3
    legend_fs = (base_legend if n_unique <= 20
                 else base_legend * 0.8 if n_unique <= 40
                 else base_legend * 0.6)
    ax.legend(
        bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False,
        fontsize=legend_fs, ncol=ncol,
        markerscale=max(1.0, 6.0 / max(s, 0.5)),
        handletextpad=0.3, labelspacing=0.35, columnspacing=0.8,
    )

    plt.tight_layout()
    return fig


def create_highlight_plot(coords, labels, highlight_values, title='',
                          embedding_name='UMAP', figsize=(12, 10),
                          point_size=3, alpha=0.8, dpi=150, rasterized=True,
                          font_sizes=None, scheme='Default', labels_on_plot=False,
                          colours=None):
    """Embedding plot with only the chosen categories coloured, the rest grey.

    Parameters
    ----------
    coords : numpy.ndarray
        (n_cells, 2) embedding coordinates.
    labels : array-like
        Categorical labels.
    highlight_values : list
        Value(s) to colour; every other cell is drawn grey in the background.
    point_size : float
        Size of the highlighted points (background uses a smaller size).
    labels_on_plot : bool
        Name the highlighted categories on the plot instead of in a legend.
    colours : dict, optional
        ``{category: colour}``, so a highlighted category keeps its colour
        from the full plot; by default the highlighted categories take the
        palette in the order given.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    fig, ax = plt.subplots(figsize=figsize)

    labels = np.asarray(labels, dtype=str)
    highlight = [str(v) for v in highlight_values]
    s = float(point_size)
    bg = np.isin(labels, highlight, invert=True)

    # Background (everything not highlighted) in light grey, drawn first.
    ax.scatter(coords[bg, 0], coords[bg, 1], c='#d9d9d9',
               s=max(s * 0.5, 0.3), alpha=0.5, rasterized=rasterized,
               edgecolors='none', linewidths=0)

    # Each highlighted category in its own colour, on top.
    palette = (_soft_palette(len(highlight), scheme) if colours is None
               else [colours[v] for v in highlight])
    for j, val in enumerate(highlight):
        m = labels == val
        ax.scatter(coords[m, 0], coords[m, 1], c=[palette[j]], s=s, alpha=alpha,
                   rasterized=rasterized, edgecolors='none', linewidths=0,
                   label=f"{val} ({m.sum():,})")

    if title:
        ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold')
    ax.set_aspect('equal')
    if labels_on_plot:
        _corner_axes(ax, embedding_name, font_sizes['legend'])
        plt.tight_layout()
        _label_centres(ax, coords, labels, highlight, font_sizes['annotation'])
        return fig

    ax.set_xlabel(f'{embedding_name} 1', fontsize=font_sizes['axis_label'])
    ax.set_ylabel(f'{embedding_name} 2', fontsize=font_sizes['axis_label'])
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False,
              fontsize=font_sizes['legend'],
              markerscale=max(1.0, 6.0 / max(s, 0.5)), handletextpad=0.3)

    plt.tight_layout()
    return fig
