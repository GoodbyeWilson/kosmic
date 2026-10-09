"""
Marker Panel Figures
====================
A marker panel ({cell_type: [gene, ...]}) drawn two ways: a dot plot of
every panel gene against the data's cell types, and one embedding per
panel coloured by the panel's module score.
"""

import textwrap

import numpy as np
import pandas as pd

from kosmic.visualisation.labels import draw_bracket
from kosmic.visualisation.style import FG, sequential_cmap

#: Dot area, in points^2, of a gene detected in every cell of a type. Dot
#: area is proportional to the percentage, so a gene in 25% of cells is a
#: quarter of this.
MAX_DOT_AREA = 140


def order_rows(mean, owner):
    """Order the data's cell types so the dot plot reads as a diagonal.

    Each type is placed under the panel whose genes it expresses most (mean
    over that panel's genes), in panel order; types under the same panel are
    ordered by that mean, highest first. Placing by expression rather than
    by name works when the annotation and the panel name types differently
    ("LYVE1+ Macrophage" against "Myeloid").
    """
    panels = list(dict.fromkeys(owner[g] for g in mean.index))
    by_panel = pd.DataFrame({p: mean.loc[[g for g in mean.index if owner[g] == p]].mean()
                             for p in panels})
    best = by_panel.idxmax(axis=1)
    rank = {p: i for i, p in enumerate(panels)}
    return sorted(mean.columns, key=lambda t: (rank[best[t]], -by_panel.loc[t, best[t]]))


def create_marker_dotplot(mean, pct, owner, *, title='Marker genes',
                          colour_label='Mean expression', figsize=None,
                          font_sizes=None):
    """Dot plot of marker genes (columns, grouped by the cell type they mark)
    against the data's cell types (rows).

    Dot area is proportional to the percentage of cells expressing the
    gene. The colour scale is capped at the 99.5th percentile of the mean
    expression shown; an arrow on the colour bar marks the cap.

    Parameters
    ----------
    mean, pct : pandas.DataFrame
        genes x cell types: mean expression and the fraction (0-1) of cells
        with non-zero expression, as from ``summarise_gene_group``.
    owner : dict
        ``{gene: panel cell type}``; the genes are drawn in ``mean``'s order
        with a bracket over each panel's genes.
    colour_label : str
        Colour bar label, naming the expression unit.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    genes = list(mean.index)
    rows = order_rows(mean, owner)
    n_g, n_r = len(genes), len(rows)
    vals = mean.loc[genes, rows].to_numpy(dtype=float)
    frac = pct.loc[genes, rows].to_numpy(dtype=float)
    vmax = float(np.nanpercentile(vals, 99.5)) or 1.0

    fig, (ax, lax) = plt.subplots(
        1, 2, figsize=figsize or (0.24 * n_g + 3.2, 0.32 * n_r + 2.4),
        gridspec_kw={'width_ratios': [max(n_g, 1) * 0.24, 1.6]})
    gx, ry = np.meshgrid(np.arange(n_g), np.arange(n_r), indexing='ij')
    pts = ax.scatter(gx.ravel(), ry.ravel(), s=frac.ravel() * MAX_DOT_AREA,
                     c=vals.ravel(), cmap=sequential_cmap(), vmin=0, vmax=vmax,
                     edgecolors='none', zorder=3)
    ax.set_xticks(range(n_g))
    ax.set_xticklabels(genes, rotation=90, fontsize=font_sizes['tick'])
    ax.set_yticks(range(n_r))
    ax.set_yticklabels(rows, fontsize=font_sizes['tick'])
    ax.set_xlim(-0.6, n_g - 0.4)
    ax.set_ylim(n_r - 0.5, -0.5)
    ax.tick_params(length=0)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)

    start = 0
    for i, g in enumerate(genes + [None]):
        if g is None or owner[g] != owner[genes[start]]:
            draw_bracket(ax, start, i - 1, owner[genes[start]],
                         fontsize=font_sizes['annotation'], rotation=90)
            start = i
    # Room above the brackets for the longest vertical panel name.
    longest = max(len(str(owner[g])) for g in genes)
    ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold',
                 pad=font_sizes['annotation'] * 0.58 * longest + 12)

    # Legend column: a row of size-key dots, then a horizontal colour bar,
    # each titled above and labelled below.
    lax.axis('off')
    steps = (20, 40, 60, 80, 100)
    sax = lax.inset_axes([0.05, 0.62, 0.9, 0.12])
    sax.scatter(range(len(steps)), [0] * len(steps),
                s=[p / 100 * MAX_DOT_AREA for p in steps], color='#7A7A7A',
                edgecolors='none', clip_on=False)
    sax.set_xticks(range(len(steps)))
    sax.set_xticklabels([str(p) for p in steps], fontsize=font_sizes['tick'])
    sax.set_xlim(-0.6, len(steps) - 0.4)
    sax.set_yticks([])
    for s in sax.spines.values():
        s.set_visible(False)
    sax.set_title('Cells expressing (%)', fontsize=font_sizes['legend'], pad=10)

    cax = lax.inset_axes([0.05, 0.22, 0.9, 0.06])
    clipped = np.nanmax(vals) > vmax
    cb = fig.colorbar(pts, cax=cax, orientation='horizontal',
                      extend='max' if clipped else 'neither')
    cax.set_title(colour_label, fontsize=font_sizes['legend'])
    cb.ax.tick_params(labelsize=font_sizes['tick'])
    cb.outline.set_edgecolor(FG)
    cb.outline.set_linewidth(0.8)
    fig.tight_layout()
    return fig


def scale_scores(scores):
    """Put a module score on a 0-1 display scale.

    The median cell is set to 0, scores below it are clipped to 0, and the
    result is divided by its 99.9th percentile and capped at 1. Anchoring on
    the top of the signal rather than on the spread of the background lets
    a type making up 0.5% of cells reach the top of the scale. It cannot
    tell a rare type from an absent one: with no cells of a type, ambient
    counts of its genes are what reach the top, and show as scattered
    speckle.
    """
    v = np.clip(np.asarray(scores, float) - np.nanmedian(scores), 0, None)
    hi = np.nanpercentile(v, 99.9)
    return np.clip(v / hi, 0, 1) if hi > 0 else np.zeros_like(v)


def panel_scores(adata, panels, seed=0):
    """Module score of each panel for every cell of *adata*.

    Uses ``sc.tl.score_genes`` (mean expression of the panel's genes minus
    that of control genes drawn from the same expression bins), so a panel
    of highly expressed genes does not score high everywhere. Panels with
    fewer than two genes in the data are skipped. Expects log-normalised
    expression in ``adata.X``.

    Returns ``{panel: (scaled scores, genes used)}``, see ``scale_scores``.
    """
    import anndata as ad
    import scanpy as sc

    work = ad.AnnData(X=adata.X, var=pd.DataFrame(index=adata.var_names.astype(str)))
    out = {}
    for name, genes in panels.items():
        present = [g for g in dict.fromkeys(genes) if g in work.var_names]
        if len(present) < 2:
            continue
        sc.tl.score_genes(work, present, score_name='_score', ctrl_size=50,
                          random_state=seed)
        out[name] = (scale_scores(work.obs['_score'].to_numpy()), present)
    return out


def create_marker_score_grid(coords, scores, *, ncols=4,
                             point_size=0.5, title=None, figsize=None,
                             font_sizes=None):
    """One embedding per panel, every cell coloured by the panel's scaled
    module score, with the panel's genes written under it.

    *scores* is ``{panel: (scaled scores, genes)}`` from ``panel_scores``,
    each score array aligned with *coords*. Cells are drawn lowest score
    first, so high-scoring cells sit on top.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    n = len(scores)
    ncols = max(1, min(ncols, n))
    nrows = -(-n // ncols)
    fig, axs = plt.subplots(nrows, ncols, figsize=figsize or (3.1 * ncols + 0.8, 3.4 * nrows),
                            squeeze=False)
    cmap = sequential_cmap()
    pts = None
    for ax, (name, (v, genes)) in zip(axs.ravel(), scores.items()):
        o = np.argsort(v, kind='stable')
        pts = ax.scatter(coords[o, 0], coords[o, 1], c=v[o], cmap=cmap, vmin=0, vmax=1,
                         s=point_size, linewidths=0, rasterized=True)
        ax.set_title(name, fontsize=font_sizes['axis_label'], fontweight='bold')
        ax.set_xlabel(textwrap.fill(' '.join(genes), 40), fontsize=font_sizes['annotation'],
                      color=FG)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')
        for s in ax.spines.values():
            s.set_visible(False)
    for ax in axs.ravel()[n:]:
        ax.axis('off')
    if title:
        fig.suptitle(title, fontsize=font_sizes['title'], fontweight='bold')
    fig.tight_layout(rect=(0, 0, 0.93, 1))
    if pts is not None:
        cax = fig.add_axes([0.945, 0.4, 0.012, 0.2])
        cb = fig.colorbar(pts, cax=cax, ticks=[0, 1])
        cb.ax.set_yticklabels(['Low', 'High'], fontsize=font_sizes['tick'])
        cb.set_label('Module score', fontsize=font_sizes['legend'])
        cb.outline.set_visible(False)
    return fig
