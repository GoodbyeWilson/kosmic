"""
Enriched Term Links
===================
Which DE genes drive the top enriched terms, as a chord diagram (genes on
the left, terms on the right, a ribbon for each gene-term membership) or
as a gene x term grid. Genes are coloured by log2 fold change.
"""

import re
import textwrap

import numpy as np
import pandas as pd

from kosmic import DEFAULT_FDR
from kosmic.de.de_analysis import excluded_genes
from kosmic.visualisation.style import CATEGORY_COLORS, FG, diverging_cmap

_GO_ID = re.compile(r'^GO:\d+\s+')


def select_term_genes(enrichment_df, de_results, n_terms=8, genes_per_term=8,
                      threshold=DEFAULT_FDR, exclude_groups=()):
    """The top terms and, for each, its DE genes with the largest |log2FC|.

    *enrichment_df* has ``Term``, ``FDR`` (the value significance is
    judged on) and ``Genes`` (a list or a comma-separated string of the
    term's genes from the query). Terms passing *threshold* are taken in
    order of significance; if none pass, the most significant are taken.
    *exclude_groups* names gene groups of
    ``kosmic.de.de_analysis.EXCLUDABLE_GENE_GROUPS`` (e.g. 'mitochondrial')
    left out of the figure; the enrichment result itself is unchanged.
    Returns ``({term: [genes]}, log2FC Series, n_passing)``; terms with no
    gene left are dropped.
    """
    df = enrichment_df.copy()
    df['FDR'] = pd.to_numeric(df['FDR'], errors='coerce')
    df = df.sort_values('FDR')
    n_passing = int((df['FDR'] < threshold).sum())
    top = df[df['FDR'] < threshold] if n_passing else df
    lfc_all = (de_results.drop_duplicates('names').set_index('names')['logfoldchanges']
               .pipe(pd.to_numeric, errors='coerce').dropna())
    lfc_all = lfc_all.drop(excluded_genes(lfc_all.index, exclude_groups))

    members = {}
    for term, genes in zip(top['Term'], top['Genes']):
        if isinstance(genes, str):
            genes = [g.strip() for g in genes.split(',')]
        genes = [g for g in dict.fromkeys(genes) if g in lfc_all.index]
        if genes:
            keep = lfc_all[genes].abs().nlargest(genes_per_term).index
            members[_GO_ID.sub('', str(term))] = [g for g in genes if g in keep]
        if len(members) == n_terms:
            break
    all_genes = list(dict.fromkeys(g for gs in members.values() for g in gs))
    return members, lfc_all[all_genes], n_passing


def order_genes(members, lfc, by='lfc'):
    """Genes from most down to most up (*by* 'lfc'), or grouped by their
    first term and then by fold change (*by* 'term')."""
    genes = list(lfc.index)
    if by == 'lfc':
        return sorted(genes, key=lfc.get)
    first = {}
    for i, gs in enumerate(members.values()):
        for g in gs:
            first.setdefault(g, i)
    return sorted(genes, key=lambda g: (first[g], lfc[g]))


def _cmap_norm(lfc):
    from matplotlib.colors import Normalize
    cmap = diverging_cmap()
    # Symmetric, capped at the 90th percentile of |log2FC| so one extreme
    # gene does not wash out the rest; the colour bar marks the cap.
    lim = max(float(np.quantile(np.abs(lfc), 0.9)), 0.1)
    return cmap, Normalize(vmin=-lim, vmax=lim, clip=True)


def _arc(r, a0, a1, n=30):
    t = np.radians(np.linspace(a0, a1, n))
    return np.column_stack([r * np.cos(t), r * np.sin(t)])


def create_term_chord(members, lfc, order='lfc', title=None, figsize=None,
                      font_sizes=None):
    """Chord diagram: genes around the left, terms around the right.

    Each gene arc is coloured by its log2 fold change and split between the
    terms it belongs to; each term arc is sized by its number of genes. A
    ribbon joins every gene to each of its terms, in the term's colour.
    Returns None when there is nothing to draw.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import PathPatch, Wedge
    from matplotlib.path import Path

    if not members:
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()
    terms = list(members)
    genes = order_genes(members, lfc, order)
    cmap, norm = _cmap_norm(lfc)

    # Genes on the left (80-280 deg, first at the top), terms on the right
    # (75 down to -75 deg, first at the top).
    gap = min(1.5, 30.0 / len(genes))
    g_w = (200.0 - gap * (len(genes) - 1)) / len(genes)
    gene_sec = {g: (80 + i * (g_w + gap), 80 + i * (g_w + gap) + g_w)
                for i, g in enumerate(genes)}
    unit = (150.0 - 4 * gap * (len(terms) - 1)) / sum(len(m) for m in members.values())
    term_sec, a = {}, 75.0
    for t in terms:
        term_sec[t] = (a - unit * len(members[t]), a)
        a -= unit * len(members[t]) + 4 * gap
    # Slots: each gene arc split among its terms, each term arc among its
    # genes in gene order, so ribbons within a term do not cross.
    gene_slots, term_slots = {}, {}
    for g in genes:
        ts = [t for t in terms if g in members[t]]
        a0, a1 = gene_sec[g]
        w = (a1 - a0) / len(ts)
        for k, t in enumerate(ts):
            gene_slots[(g, t)] = (a0 + k * w, a0 + (k + 1) * w)
    for t in terms:
        a0, a1 = term_sec[t]
        gs = sorted(members[t], key=genes.index)
        w = (a1 - a0) / len(gs)
        for k, g in enumerate(gs):
            term_slots[(g, t)] = (a0 + k * w, a0 + (k + 1) * w)

    fig, ax = plt.subplots(figsize=figsize or (10, 9))
    R, ring, n = 1.0, 0.07, 30
    codes = ([Path.MOVETO] + [Path.LINETO] * (n - 1) + [Path.CURVE3, Path.CURVE3]
             + [Path.LINETO] * (n - 1) + [Path.CURVE3, Path.CURVE3])
    for i, t in enumerate(terms):
        color = CATEGORY_COLORS[i % len(CATEGORY_COLORS)]
        for g in members[t]:
            ga, ta = _arc(R, *gene_slots[(g, t)], n), _arc(R, *term_slots[(g, t)], n)
            # Gene arc, curve through the centre, term arc, curve back.
            verts = np.vstack([ga, [[0, 0]], ta, [[0, 0]], ga[:1]])
            ax.add_patch(PathPatch(Path(verts, codes), facecolor=color,
                                   edgecolor='none', alpha=0.75))
        a0, a1 = term_sec[t]
        ax.add_patch(Wedge((0, 0), R + ring, a0, a1, width=ring, color=color))
        mid = np.radians((a0 + a1) / 2)
        r = R + ring + 0.05
        ax.text(r * np.cos(mid), r * np.sin(mid), textwrap.fill(t, 30),
                ha='left', va='center', fontsize=font_sizes['annotation'],
                fontweight='bold', color=FG)
    gene_fs = font_sizes['annotation'] * (0.8 if len(genes) > 30 else 0.9)
    for g in genes:
        a0, a1 = gene_sec[g]
        ax.add_patch(Wedge((0, 0), R + ring, a0, a1, width=ring, color=cmap(norm(lfc[g]))))
        mid = (a0 + a1) / 2
        r = R + ring + 0.03
        # Genes sit on the left: one orientation rule for every label, so
        # labels near the top and bottom do not flip.
        ax.text(r * np.cos(np.radians(mid)), r * np.sin(np.radians(mid)), g,
                rotation=mid - 180, rotation_mode='anchor', ha='right', va='center',
                fontsize=gene_fs, color=FG)

    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                      shrink=0.3, pad=0.02, location='left', extend='both')
    cb.set_label('Log$_2$ fold change', fontsize=font_sizes['legend'])
    cb.ax.tick_params(labelsize=font_sizes['tick'])
    ax.set_xlim(-1.55, 1.75)
    ax.set_ylim(-1.45, 1.45)
    ax.set_aspect('equal')
    ax.axis('off')
    if title:
        ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold')
    return fig


def create_term_grid(members, lfc, order='lfc', title=None, figsize=None,
                     font_sizes=None):
    """Gene x term grid: a row per term, a column per gene, and a dot
    coloured by log2 fold change where the gene belongs to the term.
    Returns None when there is nothing to draw."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if not members:
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()
    terms = list(members)
    genes = order_genes(members, lfc, order)
    cmap, norm = _cmap_norm(lfc)

    fig, ax = plt.subplots(figsize=figsize or (2.5 + 0.22 * len(genes),
                                               1.8 + 0.42 * len(terms)))
    for i, t in enumerate(terms):
        xs = [genes.index(g) for g in members[t]]
        ax.scatter(xs, [i] * len(xs), s=40, c=lfc[members[t]].to_numpy(), cmap=cmap,
                   norm=norm, edgecolors=FG, linewidths=0.4, zorder=3)
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels([textwrap.fill(t, 40) for t in terms], fontsize=font_sizes['tick'])
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=font_sizes['annotation'])
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(len(terms) - 0.5, -0.5)
    ax.set_xticks(np.arange(len(genes)) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(terms)) - 0.5, minor=True)
    ax.grid(which='minor', color='#E5E5E5', linewidth=0.6)
    ax.tick_params(which='both', length=0)
    for s in ax.spines.values():
        s.set_color('#BBBBBB')
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                      shrink=0.8, pad=0.01, extend='both')
    cb.set_label('Log$_2$ fold change', fontsize=font_sizes['legend'])
    cb.ax.tick_params(labelsize=font_sizes['tick'])
    if title:
        ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold')
    fig.tight_layout()
    return fig
