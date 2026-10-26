"""
Leading-Edge Heatmap
====================
The genes that drive a GSEA result, and how each changes in every dataset.

Preranked GSEA reports for each gene set its *leading edge*: the members
ranked before the running enrichment score peaks, which carry the
enrichment. The heatmap takes the strongest of them for one or more gene
sets (one block of rows per set) and shows each gene's log2 fold change in
each column -- every study of a meta-analysis and the pooled estimate, or
the one study of a DE run. A cell is marked where the gene is significant
in that column, and drawn grey where that column did not test the gene.
"""

from __future__ import annotations

import textwrap

import numpy as np
import pandas as pd

from kosmic.visualisation.style import FG, diverging_cmap

#: Colour for genes a column did not test.
UNTESTED_COLOR = '#E6E6E6'


def leading_edge_genes(leading_edge, ranking: pd.Series, nes: float,
                       n: int = 15) -> list[str]:
    """The *n* strongest leading-edge genes, in the data's spelling.

    Parameters
    ----------
    leading_edge : str or list of str
        The set's leading edge as GSEA reports it (';'-separated, or a
        list). GSEA matches genes in upper case, so these may be spelt
        differently from the data ('XBP1' for mouse 'Xbp1').
    ranking : pandas.Series
        The statistic the genes were ranked by, indexed by the data's gene
        names (the pooled z for a meta-analysis, the log2 fold change for
        a study's DE).
    nes : float
        The set's normalised enrichment score. Genes are ordered by the
        ranking in its direction, strongest first: highest first for a set
        enriched at the top (NES > 0), lowest first otherwise.
    n : int
        How many genes to keep.
    """
    if isinstance(leading_edge, str):
        leading_edge = [g for g in leading_edge.split(';') if g.strip()]
    spelling = {}
    for name in ranking.index.astype(str):
        spelling.setdefault(name.upper(), name)
    genes = [spelling[g.strip().upper()] for g in leading_edge
             if g.strip().upper() in spelling]
    scores = pd.to_numeric(ranking.reindex(genes), errors='coerce')
    order = scores.sort_values(ascending=not (nes is not None and nes >= 0),
                               na_position='last')
    return list(dict.fromkeys(order.index))[:n]


def leading_edge_table(genes: list[str], columns: list[tuple]):
    """Log2 fold changes and significance of *genes* in each column.

    Parameters
    ----------
    genes : list of str
    columns : list of (label, DataFrame, p_column)
        Each DataFrame has 'names' and 'logfoldchanges'; *p_column* names
        its adjusted p-value ('pvals_adj' for a study, 'fdr' for the
        pooled estimate).

    Returns
    -------
    (lfc, padj) : pandas.DataFrame
        genes x column labels. NaN where a column did not test the gene.
    """
    lfc, padj = {}, {}
    for label, df, p_col in columns:
        d = df.drop_duplicates('names')
        d = d.set_index(d['names'].astype(str))
        lfc[label] = pd.to_numeric(d['logfoldchanges'], errors='coerce').reindex(genes)
        padj[label] = (pd.to_numeric(d[p_col], errors='coerce').reindex(genes)
                       if p_col in d.columns else pd.Series(np.nan, index=genes))
    return (pd.DataFrame(lfc, index=genes), pd.DataFrame(padj, index=genes))


def create_leading_edge_heatmap(blocks, *, vmax: float = 2.0, alpha: float = 0.05,
                                gap_before=None, title=None,
                                colour_label='Log$_2$ fold change',
                                figsize=None, font_sizes=None):
    """Heatmap of leading-edge genes (rows) by dataset (columns).

    Parameters
    ----------
    blocks : list of dict
        One per gene set, drawn top to bottom: 'title' (the set and its
        NES / FDR), 'lfc' and 'padj' (genes x columns, from
        :func:`leading_edge_table`). All blocks share the same columns.
    vmax : float
        The colour scale runs from -vmax to +vmax; larger changes take the
        end colour.
    alpha : float
        A cell is marked '*' when its adjusted p-value is below this.
    gap_before : str, optional
        Column label to separate from those before it with a gap (the
        pooled column of a meta-analysis).
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()
    blocks = [b for b in blocks if len(b['lfc'])]
    if not blocks:
        return None

    cols = list(blocks[0]['lfc'].columns)
    n_cols = len(cols)
    heights = [len(b['lfc']) for b in blocks]
    if figsize is None:
        figsize = (max(3.6, 0.85 * n_cols + 2.6),
                   max(2.5, 0.24 * sum(heights) + 1.0 * len(blocks) + 0.6))
    # Constrained layout keeps the figure title, the wrapped set titles,
    # the tick labels and the colour bar clear of each other at any size.
    fig, axes = plt.subplots(len(blocks), 1, figsize=figsize, squeeze=False,
                             layout='constrained',
                             gridspec_kw={'height_ratios': heights})
    axes = axes[:, 0]
    cmap = diverging_cmap().with_extremes(bad=UNTESTED_COLOR)
    norm = Normalize(vmin=-vmax, vmax=vmax)
    gap = cols.index(gap_before) if gap_before in cols else None

    im = None
    for i, (ax, b) in enumerate(zip(axes, blocks)):
        values = b['lfc'][cols].to_numpy(dtype=float)
        im = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, norm=norm,
                       aspect='auto', interpolation='nearest')
        padj = b['padj'][cols].to_numpy(dtype=float)
        for r, c in zip(*np.nonzero(np.nan_to_num(padj, nan=1.0) < alpha)):
            ax.text(c, r, '*', ha='center', va='center', color=FG,
                    fontsize=font_sizes['annotation'])
        if gap is not None:
            ax.axvline(gap - 0.5, color='white', lw=4)
        ax.set_yticks(range(len(b['lfc'])))
        ax.set_yticklabels(b['lfc'].index, fontsize=font_sizes['tick'],
                           fontstyle='italic')
        ax.set_xticks(range(n_cols))
        if i == len(blocks) - 1:
            ax.set_xticklabels(cols, fontsize=font_sizes['tick'],
                               rotation=30 if n_cols > 3 else 0,
                               ha='right' if n_cols > 3 else 'center')
            ax.set_xlabel(f"* adjusted p < {alpha:g} in that column. "
                          f"Grey: not tested.", loc='left',
                          fontsize=font_sizes['annotation'], color='#666666')
        else:
            ax.set_xticklabels([])
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        # Library set names run long ('SIGNALING PATHWAYS REGULATING
        # PLURIPOTENCY OF STEM CELLS'); wrap to the heatmap's width.
        width = max(30, int(figsize[0] * 72 / font_sizes['axis_label'] * 0.9))
        ax.set_title(textwrap.fill(b['title'], width), loc='left',
                     fontsize=font_sizes['axis_label'], fontweight='bold', color=FG)

    if title:
        fig.suptitle(title, x=0.02, ha='left', fontsize=font_sizes['title'],
                     fontweight='bold', color=FG)
    cb = fig.colorbar(im, ax=list(axes), shrink=min(0.6, 6 / max(sum(heights), 1) + 0.3),
                      extend='both')
    cb.set_label(colour_label, fontsize=font_sizes['colorbar_label'], color=FG)
    cb.ax.tick_params(labelsize=font_sizes['colorbar_tick'], colors=FG)
    cb.outline.set_visible(False)
    return fig
