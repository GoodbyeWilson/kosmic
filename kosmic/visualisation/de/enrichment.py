"""
Enrichment Plots
================
The top enriched terms as a dot plot (fold enrichment, gene count and
significance) or a bar plot of -log10 significance.
"""

import re
import textwrap

import numpy as np

from kosmic.numerical import neg_log10
from kosmic import DEFAULT_FDR
from kosmic.visualisation.style import FG

_GO_ID = re.compile(r'^GO:\d+\s+')


def top_terms(enrichment_df, fdr_threshold=DEFAULT_FDR, top_n=20,
              significance_label='FDR'):
    """The first *top_n* terms passing *fdr_threshold* on the ``FDR``
    column (the input is ordered by significance), with GO IDs dropped
    from the names. If none pass, the first *top_n* rows and a title
    saying so. Returns (DataFrame, title or None)."""
    sig = enrichment_df[enrichment_df['FDR'] < fdr_threshold].head(top_n)
    title = None
    if sig.empty:
        sig = enrichment_df.head(top_n)
        title = f'No terms at {significance_label} < {fdr_threshold:g}: top {len(sig)} shown'
    sig = sig.copy()
    sig['Term'] = sig['Term'].astype(str).str.replace(_GO_ID, '', regex=True)
    return sig, title


def create_enrichment_dot_plot(enrichment_df, *, fdr_threshold=DEFAULT_FDR, top_n=20,
                               term_wrap=45, figsize=None, font_sizes=None,
                               title='Top Enriched Terms', significance_label='FDR'):
    """Dot plot of the top enriched terms, most significant at the top.

    x is the fold enrichment (the share of the query's genes in the term
    over the share expected from the background), dot area the number of
    query genes in the term, and colour -log10 of the ``FDR`` column, named
    by *significance_label*. Needs ``Term``, ``FDR``, ``Fold_Enrichment``
    and ``Gene_Count``. Returns None when there is nothing to plot.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    if enrichment_df is None or enrichment_df.empty:
        return None
    df, note = top_terms(enrichment_df, fdr_threshold, top_n, significance_label)
    if df.empty:
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    n = len(df)
    y = np.arange(n)[::-1]
    x = df['Fold_Enrichment'].to_numpy(dtype=float)
    count = df['Gene_Count'].to_numpy(dtype=float)
    sig = neg_log10(df['FDR'].to_numpy(dtype=float))

    def area(c):
        return 30 + 220 * np.sqrt(np.asarray(c, float) / max(count.max(), 1))

    fig, ax = plt.subplots(figsize=figsize or (8, 1.5 + 0.38 * n))
    cmap = LinearSegmentedColormap.from_list('sig', ['#D9D9D9', FG])
    pts = ax.scatter(x, y, s=area(count), c=sig, cmap=cmap, edgecolors=FG,
                     linewidths=0.5, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([textwrap.fill(t, term_wrap) for t in df['Term']],
                       fontsize=font_sizes['tick'])
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(left=min(1.0, float(np.nanmin(x)) * 0.9) if n else 0)
    ax.set_xlabel('Fold enrichment', fontsize=font_sizes['axis_label'])
    ax.set_title(note or title, fontsize=font_sizes['title'], fontweight='bold')
    ax.tick_params(labelsize=font_sizes['tick'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    cb = fig.colorbar(pts, ax=ax, shrink=0.45, pad=0.02, anchor=(0, 1))
    cb.set_label(f'$-$log$_{{10}}$ {significance_label}', fontsize=font_sizes['legend'])
    cb.ax.tick_params(labelsize=font_sizes['tick'])
    steps = sorted({int(v) for v in np.quantile(count, [0, 0.5, 1])})
    handles = [ax.scatter([], [], s=area(v), color='white', edgecolors=FG,
                          linewidths=0.5, label=f'{v}') for v in steps]
    ax.legend(handles=handles, title='Genes', frameon=False, loc='lower left',
              bbox_to_anchor=(1.02, 0), fontsize=font_sizes['legend'],
              title_fontsize=font_sizes['legend'], labelspacing=1.2, borderaxespad=0)
    fig.tight_layout()
    return fig


def create_enrichment_bar_plot(
    enrichment_df,
    *,
    fdr_threshold: float = DEFAULT_FDR,
    top_n: int = 25,
    term_truncate: int = 60,
    figsize=None,
    font_sizes=None,
    bg_color: str = 'white',
    fg_color: str = FG,
    bar_color: str = '#7A7A7A',
    sig_line_color: str = FG,
    title: str = 'Top Enriched Terms',
    significance_label: str = 'FDR',
):
    """Horizontal -log10 significance bar plot for the top enriched terms.

    Bars and the threshold line use the ``FDR`` column, which holds the
    value significance is judged on; *significance_label* names it on the
    axis (e.g. 'elim P (uncorrected)' for GO elim results, which are not
    adjusted). Keeps the first ``top_n`` rows passing the threshold
    (assumes the input is ordered by significance). If none pass, the top
    rows are shown and the title says so. Returns ``None`` when the input
    is empty.

    Parameters
    ----------
    enrichment_df : pandas.DataFrame
        Must contain ``Term`` and ``FDR`` columns.
    bg_color, fg_color, bar_color, sig_line_color : str
        Override these to theme the plot for an in-app dark preview;
        defaults are the house publication style.

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if enrichment_df is None or enrichment_df.empty:
        return None

    sig_df, note = top_terms(enrichment_df, fdr_threshold, top_n, significance_label)
    title = note or title
    if sig_df.empty:
        return None

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    neg_log_p = neg_log10(sig_df['FDR'])
    terms = [
        t[:term_truncate] + '...' if len(t) > term_truncate else t
        for t in sig_df['Term'].values
    ]

    if figsize is None:
        height = max(3, len(terms) * 0.35)
        figsize = (8, height)
    else:
        # Auto-scale height to the term count, keep caller's width.
        figsize = (figsize[0], max(3, len(terms) * 0.35))

    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor(bg_color)
    ax.set_facecolor(bg_color)

    ax.barh(range(len(terms)), neg_log_p, color=bar_color, edgecolor='none')
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels(terms, fontsize=font_sizes.get('tick', 9))
    ax.invert_yaxis()
    ax.set_xlabel(f'-log10 {significance_label}', fontsize=font_sizes.get('axis_label', 10))
    ax.set_title(title, fontsize=font_sizes.get('title', 12), fontweight='bold')

    ax.tick_params(colors=fg_color)
    for spine in ax.spines.values():
        spine.set_edgecolor(fg_color)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.label.set_color(fg_color)
    ax.yaxis.label.set_color(fg_color)
    ax.title.set_color(fg_color)

    ax.axvline(-np.log10(fdr_threshold), color=sig_line_color, linestyle=':',
               alpha=0.7, linewidth=1)
    plt.tight_layout()
    return fig
