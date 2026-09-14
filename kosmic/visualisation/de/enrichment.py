"""
Enrichment Bar Plot
====================
Horizontal bar plot of -log10(P) for the top enriched terms.
"""

import numpy as np

from kosmic.numerical import neg_log10
from kosmic import DEFAULT_FDR


def create_enrichment_bar_plot(
    enrichment_df,
    *,
    fdr_threshold: float = DEFAULT_FDR,
    top_n: int = 25,
    term_truncate: int = 60,
    figsize=None,
    font_sizes=None,
    bg_color: str = 'white',
    fg_color: str = 'black',
    bar_color: str = '#4682B4',  # steelblue
    sig_line_color: str = '#D62728',  # red
    title: str = 'Top Enriched Terms',
):
    """Horizontal -log10(P) bar plot for the top enriched terms.

    Filters to ``FDR < fdr_threshold`` and keeps the first ``top_n`` rows
    (assumes the input is already ordered by significance). Falls back to
    the unfiltered head if no rows pass the FDR cut. Returns ``None`` when
    the input is empty.

    Parameters
    ----------
    enrichment_df : pandas.DataFrame
        Must contain ``Term``, ``P_value``, and ``FDR`` columns.
    bg_color, fg_color, bar_color, sig_line_color : str
        Override these to theme the plot for an in-app dark preview;
        defaults are publication-style (white/black/steelblue/red).

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if enrichment_df is None or enrichment_df.empty:
        return None

    sig_df = enrichment_df[enrichment_df['FDR'] < fdr_threshold].head(top_n)
    if sig_df.empty:
        sig_df = enrichment_df.head(top_n)
    if sig_df.empty:
        return None

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    neg_log_p = neg_log10(sig_df['P_value'])
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
    ax.set_xlabel('-log10(P-value)', fontsize=font_sizes.get('axis_label', 10))
    ax.set_title(title, fontsize=font_sizes.get('title', 12), fontweight='bold')

    ax.tick_params(colors=fg_color)
    for spine in ax.spines.values():
        spine.set_edgecolor(fg_color)
    ax.xaxis.label.set_color(fg_color)
    ax.yaxis.label.set_color(fg_color)
    ax.title.set_color(fg_color)

    ax.axvline(-np.log10(fdr_threshold), color=sig_line_color, linestyle='--',
               alpha=0.7, linewidth=1)
    plt.tight_layout()
    return fig
