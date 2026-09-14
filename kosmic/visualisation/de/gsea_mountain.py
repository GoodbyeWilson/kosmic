# GSEA running-enrichment ("mountain") plot: the running enrichment score
# across the ranked gene list for one gene set, with a hit rug beneath.
# Matplotlib mirror of the DE window's in-app fgsea enrichment curve, for
# publication export.

import numpy as np


def create_gsea_mountain_plot(res_curve, hits, term, nes=None, fdr=None,
                              disease_label=None, control_label=None,
                              figsize=None, font_sizes=None, theme_colors=None):
    """Running-enrichment curve + hit rug for a single gene set.

    'res_curve' is the running enrichment score across the ranked genes,
    'hits' the rank positions of the set's genes. The peak (the ES) is
    marked. Returns a Figure, or None if there is nothing to plot.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    res_curve = np.asarray(res_curve, dtype=float)
    hits = np.asarray(hits, dtype=int)
    n = len(res_curve)
    if n == 0:
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    tc = theme_colors or {}
    curve_color = tc.get('plot_disease', '#2E8B57')
    fg = tc.get('fg', '#333333')

    if figsize is None:
        figsize = (6, 3.2)

    x = np.arange(n)
    fig, (ax, ax_rug) = plt.subplots(
        2, 1, figsize=figsize, sharex=True,
        gridspec_kw={'height_ratios': [6, 1], 'hspace': 0.06})

    ax.fill_between(x, res_curve, 0, color=curve_color, alpha=0.15, linewidth=0)
    ax.plot(x, res_curve, color=curve_color, linewidth=1.6, zorder=3)
    ax.axhline(0, color=fg, linewidth=0.8)
    peak = int(np.argmax(np.abs(res_curve)))
    ax.plot(peak, res_curve[peak], 'o', color=curve_color, markersize=6,
            markeredgecolor=fg, markeredgewidth=0.5, zorder=4)
    ax.set_ylabel('Running enrichment score',
                  fontsize=font_sizes['axis_label'])

    title = str(term).replace('_', ' ')
    bits = []
    if nes is not None and np.isfinite(nes):
        bits.append(f'NES = {nes:+.2f}')
    if fdr is not None and np.isfinite(fdr):
        bits.append(f'FDR = {fdr:.1e}')
    if bits:
        title += '  (' + ', '.join(bits) + ')'
    ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=font_sizes['tick'])

    ax_rug.vlines(hits, 0, 1, color=fg, linewidth=0.5)
    ax_rug.set_yticks([])
    ax_rug.set_xlim(0, n)
    xlab = 'Gene rank'
    if disease_label and control_label:
        xlab = f'Gene rank  ({disease_label}-up  to  {control_label}-up)'
    ax_rug.set_xlabel(xlab, fontsize=font_sizes['axis_label'])
    ax_rug.tick_params(labelsize=font_sizes['tick'])
    for s in ('top', 'right', 'left'):
        ax_rug.spines[s].set_visible(False)

    fig.tight_layout()
    return fig
