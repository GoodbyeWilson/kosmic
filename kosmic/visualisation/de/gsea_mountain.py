# GSEA running-enrichment ("mountain") plot: the running enrichment score
# across the ranked gene list for one gene set, with a hit rug beneath.
# Matplotlib mirror of the DE window's in-app fgsea enrichment curve, for
# publication export.

import textwrap

import numpy as np

from kosmic.de.fgsea import PERMUTATIONS
from kosmic.visualisation.style import DOWN_COLOR, FG, UP_COLOR


def _format_fdr(fdr):
    # 0 is a permutation floor, not a value: no permutation beat the score.
    return f'FDR < {1 / PERMUTATIONS:g}' if fdr == 0 else f'FDR = {fdr:.2g}'


def create_gsea_mountain_plot(res_curve, hits, term, nes=None, fdr=None,
                              disease_label=None, control_label=None,
                              figsize=None, font_sizes=None):
    """Running-enrichment curve + hit rug for a single gene set.

    'res_curve' is the running enrichment score across the ranked genes,
    'hits' the rank positions of the set's genes. The peak (the ES) is
    marked. The curve is red when the set is enriched toward the
    disease-up end (ES > 0) and blue toward the control-up end. Returns a
    Figure, or None if there is nothing to plot.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    res_curve = np.asarray(res_curve, dtype=float)
    hits = np.asarray(hits, dtype=int)
    n = len(res_curve)
    if n == 0:
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    peak = int(np.argmax(np.abs(res_curve)))
    curve_color = UP_COLOR if res_curve[peak] > 0 else DOWN_COLOR
    fg = FG

    if figsize is None:
        figsize = (7, 4.5)

    x = np.arange(n)
    fig, (ax, ax_rug) = plt.subplots(
        2, 1, figsize=figsize, sharex=True,
        gridspec_kw={'height_ratios': [6, 1], 'hspace': 0.06})

    ax.fill_between(x, res_curve, 0, color=curve_color, alpha=0.15, linewidth=0)
    ax.plot(x, res_curve, color=curve_color, linewidth=1.6, zorder=3)
    ax.axhline(0, color=fg, linewidth=0.8)
    ax.plot(peak, res_curve[peak], 'o', color=curve_color, markersize=6,
            markeredgecolor=fg, markeredgewidth=0.5, zorder=4)
    ax.set_ylabel('Running enrichment score',
                  fontsize=font_sizes['axis_label'])

    title = textwrap.fill(str(term).replace('_', ' '), 45)
    ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold')
    stats = []
    if nes is not None and np.isfinite(nes):
        stats.append(f'NES = {nes:+.2f}')
    if fdr is not None and np.isfinite(fdr):
        stats.append(_format_fdr(fdr))
    if stats:
        # A falling curve leaves the lower left empty, a rising one the
        # upper right.
        down = res_curve[peak] < 0
        ax.text(0.02 if down else 0.98, 0.05 if down else 0.95,
                '\n'.join(stats), transform=ax.transAxes,
                ha='left' if down else 'right', va='bottom' if down else 'top',
                fontsize=font_sizes['annotation'], color=fg)
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
