# Preranked GSEA (fgsea) NES bar chart: normalised enrichment score per pathway
# against the genome-wide gene ranking, coloured by direction where significant. Matplotlib
# mirror of the in-app NES bar, for publication export.

import numpy as np
import pandas as pd

from kosmic import DEFAULT_FDR
from kosmic.visualisation.style import DOWN_COLOR, FG, NS_COLOR, UP_COLOR


def create_fgsea_nes_plot(res, disease_label=None, control_label=None,
                          fdr_threshold=DEFAULT_FDR, figsize=None,
                          font_sizes=None, title=None, x_limit=None):
    """Horizontal NES bar chart across pathways, sorted.

    'res' is the fgsea summary with columns 'names', 'nes', 'pvals_adj'. Bars
    extend from zero; NES < 0 means the pathway is enriched toward the
    control-up (disease-down) end of the ranking. Pathways passing
    *fdr_threshold* are red (NES > 0) or blue (NES < 0), the rest grey.
    The NES axis runs from -*x_limit* to +*x_limit*; without one, from
    minus to plus the largest |NES| shown, rounded up to the next 0.5. The
    plot is drawn in a full box.
    Returns a Figure, or None if there is nothing to plot.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    if res is None or len(res) == 0:
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    from kosmic.de.fgsea import select_top_nes
    df = select_top_nes(res)
    if df is None or df.empty:
        return None
    df = df.copy()
    df['pvals_adj'] = pd.to_numeric(df['pvals_adj'], errors='coerce')

    from kosmic.reference.pathways.libraries import term_label
    names = [term_label(n) for n in df['names']]
    nes = df['nes'].to_numpy(dtype=float)
    padj = df['pvals_adj'].to_numpy(dtype=float)
    sig = np.isfinite(padj) & (padj < fdr_threshold)
    colors = np.where(~sig, NS_COLOR, np.where(nes > 0, UP_COLOR, DOWN_COLOR))

    n = len(names)
    if figsize is None:
        figsize = (7, max(2.5, 0.45 * n + 1.2))
    fig, ax = plt.subplots(figsize=figsize)
    y = np.arange(n)
    ax.barh(y, nes, color=colors, edgecolor='none', height=0.72, zorder=3)
    ax.axvline(0, color=FG, linewidth=0.8, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=font_sizes['tick'])
    ax.set_ylim(-0.6, n - 0.4)
    if x_limit and x_limit > 0:
        lim = float(x_limit)
    else:
        finite = np.abs(nes[np.isfinite(nes)])
        lim = max(0.5, float(np.ceil(finite.max() * 2) / 2)) if len(finite) else 1.0
    ax.set_xlim(-lim, lim)
    ax.set_xlabel('Normalised enrichment score (NES)',
                  fontsize=font_sizes['axis_label'], fontweight='bold')
    if title is None and disease_label and control_label:
        title = f'Pathway enrichment: {disease_label} vs {control_label}'
    if title:
        # Anchored to the figure, not the axes: long set names push the
        # axes right, and a title centred over them runs off the page.
        fig.suptitle(title, x=0.01, ha='left', fontsize=font_sizes['title'],
                     fontweight='bold', color=FG)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(FG)
        spine.set_linewidth(0.8)

    ax.legend(
        handles=[Patch(facecolor=UP_COLOR, label=f'Up, FDR < {fdr_threshold:g}'),
                 Patch(facecolor=DOWN_COLOR, label=f'Down, FDR < {fdr_threshold:g}'),
                 Patch(facecolor=NS_COLOR, label='Not significant')],
        fontsize=font_sizes['legend'], loc='best', frameon=False)
    fig.tight_layout()
    return fig
