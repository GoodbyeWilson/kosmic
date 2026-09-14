# Preranked GSEA (fgsea) NES bar chart: normalised enrichment score per pathway
# against the genome-wide gene ranking, coloured by FDR significance. Matplotlib
# mirror of the in-app NES bar, for publication export.

import numpy as np
import pandas as pd

from kosmic import DEFAULT_FDR


def create_fgsea_nes_plot(res, disease_label=None, control_label=None,
                          fdr_threshold=DEFAULT_FDR, figsize=None,
                          font_sizes=None, theme_colors=None, title=None):
    """Horizontal NES bar chart across pathways, sorted, coloured by FDR.

    'res' is the fgsea summary with columns 'names', 'nes', 'pvals_adj'. Bars
    extend from zero; NES < 0 means the pathway is enriched toward the
    control-up (disease-down) end of the ranking. Returns a Figure, or None if
    there is nothing to plot.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    if res is None or len(res) == 0:
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    tc = theme_colors or {}
    sig_color = tc.get('plot_disease', '#8E6FB0')
    ns_color = tc.get('plot_ns', '#B0B0B0')

    from kosmic.de.fgsea import select_top_nes
    df = select_top_nes(res)
    if df is None or df.empty:
        return None
    df = df.copy()
    df['pvals_adj'] = pd.to_numeric(df['pvals_adj'], errors='coerce')

    names = [str(n).replace('_', ' ') for n in df['names']]
    nes = df['nes'].to_numpy(dtype=float)
    padj = df['pvals_adj'].to_numpy(dtype=float)
    colors = [sig_color if (np.isfinite(p) and p < fdr_threshold) else ns_color
              for p in padj]

    n = len(names)
    if figsize is None:
        figsize = (7, max(2.5, 0.45 * n + 1.2))
    fig, ax = plt.subplots(figsize=figsize)
    y = np.arange(n)
    ax.barh(y, nes, color=colors, edgecolor='none', height=0.72, zorder=3)
    ax.axvline(0, color='black', linewidth=0.8, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=font_sizes['tick'])
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel('Normalised enrichment score (NES)',
                  fontsize=font_sizes['axis_label'], fontweight='bold')
    if title is None and disease_label and control_label:
        title = f'Pathway enrichment: {disease_label} vs {control_label}'
    if title:
        ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.legend(
        handles=[Patch(facecolor=sig_color, label=f'FDR < {fdr_threshold}'),
                 Patch(facecolor=ns_color, label='ns')],
        fontsize=font_sizes['legend'], loc='best', frameon=False)
    fig.tight_layout()
    return fig
