"""
Bubble Plot
============
Pathway bubble plot: x = effect size, bubble size = significance, colour =
direction for significant pathways (grey otherwise).
"""

import numpy as np

from kosmic import DEFAULT_FDR
from kosmic.numerical import neg_log10
from kosmic.visualisation.style import DOWN_COLOR, FG, NS_COLOR, UP_COLOR


def create_bubble_plot(stats_df, disease_label, control_label,
                       plot_title=None, y_label=None, figsize=(8, 6),
                       font_sizes=None, fdr_threshold=DEFAULT_FDR):
    """Create a pathway bubble plot from scoring statistics.

    Pathways passing *fdr_threshold* are red (up in disease) or blue
    (down); the rest are grey. Bubble area scales with -log10 FDR; the
    size legend gives the scale.

    Parameters
    ----------
    stats_df : pandas.DataFrame
        Must have columns: Pathway, Effect_Size, Log10_P_Corrected (or the
        DE-format names, logfoldchanges, pvals_adj).
    disease_label, control_label : str
        Condition labels for axis annotation.
    plot_title : str, optional
        Custom plot title.
    y_label : str, optional
        Custom y-axis label.
    figsize : tuple
        Figure size.

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

    df = stats_df.copy()

    # Accept standard DE format (names, logfoldchanges, pvals_adj)
    if 'names' in df.columns and 'Pathway' not in df.columns:
        df = df.rename(columns={'names': 'Pathway'})
    if 'logfoldchanges' in df.columns and 'Effect_Size' not in df.columns:
        df = df.rename(columns={'logfoldchanges': 'Effect_Size'})
    if 'Log10_P_Corrected' not in df.columns and 'pvals_adj' in df.columns:
        df['Log10_P_Corrected'] = neg_log10(df['pvals_adj'])

    df = df.sort_values('Effect_Size', ascending=True)
    labels = df['Pathway'].str.replace('_', ' ').values
    x = df['Effect_Size'].to_numpy(dtype=float)
    log_p = df['Log10_P_Corrected'].to_numpy(dtype=float)
    sig = log_p > -np.log10(fdr_threshold)
    colors = np.where(~sig, NS_COLOR, np.where(x > 0, UP_COLOR, DOWN_COLOR))

    def area(v):
        return 40 + 120 * np.asarray(v, dtype=float)

    fig, ax = plt.subplots(figsize=figsize)
    y = np.arange(len(df))
    ax.axvline(0, color=FG, linewidth=0.8)
    ax.scatter(x, y, s=area(log_p), c=colors, edgecolors=FG, linewidth=0.5,
               zorder=3)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=font_sizes['tick'])
    lim = max(float(np.abs(x).max()) * 1.15, 0.1)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-0.7, len(df) - 0.3)
    ax.set_xlabel(f'Effect size ({disease_label} vs {control_label})',
                  fontsize=font_sizes['axis_label'])
    if y_label:
        ax.set_ylabel(y_label, fontsize=font_sizes['axis_label'])
    ax.set_title(plot_title or f'Gene Set Activity: {disease_label} vs {control_label}',
                 fontsize=font_sizes['title'], fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=font_sizes['tick'])

    # Size legend from round values within the data range, plus the colours.
    top = max(1, int(np.ceil(np.nanmax(log_p))))
    steps = sorted({1, max(1, top // 2), top})
    handles = [ax.scatter([], [], s=area(v), color='white', edgecolors=FG,
                          linewidth=0.5, label=f'{v}') for v in steps]
    handles += [ax.scatter([], [], s=area(1), color=c, edgecolors=FG,
                           linewidth=0.5, label=lab)
                for c, lab in ((UP_COLOR, f'Up, FDR < {fdr_threshold:g}'),
                               (DOWN_COLOR, f'Down, FDR < {fdr_threshold:g}'),
                               (NS_COLOR, 'Not significant'))]
    ax.legend(handles=handles, title='$-$log$_{10}$ FDR', frameon=False,
              fontsize=font_sizes['legend'], title_fontsize=font_sizes['legend'],
              loc='upper left', bbox_to_anchor=(1.02, 1), labelspacing=1.2,
              borderaxespad=0)

    fig.tight_layout()
    return fig
