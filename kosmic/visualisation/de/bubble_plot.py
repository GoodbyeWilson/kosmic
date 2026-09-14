"""
Bubble Plot
============
Pathway bubble plot: bubble size = significance, color = effect size direction.
"""

import numpy as np
from matplotlib.patches import Patch

from kosmic.numerical import neg_log10


def create_bubble_plot(stats_df, disease_label, control_label,
                       plot_title=None, y_label=None, figsize=(12, 10),
                       font_sizes=None):
    """Create a pathway bubble plot from scoring statistics.

    Parameters
    ----------
    stats_df : pandas.DataFrame
        Must have columns: Pathway, Effect_Size, Log10_P_Corrected.
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

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

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
    df['Pathway_Clean'] = df['Pathway'].str.replace('_', ' ')

    fig, ax = plt.subplots(figsize=figsize)

    y_positions = np.arange(len(df))
    x = df['Effect_Size'].values

    bubble_sizes = df['Log10_P_Corrected'].values * 500
    bubble_sizes = np.maximum(bubble_sizes, 300)

    colors = []
    for es in df['Effect_Size']:
        if es > 0.3:
            colors.append('#d62728')
        elif es > 0.1:
            colors.append('#ff7f0e')
        elif es > -0.1:
            colors.append('#2ca02c')
        else:
            colors.append('#1f77b4')

    ax.scatter(x, y_positions, s=bubble_sizes, c=colors, alpha=0.7,
               edgecolors='black', linewidth=0.8)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(df['Pathway_Clean'].values, fontsize=font_sizes['tick'])
    ax.grid(False)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

    ax.set_xlabel(f'Effect Size ({disease_label} vs {control_label})',
                  fontsize=font_sizes['axis_label'], fontweight='bold')
    ax.set_ylabel(y_label or 'Gene Sets', fontsize=font_sizes['axis_label'], fontweight='bold')
    title = plot_title or f'Gene Set Activity: {disease_label} vs {control_label}'
    ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold', pad=20)

    for i, (idx, row) in enumerate(df.iterrows()):
        ax.text(row['Effect_Size'] + 0.02, i, f'{row["Effect_Size"]:.3f}',
                va='center', fontsize=font_sizes['annotation'], fontweight='bold')

    legend_elements = [
        Patch(facecolor='#d62728', alpha=0.7, label='Strong UP (>0.3)'),
        Patch(facecolor='#ff7f0e', alpha=0.7, label='Moderate UP (0.1-0.3)'),
        Patch(facecolor='#2ca02c', alpha=0.7, label='Neutral (-0.1-0.1)'),
        Patch(facecolor='#1f77b4', alpha=0.7, label='DOWN (<-0.1)'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=font_sizes['legend'], title='Effect Size')

    ax.text(0.02, 0.98, 'Bubble Size = Statistical Significance\n(Larger = Lower p-value)',
            transform=ax.transAxes, fontsize=font_sizes['axis_label'], va='top', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='gray', alpha=0.8))

    plt.tight_layout()
    plt.subplots_adjust(left=0.3)

    return fig
