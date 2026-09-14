"""
Pathway Summary Plots
======================
Bar plots showing significant gene counts by pathway.
"""

import numpy as np
import pandas as pd
from kosmic import DEFAULT_FDR, DEFAULT_LFC_THRESHOLD


def create_pathway_bar_plot(de_results, pathway_coverage, disease_cond, control_cond,
                            pval_threshold=DEFAULT_FDR, logfc_threshold=DEFAULT_LFC_THRESHOLD,
                            figsize=None, font_sizes=None):
    """Create horizontal bar plot of significant genes per pathway.

    Parameters
    ----------
    de_results : pandas.DataFrame
        DE results with columns: names, logfoldchanges, pvals_adj.
    pathway_coverage : dict
        {pathway: {genes, total_genes, available_genes, coverage_pct}}.
    disease_cond, control_cond : str
        Condition labels.
    pval_threshold : float
        Significance threshold.
    logfc_threshold : float
        Fold change threshold.
    figsize : tuple, optional
        Auto-sized if None.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    significant = de_results[
        (de_results['pvals_adj'] < pval_threshold) &
        (de_results['logfoldchanges'].abs() > logfc_threshold)
    ]

    pathway_summary = []
    for pathway_name, info in pathway_coverage.items():
        pathway_genes = set(info['genes'])
        sig_up = sum(1 for _, row in significant.iterrows()
                     if row['names'] in pathway_genes and row['logfoldchanges'] > 0)
        sig_down = sum(1 for _, row in significant.iterrows()
                       if row['names'] in pathway_genes and row['logfoldchanges'] < 0)
        pathway_summary.append({
            'pathway': pathway_name.replace('_', ' '),
            'total_genes': info['total_genes'],
            'available_genes': info['available_genes'],
            'sig_up': sig_up,
            'sig_down': sig_down,
            'total_sig': sig_up + sig_down,
        })

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    pathway_df = pd.DataFrame(pathway_summary).sort_values('total_sig', ascending=True)

    if figsize is None:
        figsize = (14, max(8, len(pathway_df) * 0.5))

    fig, ax = plt.subplots(figsize=figsize)
    y_pos = np.arange(len(pathway_df))

    ax.barh(y_pos, pathway_df['sig_up'], left=0, color='red', alpha=0.8,
            label=f'Upregulated in {disease_cond}')
    ax.barh(y_pos, -pathway_df['sig_down'], left=0, color='blue', alpha=0.8,
            label=f'Downregulated in {disease_cond}')

    y_labels = [f"{row['pathway']} ({row['available_genes']}/{row['total_genes']} genes)"
                for _, row in pathway_df.iterrows()]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=font_sizes['tick'])
    ax.set_xlabel('Number of Significant Genes', fontsize=font_sizes['axis_label'], fontweight='bold')
    ax.set_title(f'Significant Metabolic Genes by Pathway\n({disease_cond} vs {control_cond})',
                 fontsize=font_sizes['title'], fontweight='bold')
    ax.axvline(x=0, color='black', linewidth=0.8)
    ax.legend(loc='lower right', fontsize=font_sizes['legend'])

    for i, (_, row) in enumerate(pathway_df.iterrows()):
        if row['sig_up'] > 0:
            ax.text(row['sig_up'] + 0.1, i, str(int(row['sig_up'])),
                    va='center', fontsize=font_sizes['annotation'], fontweight='bold')
        if row['sig_down'] > 0:
            ax.text(-row['sig_down'] - 0.1, i, str(int(row['sig_down'])),
                    va='center', ha='right', fontsize=font_sizes['annotation'], fontweight='bold')

    plt.tight_layout()
    return fig
