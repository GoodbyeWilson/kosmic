"""
Top DE Plots
=============
Pathway summary split bar chart and per-pathway gene bar charts.
"""

import numpy as np
import pandas as pd


def create_pathway_summary_chart(de_results, pathway_gene_sets, disease_label, control_label,
                                 label='Significant', figsize=None, font_sizes=None):
    """Create a split bar chart of up/down DE genes per pathway.

    Parameters
    ----------
    de_results : pandas.DataFrame
        DE results with 'names', 'logfoldchanges' columns.
    pathway_gene_sets : dict
        {pathway_name: [gene_list]}.
    disease_label, control_label : str
        Condition labels.
    label : str
        Label for plot (e.g. 'Significant', 'All Genes').
    figsize : tuple, optional
        Figure size. Auto-computed if None.

    Returns
    -------
    matplotlib.figure.Figure or None
        None if no pathways have DE genes.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    df = de_results
    pathway_stats = []
    for pathway_name, genes in pathway_gene_sets.items():
        genes_upper = {g.upper() for g in genes}
        pw_genes = df[df['names'].str.upper().isin(genes_upper)]
        n_up = int((pw_genes['logfoldchanges'] > 0).sum())
        n_down = int((pw_genes['logfoldchanges'] < 0).sum())
        if n_up + n_down > 0:
            pathway_stats.append({
                'pathway': pathway_name.replace('_', ' '),
                'n_up': n_up,
                'n_down': n_down,
                'total': n_up + n_down,
            })

    if not pathway_stats:
        return None

    pw_df = pd.DataFrame(pathway_stats).sort_values('total', ascending=True)
    n_pw = len(pw_df)
    fig_h = max(4, n_pw * 0.5 + 2)
    if figsize is None:
        figsize = (10, fig_h)

    fig, ax = plt.subplots(figsize=figsize)
    y = np.arange(n_pw)

    ax.barh(y, pw_df['n_up'].values, height=0.6, color='#d62728', alpha=0.85, label='Upregulated')
    ax.barh(y, -pw_df['n_down'].values, height=0.6, color='#1f77b4', alpha=0.85, label='Downregulated')

    ax.set_yticks(y)
    ax.set_yticklabels(pw_df['pathway'].values, fontsize=font_sizes['tick'])
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel(f'Number of DE Genes ({label})', fontsize=font_sizes['axis_label'], fontweight='bold')
    ax.set_title(f'Top DE Genes by Pathway ({label})\n{disease_label} vs {control_label}',
                 fontsize=font_sizes['title'], fontweight='bold')
    ax.legend(loc='lower right', fontsize=font_sizes['legend'])

    for i, row in enumerate(pw_df.itertuples()):
        if row.n_up > 0:
            ax.text(row.n_up + 0.1, i, str(row.n_up), va='center',
                    fontsize=font_sizes['annotation'], fontweight='bold')
        if row.n_down > 0:
            ax.text(-row.n_down - 0.1, i, str(row.n_down), va='center', ha='right',
                    fontsize=font_sizes['annotation'], fontweight='bold')

    plt.tight_layout()
    plt.subplots_adjust(left=0.3)

    return fig


def create_gene_bar_chart(de_results, pathway_name, pathway_genes,
                          disease_label, control_label, label='Significant',
                          top_n=15, figsize=None, font_sizes=None):
    """Create a per-pathway top DE gene horizontal bar chart.

    Parameters
    ----------
    de_results : pandas.DataFrame
        Must have 'names', 'logfoldchanges', 'pvals_adj', 'abs_logfoldchange'.
    pathway_name : str
        Pathway name for title.
    pathway_genes : list of str
        Genes in this pathway.
    disease_label, control_label : str
        Condition labels.
    label : str
        Label (e.g. 'Significant', 'All Genes').
    top_n : int
        Max genes to show.
    figsize : tuple, optional

    Returns
    -------
    matplotlib.figure.Figure or None
        None if no genes found.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    genes_upper = {g.upper() for g in pathway_genes}
    pw_df = de_results[de_results['names'].str.upper().isin(genes_upper)].copy()
    if len(pw_df) == 0:
        return None

    pw_df = pw_df.sort_values('abs_logfoldchange', ascending=False).head(top_n)
    pw_df = pw_df.sort_values('abs_logfoldchange', ascending=True)

    n_genes = len(pw_df)
    fig_h = max(4, n_genes * 0.4 + 2)
    if figsize is None:
        figsize = (10, fig_h)

    fig, ax = plt.subplots(figsize=figsize)
    y = np.arange(n_genes)
    lfc_vals = pw_df['logfoldchanges'].values
    colors = ['#d62728' if lfc > 0 else '#1f77b4' for lfc in lfc_vals]

    ax.barh(y, np.abs(lfc_vals), height=0.6, color=colors, alpha=0.85,
            edgecolor='black', linewidth=0.5)

    ax.set_yticks(y)
    ax.set_yticklabels(pw_df['names'].values, fontsize=font_sizes['tick'])
    ax.set_xlabel('|log2 Fold Change|', fontsize=font_sizes['axis_label'], fontweight='bold')
    ax.set_title(f'{pathway_name.replace("_", " ")} — Top DE Genes ({label})\n'
                 f'{disease_label} vs {control_label}',
                 fontsize=font_sizes['title'], fontweight='bold')

    for i, (_, row) in enumerate(pw_df.iterrows()):
        lfc = row['logfoldchanges']
        pval_text = ''
        if 'pvals_adj' in row.index and pd.notna(row['pvals_adj']):
            pval_text = f'  (padj={row["pvals_adj"]:.2e})'
        ax.text(np.abs(lfc) + 0.02, i, f'log2FC={lfc:+.2f}{pval_text}',
                va='center', fontsize=font_sizes['annotation'], fontweight='bold',
                color='#d62728' if lfc > 0 else '#1f77b4')

    legend_elements = [
        Patch(facecolor='#d62728', alpha=0.85, label=f'Up in {disease_label}'),
        Patch(facecolor='#1f77b4', alpha=0.85, label=f'Down in {disease_label}'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=font_sizes['annotation'])
    plt.tight_layout()

    return fig
