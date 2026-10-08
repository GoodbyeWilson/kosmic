"""
Top DE Plots
=============
Pathway summary split bar chart and per-pathway gene bar charts.
"""

import numpy as np
import pandas as pd

from kosmic.visualisation.style import DOWN_COLOR, FG, UP_COLOR


def draw_split_bars(ax, labels, n_up, n_down, disease_label, xlabel, font_sizes):
    """Up counts to the right of zero, down counts to the left, one row
    per label (first row at the bottom). The x axis shows counts on both
    sides."""
    from matplotlib.ticker import FuncFormatter

    y = np.arange(len(labels))
    ax.barh(y, n_up, height=0.65, color=UP_COLOR, label=f'Up in {disease_label}')
    ax.barh(y, -np.asarray(n_down), height=0.65, color=DOWN_COLOR,
            label=f'Down in {disease_label}')
    ax.axvline(0, color=FG, linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=font_sizes['tick'])
    ax.set_ylim(-0.6, len(labels) - 0.4)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{abs(v):g}'))
    ax.set_xlabel(xlabel, fontsize=font_sizes['axis_label'])
    ax.tick_params(labelsize=font_sizes['tick'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False,
              fontsize=font_sizes['legend'])


def create_pathway_summary_chart(de_results, pathway_gene_sets, disease_label, control_label,
                                 label='Significant', figsize=None, font_sizes=None,
                                 pathway_coverage=None):
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
    pathway_coverage : dict, optional
        {pathway: {available_genes, total_genes, ...}}. When given, each
        label ends with the gene set's coverage: genes measured in the
        dataset / genes in the set.
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
            name = pathway_name.replace('_', ' ')
            cov = (pathway_coverage or {}).get(pathway_name)
            if cov:
                name += f" ({cov['available_genes']}/{cov['total_genes']})"
            pathway_stats.append({
                'pathway': name,
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
    draw_split_bars(ax, pw_df['pathway'].values, pw_df['n_up'].values,
                    pw_df['n_down'].values, disease_label,
                    f'DE genes ({label.lower()})', font_sizes)
    ax.set_title(f'DE genes by pathway ({label.lower()})\n{disease_label} vs {control_label}',
                 fontsize=font_sizes['title'], fontweight='bold')
    fig.tight_layout()
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
    colors = [UP_COLOR if lfc > 0 else DOWN_COLOR for lfc in lfc_vals]

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
                color=UP_COLOR if lfc > 0 else DOWN_COLOR)

    legend_elements = [
        Patch(facecolor=UP_COLOR, label=f'Up in {disease_label}'),
        Patch(facecolor=DOWN_COLOR, label=f'Down in {disease_label}'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=font_sizes['annotation'])
    plt.tight_layout()

    return fig
