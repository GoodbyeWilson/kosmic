"""
Variable Gene Plots
====================
Two-panel plot: coefficient of variation bars + log2FC dots for top variable genes.
"""

import numpy as np


def compute_variable_genes(gene_data_df, disease_label, control_label,
                           gene_to_subfamily=None, top_n=20):
    """Compute coefficient of variation and log2FC for genes.

    Parameters
    ----------
    gene_data_df : pandas.DataFrame
        Must have columns: Gene, Sample, Condition, Mean_Expression.
    disease_label, control_label : str
        Condition labels.
    gene_to_subfamily : dict, optional
        {gene: subfamily_name} for coloring.
    top_n : int
        Number of top variable genes to return.

    Returns
    -------
    pandas.DataFrame or None
        Columns: Gene, Mean, Std, CV, Disease_Mean, Control_Mean, Log2FC, Subfamily.
        None if fewer than 2 genes with expression.
    """
    gene_cv = gene_data_df.groupby('Gene')['Mean_Expression'].agg(['mean', 'std']).reset_index()
    gene_cv.columns = ['Gene', 'Mean', 'Std']
    gene_cv['CV'] = gene_cv['Std'] / (gene_cv['Mean'] + 1e-10)

    # Per-condition means
    cond_means = gene_data_df.groupby(['Gene', 'Condition'])['Mean_Expression'].mean().unstack(fill_value=0)
    if disease_label in cond_means.columns and control_label in cond_means.columns:
        gene_cv['Disease_Mean'] = gene_cv['Gene'].map(cond_means[disease_label].to_dict()).fillna(0)
        gene_cv['Control_Mean'] = gene_cv['Gene'].map(cond_means[control_label].to_dict()).fillna(0)
        gene_cv['Log2FC'] = np.log2((gene_cv['Disease_Mean'] + 0.01) / (gene_cv['Control_Mean'] + 0.01))
    else:
        gene_cv['Log2FC'] = 0.0

    if gene_to_subfamily:
        gene_cv['Subfamily'] = gene_cv['Gene'].map(gene_to_subfamily).fillna('')
    else:
        gene_cv['Subfamily'] = ''

    # Filter and select top N
    gene_cv = gene_cv[gene_cv['Mean'] > 0.01]
    top_genes = gene_cv.nlargest(min(top_n, len(gene_cv)), 'CV')

    if len(top_genes) < 2:
        return None

    return top_genes.sort_values('CV', ascending=True)


def create_variable_gene_plot(top_genes_df, parent_name, disease_label, control_label,
                              figsize=None, font_sizes=None):
    """Create a 2-panel variable gene plot: CV bars + log2FC dots.

    Parameters
    ----------
    top_genes_df : pandas.DataFrame
        From compute_variable_genes. Must have: Gene, CV, Log2FC, Subfamily.
    parent_name : str
        Parent category name for title.
    disease_label, control_label : str
        Condition labels.
    figsize : tuple, optional

    Returns
    -------
    matplotlib.figure.Figure
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

    n = len(top_genes_df)
    if figsize is None:
        figsize = (16, max(6, n * 0.35))

    fig, axes = plt.subplots(1, 2, figsize=figsize, gridspec_kw={'width_ratios': [2, 1]})

    # Left panel: CV bar plot
    ax1 = axes[0]
    subfamilies = top_genes_df['Subfamily'].unique()
    cmap = plt.cm.get_cmap('tab20', max(len(subfamilies), 1))
    subfamily_colors = {sf: cmap(i) for i, sf in enumerate(subfamilies)}

    ax1.barh(range(n), top_genes_df['CV'].values,
             color=[subfamily_colors[sf] for sf in top_genes_df['Subfamily']],
             alpha=0.8, edgecolor='black', linewidth=0.5)
    ax1.set_yticks(range(n))
    ax1.set_yticklabels(top_genes_df['Gene'].values, fontsize=font_sizes['tick'])
    ax1.set_xlabel('Coefficient of Variation', fontsize=font_sizes['axis_label'], fontweight='bold')
    ax1.set_title(f'Top Variable Genes — {parent_name.replace("_", " ")}',
                  fontsize=font_sizes['title'], fontweight='bold')

    if len(subfamilies) > 1 or (len(subfamilies) == 1 and subfamilies[0]):
        legend_patches = [Patch(facecolor=subfamily_colors[sf], label=sf) for sf in subfamilies if sf]
        if legend_patches:
            ax1.legend(handles=legend_patches, loc='lower right',
                       fontsize=font_sizes['annotation'],
                       title='Subfamily', title_fontsize=font_sizes['tick'])

    # Right panel: Log2FC dot plot
    ax2 = axes[1]
    colors_fc = ['#d62728' if fc > 0 else '#1f77b4' for fc in top_genes_df['Log2FC']]
    ax2.scatter(top_genes_df['Log2FC'].values, range(n),
                c=colors_fc, s=60, alpha=0.8, edgecolors='black', linewidth=0.5)
    ax2.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax2.set_yticks(range(n))
    ax2.set_yticklabels([])
    ax2.set_xlabel(f'Log2 FC ({disease_label}/{control_label})',
                   fontsize=font_sizes['axis_label'], fontweight='bold')
    ax2.set_title('Direction of Change', fontsize=font_sizes['title'], fontweight='bold')

    plt.tight_layout()
    return fig
