"""
Gene Heatmap
=============
Z-score heatmap of gene expression per pathway with condition annotation.
"""

import numpy as np
import pandas as pd

from kosmic.scrna.inspect.detection import detect_species, format_gene_for_species
from kosmic import DE_MIN_CELLS


def create_gene_heatmap(sample_gene_df, pathway_genes, pathway_name,
                        control_label, disease_label, figsize=None,
                        cmap='viridis', vmin_pct=5, vmax_pct=95,
                        theme_colors=None, font_sizes=None):
    """Create a z-score heatmap for a single pathway.

    Parameters
    ----------
    sample_gene_df : pandas.DataFrame
        Must have 'sample', 'condition' columns plus one column per gene with
        expression values (already z-scored or raw — z-scoring applied if no negatives).
    pathway_genes : list of str
        Genes to include (must be columns in sample_gene_df).
    pathway_name : str
        Pathway name for title.
    control_label, disease_label : str
        Condition labels.
    figsize : tuple, optional
        Figure size. Auto-computed if None.
    cmap : str
        Matplotlib colormap name.
    vmin_pct, vmax_pct : float
        Percentile clipping for the color scale (0–100).
    theme_colors : dict, optional
        Theme color dict with keys like 'bg_primary', 'fg_primary', etc.
        If None, uses default light styling.

    Returns
    -------
    matplotlib.figure.Figure or None
        None if fewer than 3 genes available.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import seaborn as sns

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    available_genes = [g for g in pathway_genes if g in sample_gene_df.columns]
    if len(available_genes) < 3:
        return None

    # Theme colors
    if theme_colors:
        bg = theme_colors.get('bg_primary', '#ffffff')
        fg = theme_colors.get('fg_primary', '#000000')
        fg2 = theme_colors.get('fg_secondary', '#666666')
    else:
        bg, fg, fg2 = '#ffffff', '#000000', '#666666'

    df = sample_gene_df.copy()

    # Z-score if data is all non-negative
    if not np.any(df[available_genes].values < 0):
        for gene in available_genes:
            mean_v = df[gene].mean()
            std_v = df[gene].std()
            if std_v > 0:
                df[gene] = (df[gene] - mean_v) / std_v
            else:
                df[gene] = 0

    # Build heatmap matrix (genes x samples)
    heatmap_data = df.set_index('sample')[available_genes].T

    # Sort samples by condition
    control_samples = df[df['condition'] == control_label]['sample'].tolist()
    disease_samples = df[df['condition'] == disease_label]['sample'].tolist()
    sample_order = control_samples + disease_samples
    heatmap_data = heatmap_data[[s for s in sample_order if s in heatmap_data.columns]]

    if figsize is None:
        fig_height = max(6, len(available_genes) * 0.4)
        fig_width = max(12, len(sample_order) * 0.3)
        figsize = (fig_width, fig_height)

    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor(bg)
    ax.set_facecolor(bg)

    vmin = np.percentile(heatmap_data.values, vmin_pct)
    vmax = np.percentile(heatmap_data.values, vmax_pct)

    # Default font sizes if not provided
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    cbar_kws = {'label': 'Z-score Expression'}
    sns.heatmap(heatmap_data, cmap=cmap, vmin=vmin, vmax=vmax,
                cbar_kws=cbar_kws, xticklabels=True, yticklabels=True, ax=ax)

    # Style colorbar
    cbar = ax.collections[0].colorbar
    if cbar is not None:
        cbar.ax.yaxis.label.set_color(fg)
        cbar.ax.yaxis.label.set_fontsize(font_sizes['colorbar_label'])
        cbar.ax.tick_params(colors=fg2, labelsize=font_sizes['colorbar_tick'])
        for spine in cbar.ax.spines.values():
            spine.set_edgecolor(fg2)

    # Condition color bar at top
    control_color = (theme_colors.get('plot_control', '#3498DB')
                     if theme_colors else '#3498DB')
    disease_color = (theme_colors.get('plot_disease', '#E74C3C')
                     if theme_colors else '#E74C3C')
    condition_map = df.set_index('sample')['condition'].to_dict()
    color_map = {disease_label: disease_color, control_label: control_color}

    for i, sample in enumerate(heatmap_data.columns):
        cond = condition_map.get(sample, '')
        color = color_map.get(cond, '#999999')
        ax.add_patch(plt.Rectangle((i, len(available_genes)), 1, 0.3,
                                   color=color, clip_on=False))

    legend_elements = [
        Patch(facecolor=control_color, label=control_label),
        Patch(facecolor=disease_color, label=disease_label),
    ]
    legend = ax.legend(handles=legend_elements, loc='upper left',
                       bbox_to_anchor=(0, 1.08), ncol=2, frameon=True,
                       fontsize=font_sizes['legend'])
    legend.get_frame().set_facecolor(bg)
    legend.get_frame().set_edgecolor(fg2)
    for text in legend.get_texts():
        text.set_color(fg)

    ax.set_title(
        f'{pathway_name.replace("_", " ")} Gene Expression\n'
        f'({control_label} → {disease_label})',
        fontsize=font_sizes['title'], fontweight='bold', color=fg, pad=20)
    ax.set_xlabel('Samples', fontsize=font_sizes['axis_label'], color=fg)
    ax.set_ylabel('Genes', fontsize=font_sizes['axis_label'], color=fg)
    ax.tick_params(axis='x', rotation=45, colors=fg2, labelsize=font_sizes['tick'])
    ax.tick_params(axis='y', colors=fg2, labelsize=font_sizes['tick'])
    plt.setp(ax.get_xticklabels(), ha='right')

    for spine in ax.spines.values():
        spine.set_edgecolor(fg2)

    fig.tight_layout()

    return fig


def prepare_heatmap_data(adata, pathway_gene_sets, sample_col, condition_col,
                         min_cells=DE_MIN_CELLS, allowed_genes=None):
    """Prepare sample-level gene expression data for heatmaps.

    Parameters
    ----------
    adata : anndata.AnnData
        Input data.
    pathway_gene_sets : dict
        {pathway_name: [gene_list]}.
    sample_col : str
        Sample grouping column.
    condition_col : str
        Condition column.
    min_cells : int
        Minimum cells per sample.
    allowed_genes : set of str, optional
        When given, only these genes are kept per pathway. Callers pass the
        DE-tested gene set so the heatmap shows only genes the DE window kept
        (e.g. drops a one-donor gene) -- applied here so the in-app and exported
        heatmaps stay identical.

    Returns
    -------
    tuple of (sample_df, available_pathways)
        sample_df: DataFrame with 'sample', 'condition', + gene columns.
        available_pathways: dict of {pathway: [available_genes]}.
    """
    import scanpy as sc

    species = detect_species(list(
        adata.raw.var_names if adata.raw is not None else adata.var_names
    ))
    coverage_var = set(adata.raw.var_names if adata.raw is not None else adata.var_names)

    # Find available pathways
    available_pathways = {}
    for pathway_name, genes in pathway_gene_sets.items():
        genes_fmt = [format_gene_for_species(g, species) for g in genes]
        available = [g for g in genes_fmt if g in coverage_var]
        if allowed_genes is not None:
            available = [g for g in available if g in allowed_genes]
        if len(available) >= 3:
            available_pathways[pathway_name] = available

    # Collect all genes
    all_genes = list(dict.fromkeys(g for gs in available_pathways.values() for g in gs))

    # Use raw data, normalize for consistent expression values
    use_raw = adata.raw is not None
    if use_raw:
        valid_genes = [g for g in all_genes if g in adata.raw.var_names]
        adata_expr = adata.raw[:, valid_genes].to_adata()
        adata_expr.obs = adata.obs.copy()
        sc.pp.normalize_total(adata_expr, target_sum=1e4)
        sc.pp.log1p(adata_expr)
        all_genes = list(adata_expr.var_names)
        var_names = adata_expr.var_names
    else:
        adata_expr = adata
        all_genes = [g for g in all_genes if g in adata.var_names]
        var_names = adata.var_names

    # Build sample-level data
    sample_data = []
    for sample_id in adata.obs[sample_col].unique():
        sample_mask = adata.obs[sample_col] == sample_id
        sample_cells = adata_expr[sample_mask]
        X = sample_cells.X.toarray() if hasattr(sample_cells.X, 'toarray') else np.array(sample_cells.X)

        if X.shape[0] < min_cells:
            continue

        condition = adata.obs.loc[sample_mask, condition_col].iloc[0]
        row = {'sample': sample_id, 'condition': condition}

        for gene in all_genes:
            if gene in var_names:
                gene_idx = var_names.get_loc(gene)
                row[gene] = np.mean(X[:, gene_idx])

        sample_data.append(row)

    sample_df = pd.DataFrame(sample_data)
    # Filter available_pathways to genes actually in DataFrame
    for pw in list(available_pathways.keys()):
        available_pathways[pw] = [g for g in available_pathways[pw] if g in sample_df.columns]
        if len(available_pathways[pw]) < 3:
            del available_pathways[pw]

    return sample_df, available_pathways
