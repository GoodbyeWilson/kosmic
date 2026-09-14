"""
Heatmap + Forest Plot Meta-Analysis Visualization
==================================================
Creates publication-ready figures combining a heatmap of log fold changes
across multiple datasets with a forest plot showing pooled effect sizes
and confidence intervals from meta-analysis.

Inspired by Figure 2J from PMC10521807.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from pathlib import Path
from typing import List
import warnings
from kosmic import DEFAULT_FDR
warnings.filterwarnings('ignore')










def create_multi_dataset_figure(meta_df: pd.DataFrame,
                                 datasets: List[pd.DataFrame],
                                 dataset_labels: List[str],
                                 output_path: str = None,
                                 title: str = "Meta-Analysis of Differential Expression",
                                 condition_labels: tuple = ("Control", "Disease"),
                                 group_by_pathway: bool = True,
                                 figsize: tuple = None,
                                 cmap: str = 'RdBu_r',
                                 show_individual_studies: bool = True):
    """
    Create heatmap + forest plot for multiple datasets.

    Layout:
    [Pathway labels | Heatmap (datasets as columns) | Forest plot with pooled + individual effects]
    """

    # Sort and group
    if group_by_pathway and 'pathways' in meta_df.columns:
        meta_df['pathway_clean'] = meta_df['pathways'].fillna('Other')
        meta_df = meta_df.sort_values(['pathway_clean', 'logfoldchanges'], ascending=[True, False])
    else:
        meta_df = meta_df.sort_values('logfoldchanges', ascending=False)

    n_genes = len(meta_df)
    n_datasets = len(datasets)

    # Build heatmap matrix
    heatmap_matrix = np.full((n_genes, n_datasets), np.nan)
    gene_order = meta_df['names'].tolist()

    for j, df in enumerate(datasets):
        for i, gene in enumerate(gene_order):
            gene_rows = df[df['names'] == gene]
            if len(gene_rows) > 0:
                heatmap_matrix[i, j] = gene_rows.iloc[0]['logfoldchanges']

    # Auto-size figure
    if figsize is None:
        height = max(8, n_genes * 0.35 + 3)
        width = 6 + n_datasets * 0.8 + (4 if show_individual_studies else 2)
        figsize = (width, height)

    fig = plt.figure(figsize=figsize)

    # Layout ratios
    pathway_ratio = 0.6 if group_by_pathway and 'pathways' in meta_df.columns else 0
    heatmap_ratio = max(1.5, n_datasets * 0.5)
    forest_ratio = 3.5 if show_individual_studies else 2.5
    cbar_ratio = 0.2

    if pathway_ratio > 0:
        gs = GridSpec(1, 4, width_ratios=[pathway_ratio, heatmap_ratio, forest_ratio, cbar_ratio], wspace=0.08)
        ax_pathway = fig.add_subplot(gs[0])
        ax_heatmap = fig.add_subplot(gs[1])
        ax_forest = fig.add_subplot(gs[2])
        ax_cbar = fig.add_subplot(gs[3])
    else:
        gs = GridSpec(1, 3, width_ratios=[heatmap_ratio, forest_ratio, cbar_ratio], wspace=0.08)
        ax_pathway = None
        ax_heatmap = fig.add_subplot(gs[0])
        ax_forest = fig.add_subplot(gs[1])
        ax_cbar = fig.add_subplot(gs[2])

    genes = meta_df['names'].values
    pooled_logfc = meta_df['logfoldchanges'].values
    ci_lower = meta_df['ci_lower'].values
    ci_upper = meta_df['ci_upper'].values

    # Significance
    if 'pvals_pooled' in meta_df.columns:
        significant = meta_df['pvals_pooled'].values < DEFAULT_FDR
    elif 'pvals_adj' in meta_df.columns:
        significant = meta_df['pvals_adj'].values < DEFAULT_FDR
    else:
        significant = np.ones(n_genes, dtype=bool)

    # Color scale
    vmax = np.nanmax(np.abs(heatmap_matrix))
    vmax = max(vmax, np.max(np.abs(pooled_logfc)), 1)
    vmin = -vmax

    # ===== PATHWAY LABELS =====
    if ax_pathway is not None:
        pathway_labels = meta_df['pathway_clean'].values
        unique_pathways = []
        pathway_positions = []
        current_pathway = None
        start_pos = 0

        for i, pw in enumerate(pathway_labels):
            if pw != current_pathway:
                if current_pathway is not None:
                    unique_pathways.append(current_pathway)
                    pathway_positions.append((start_pos, i - 1))
                current_pathway = pw
                start_pos = i
        unique_pathways.append(current_pathway)
        pathway_positions.append((start_pos, n_genes - 1))

        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_pathways)))
        for idx, ((start, end), pw) in enumerate(zip(pathway_positions, unique_pathways)):
            rect = mpatches.Rectangle((0, start - 0.5), 1, end - start + 1,
                                       facecolor=colors[idx], edgecolor='white', linewidth=0.5)
            ax_pathway.add_patch(rect)
            mid = (start + end) / 2
            pw_display = pw[:20] + '...' if len(pw) > 20 else pw
            ax_pathway.text(0.5, mid, pw_display, ha='center', va='center',
                           fontsize=7, rotation=90, fontweight='bold')

        ax_pathway.set_xlim(0, 1)
        ax_pathway.set_ylim(-0.5, n_genes - 0.5)
        ax_pathway.set_xticks([])
        ax_pathway.set_yticks([])
        ax_pathway.invert_yaxis()
        ax_pathway.set_title('Pathway', fontsize=9, fontweight='bold')
        for spine in ax_pathway.spines.values():
            spine.set_visible(False)

    # ===== HEATMAP =====
    im = ax_heatmap.imshow(heatmap_matrix, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)

    ax_heatmap.set_yticks(range(n_genes))
    ax_heatmap.set_yticklabels(genes, fontsize=7)
    ax_heatmap.set_xticks(range(n_datasets))

    xlabels = [lbl[:12] + '..' if len(lbl) > 12 else lbl for lbl in dataset_labels]
    ax_heatmap.set_xticklabels(xlabels, fontsize=8, rotation=45, ha='right')
    ax_heatmap.set_title('Log2FC by Dataset', fontsize=9, fontweight='bold')

    for i in range(n_genes):
        if significant[i]:
            for j in range(n_datasets):
                if not np.isnan(heatmap_matrix[i, j]):
                    val = heatmap_matrix[i, j]
                    color = 'white' if abs(val) > vmax * 0.5 else 'black'
                    ax_heatmap.text(j, i, '*', ha='center', va='center',
                                   fontsize=8, fontweight='bold', color=color)

    # ===== FOREST PLOT =====
    study_colors = plt.cm.tab10(np.linspace(0, 1, n_datasets))

    for i in range(n_genes):
        if show_individual_studies and 'study_effects' in meta_df.columns:
            study_effects = meta_df.iloc[i]['study_effects']
            if isinstance(study_effects, list):
                for k, effect in enumerate(study_effects):
                    x = effect['logfc']
                    se = effect.get('se', 0.3)
                    if np.isfinite(x) and np.isfinite(se):
                        ax_forest.hlines(i, x - 1.96*se, x + 1.96*se,
                                        color=study_colors[k % len(study_colors)],
                                        alpha=0.5, linewidth=1)
                        ax_forest.scatter(x, i, color=study_colors[k % len(study_colors)],
                                         s=20, marker='o', alpha=0.6, zorder=2)

        color = '#d62728' if pooled_logfc[i] > 0 else '#1f77b4'
        alpha = 1.0 if significant[i] else 0.4

        ax_forest.hlines(i, ci_lower[i], ci_upper[i], color=color, alpha=alpha, linewidth=2.5)

        cap_height = 0.2
        ax_forest.vlines(ci_lower[i], i - cap_height, i + cap_height, color=color, alpha=alpha, linewidth=1.5)
        ax_forest.vlines(ci_upper[i], i - cap_height, i + cap_height, color=color, alpha=alpha, linewidth=1.5)

        marker = 'D' if significant[i] else 'o'
        ax_forest.scatter(pooled_logfc[i], i, color=color, s=80, marker=marker,
                         zorder=4, alpha=alpha, edgecolors='black', linewidth=0.8)

    ax_forest.axvline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)

    ax_forest.set_ylim(-0.5, n_genes - 0.5)
    ax_forest.invert_yaxis()
    ax_forest.set_yticks([])
    ax_forest.set_xlabel(f'Log2 Fold Change\n({condition_labels[1]} vs {condition_labels[0]})', fontsize=9)
    ax_forest.set_title('Pooled Effect Size (95% CI)', fontsize=9, fontweight='bold')
    ax_forest.spines['top'].set_visible(False)
    ax_forest.spines['right'].set_visible(False)
    ax_forest.spines['left'].set_visible(False)
    ax_forest.grid(axis='x', alpha=0.3, linestyle=':')

    # Legend
    legend_elements = [
        plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='#d62728',
                   markersize=8, markeredgecolor='black', label='Upregulated (pooled, p<0.05)'),
        plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='#1f77b4',
                   markersize=8, markeredgecolor='black', label='Downregulated (pooled, p<0.05)'),
    ]

    if show_individual_studies and n_datasets > 1:
        for k, label in enumerate(dataset_labels[:min(5, n_datasets)]):
            legend_elements.append(
                plt.Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=study_colors[k], markersize=6, alpha=0.6,
                          label=label[:15])
            )

    ax_forest.legend(handles=legend_elements, loc='upper right', fontsize=7, framealpha=0.9)

    # ===== COLORBAR =====
    cbar = plt.colorbar(im, cax=ax_cbar)
    cbar.set_label('Log2 FC', fontsize=8)

    # Title
    n_sig = np.sum(significant)
    subtitle = f'{n_genes} genes, {n_sig} significant (p<0.05)'
    if n_datasets > 1:
        subtitle += f', {n_datasets} datasets'

    fig.suptitle(f'{title}\n{subtitle}', fontsize=11, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Save
    if output_path:
        output_path = Path(output_path)
        plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')

        if output_path.suffix.lower() == '.png':
            pdf_path = output_path.with_suffix('.pdf')
            plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    else:
        plt.show()

    plt.close()
    return fig
