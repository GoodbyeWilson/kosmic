"""
Dataset Summary
================
Compute pathway-level summary stats and create a colored table figure.
"""

import numpy as np
import pandas as pd
from kosmic import DEFAULT_FDR


def compute_dataset_summary(de_results, pathway_gene_sets, pathway_coverage):
    """Compute pathway-level summary statistics from DE results.

    Parameters
    ----------
    de_results : pandas.DataFrame
        DE results with 'names', 'logfoldchanges', 'pvals_adj' columns.
    pathway_gene_sets : dict
        {pathway_name: [gene_list]}.
    pathway_coverage : dict
        From prepare_gene_coverage; {pathway: {found, total_genes, ...}}.

    Returns
    -------
    pandas.DataFrame
        Columns: pathway, direction, mean_log2fc, mean_padj, neg_log10_padj,
        sig_genes_up, sig_genes_down, total_sig_genes, total_genes_tested,
        total_genes_defined, coverage_pct.
    """
    rows = []
    for pathway_name, genes in pathway_gene_sets.items():
        cov_info = pathway_coverage.get(pathway_name, {})
        genes_found = cov_info.get('found', 0) if isinstance(cov_info, dict) else 0
        genes_defined = len(genes)
        coverage_pct = (genes_found / genes_defined * 100) if genes_defined > 0 else 0

        pathway_genes_upper = [g.upper() for g in genes]
        if de_results is not None and len(de_results) > 0:
            pw_de = de_results[de_results['names'].str.upper().isin(pathway_genes_upper)]
            sig_mask = pw_de['pvals_adj'] < DEFAULT_FDR
            sig_up = int((sig_mask & (pw_de['logfoldchanges'] > 0)).sum())
            sig_down = int((sig_mask & (pw_de['logfoldchanges'] < 0)).sum())
            total_tested = len(pw_de)
        else:
            sig_up = sig_down = total_tested = 0

        total_sig = sig_up + sig_down

        mean_log2fc = np.nan
        mean_padj = np.nan
        if de_results is not None and len(de_results) > 0:
            pw_sig = pw_de[pw_de['pvals_adj'] < DEFAULT_FDR] if total_tested > 0 else pd.DataFrame()
            if len(pw_sig) > 0:
                mean_log2fc = pw_sig['logfoldchanges'].mean()
                mean_padj = pw_sig['pvals_adj'].mean()

        if total_sig > 0 and pd.notna(mean_log2fc):
            direction = "UP" if mean_log2fc > 0 else "DOWN"
        elif total_sig > 0:
            direction = "UP" if sig_up > sig_down else ("DOWN" if sig_down > sig_up else "NS")
        else:
            direction = "NS"

        neg_log10_padj = -np.log10(mean_padj) if pd.notna(mean_padj) and mean_padj > 0 else 0.0

        rows.append({
            'pathway': pathway_name,
            'direction': direction,
            'mean_log2fc': round(mean_log2fc, 4) if pd.notna(mean_log2fc) else np.nan,
            'mean_padj': round(mean_padj, 6) if pd.notna(mean_padj) else np.nan,
            'neg_log10_padj': round(neg_log10_padj, 4),
            'sig_genes_up': sig_up,
            'sig_genes_down': sig_down,
            'total_sig_genes': total_sig,
            'total_genes_tested': total_tested,
            'total_genes_defined': genes_defined,
            'coverage_pct': round(coverage_pct, 1),
        })

    summary_df = pd.DataFrame(rows)
    summary_df['_abs_lfc'] = summary_df['mean_log2fc'].abs()
    summary_df = summary_df.sort_values('_abs_lfc', ascending=False, na_position='last')
    summary_df = summary_df.drop(columns=['_abs_lfc'])
    return summary_df


def create_summary_figure(summary_df, dataset_name, disease_label, control_label,
                          geneset_label, figsize=None, font_sizes=None):
    """Create a colored table figure from summary data.

    Parameters
    ----------
    summary_df : pandas.DataFrame
        From compute_dataset_summary.
    dataset_name : str
        Dataset identifier.
    disease_label, control_label : str
        Condition labels.
    geneset_label : str
        Gene set label (e.g. 'Metabolic', 'GPCR').
    figsize : tuple, optional

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    label_lower = geneset_label.lower()
    is_family = 'gpcr' in label_lower or 'ion' in label_lower
    group_noun = 'families' if is_family else 'pathways'
    group_singular = 'Family' if is_family else 'Pathway'

    sig_df = summary_df[summary_df['direction'].isin(['UP', 'DOWN'])].copy()

    if len(sig_df) == 0:
        fig, ax = plt.subplots(figsize=figsize or (8, 2))
        ax.axis('off')
        title = f'Dataset Summary: {disease_label} vs {control_label}'
        subtitle = f'{dataset_name} | {geneset_label}'
        ax.set_title(f'{title}\n{subtitle}', fontsize=font_sizes['tick'], fontweight='bold', pad=10)
        ax.text(0.5, 0.4, 'No significant changes detected',
                ha='center', va='center', fontsize=font_sizes['axis_label'], color='#999999',
                transform=ax.transAxes)
        plt.tight_layout()
        return fig

    up_df = sig_df[sig_df['direction'] == 'UP'].sort_values('total_sig_genes', ascending=False)
    down_df = sig_df[sig_df['direction'] == 'DOWN'].sort_values('total_sig_genes', ascending=False)

    display_rows = []
    if len(up_df) > 0:
        display_rows.append(('header', f'UPREGULATED ({len(up_df)})'))
        for _, r in up_df.iterrows():
            display_rows.append(('row', r))
    if len(down_df) > 0:
        if len(up_df) > 0:
            display_rows.append(('spacer', None))
        display_rows.append(('header', f'DOWNREGULATED ({len(down_df)})'))
        for _, r in down_df.iterrows():
            display_rows.append(('row', r))

    n_display = len(display_rows)
    fig_height = max(3, 1.0 + n_display * 0.4)
    if figsize is None:
        figsize = (10, fig_height)

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.5, n_display - 0.5)
    ax.invert_yaxis()
    ax.axis('off')

    col_x = [0.3, 4.2, 5.5, 6.8, 8.2, 9.2]
    col_labels = [group_singular, 'Mean log2FC', '-log10(padj)', 'Sig Up', 'Sig Down', 'Total']

    for x, label in zip(col_x, col_labels):
        ax.text(x, -0.7, label, fontsize=font_sizes['annotation'], fontweight='bold',
                ha='left' if x < 4 else 'center', va='center')

    rdbu = plt.cm.RdBu_r
    oranges = plt.cm.Oranges

    lfc_vals = sig_df['mean_log2fc'].dropna()
    lfc_max = max(abs(lfc_vals.max()), abs(lfc_vals.min())) if len(lfc_vals) > 0 else 1
    lfc_norm = mcolors.TwoSlopeNorm(vmin=-lfc_max, vcenter=0, vmax=lfc_max) if lfc_max > 0 else None

    nlp_vals = sig_df['neg_log10_padj']
    nlp_max = nlp_vals.max() if nlp_vals.max() > 0 else 1
    nlp_norm = mcolors.Normalize(vmin=0, vmax=nlp_max)

    for i, (row_type, data) in enumerate(display_rows):
        y = i
        if row_type == 'header':
            color = '#d62728' if 'UP' in data else '#1f77b4'
            rect = plt.Rectangle((0.0, y - 0.35), 10, 0.7,
                                 facecolor=color, alpha=0.12, edgecolor='none')
            ax.add_patch(rect)
            ax.text(0.15, y, data, fontsize=font_sizes['annotation'], fontweight='bold',
                    va='center', ha='left', color=color)
            continue
        if row_type == 'spacer':
            continue

        row = data
        name = row['pathway'].replace('_', ' ')
        if len(name) > 40:
            name = name[:37] + '...'
        ax.text(col_x[0], y, name, fontsize=font_sizes['colorbar_tick'], va='center', ha='left')

        lfc = row['mean_log2fc']
        if pd.notna(lfc) and lfc_norm is not None:
            bg = rdbu(lfc_norm(lfc))
            rect = plt.Rectangle((col_x[1] - 0.4, y - 0.35), 0.8, 0.7,
                                 facecolor=bg, edgecolor='#cccccc', linewidth=0.5)
            ax.add_patch(rect)
            text_color = 'white' if abs(lfc) > lfc_max * 0.6 else 'black'
            ax.text(col_x[1], y, f'{lfc:.2f}', fontsize=font_sizes['colorbar_tick'], va='center',
                    ha='center', color=text_color)
        else:
            ax.text(col_x[1], y, 'N/A', fontsize=font_sizes['colorbar_tick'], va='center',
                    ha='center', color='#999999')

        nlp = row['neg_log10_padj']
        if nlp > 0:
            bg = oranges(nlp_norm(nlp))
            rect = plt.Rectangle((col_x[2] - 0.4, y - 0.35), 0.8, 0.7,
                                 facecolor=bg, edgecolor='#cccccc', linewidth=0.5)
            ax.add_patch(rect)
            text_color = 'white' if nlp > nlp_max * 0.6 else 'black'
            ax.text(col_x[2], y, f'{nlp:.1f}', fontsize=font_sizes['colorbar_tick'], va='center',
                    ha='center', color=text_color)
        else:
            ax.text(col_x[2], y, '0.0', fontsize=font_sizes['colorbar_tick'], va='center',
                    ha='center', color='#999999')

        ax.text(col_x[3], y, str(row['sig_genes_up']), fontsize=font_sizes['colorbar_tick'],
                va='center', ha='center',
                color='#d62728' if row['sig_genes_up'] > 0 else '#999999')
        ax.text(col_x[4], y, str(row['sig_genes_down']), fontsize=font_sizes['colorbar_tick'],
                va='center', ha='center',
                color='#1f77b4' if row['sig_genes_down'] > 0 else '#999999')
        ax.text(col_x[5], y, str(row['total_sig_genes']), fontsize=font_sizes['colorbar_tick'],
                va='center', ha='center', fontweight='bold')

    n_total = len(summary_df)
    n_sig = len(sig_df)
    title = f'Dataset Summary: {disease_label} vs {control_label}'
    subtitle = f'{dataset_name} | {geneset_label} | {n_sig}/{n_total} {group_noun} with significant changes'
    ax.set_title(f'{title}\n{subtitle}', fontsize=font_sizes['tick'], fontweight='bold', pad=15)

    plt.tight_layout()
    return fig
