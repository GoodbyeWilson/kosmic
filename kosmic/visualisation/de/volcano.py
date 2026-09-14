"""
Volcano Plot Rendering
=======================
Publication-quality volcano plots for DE results.
"""

import numpy as np
import pandas as pd

from kosmic.numerical import neg_log10
from kosmic import DEFAULT_FDR


def create_volcano_plot(de_results, disease_cond, control_cond, dataset_name='',
                        pval_col='pvals_adj', pval_threshold=DEFAULT_FDR,
                        logfc_threshold=0.25, figsize=(8, 8),
                        genes_of_interest=None, max_labels=10,
                        font_sizes=None, point_size=80,
                        color_scheme=0, x_limit=None, y_limit=None, title=None):
    """Create a volcano plot from DE results.

    Parameters
    ----------
    de_results : pandas.DataFrame
        Must have columns: names, logfoldchanges, and the pval_col.
    disease_cond, control_cond : str
        Condition labels for axis labels.
    dataset_name : str
        Dataset name for title.
    pval_col : str
        'pvals_adj' for adjusted, 'pvals' for nominal.
    pval_threshold : float
        Significance threshold.
    logfc_threshold : float
        Fold change threshold.
    figsize : tuple
        Figure size.
    genes_of_interest : list, optional
        Gene names to prioritize labeling.
    max_labels : int
        Max gene labels on plot.
    x_limit : float, optional
        Symmetric x-axis limit (+/- log2FC). None or 0 uses the dynamic range.
    y_limit : float, optional
        Top of the y-axis (-log10 p). None or 0 leaves it automatic (with a
        small headroom margin so top labels clear the legend box).
    title : str, optional
        Headline for the title's first line. Defaults to a generic
        "Differential Expression" -- callers running a pathway-restricted
        or otherwise non-genome-wide analysis should pass their own so the
        title doesn't claim "Metabolic" for genes it didn't test.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    # Guard against a non-unique index or duplicate column labels: either turns
    # a scalar lookup into a Series/DataFrame and breaks the boolean tests below.
    de_results = de_results.loc[:, ~de_results.columns.duplicated()].reset_index(drop=True)
    # Drop DataFrame-valued attrs (e.g. 'genome_wide_ranking'): pandas propagates
    # attrs through nsmallest/concat and compares them with '==', which raises
    # "truth value of a DataFrame is ambiguous" when an attr value is a frame.
    de_results.attrs = {}

    # Fall back if the requested p-value column is absent.
    if pval_col not in de_results.columns:
        pval_col = 'pvals_adj' if 'pvals_adj' in de_results.columns else 'pvals'

    fig, ax = plt.subplots(figsize=figsize)

    is_nominal = pval_col == 'pvals'
    pvals = pd.to_numeric(de_results[pval_col], errors='coerce')
    lfc = pd.to_numeric(de_results['logfoldchanges'], errors='coerce')
    neg_log10_pvals = pd.Series(neg_log10(pvals), index=de_results.index)

    # Color schemes: (up_color, down_color)
    _COLOR_SCHEMES = [
        ('#d62728', '#1f77b4'),   # red / blue
        ('#e69620', '#1eaa9e'),   # orange / teal
        ('#c832c8', '#32c8dc'),   # magenta / cyan
    ]
    up_color, down_color = _COLOR_SCHEMES[min(color_scheme, len(_COLOR_SCHEMES) - 1)]

    # Colour points: significant + up / down / not significant. Vectorised so it
    # never depends on per-row index lookups.
    sig = (pvals < pval_threshold) & (lfc.abs() > logfc_threshold)
    colors = list(np.where(sig & (lfc > 0), up_color,
                           np.where(sig & (lfc < 0), down_color, 'gray')))

    ax.scatter(
        lfc, neg_log10_pvals,
        c=colors, alpha=1.0, s=point_size, edgecolors='none', linewidth=0,
    )

    # x-axis: user-set symmetric limit, else a dynamic range from the data.
    if x_limit and x_limit > 0:
        fc_lim = float(x_limit)
    else:
        fc_max = lfc.abs().max()
        fc_lim = max(fc_max * 1.2, logfc_threshold * 1.5, 0.5)
    ax.set_xlim(-fc_lim, fc_lim)
    if y_limit and y_limit > 0:
        ax.set_ylim(top=float(y_limit))
    else:
        # Headroom so a point/label sitting at the data max doesn't collide
        # with the legend box pinned to the top-right corner.
        y_max = float(neg_log10_pvals.max()) if len(neg_log10_pvals) else 1.0
        if np.isfinite(y_max) and y_max > 0:
            ax.set_ylim(top=y_max * 1.15)

    # Threshold lines
    ax.axhline(y=-np.log10(pval_threshold), color='black', linestyle='--', alpha=0.8, linewidth=1)
    ax.axvline(x=logfc_threshold, color='black', linestyle='--', alpha=0.8, linewidth=1)
    ax.axvline(x=-logfc_threshold, color='black', linestyle='--', alpha=0.8, linewidth=1)

    # Label genes
    genes_to_label = _select_genes_to_label(
        de_results, pval_col, pval_threshold, logfc_threshold,
        genes_of_interest, max_labels,
    )
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    for idx in genes_to_label:
        row = de_results.loc[idx]
        x, y = row['logfoldchanges'], neg_log10_pvals.loc[idx]
        # Always offsetting up-and-right ran labels for points near the
        # right/top edge (e.g. the single strongest hit) past the axes
        # border and under the legend pinned at 'upper right'. Flip the
        # offset toward whichever quadrant has room instead.
        near_right = x > x0 + 0.75 * (x1 - x0)
        near_top = y > y0 + 0.8 * (y1 - y0)
        dx, ha = (-5, 'right') if near_right else (5, 'left')
        dy, va = (-5, 'top') if near_top else (5, 'bottom')
        ax.annotate(
            row['names'], (x, y),
            xytext=(dx, dy), textcoords='offset points',
            fontsize=font_sizes['tick'], fontweight='bold',
            ha=ha, va=va,
        )

    pval_label = 'Nominal P-value' if is_nominal else 'Adjusted P-value'
    ax.set_xlabel(f'Log2 Fold Change ({disease_cond} vs {control_cond})',
                  fontsize=font_sizes['axis_label'], fontweight='bold')
    ax.set_ylabel(f'-Log10 {pval_label}',
                  fontsize=font_sizes['axis_label'], fontweight='bold')

    title_suffix = ' (Nominal P-values)' if is_nominal else ''
    headline = title or 'Differential Expression'
    ax.set_title(
        f'{headline}{title_suffix}\n'
        f'{disease_cond} vs {control_cond} ({dataset_name})',
        fontsize=font_sizes['title'], fontweight='bold',
    )

    n_up = colors.count(up_color)
    n_down = colors.count(down_color)
    legend_elements = [
        Patch(facecolor=up_color, alpha=1.0, label=f'Upregulated in {disease_cond} (n={n_up})'),
        Patch(facecolor=down_color, alpha=1.0, label=f'Downregulated in {disease_cond} (n={n_down})'),
        Patch(facecolor='gray', alpha=1.0, label='Not significant'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=font_sizes['legend'])

    if is_nominal:
        info = f'Nominal p < {pval_threshold}\n|log2FC| > {logfc_threshold}\nTotal significant: {n_up + n_down}'
    else:
        info = (f'Patient-level analysis\n(Pseudobulk aggregation)\n\n'
                f'padj < {pval_threshold}\n|log2FC| > {logfc_threshold}\nTotal significant: {n_up + n_down}')
    ax.text(0.02, 0.98, info, transform=ax.transAxes, fontsize=font_sizes['tick'],
            verticalalignment='top', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.9, edgecolor='gray'))

    plt.tight_layout()
    return fig


def _select_genes_to_label(de_results, pval_col, pval_threshold, logfc_threshold,
                           genes_of_interest, max_labels):
    """Select gene indices to label on volcano plot."""
    if genes_of_interest is None:
        genes_of_interest = [
            'ATP5F1A', 'COX4I1', 'GAPDH', 'HK1', 'HK2', 'LDHA', 'LDHB', 'PKM',
            'CS', 'IDH1', 'IDH2', 'SDHA', 'SDHB',
            'Atp5f1a', 'Cox4i1', 'Gapdh', 'Hk1', 'Hk2', 'Ldha', 'Ldhb', 'Pkm',
            'Cs', 'Idh1', 'Idh2', 'Sdha', 'Sdhb',
        ]

    indices = []

    # Genes of interest that are significant
    for gene in genes_of_interest:
        matches = de_results[de_results['names'] == gene]
        if len(matches) > 0:
            row = matches.iloc[0]
            if row[pval_col] < pval_threshold and abs(row['logfoldchanges']) > logfc_threshold:
                indices.append(matches.index[0])

    # Top significant by p-value
    sig = de_results[
        (de_results[pval_col] < pval_threshold) &
        (de_results['logfoldchanges'].abs() > logfc_threshold)
    ]
    if len(sig) > 0:
        for subset in [sig[sig['logfoldchanges'] > 0], sig[sig['logfoldchanges'] < 0]]:
            if len(subset) > 0:
                top = subset.nsmallest(5, pval_col)
                indices.extend(top.index.tolist())

    # Deduplicate preserving order
    seen = set()
    unique = []
    for idx in indices:
        if idx not in seen:
            unique.append(idx)
            seen.add(idx)
    return unique[:max_labels]
