"""
Volcano Plot Rendering
=======================
Publication-quality volcano plots for DE results.
"""

import numpy as np
import pandas as pd

from kosmic.numerical import neg_log10
from kosmic import DEFAULT_FDR
from kosmic.visualisation.labels import place_labels
from kosmic.visualisation.style import DOWN_COLOR, FG, NS_COLOR, UP_COLOR


def create_volcano_plot(de_results, disease_cond, control_cond, dataset_name='',
                        pval_col='pvals_adj', pval_threshold=DEFAULT_FDR,
                        logfc_threshold=0.25, figsize=(8, 8),
                        genes_of_interest=None, max_labels=10,
                        font_sizes=None, point_size=9,
                        color_scheme=0, x_limit=None, y_limit=None, title=None,
                        xlabel=None, subtitle=None):
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
        Significant genes to label before the top hits.
    max_labels : int
        Max gene labels on plot. Without *genes_of_interest*, the most
        significant genes are labelled, split between up and down.
    x_limit : float, optional
        Symmetric x-axis limit (+/- log2FC). None or 0 uses the dynamic range.
    y_limit : float, optional
        Top of the y-axis (-log10 p). None or 0 leaves it automatic.
    title : str, optional
        Headline for the title's first line. Defaults to a generic
        "Differential Expression" -- callers running a pathway-restricted
        or otherwise non-genome-wide analysis should pass their own so the
        title doesn't claim "Metabolic" for genes it didn't test.
    xlabel, subtitle : str, optional
        Replace the default x-axis label and the title's second line (e.g.
        for pooled meta-analysis results).

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
        (UP_COLOR, DOWN_COLOR),   # red / blue
        ('#e69620', '#1eaa9e'),   # orange / teal
        ('#c832c8', '#32c8dc'),   # magenta / cyan
    ]
    up_color, down_color = _COLOR_SCHEMES[min(color_scheme, len(_COLOR_SCHEMES) - 1)]

    # Significant up / down / not significant. Vectorised so it never
    # depends on per-row index lookups.
    sig = (pvals < pval_threshold) & (lfc.abs() > logfc_threshold)
    up, down = sig & (lfc > 0), sig & (lfc < 0)
    # Not-significant genes first, so coloured points sit on top.
    for mask, color in ((~sig, NS_COLOR), (down, down_color), (up, up_color)):
        ax.scatter(lfc[mask], neg_log10_pvals[mask], color=color,
                   s=point_size, edgecolors='none', linewidth=0)

    # x-axis: user-set symmetric limit, else the data range.
    if x_limit and x_limit > 0:
        fc_lim = float(x_limit)
    else:
        fc_lim = max(float(lfc.abs().max()) * 1.05, logfc_threshold * 1.5, 0.5)
    ax.set_xlim(-fc_lim, fc_lim)
    if y_limit and y_limit > 0:
        ax.set_ylim(top=float(y_limit))
    else:
        # Headroom for labels at the top.
        y_max = float(neg_log10_pvals.max()) if len(neg_log10_pvals) else 1.0
        if np.isfinite(y_max) and y_max > 0:
            ax.set_ylim(bottom=0, top=y_max * 1.1)

    ax.axhline(y=-np.log10(pval_threshold), color=FG, linestyle=':', linewidth=0.8)

    genes_to_label = _select_genes_to_label(
        de_results, pval_col, pval_threshold, logfc_threshold,
        genes_of_interest, max_labels,
    )
    labels = [(de_results.loc[i, 'names'], de_results.loc[i, 'logfoldchanges'],
               neg_log10_pvals.loc[i]) for i in genes_to_label]

    stat = 'P' if is_nominal else 'FDR'
    ax.set_xlabel(xlabel or f'Log$_2$ fold change ({disease_cond} vs {control_cond})',
                  fontsize=font_sizes['axis_label'])
    ax.set_ylabel(f'$-$Log$_{{10}}$ {stat}', fontsize=font_sizes['axis_label'])

    headline = title or 'Differential Expression'
    if is_nominal:
        headline += ' (nominal P)'
    if subtitle is None:
        subtitle = f'{disease_cond} vs {control_cond}'
        if dataset_name:
            subtitle += f' ({dataset_name})'
    ax.set_title(f'{headline}\n{subtitle}', fontsize=font_sizes['title'],
                 fontweight='bold')
    ax.tick_params(direction='out', labelsize=font_sizes['tick'])

    fig.tight_layout()
    # After the layout: placement works in page coordinates, which the
    # layout changes.
    place_labels(ax, labels, font_sizes['annotation'])
    return fig


def _select_genes_to_label(de_results, pval_col, pval_threshold, logfc_threshold,
                           genes_of_interest, max_labels):
    """Select gene indices to label on volcano plot."""
    genes_of_interest = genes_of_interest or []

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
                top = subset.nsmallest(max(1, max_labels // 2), pval_col)
                indices.extend(top.index.tolist())

    # Deduplicate preserving order
    seen = set()
    unique = []
    for idx in indices:
        if idx not in seen:
            unique.append(idx)
            seen.add(idx)
    return unique[:max_labels]
