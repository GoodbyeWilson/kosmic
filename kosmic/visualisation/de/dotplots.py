"""
Dotplot Visualizations
=======================
Patient-level pathway dotplots, per-gene dotplots, and expression dotplots.
"""

import numpy as np
import pandas as pd

from kosmic import DEFAULT_FDR


def pathway_scores_to_sample_df(pathway_score_data, standardize=False):
    """Reshape the DE window's cached per-donor pathway scores for the dotplot.

    'pathway_score_data' is the dict the pathway-scoring run caches
    ('names', 'matrix' (samples x pathways), 'sample_df' with a 'condition'
    column) using whichever scoring method the user selected. When 'standardize'
    is set, each pathway's per-donor scores are z-scored across samples (mean 0,
    SD 1) so pathways sit on a comparable scale. Returns (sample_df,
    pathway_names) in the shape 'create_patient_dotplot' expects.
    """
    names = list(pathway_score_data.get('names') or [])
    matrix = np.asarray(pathway_score_data.get('matrix'), dtype=float)
    base = pathway_score_data.get('sample_df')
    out = pd.DataFrame({'condition': np.asarray(base['condition'])})
    if 'sample' in getattr(base, 'columns', []):
        out['sample'] = np.asarray(base['sample'])
    for col, pw in enumerate(names):
        vals = matrix[:, col].astype(float)
        if standardize:
            sd = np.nanstd(vals)
            vals = (vals - np.nanmean(vals)) / sd if sd > 0 else vals - np.nanmean(vals)
        out[f'{pw}_raw'] = vals
    return out, names


def create_patient_dotplot(sample_df, pathway_names, disease_label, control_label,
                           figsize=None, font_sizes=None,
                           ylabel='Pathway score (per donor)', suptitle=None,
                           shared_y=False, theme_colors=None):
    """Create a grid of pathway dotplots with MWU statistics.

    Parameters
    ----------
    sample_df : pandas.DataFrame
        From compute_patient_pathway_scores.
    pathway_names : list of str
        Pathways to plot.
    disease_label, control_label : str
        Condition labels.
    figsize : tuple, optional

    Returns
    -------
    tuple of (fig, pathway_stats)
        fig: matplotlib Figure
        pathway_stats: dict of {pathway: {p_value, p_corrected, disease_mean, ...}}
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.stats import mannwhitneyu
    from kosmic.numerical import bh_fdr

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    tc = theme_colors or {}
    control_color = tc.get('plot_control', '#3498DB')
    disease_color = tc.get('plot_disease', '#E74C3C')

    # Stats
    pathway_stats = {}
    p_values = []
    pw_for_correction = []

    for pw in pathway_names:
        raw_col = f'{pw}_raw'
        d_scores = sample_df[sample_df['condition'] == disease_label][raw_col].values
        c_scores = sample_df[sample_df['condition'] == control_label][raw_col].values
        if len(d_scores) >= 2 and len(c_scores) >= 2:
            try:
                stat, pval = mannwhitneyu(d_scores, c_scores, alternative='two-sided')
                pathway_stats[pw] = {
                    'p_value': pval,
                    'disease_mean': np.mean(d_scores),
                    'control_mean': np.mean(c_scores),
                    'disease_sem': np.std(d_scores) / np.sqrt(len(d_scores)),
                    'control_sem': np.std(c_scores) / np.sqrt(len(c_scores)),
                    'disease_var': np.var(d_scores, ddof=1),
                    'control_var': np.var(c_scores, ddof=1),
                    'n_disease': len(d_scores),
                    'n_control': len(c_scores),
                }
                p_values.append(pval)
                pw_for_correction.append(pw)
            except Exception:
                pathway_stats[pw] = {'p_value': np.nan}
        else:
            pathway_stats[pw] = {'p_value': np.nan}

    if p_values:
        p_corrected = bh_fdr(p_values)
        for i, pw in enumerate(pw_for_correction):
            pathway_stats[pw]['p_corrected'] = p_corrected[i]

    # Plot
    n_pathways = len(pathway_names)
    n_cols = 3
    n_rows = int(np.ceil(n_pathways / n_cols))
    if figsize is None:
        figsize = (15, 4 * n_rows)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    # A centred score (z-scores, or per-donor log2FC) has negatives -> draw a
    # zero reference and don't floor the axis at 0.
    all_vals = sample_df[[f'{pw}_raw' for pw in pathway_names]].to_numpy(dtype=float)
    centred = bool(np.any(all_vals < 0)) if all_vals.size else False

    def _box(ax, x, vals, color):
        """Tukey box + median + 1.5-IQR whiskers, matching the DE-window plot."""
        vals = np.asarray(vals, dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            return
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        iqr = q3 - q1
        if iqr > 0:
            lo = float(vals[vals >= q1 - 1.5 * iqr].min())
            hi = float(vals[vals <= q3 + 1.5 * iqr].max())
        else:
            lo, hi = float(vals.min()), float(vals.max())
        w, cap = 0.5, 0.125
        xl, xr = x - w / 2, x + w / 2
        line_kw = dict(color=color, solid_capstyle='butt', clip_on=False, zorder=2)
        ax.plot([xl, xr, xr, xl, xl], [q1, q1, q3, q3, q1], linewidth=2, **line_kw)
        ax.plot([xl, xr], [med, med], linewidth=3, **line_kw)
        ax.plot([x, x], [q3, hi], linewidth=2, **line_kw)
        ax.plot([x, x], [q1, lo], linewidth=2, **line_kw)
        ax.plot([x - cap, x + cap], [hi, hi], linewidth=2, **line_kw)
        ax.plot([x - cap, x + cap], [lo, lo], linewidth=2, **line_kw)

    rng = np.random.default_rng(42)
    ranges = []
    for i, pw in enumerate(pathway_names):
        ax = axes[i // n_cols, i % n_cols]
        raw_col = f'{pw}_raw'
        d_scores = sample_df[sample_df['condition'] == disease_label][raw_col].to_numpy(dtype=float)
        c_scores = sample_df[sample_df['condition'] == control_label][raw_col].to_numpy(dtype=float)

        if centred:
            ax.axhline(0, color='grey', linewidth=0.8, zorder=0, clip_on=False)

        _box(ax, 0, c_scores, control_color)
        _box(ax, 1, d_scores, disease_color)
        # Single filled circles: no edge stroke, no clip mask -> clean vector.
        ax.scatter(rng.uniform(-0.12, 0.12, len(c_scores)), c_scores,
                   color=control_color, s=55, linewidths=0, edgecolors='none',
                   zorder=3, clip_on=False)
        ax.scatter(1 + rng.uniform(-0.12, 0.12, len(d_scores)), d_scores,
                   color=disease_color, s=55, linewidths=0, edgecolors='none',
                   zorder=3, clip_on=False)

        alls = np.concatenate([c_scores, d_scores])
        alls = alls[np.isfinite(alls)]
        ranges.append((float(alls.min()), float(alls.max())) if len(alls) else (0.0, 1.0))

        ax.set_xlim(-0.5, 1.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([control_label, disease_label], ha='center',
                           fontsize=font_sizes['title'], fontweight='bold')
        ax.set_title(pw.replace('_', ' '), fontsize=font_sizes['title'], fontweight='bold')
        ax.set_ylabel(ylabel, fontsize=font_sizes['axis_label'], fontweight='bold')
        ax.grid(False)

    # Y-limits: one shared scale across the grid, or per-subplot.
    if shared_y and ranges:
        gmin = min(r[0] for r in ranges)
        gmax = max(r[1] for r in ranges)
        span = (gmax - gmin) or 1.0
        shared_lim = ((0 if (not centred and gmin >= 0) else gmin - 0.1 * span),
                      gmax + 0.25 * span)  # headroom for the significance label

    for i, pw in enumerate(pathway_names):
        ax = axes[i // n_cols, i % n_cols]
        if shared_y and ranges:
            ax.set_ylim(*shared_lim)
        else:
            lo, hi = ranges[i]
            span = (hi - lo) or 1.0
            ax.set_ylim((0 if (not centred and lo >= 0) else lo - 0.1 * span),
                        hi + 0.25 * span)
        st = pathway_stats.get(pw, {})
        if 'p_corrected' in st:
            pval = st['p_corrected']
            label = ('*' if pval < DEFAULT_FDR
                     else 'p < 0.001' if pval < 0.001
                     else f'p = {pval:.3f}' if pval < 0.01 else f'p = {pval:.2f}')
            ax.plot([0.2, 0.8], [0.9, 0.9], transform=ax.transAxes,
                    color='black', linewidth=1, solid_capstyle='butt',
                    clip_on=False)
            ax.text(0.5, 0.92, label, transform=ax.transAxes, ha='center',
                    va='bottom', fontsize=font_sizes['tick'], fontweight='bold')

    for i in range(n_pathways, n_rows * n_cols):
        axes[i // n_cols, i % n_cols].set_visible(False)

    if suptitle is None:
        suptitle = f'Pathway scores: {disease_label} vs {control_label}'
    plt.suptitle(suptitle,
                 fontsize=font_sizes['suptitle'], fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)

    return fig, pathway_stats


def create_gene_dotplot(sample_df, pathway_genes, pathway_name,
                        disease_label, control_label, figsize=None,
                        font_sizes=None, theme_colors=None):
    """Create a per-gene dotplot grid for a single pathway.

    Parameters
    ----------
    sample_df : pandas.DataFrame
        Must have 'condition' column + gene columns with mean expression values.
    pathway_genes : list of str
        Genes to plot (must be columns in sample_df).
    pathway_name : str
        Pathway name for title.
    disease_label, control_label : str
        Condition labels.
    figsize : tuple, optional

    Returns
    -------
    matplotlib.figure.Figure or None
        None if fewer than 3 genes available.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    control_color = (theme_colors.get('plot_control', '#3498DB')
                     if theme_colors else '#3498DB')
    disease_color = (theme_colors.get('plot_disease', '#E74C3C')
                     if theme_colors else '#E74C3C')

    genes = [g for g in pathway_genes if g in sample_df.columns]
    if len(genes) < 3:
        return None

    n_genes = len(genes)
    n_cols = min(6, n_genes)
    n_rows = int(np.ceil(n_genes / n_cols))
    if figsize is None:
        figsize = (4 * n_cols, 4 * n_rows)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    # Compute y-axis range
    all_vals = []
    for g in genes:
        all_vals.extend(sample_df[g].values)
    all_vals = np.array(all_vals)
    if np.any(all_vals < 0):
        y_min = np.min(all_vals) * 1.15
        y_max = np.max(all_vals) * 1.15
    else:
        y_min = 0
        y_max = np.max(all_vals) * 1.1 + 0.1

    for i, gene in enumerate(genes):
        r, c = i // n_cols, i % n_cols
        ax = axes[r, c]

        d_vals = sample_df[sample_df['condition'] == disease_label][gene].values
        c_vals = sample_df[sample_df['condition'] == control_label][gene].values

        x_c = np.random.normal(0, 0.04, len(c_vals))
        x_d = np.random.normal(1, 0.04, len(d_vals))

        ax.scatter(x_c, c_vals, color=control_color, alpha=0.7, s=50, edgecolors='darkblue', linewidth=0.5)
        ax.scatter(x_d, d_vals, color=disease_color, alpha=0.7, s=50, edgecolors='darkred', linewidth=0.5)

        c_mean = np.mean(c_vals) if len(c_vals) > 0 else 0
        d_mean = np.mean(d_vals) if len(d_vals) > 0 else 0

        ax.bar(0, c_mean, width=0.3, color=control_color, alpha=0.3, edgecolor='darkblue')
        ax.bar(1, d_mean, width=0.3, color=disease_color, alpha=0.3, edgecolor='darkred')
        ax.hlines(c_mean, -0.15, 0.15, colors='darkblue', linewidth=2)
        ax.hlines(d_mean, 0.85, 1.15, colors='darkred', linewidth=2)

        ax.set_xlim(-0.4, 1.4)
        ax.set_ylim(y_min, y_max)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([control_label, disease_label], fontsize=font_sizes['tick'])
        ax.set_title(gene, fontsize=font_sizes['tick'], fontweight='bold')
        ax.set_ylabel('Mean Expression', fontsize=font_sizes['annotation'])
        ax.grid(True, alpha=0.3, axis='y')

        if c_mean > 0 and d_mean > 0:
            fc_text = f'FC: {d_mean / c_mean:.2f}'
        elif np.any(d_vals < 0) or np.any(c_vals < 0):
            fc_text = f'Diff: {d_mean - c_mean:.2f}'
        else:
            fc_text = 'FC: 1.00'
        ax.text(0.5, y_min + 0.05, fc_text, ha='center', fontsize=font_sizes['annotation'],
                style='italic',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

    for i in range(n_genes, n_rows * n_cols):
        axes[i // n_cols, i % n_cols].set_visible(False)

    plt.suptitle(f'{pathway_name.replace("_", " ")} Genes\n{disease_label} vs {control_label}',
                 fontsize=font_sizes['title'], fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)

    return fig


def create_expression_dotplot(stats_df, resolved_genes, conditions,
                              expressing_only=False, figsize=None,
                              dark_mode=True, theme_colors=None,
                              font_sizes=None, dotplot_cmap='Reds'):
    """Create a scanpy-style multi-gene dot plot.

    Parameters
    ----------
    stats_df : pandas.DataFrame
        From compute_dotplot_stats. Must have: gene, condition, pct_expressing,
        mean_expression, mean_all_cells, n_cells.
    resolved_genes : list of str
        Ordered gene names.
    conditions : list of str
        Ordered condition names.
    expressing_only : bool
        Label for expressing-only mode.
    figsize : tuple, optional
    dark_mode : bool
        If True, use dark background with light text.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import Normalize
    import matplotlib.cm as cm

    # Theme colours — use passed-in colours if available, else sensible defaults
    if theme_colors:
        bg_color = theme_colors.get('bg_primary', '#1e1e1e')
        fg_color = theme_colors.get('fg_primary', '#d4d4d4')
        fg2_color = theme_colors.get('fg_secondary', '#858585')
        edge_color = theme_colors.get('border', '#3c3c3c')
    elif dark_mode:
        bg_color = '#1e1e1e'
        fg_color = '#d4d4d4'
        fg2_color = '#858585'
        edge_color = '#3c3c3c'
    else:
        bg_color = 'white'
        fg_color = 'black'
        fg2_color = '#555555'
        edge_color = 'black'

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    n_genes = len(resolved_genes)
    n_conditions = len(conditions)

    x_positions = {gene: i for i, gene in enumerate(resolved_genes)}
    y_positions = {cond: i for i, cond in enumerate(conditions)}

    all_pcts = stats_df['pct_expressing'].values
    all_means = stats_df['mean_expression'].values

    mean_min, mean_max = all_means.min(), all_means.max()
    norm = Normalize(vmin=mean_min, vmax=mean_max)
    max_dot_size = 350

    if figsize is None:
        fig_width = max(4, 1.0 * n_genes + 3.5)
        fig_height = max(2.5, 0.7 * n_conditions + 2.0)
        figsize = (fig_width, fig_height)

    fig = plt.figure(figsize=figsize, facecolor=bg_color)
    gs = gridspec.GridSpec(1, 2, width_ratios=[n_genes, 1.8],
                           wspace=0.4, left=0.15, right=0.95,
                           top=0.92, bottom=0.25)
    ax = fig.add_subplot(gs[0, 0])
    ax.set_facecolor(bg_color)
    legend_ax = fig.add_subplot(gs[0, 1])
    legend_ax.axis('off')
    legend_ax.set_facecolor(bg_color)

    cmap = cm.get_cmap(dotplot_cmap)
    for _, row in stats_df.iterrows():
        x = x_positions.get(row['gene'])
        y = y_positions.get(row['condition'])
        if x is None or y is None:
            continue
        size = (row['pct_expressing'] / 100.0) * max_dot_size
        color = cmap(norm(row['mean_expression']))
        ax.scatter(x, y, s=size, c=[color], edgecolors=edge_color, linewidths=0.5, zorder=3)

    ax.set_xticks(range(n_genes))
    ax.set_xticklabels(resolved_genes, rotation=45, ha='right',
                       fontsize=font_sizes['tick'], fontstyle='italic', color=fg_color)
    ax.set_yticks(range(n_conditions))
    ax.set_yticklabels(conditions, fontsize=font_sizes['tick'], color=fg_color)
    ax.set_xlim(-0.5, n_genes - 0.5)
    ax.set_ylim(-0.5, n_conditions - 0.5)
    ax.invert_yaxis()
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['left'].set_color(edge_color)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['bottom'].set_color(edge_color)
    ax.tick_params(axis='both', which='both', length=3, width=0.5, colors=fg2_color)

    # Size legend
    legend_pcts = [p for p in [10, 30, 50, 70, 90] if p <= max(all_pcts.max(), 10) + 10]
    if not legend_pcts:
        legend_pcts = [10]

    legend_y_top = 0.95
    legend_x = 0.35

    legend_ax.text(legend_x, legend_y_top + 0.05, 'Fraction of cells\nin group (%)',
                   fontsize=font_sizes['annotation'], fontweight='bold', ha='center', va='top',
                   transform=legend_ax.transAxes, color=fg_color)

    for i, pct in enumerate(legend_pcts):
        y_pos = legend_y_top - 0.12 - i * 0.12
        size = (pct / 100.0) * max_dot_size
        legend_ax.scatter(legend_x - 0.15, y_pos, s=size,
                          c='lightgray', edgecolors=edge_color, linewidths=0.5,
                          transform=legend_ax.transAxes, clip_on=False, zorder=5)
        legend_ax.text(legend_x + 0.15, y_pos, f'{pct}', fontsize=font_sizes['annotation'],
                       va='center', ha='left', transform=legend_ax.transAxes, color=fg2_color)

    # Colorbar
    legend_bbox = legend_ax.get_position()
    cbar_frac = max(0.05, legend_y_top - 0.12 - len(legend_pcts) * 0.12 - 0.10)
    cbar_bottom = legend_bbox.y0 + cbar_frac * legend_bbox.height
    cbar_ax = fig.add_axes([
        legend_bbox.x0 + 0.02, cbar_bottom,
        legend_bbox.width * 0.75, 0.015,
    ])
    cbar_ax.set_facecolor(bg_color)
    sm = cm.ScalarMappable(cmap=dotplot_cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar.outline.set_linewidth(0.5)
    cbar.outline.set_edgecolor(edge_color)
    cbar.ax.tick_params(labelsize=font_sizes['colorbar_tick'], width=0.5, length=2, colors=fg2_color)
    mean_label = "Mean expression\n(expressing only)" if expressing_only else "Mean expression\nin group"
    cbar.ax.set_title(mean_label, fontsize=font_sizes['annotation'], fontweight='bold', pad=4, color=fg_color)

    return fig
