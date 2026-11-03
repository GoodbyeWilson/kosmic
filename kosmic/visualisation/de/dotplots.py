"""
Dotplot Visualizations
=======================
Patient-level pathway dotplots.
"""

import numpy as np
import pandas as pd


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
                           shared_y=False, theme_colors=None, fdr_by_pathway=None):
    """Create a grid of per-donor pathway score dotplots.

    Parameters
    ----------
    sample_df : pandas.DataFrame
        From compute_patient_pathway_scores.
    pathway_names : list of str
        Pathways to plot.
    disease_label, control_label : str
        Condition labels.
    figsize : tuple, optional
    fdr_by_pathway : dict, optional
        {pathway: FDR} from Pathway DE, shown above each panel. The figure
        computes no test of its own; a pathway without an FDR gets no label.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fdr_by_pathway = fdr_by_pathway or {}
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    from kosmic.visualisation.style import CONTROL_COLOR, DISEASE_COLOR, FG
    tc = theme_colors or {}
    control_color = tc.get('plot_control', CONTROL_COLOR)
    disease_color = tc.get('plot_disease', DISEASE_COLOR)

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
                           fontsize=font_sizes['axis_label'])
        ax.set_title(pw.replace('_', ' '), fontsize=font_sizes['title'], fontweight='bold')
        if i % n_cols == 0:   # one y label per row
            ax.set_ylabel(ylabel, fontsize=font_sizes['axis_label'])
        ax.tick_params(labelsize=font_sizes['tick'])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
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
        pval = fdr_by_pathway.get(pw)
        if pval is not None and np.isfinite(pval):
            # The Pathway DE result for this pathway.
            label = ('FDR < 0.001' if pval < 0.001
                     else f'FDR = {pval:.3f}' if pval < 0.01 else f'FDR = {pval:.2f}')
            ax.plot([0.2, 0.8], [0.9, 0.9], transform=ax.transAxes,
                    color=FG, linewidth=0.8, solid_capstyle='butt',
                    clip_on=False)
            ax.text(0.5, 0.92, label, transform=ax.transAxes, ha='center',
                    va='bottom', fontsize=font_sizes['tick'], color=FG)

    for i in range(n_pathways, n_rows * n_cols):
        axes[i // n_cols, i % n_cols].set_visible(False)

    if suptitle is None:
        suptitle = f'Pathway scores: {disease_label} vs {control_label}'
    plt.suptitle(suptitle,
                 fontsize=font_sizes['suptitle'], fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)

    return fig
