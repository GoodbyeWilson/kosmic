"""
Second-Level (Family) Forest Plot
===================================
Publication-quality forest plot for family-level pooled estimates.
"""

import numpy as np
import pandas as pd


def create_family_forest_plot(
    family_names,
    effect_sizes,
    standard_errors,
    overall_ma=None,
    title="Second-Level Meta-Analysis Forest Plot",
    disease_label="HF",
    control_label="NF",
    figsize=(10, 8),
    save_path=None,
    font_sizes=None,
):
    """Create a forest plot where each row is a family-level pooled estimate.

    Parameters
    ----------
    family_names : list of str
        Labels for each family (y-axis).
    effect_sizes : array-like
        Pooled log2FC for each family.
    standard_errors : array-like
        SE of pooled log2FC for each family.
    overall_ma : dict or None
        Optional override dict with keys ``logfoldchanges``, ``ci_lower``,
        ``ci_upper``, ``heterogeneity_i2`` (0–1), ``tau_squared``,
        ``n_studies``. When ``None``, computed via ``dl_fast`` on the input.
    title : str
        Plot title.
    disease_label, control_label : str
        Labels for direction annotations.
    figsize : tuple
        Figure size in inches.
    save_path : str or None
        If provided, save PNG (300 DPI) and PDF to this path (without extension).

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

    effects = np.array(effect_sizes, dtype=float)
    ses = np.array(standard_errors, dtype=float)
    variances = ses ** 2
    n = len(effects)

    if overall_ma is None:
        from kosmic.meta_analysis.pooling.dl import dl_fast
        # dl_fast pools rows that share a `names` value across a list of
        # DataFrames. Build one "overall" gene with K rows (one per family).
        per_family = [
            pd.DataFrame({
                'names': ['overall'],
                'logfoldchanges': [effects[i]],
                'se': [ses[i]],
                'pvals': [1.0],  # unused by dl_fast beyond matrix build
                'dataset': [f'Family_{i}'],
            })
            for i in range(n)
        ]
        # min_studies=1: the caller chose how many families to pass; pool
        # whatever we have (even a single-family edge case plots sanely).
        overall_ma = dl_fast(per_family, min_studies=1).iloc[0].to_dict()

    # Compute per-family CIs
    ci_lo = effects - 1.96 * ses
    ci_hi = effects + 1.96 * ses

    # Sort by effect size (most negative at top)
    order = np.argsort(effects)
    effects_s = effects[order]
    ci_lo_s = ci_lo[order]
    ci_hi_s = ci_hi[order]
    labels_s = [family_names[i] for i in order]

    # Weight percentages (random-effects weights on the family estimates)
    weights_re = 1.0 / (variances[order] + overall_ma['tau_squared'])
    weight_pcts = 100.0 * weights_re / weights_re.sum()

    fig, ax = plt.subplots(figsize=figsize)

    # --- Individual family rows ---
    for i in range(n):
        color = '#d62728' if effects_s[i] < 0 else '#1f77b4'

        # CI line
        ax.plot([ci_lo_s[i], ci_hi_s[i]], [i, i],
                color=color, linewidth=2, alpha=0.7, zorder=2)

        # Diamond marker for the pooled family estimate
        ax.scatter([effects_s[i]], [i], s=120, c=color, marker='D',
                   edgecolors='black', linewidth=0.8, alpha=0.85, zorder=3)

        # Annotate effect size + CI on the right
        ax.text(ci_hi_s[i] + 0.02, i,
                f" {effects_s[i]:.3f} [{ci_lo_s[i]:.3f}, {ci_hi_s[i]:.3f}]  ({weight_pcts[i]:.1f}%)",
                va='center', ha='left', fontsize=font_sizes['annotation'], color='#333333')

    # --- Overall pooled diamond ---
    pooled_y = n + 1.0
    diamond_hw = 0.35  # half-height of diamond

    # Draw diamond as a filled polygon. The pooled-effect column is
    # ``logfoldchanges`` from dl_fast; an override dict may key it as
    # ``mean``.
    pooled_effect = overall_ma.get('logfoldchanges', overall_ma.get('mean'))
    diamond_x = [overall_ma['ci_lower'], pooled_effect,
                 overall_ma['ci_upper'], pooled_effect]
    diamond_y = [pooled_y, pooled_y - diamond_hw,
                 pooled_y, pooled_y + diamond_hw]
    ax.fill(diamond_x, diamond_y, color='red', alpha=0.7, zorder=4)
    ax.plot(diamond_x + [diamond_x[0]], diamond_y + [diamond_y[0]],
            color='darkred', linewidth=1.2, zorder=5)

    # Shaded CI band
    ax.axvspan(overall_ma['ci_lower'], overall_ma['ci_upper'],
               alpha=0.08, color='red', zorder=0)

    # Reference line at 0
    ax.axvline(0, color='black', linestyle='--', alpha=0.5, linewidth=0.8, zorder=1)

    # Separator line above overall diamond
    ax.axhline(n + 0.3, color='grey', linewidth=0.5, linestyle='-', alpha=0.5)

    # Y-axis labels
    yticks = list(range(n)) + [pooled_y]
    ylabels = labels_s + ['OVERALL']
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=font_sizes['tick'])
    ax.get_yticklabels()[-1].set_fontweight('bold')
    ax.get_yticklabels()[-1].set_color('red')
    ax.invert_yaxis()

    # X-axis
    ax.set_xlabel('Pooled log2FC (95% CI)', fontsize=font_sizes['axis_label'])
    ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold')
    ax.grid(True, alpha=0.2, axis='x')

    # Direction labels
    ax.text(0.01, 0.99, f'Reduced in {disease_label}',
            transform=ax.transAxes, fontsize=font_sizes['tick'], va='top', ha='left',
            color='#d62728', fontstyle='italic')
    ax.text(0.99, 0.99, f'Increased in {disease_label}',
            transform=ax.transAxes, fontsize=font_sizes['tick'], va='top', ha='right',
            color='#1f77b4', fontstyle='italic')

    # Summary annotation box.
    # dl_fast stores I^2 fractional [0, 1] under 'heterogeneity_i2';
    # an override dict may use 'I_squared' as a percentage [0, 100].
    fc = 2 ** pooled_effect
    fc_lo = 2 ** overall_ma['ci_lower']
    fc_hi = 2 ** overall_ma['ci_upper']
    if 'heterogeneity_i2' in overall_ma:
        i_squared_pct = overall_ma['heterogeneity_i2'] * 100
    else:
        i_squared_pct = overall_ma['I_squared']
    ann = (f"Pooled log2FC = {pooled_effect:.4f} "
           f"[{overall_ma['ci_lower']:.4f}, {overall_ma['ci_upper']:.4f}]\n"
           f"FC ({disease_label}/{control_label}) = {fc:.3f} [{fc_lo:.3f}, {fc_hi:.3f}]\n"
           f"I\u00B2 = {i_squared_pct:.1f}%,  "
           f"\u03C4\u00B2 = {overall_ma['tau_squared']:.4f}\n"
           f"k = {overall_ma['n_studies']} families")
    ax.text(0.98, 0.02, ann, transform=ax.transAxes, fontsize=font_sizes['annotation'],
            va='bottom', ha='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow',
                      edgecolor='grey', alpha=0.9))

    fig.tight_layout()

    # Save if requested
    if save_path is not None:
        fig.savefig(f"{save_path}.png", dpi=300, bbox_inches='tight', facecolor='white')
        fig.savefig(f"{save_path}.pdf", dpi=300, bbox_inches='tight', facecolor='white')

    return fig
