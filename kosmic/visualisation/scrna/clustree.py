"""
Clustree Stability Graph
==========================
Matplotlib rendering of Leiden clustering stability across resolutions.
"""

import numpy as np
from collections import Counter


def create_clustree_figure(leiden_per_level, level_positions, figsize=None,
                           font_sizes=None):
    """Render a clustree-style stability graph.

    Parameters
    ----------
    leiden_per_level : list of (resolution, labels_array)
        Leiden results at each resolution level.
    level_positions : list of dict
        {cluster_id: x_position} for each level (from compute_clustree_data).
    figsize : tuple, optional
        Figure size. Auto-calculated if None.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    n_levels = len(leiden_per_level)
    max_clusters = max(len(set(lbl)) for _, lbl in leiden_per_level)

    if figsize is None:
        fig_w = max(8, max_clusters * 0.9)
        fig_h = max(6, n_levels * 2.4)
        figsize = (fig_w, fig_h)

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')

    # Precompute cell counts per (level, cluster)
    all_counts = {}
    for li, (_, lbl) in enumerate(leiden_per_level):
        for k in set(lbl):
            all_counts[(li, k)] = int(np.sum(lbl == k))
    max_count = max(all_counts.values())

    level_cmap = plt.get_cmap('tab10', n_levels)

    # Draw edges (behind nodes)
    for li in range(n_levels - 1):
        _, lbl_curr = leiden_per_level[li]
        _, lbl_next = leiden_per_level[li + 1]
        pos_curr = level_positions[li]
        pos_next = level_positions[li + 1]
        transitions = Counter(zip(lbl_curr, lbl_next))

        for (k_c, k_n), count in transitions.items():
            if k_c not in pos_curr or k_n not in pos_next:
                continue
            proportion = count / all_counts[(li, k_c)]
            lw = 0.4 + proportion * 7
            alpha = 0.25 + proportion * 0.55
            color = '#2166ac' if proportion > 0.5 else '#fc8d59'
            ax.plot(
                [pos_curr[k_c], pos_next[k_n]], [li, li + 1], '-',
                color=color, lw=lw, alpha=alpha, solid_capstyle='round',
            )

    # Draw nodes (in front of edges)
    for li, (_, lbl) in enumerate(leiden_per_level):
        color = level_cmap(li)
        pos = level_positions[li]
        node_sizes = [200 + (all_counts[(li, k)] / max_count) * 900 for k in pos]
        ax.scatter(
            list(pos.values()), [li] * len(pos),
            s=node_sizes, c=[color] * len(pos),
            zorder=5, edgecolors='white', linewidths=1.5,
        )
        fsize = max(5, font_sizes['tick'] - n_levels // 4)
        for k, x in pos.items():
            ax.text(x, li, str(k), ha='center', va='center',
                    fontsize=fsize, fontweight='bold', color='white', zorder=6)

    # Axes
    ax.set_yticks(range(n_levels))
    ax.set_yticklabels(
        [f"res={r:.2f}  ({len(set(lbl))} clusters)"
         for r, lbl in leiden_per_level],
        fontsize=font_sizes['tick'],
    )
    ax.invert_yaxis()
    ax.set_xlim(-0.12, 1.12)
    ax.set_xticks([])
    for spine in ['top', 'right', 'bottom']:
        ax.spines[spine].set_visible(False)

    legend_elements = [
        mpatches.Patch(color='#2166ac', alpha=0.85,
                       label='Dominant flow (>50% of source cluster)'),
        mpatches.Patch(color='#fc8d59', alpha=0.85,
                       label='Secondary flow (≤50% — unstable split)'),
    ]
    ax.legend(handles=legend_elements, loc='lower right',
              fontsize=font_sizes['annotation'], framealpha=0.85)
    ax.set_title(
        "Leiden Clustering Stability Analysis\n"
        "Edge width \u221d cell proportion  \u00b7  "
        "blue = stable dominant flow  \u00b7  orange = unstable / splitting",
        fontsize=font_sizes['tick'], pad=14,
    )

    plt.tight_layout()
    return fig
