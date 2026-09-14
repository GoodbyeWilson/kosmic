"""
PCA visualisation helpers (no GUI imports).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def find_elbow(variance_ratio: np.ndarray, min_pcs: int = 10,
               cumulative_threshold: float = 0.80) -> int:
    """Detect the suggested number of PCs for scRNA-seq analysis.

    Hybrid approach: takes the maximum of three estimates to avoid
    under-estimation, which is the most common failure mode in scRNA.

    1. **Max-distance-to-line** (classic elbow/Kneedle method)
    2. **Cumulative variance threshold** (default 80%)
    3. **Hard floor** (default 10 PCs)

    Parameters
    ----------
    variance_ratio : array-like
        Per-PC explained variance ratio.
    min_pcs : int
        Absolute minimum PCs to suggest (default 10).
    cumulative_threshold : float
        Cumulative variance fraction target (default 0.80).

    Returns
    -------
    int
        1-based PC number at the suggested cutoff.
    """
    y = np.asarray(variance_ratio)
    n = len(y)
    if n < 3:
        return n

    # Method 1: max-distance-to-line (classic elbow)
    x = np.linspace(0, 1, n)
    p1 = np.array([x[0], y[0]])
    p2 = np.array([x[-1], y[-1]])
    line_vec = p2 - p1
    line_len = np.linalg.norm(line_vec)
    if line_len == 0:
        return min_pcs

    line_unit = line_vec / line_len
    dists = np.abs(np.cross(line_unit, p1 - np.column_stack([x, y])))
    elbow_pc = int(np.argmax(dists)) + 1  # 1-based

    # Method 2: cumulative variance threshold
    cumulative = np.cumsum(y)
    above = np.where(cumulative >= cumulative_threshold)[0]
    cum_pc = int(above[0]) + 1 if len(above) > 0 else n

    # Take the max of all three
    suggested = max(elbow_pc, cum_pc, min_pcs)
    # Cap at available PCs
    return min(suggested, n)


def create_elbow_plot(variance_ratio: np.ndarray, *,
                      suggested_pcs: int | None = None,
                      dark_mode: bool = True):
    """Elbow plot showing per-PC variance explained and the cumulative curve.

    Parameters
    ----------
    variance_ratio : array-like
        Per-PC explained variance ratio (from ``adata.uns['pca']['variance_ratio']``).
    suggested_pcs : int or None
        If given, mark this PC on the plot as the suggested elbow point.
    dark_mode : bool
        If True, use a dark background suitable for the app theme.

    Returns
    -------
    matplotlib.figure.Figure
    """
    variance_ratio = np.asarray(variance_ratio)
    n_pcs = len(variance_ratio)
    x = np.arange(1, n_pcs + 1)
    cumulative = np.cumsum(variance_ratio)

    if dark_mode:
        plt.style.use('dark_background')
    else:
        plt.style.use('default')

    fig, ax1 = plt.subplots(figsize=(8, 4.5), dpi=120)

    # Per-PC variance
    ax1.bar(x, variance_ratio * 100, color='#5dade2', alpha=0.7, width=0.8,
            label='Per PC')
    ax1.set_xlabel('Principal Component')
    ax1.set_ylabel('Variance Explained (%)')
    ax1.set_xlim(0.5, n_pcs + 0.5)

    # Cumulative on secondary axis
    ax2 = ax1.twinx()
    ax2.plot(x, cumulative * 100, color='#e74c3c', linewidth=2,
             marker='o', markersize=3, label='Cumulative')
    ax2.set_ylabel('Cumulative Variance (%)')
    ax2.set_ylim(0, min(100, cumulative[-1] * 100 + 5))

    # Reference lines at common thresholds
    for pct in (80, 90):
        if cumulative[-1] * 100 >= pct:
            idx = int(np.searchsorted(cumulative * 100, pct))
            ax2.axhline(pct, color='#aaa', linewidth=0.8, linestyle='--', alpha=0.5)
            ax2.axvline(idx + 1, color='#f39c12', linewidth=1, linestyle='--', alpha=0.7)
            ax2.annotate(f'{pct}% at PC {idx + 1}',
                         xy=(idx + 1, pct), xytext=(idx + 4, pct - 3),
                         fontsize=9, color='#f39c12',
                         arrowprops=dict(arrowstyle='->', color='#f39c12', lw=1))

    # Mark suggested elbow
    if suggested_pcs is not None and 1 <= suggested_pcs <= n_pcs:
        cum_val = cumulative[suggested_pcs - 1] * 100
        ax1.axvline(suggested_pcs, color='#2ecc71', linewidth=2.5,
                    linestyle='-', alpha=0.9, zorder=5)
        ax1.annotate(f'Elbow: PC {suggested_pcs}\n({cum_val:.0f}% var)',
                     xy=(suggested_pcs, variance_ratio[suggested_pcs - 1] * 100),
                     xytext=(suggested_pcs + 3, ax1.get_ylim()[1] * 0.8),
                     fontsize=10, fontweight='bold', color='#2ecc71',
                     arrowprops=dict(arrowstyle='->', color='#2ecc71', lw=1.5))

    fig.legend(loc='upper right', bbox_to_anchor=(0.88, 0.95), fontsize=9)
    fig.tight_layout()
    return fig
