"""
PCA visualisation helpers (no GUI imports).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def find_elbow(variance_ratio: np.ndarray, min_pcs: int = 10) -> int:
    """The suggested number of PCs: the elbow of the variance curve.

    The elbow is the point of maximum distance from the straight line
    joining the first and last PC (the Kneedle method), taken on the
    log of the variance so the knee marks where the decay slows rather
    than the corner of the first few PCs' steep drop (on linear scale a
    heart snRNA-seq scree "knees" at PC 6-7; on log scale at PC 18),
    floored at *min_pcs*. See :func:`elbow_report` for whether an elbow
    was found.

    Earlier versions also took the PC at which cumulative variance
    reached 80% and returned the largest of the estimates. Single-cell
    data never gets near 80% of total variance within the PCs computed
    (a heart snRNA-seq study sits at 50-65% by PC 80), so that estimate
    was always the last PC and the "elbow" was reported at the ceiling.

    Parameters
    ----------
    variance_ratio : array-like
        Per-PC explained variance ratio.
    min_pcs : int
        Minimum PCs to suggest (default 10).

    Returns
    -------
    int
        1-based PC number at the suggested cutoff.
    """
    return elbow_report(variance_ratio, min_pcs)['suggested']


def elbow_report(variance_ratio: np.ndarray, min_pcs: int = 10) -> dict:
    """Elbow of the variance curve, with whether it is a real one.

    Returns a dict: ``elbow`` (1-based PC of maximum distance from the
    end-to-end line), ``found`` (False when the elbow lands on the first
    or last PC, i.e. the curve has no knee within the computed PCs),
    ``suggested`` (``elbow`` floored at *min_pcs*; *min_pcs* when not
    found), ``cumulative`` (cumulative variance fraction per PC).
    """
    ratio = np.asarray(variance_ratio, dtype=float)
    n = len(ratio)
    if n < 3:
        return {'elbow': n, 'found': False, 'suggested': n,
                'cumulative': np.cumsum(ratio) if n else np.array([])}
    y = np.log(np.clip(ratio, 1e-12, None))

    # Method 1: max-distance-to-line (classic elbow)
    x = np.linspace(0, 1, n)
    p1 = np.array([x[0], y[0]])
    p2 = np.array([x[-1], y[-1]])
    line_vec = p2 - p1
    line_len = np.linalg.norm(line_vec)
    if line_len == 0:
        return {'elbow': n, 'found': False, 'suggested': min(min_pcs, n),
                'cumulative': np.cumsum(ratio)}

    line_unit = line_vec / line_len
    dists = np.abs(np.cross(line_unit, p1 - np.column_stack([x, y])))
    elbow_pc = int(np.argmax(dists)) + 1  # 1-based
    # A knee is real when the curve bows away from its end-to-end line by
    # more than 1% of that line's length (a scree with a knee sits at
    # 3-40%; a straight line at ~0) and the knee is not an end point.
    found = bool(1 < elbow_pc < n and dists.max() / line_len > 0.01)
    suggested = min(max(elbow_pc, min_pcs), n) if found else min(min_pcs, n)
    return {'elbow': elbow_pc, 'found': found, 'suggested': suggested,
            'cumulative': np.cumsum(ratio)}


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
