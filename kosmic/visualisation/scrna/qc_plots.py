"""
QC Plots
========
Pure matplotlib functions for quality-control visualisations. No GUI dependency.
All functions return matplotlib Figures (or None when the input cannot produce one).
"""

from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _apply_theme(fig, axes, dark_mode: bool):
    """Apply dark/light theme colors to figure and axes."""
    if dark_mode:
        bg = "#1e1e1e"
        fg = "#cccccc"
        grid = "#333333"
    else:
        bg = "#ffffff"
        fg = "#333333"
        grid = "#dddddd"

    fig.patch.set_facecolor(bg)
    for ax in axes:
        ax.set_facecolor(bg)
        ax.tick_params(colors=fg, labelsize=8)
        ax.xaxis.label.set_color(fg)
        ax.yaxis.label.set_color(fg)
        ax.title.set_color(fg)
        for spine in ax.spines.values():
            spine.set_color(grid)


def create_qc_violin_plots(
    obs_df: pd.DataFrame,
    metrics: Optional[List[str]] = None,
    thresholds: Optional[Dict[str, Dict[str, float]]] = None,
    dark_mode: bool = True,
    figsize: Optional[tuple] = None,
    font_sizes: Optional[dict] = None,
):
    """Side-by-side violin plots for QC metrics.

    Parameters
    ----------
    obs_df : pd.DataFrame
        adata.obs with QC metric columns.
    metrics : list of str, optional
        Columns to plot. Defaults to n_genes_by_counts, total_counts, pct_counts_mt.
    thresholds : dict, optional
        Per-metric threshold lines, e.g.
        {"n_genes_by_counts": {"min": 200, "max": 6000}, "pct_counts_mt": {"max": 20}}
    dark_mode : bool
        Use dark theme colors.

    Returns
    -------
    matplotlib.figure.Figure or None
        ``None`` when none of the requested metrics are present.
    """
    if metrics is None:
        metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]

    available = [m for m in metrics if m in obs_df.columns]
    if not available:
        return None

    n = len(available)
    fig_size = figsize if figsize else (3.2 * n, 3.0)
    fig, axes = plt.subplots(1, n, figsize=fig_size)
    if n == 1:
        axes = [axes]

    _apply_theme(fig, axes, dark_mode)
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()
    from kosmic.visualisation.style import FG, NS_COLOR

    labels = {
        "n_genes_by_counts": "Genes per cell",
        "total_counts": "Total counts",
        "pct_counts_mt": "MT %",
    }

    for ax, metric in zip(axes, available):
        data = obs_df[metric].dropna().values
        title_kw = dict(fontsize=font_sizes['title'], fontweight='bold')
        if len(data) == 0:
            ax.set_title(labels.get(metric, metric), **title_kw)
            continue

        vp = ax.violinplot(data, showmedians=True, showextrema=False)
        for body in vp["bodies"]:
            body.set_facecolor("#4fc1ff" if dark_mode else NS_COLOR)
            body.set_edgecolor("none" if dark_mode else FG)
            body.set_linewidth(0.6)
            body.set_alpha(0.7 if dark_mode else 1.0)
        vp["cmedians"].set_color("#ff6b6b" if dark_mode else FG)
        vp["cmedians"].set_linewidth(1.5)

        ax.set_title(labels.get(metric, metric), **title_kw)
        ax.set_xticks([])
        ax.tick_params(labelsize=font_sizes['tick'])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        if not dark_mode:
            ax.spines['left'].set_color(FG)

        if thresholds and metric in thresholds:
            t = thresholds[metric]
            for key in ("min", "max"):
                if key in t:
                    ax.axhline(t[key], color=FG, linestyle=":", linewidth=1)

    fig.tight_layout(pad=1.0)
    return fig


def create_scrublet_histogram(
    doublet_scores: np.ndarray,
    threshold: Optional[float],
    dark_mode: bool = True,
):
    """Histogram of Scrublet doublet scores, with the threshold line when
    there is a single one (per-sample runs have one threshold each).

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, ax = plt.subplots(figsize=(4.5, 3.0))
    _apply_theme(fig, [ax], dark_mode)

    color = "#4fc1ff" if dark_mode else "#1976d2"
    ax.hist(doublet_scores, bins=50, color=color, alpha=0.7, edgecolor="none")
    if threshold is not None:
        ax.axvline(threshold, color="#ff4444", linestyle="--", linewidth=1.5,
                   label=f"Threshold: {threshold:.3f}")
        ax.legend(fontsize=8, framealpha=0.6)
    ax.set_xlabel("Doublet score", fontsize=9)
    ax.set_ylabel("Count", fontsize=9)
    ax.set_title("Scrublet Doublet Scores", fontsize=10)

    fig.tight_layout(pad=1.0)
    return fig
