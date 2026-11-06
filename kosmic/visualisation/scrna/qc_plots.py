"""
QC Plots
========
Pure matplotlib functions for quality-control visualisations. No GUI dependency.
All functions return matplotlib Figures (or None when the input cannot produce one).
"""

import textwrap
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


#: The QC metrics drawn by default, and the optional ones drawn when present.
QC_METRICS = ("n_genes_by_counts", "total_counts", "pct_counts_mt")
EXTRA_QC_METRICS = ("decontX_contamination", "doublet_score")

QC_METRIC_LABELS = {
    "n_genes_by_counts": "Genes per cell",
    "total_counts": "Counts per cell",
    "pct_counts_mt": "Mitochondrial counts (%)",
    "decontX_contamination": "Ambient RNA, DecontX (%)",
    "doublet_score": "Doublet score",
}
# Stored as a fraction, drawn as a percentage.
_PERCENT_OF_FRACTION = {"decontX_contamination"}
# Drawn on a log axis when 'log_counts' is set.
_LOG_METRICS = {"n_genes_by_counts", "total_counts"}
# A violin's outline is estimated from at most this many cells (fewer when
# there are many violins, as by sample); the median line always uses every
# cell.
_MAX_CELLS_PER_VIOLIN = 20_000
_MAX_CELLS_PER_VIOLIN_MANY = 5_000


def _log_ticks(ax) -> None:
    """Ticks on an axis drawn in log10 units, labelled in the original
    units: 1, 2 and 5 times powers of ten when that gives at most six, else
    1 and 3 times (nearly even steps on a log axis), else powers of ten."""
    lo, hi = ax.get_ylim()
    decades = range(int(np.floor(lo)), int(np.ceil(hi)) + 1)
    for steps in ((1, 2, 5), (1, 3), (1,)):
        ticks = [np.log10(m * 10.0 ** e) for e in decades for m in steps
                 if lo <= np.log10(m * 10.0 ** e) <= hi]
        if len(ticks) <= 6:
            break
    ax.set_yticks(ticks, [f'{10 ** t:,.0f}' for t in ticks])


def _violin_values(values: np.ndarray, rng, cap: int = _MAX_CELLS_PER_VIOLIN) -> np.ndarray:
    if len(values) > cap:
        return rng.choice(values, cap, replace=False)
    return values


def create_qc_violin_plots(
    obs_df: pd.DataFrame,
    metrics: Optional[List[str]] = None,
    thresholds: Optional[Dict[str, Dict[str, float]]] = None,
    groups: Optional[pd.Series] = None,
    group_order: Optional[List[str]] = None,
    group_colours: Optional[Dict[str, str]] = None,
    log_counts: bool = False,
    y_limits: Optional[Dict[str, tuple]] = None,
    figsize: Optional[tuple] = None,
    font_sizes: Optional[dict] = None,
):
    """Violin plots of cell-level QC metrics, pooled or one violin per group.

    Pooled, or with a few groups (such as condition), the metrics sit side by
    side with the metric as each panel's title. With many groups (such as
    samples) the metrics are stacked, one row each, with the groups along a
    shared x axis.

    Parameters
    ----------
    obs_df : pd.DataFrame
        adata.obs with QC metric columns.
    metrics : list of str, optional
        Columns to plot, in order. Defaults to ``QC_METRICS``. Columns not in
        ``obs_df`` are left out.
    thresholds : dict, optional
        Per-metric limits drawn as dotted lines, e.g.
        ``{"n_genes_by_counts": {"min": 300, "max": 5000}}`` (see
        ``kosmic.provenance.qc_limits``).
    groups : pd.Series, optional
        A group label per cell, aligned with ``obs_df``. None pools all cells.
    group_order : list of str, optional
        Groups to draw, in order. Defaults to the sorted labels.
    group_colours : dict, optional
        Fill colour per group. Groups without one, and the pooled violin, are
        grey.
    log_counts : bool
        Draw genes and counts per cell on a log10 axis.
    y_limits : dict, optional
        ``{metric: (low, high)}`` in the metric's units (a percentage for
        DecontX contamination); None for an end set from the data. Typing
        the same ranges for several studies gives their figures the same
        scales.
    figsize, font_sizes
        Figure size in inches and the font-size dict
        (``default_font_sizes()`` when None).

    Returns
    -------
    matplotlib.figure.Figure or None
        ``None`` when none of the requested metrics are present.
    """
    from matplotlib.ticker import FuncFormatter

    from kosmic.visualisation.style import FG, NS_COLOR

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()
    metrics = [m for m in (metrics or QC_METRICS) if m in obs_df.columns]
    if not metrics:
        return None

    if groups is None:
        order = ['All cells']
        labels = pd.Series('All cells', index=obs_df.index)
    else:
        labels = pd.Series(np.asarray(groups).astype(str), index=obs_df.index)
        order = [g for g in (group_order or sorted(labels.unique())) if (labels == g).any()]
    colours = group_colours or {}
    stacked = len(order) > 4
    n = len(metrics)
    if stacked:
        size = figsize or (max(6.0, 0.25 * len(order) + 1.5), 1.9 * n + 0.8)
        fig, axes = plt.subplots(n, 1, figsize=size, sharex=True, squeeze=False)
    else:
        size = figsize or (3.2 * n, 3.2)
        fig, axes = plt.subplots(1, n, figsize=size, squeeze=False)
    axes = axes.ravel()
    fig.patch.set_facecolor('white')
    rng = np.random.default_rng(0)
    positions = np.arange(1, len(order) + 1)

    for ax, metric in zip(axes, metrics):
        values = pd.to_numeric(obs_df[metric], errors='coerce')
        if metric in _PERCENT_OF_FRACTION:
            values = values * 100
        log = log_counts and metric in _LOG_METRICS
        if log:
            values = np.log10(values.where(values > 0))
        by_group = values.groupby(labels.to_numpy(), observed=True)
        per_group = [by_group.get_group(g).dropna().to_numpy() if g in by_group.groups
                     else np.array([]) for g in order]
        drawn = [(p, v) for p, v in zip(positions, per_group) if len(v) > 1]
        if drawn:
            cap = _MAX_CELLS_PER_VIOLIN_MANY if stacked else _MAX_CELLS_PER_VIOLIN
            vp = ax.violinplot([_violin_values(v, rng, cap) for _, v in drawn],
                               positions=[p for p, _ in drawn], widths=0.8,
                               showmedians=False, showextrema=False)
            for body, (p, _v) in zip(vp['bodies'], drawn):
                body.set_facecolor(colours.get(order[p - 1], NS_COLOR))
                body.set_edgecolor(FG)
                body.set_linewidth(0.4 if stacked else 0.6)
                body.set_alpha(1.0)
            for p, v in drawn:
                ax.hlines(np.median(v), p - 0.3, p + 0.3, color=FG,
                          linewidth=1.0 if stacked else 1.5)
        for key, value in (thresholds or {}).get(metric, {}).items():
            if log and value <= 0:
                continue
            ax.axhline(np.log10(value) if log else value, color=FG, linestyle=':',
                       linewidth=0.8)
        low, high = (y_limits or {}).get(metric, (None, None))
        if log:
            low = np.log10(low) if low and low > 0 else None
            high = np.log10(high) if high and high > 0 else None
        if low is not None or high is not None:
            ax.set_ylim(low, high)
        if log:
            _log_ticks(ax)
        else:
            ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _pos: f'{y:,g}'))
        name = QC_METRIC_LABELS.get(metric, metric)
        if stacked:
            ax.set_ylabel(textwrap.fill(name, 14), fontsize=font_sizes['axis_label'], color=FG)
        else:
            # Wrapped so neighbouring titles do not run into each other.
            ax.set_title(textwrap.fill(name, 16), fontsize=font_sizes['title'],
                         fontweight='bold', color=FG)
        ax.tick_params(labelsize=font_sizes['tick'], colors=FG)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        ax.spines['left'].set_color(FG)
        ax.spines['bottom'].set_color(FG)
        ax.set_xlim(0.4, len(order) + 0.6)
        if groups is None:
            ax.set_xticks([])
            ax.spines['bottom'].set_visible(False)
        elif not stacked:
            ax.set_xticks(positions, order, fontsize=font_sizes['tick'])
            ax.tick_params(axis='x', length=0)
    if stacked:
        axes[-1].set_xticks(positions, order, rotation=90, fontsize=font_sizes['tick'] * 0.8)
        axes[-1].tick_params(axis='x', length=0)
    fig.tight_layout(pad=1.0)
    if stacked:
        fig.align_ylabels(axes)
    return fig


def create_cell_filtering_plot(
    steps: List[dict],
    figsize: Optional[tuple] = None,
    font_sizes: Optional[dict] = None,
):
    """Flow chart of the cells left after each filtering step.

    The cells remaining after each step are in a box down the left; between
    two boxes an arrow leads down, and a box beside it names the step and
    the cells it removed (a CONSORT-style diagram).

    Parameters
    ----------
    steps : list of dict
        ``{'step', 'n_cells', 'n_removed'}`` per step, the starting cells
        first, as ``kosmic.provenance.cell_filtering_steps`` returns them.
    figsize, font_sizes
        Figure size in inches and the font-size dict
        (``default_font_sizes()`` when None).

    Returns
    -------
    matplotlib.figure.Figure or None
        ``None`` when there are no steps.
    """
    from matplotlib.patches import FancyBboxPatch

    from kosmic.visualisation.style import FG, FG_SECONDARY, MISSING_COLOR

    if not steps:
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()
    n = len(steps)
    fig, ax = plt.subplots(figsize=figsize or (6.0, 1.05 * n + 0.3))
    fig.patch.set_facecolor('white')
    ax.axis('off')
    # Data units: x 0-10, each step 2 units tall, the first at the top.
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.7, 2 * (n - 1) + 0.7)
    main_x, main_w, side_x, side_w, box_h = 0.2, 3.8, 5.4, 4.4, 0.9

    def box(x, y, w, text, fill, colour, weight='normal'):
        ax.add_patch(FancyBboxPatch((x, y - box_h / 2), w, box_h,
                                    boxstyle='round,pad=0,rounding_size=0.12',
                                    facecolor=fill, edgecolor=colour, linewidth=0.8))
        ax.text(x + w / 2, y, text, ha='center', va='center', color=colour,
                fontsize=font_sizes['annotation'], fontweight=weight)

    tops = [2 * (n - 1 - i) for i in range(n)]
    for i, (s, y) in enumerate(zip(steps, tops)):
        label = (f"{s['step']}\n{s['n_cells']:,} cells" if i == 0
                 else f"{s['n_cells']:,} cells")
        box(main_x, y, main_w, label, 'white', FG, 'bold' if i in (0, n - 1) else 'normal')
        if i == 0:
            continue
        cx, mid = main_x + main_w / 2, (tops[i - 1] + y) / 2
        ax.annotate('', xy=(cx, y + box_h / 2), xytext=(cx, tops[i - 1] - box_h / 2),
                    arrowprops={'arrowstyle': '-|>', 'color': FG, 'lw': 0.8,
                                'shrinkA': 0, 'shrinkB': 0})
        ax.annotate('', xy=(side_x, mid), xytext=(cx, mid),
                    arrowprops={'arrowstyle': '-|>', 'color': FG_SECONDARY, 'lw': 0.8,
                                'shrinkA': 0, 'shrinkB': 0})
        box(side_x, mid, side_w, f"{s['step']}\n{s['n_removed']:,} removed",
            MISSING_COLOR, FG_SECONDARY)
    fig.tight_layout(pad=0.5)
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
