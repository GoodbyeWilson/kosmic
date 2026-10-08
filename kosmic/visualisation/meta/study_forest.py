"""
Study Forest and Leave-One-Out Plots
====================================
One gene or pathway across the studies of a meta-analysis: each study's
effect with its 95% CI and the pooled estimate (forest), or the pooled
estimate re-computed with each study left out in turn (leave-one-out).
Also the text reports exported beside these figures.
"""

import numpy as np
import pandas as pd

from kosmic.visualisation.style import DOWN_COLOR, FG, NS_COLOR, UP_COLOR

Z95 = 1.959964
WHISKER_COLOR = '#9AA0A6'
POOLED_COLOR = '#7A7A7A'


def study_effects(item, studies):
    """Each study's log2 fold change, SE and 95% CI for *item*.

    *studies* is a list of (study name, DE table) with 'names',
    'logfoldchanges' and 'se'. Studies without a finite effect and positive
    SE for the item are left out. Returns a DataFrame in study order.
    """
    rows = []
    for name, df in studies:
        hit = df[df['names'].astype(str) == str(item)]
        if hit.empty:
            continue
        lfc = float(pd.to_numeric(hit['logfoldchanges'].iloc[0], errors='coerce'))
        se = float(pd.to_numeric(hit['se'].iloc[0], errors='coerce'))
        if np.isfinite(lfc) and np.isfinite(se) and se > 0:
            rows.append({'study': name, 'logfoldchanges': lfc, 'se': se,
                         'ci_lower': lfc - Z95 * se, 'ci_upper': lfc + Z95 * se})
    return pd.DataFrame(rows, columns=['study', 'logfoldchanges', 'se',
                                       'ci_lower', 'ci_upper'])


def _as_dict(pooled):
    """The item's meta-analysis row (a Series or dict) as a dict."""
    return {} if pooled is None else dict(pooled)


def loo_pooling_key(pooled):
    """The pooling method leave-one-out re-pools with: the one recorded in
    the meta-analysis row when it is a single analytical method, otherwise
    REML (e.g. for a consensus of several methods)."""
    from kosmic.meta_analysis.pooling_dispatch import get_analytical_pool_fn
    key = str(_as_dict(pooled).get('pooling_method', 'reml')).lower()
    try:
        get_analytical_pool_fn(key, min_studies=2)
        return key
    except KeyError:
        return 'reml'


def loo_estimates(item, studies, pooled=None):
    """The pooled estimate for *item* with each study left out in turn.

    Re-pools with the pooling function named in *pooled* (the item's row
    of the meta-analysis table). Needs at least three studies with an
    effect. Returns a DataFrame: omitted, logfoldchanges, ci_lower,
    ci_upper, i2_pct, n_studies.
    """
    eff = study_effects(item, studies)
    cols = ['omitted', 'logfoldchanges', 'ci_lower', 'ci_upper',
            'i2_pct', 'n_studies']
    if len(eff) < 3:
        return pd.DataFrame(columns=cols)
    from kosmic.meta_analysis.pooling_dispatch import get_analytical_pool_fn
    pool = get_analytical_pool_fn(loo_pooling_key(pooled), min_studies=2)
    frames = [pd.DataFrame({'names': [item], 'logfoldchanges': [r.logfoldchanges],
                            'se': [r.se]}) for r in eff.itertuples()]
    rows = []
    for i, left_out in enumerate(eff['study']):
        res = pool([f for j, f in enumerate(frames) if j != i])
        r = res.iloc[0] if len(res) else {}
        rows.append({'omitted': left_out,
                     'logfoldchanges': float(r.get('logfoldchanges', np.nan)),
                     'ci_lower': float(r.get('ci_lower', np.nan)),
                     'ci_upper': float(r.get('ci_upper', np.nan)),
                     # I2 as a percentage (the tables hold a fraction).
                     'i2_pct': 100 * float(r.get('heterogeneity_i2', np.nan)),
                     'n_studies': len(frames) - 1})
    return pd.DataFrame(rows, columns=cols)


def _band_color(lo, hi):
    """Shading for the pooled CI: red or blue when it excludes zero, grey
    when it crosses zero."""
    if lo > 0:
        return UP_COLOR
    if hi < 0:
        return DOWN_COLOR
    return NS_COLOR


def _sym_limit(values):
    m = float(np.nanmax(np.abs(values))) if len(values) else 1.0
    step = 0.5 if m <= 2 else 1.0
    return max(step, np.ceil(m * 1.1 / step) * step)


def _style_axes(ax, title, xlabel, ylabel, font_sizes):
    ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold')
    ax.set_xlabel(xlabel, fontsize=font_sizes['axis_label'], fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=font_sizes['axis_label'], fontweight='bold')
    ax.tick_params(direction='out', labelsize=font_sizes['tick'])
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
        spine.set_color(FG)


def create_study_forest(name, effects, pooled=None, title=None, figsize=None,
                        font_sizes=None, x_limit=None):
    """Forest plot: one row per study, the pooled estimate below.

    *effects* is from :func:`study_effects`; *pooled* is the item's row of
    the meta-analysis table (logfoldchanges, ci_lower, ci_upper,
    tau_squared). Study diamonds are red (up) or blue (down) and sized by
    their random-effects weight, 1 / (SE^2 + tau^2). The pooled 95% CI is
    shaded behind every row. Returns None when no study has an effect.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if effects is None or effects.empty:
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    n = len(effects)
    fig, ax = plt.subplots(figsize=figsize or (7, 1.6 + 0.55 * (n + 1)))
    y = np.arange(n)[::-1] + 1          # first study at the top, pooled at 0
    lfc = effects['logfoldchanges'].to_numpy(float)
    lo, hi = effects['ci_lower'].to_numpy(float), effects['ci_upper'].to_numpy(float)

    tau2 = float(_as_dict(pooled).get('tau_squared', 0) or 0)
    w = 1.0 / (effects['se'].to_numpy(float) ** 2 + max(tau2, 0.0))
    sizes = (6 + 6 * np.sqrt(w / w.max())) ** 2

    p_lfc = float(_as_dict(pooled).get('logfoldchanges', np.nan))
    p_lo = float(_as_dict(pooled).get('ci_lower', np.nan))
    p_hi = float(_as_dict(pooled).get('ci_upper', np.nan))
    has_pooled = np.isfinite([p_lfc, p_lo, p_hi]).all()
    if has_pooled:
        ax.axvspan(p_lo, p_hi, color=_band_color(p_lo, p_hi), alpha=0.2,
                   linewidth=0, zorder=0)
    ax.axvline(0, color=FG, linestyle='--', linewidth=0.9, zorder=1)
    ax.hlines(y, lo, hi, color=WHISKER_COLOR, linewidth=1.4, zorder=2)
    ax.scatter(lfc, y, s=sizes, marker='D', zorder=3, edgecolors=FG, linewidths=0.5,
               c=[UP_COLOR if v > 0 else DOWN_COLOR for v in lfc])
    labels = list(effects['study'].astype(str).str.replace('_', ' '))
    ticks = list(y)
    if has_pooled:
        hh = 0.3
        ax.fill([p_lo, p_lfc, p_hi, p_lfc], [0, hh, 0, -hh], facecolor=POOLED_COLOR,
                edgecolor='black', linewidth=1.0, zorder=4)
        labels.append('Pooled')
        ticks.append(0)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=font_sizes['tick'])
    ax.set_ylim(-0.7, n + 0.7)
    lim = float(x_limit) if x_limit else _sym_limit(np.r_[lo, hi, [p_lo, p_hi] if has_pooled else []])
    ax.set_xlim(-lim, lim)
    _style_axes(ax, title or f'{name} Expression Across Studies',
                'Log$_2$ Fold Change', 'Studies', font_sizes)
    fig.tight_layout()
    return fig


def create_loo_plot(name, loo, pooled=None, title=None, figsize=None,
                    font_sizes=None, x_limit=None):
    """Leave-one-out plot: one row per study left out, showing the pooled
    estimate and 95% CI of the remaining studies. The all-studies pooled
    estimate is a solid line with its CI shaded. Returns None when there
    are no folds (fewer than three studies)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if loo is None or loo.empty:
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    n = len(loo)
    fig, ax = plt.subplots(figsize=figsize or (7, 1.6 + 0.55 * n))
    y = np.arange(n)[::-1]
    lfc = loo['logfoldchanges'].to_numpy(float)
    lo, hi = loo['ci_lower'].to_numpy(float), loo['ci_upper'].to_numpy(float)

    p_lfc = float(_as_dict(pooled).get('logfoldchanges', np.nan))
    p_lo = float(_as_dict(pooled).get('ci_lower', np.nan))
    p_hi = float(_as_dict(pooled).get('ci_upper', np.nan))
    if np.isfinite([p_lfc, p_lo, p_hi]).all():
        ax.axvspan(p_lo, p_hi, color=_band_color(p_lo, p_hi), alpha=0.2,
                   linewidth=0, zorder=0)
        ax.axvline(p_lfc, color=FG, linewidth=1.0, zorder=1, label='All studies')
        ax.legend(loc='lower right', frameon=False, fontsize=font_sizes['legend'])
    ax.axvline(0, color=FG, linestyle='--', linewidth=0.9, zorder=1)
    ax.hlines(y, lo, hi, color=WHISKER_COLOR, linewidth=1.4, zorder=2)
    ax.scatter(lfc, y, s=90, marker='D', zorder=3, edgecolors=FG, linewidths=0.5,
               c=[UP_COLOR if v > 0 else DOWN_COLOR for v in lfc])
    ax.set_yticks(y)
    ax.set_yticklabels(loo['omitted'].astype(str).str.replace('_', ' '),
                       fontsize=font_sizes['tick'])
    ax.set_ylim(-0.7, n - 0.3)
    lim = float(x_limit) if x_limit else _sym_limit(np.r_[lo, hi, p_lo, p_hi])
    ax.set_xlim(-lim, lim)
    _style_axes(ax, title or f'{name}: Leave-One-Out Sensitivity',
                'Pooled Log$_2$ Fold Change', 'Study Removed', font_sizes)
    fig.tight_layout()
    return fig


def _pooled_lines(name, pooled, n_studies):
    p = _as_dict(pooled)

    def num(key, fmt):
        v = p.get(key)
        try:
            return format(float(v), fmt)
        except (TypeError, ValueError):
            return 'NA'

    def pct(key):   # the tables hold I2 as a fraction
        try:
            return f'{100 * float(p.get(key)):.1f}%'
        except (TypeError, ValueError):
            return 'NA'
    return [
        name,
        f"Pooling: {p.get('pooling_method', 'NA')}, {n_studies} studies",
        f"Pooled log2FC: {num('logfoldchanges', '.4f')} "
        f"(95% CI {num('ci_lower', '.4f')} to {num('ci_upper', '.4f')})",
        f"FDR: {num('fdr', '.3g')}",
        f"Heterogeneity: I2 = {pct('heterogeneity_i2')}, "
        f"tau2 = {num('tau_squared', '.4g')}",
        '',
    ]


def forest_report(name, effects, pooled=None):
    """Text of the numbers behind a forest plot: pooled estimate,
    heterogeneity, and each study's log2FC, SE, 95% CI and weight."""
    tau2 = float(_as_dict(pooled).get('tau_squared', 0) or 0)
    w = 1.0 / (effects['se'] ** 2 + max(tau2, 0.0))
    table = effects.assign(weight_pct=100 * w / w.sum())
    return '\n'.join(_pooled_lines(name, pooled, len(effects))) + \
        table.to_string(index=False, float_format=lambda v: f'{v:.4f}') + '\n'


def loo_report(name, loo, pooled=None):
    """Text of the numbers behind a leave-one-out plot."""
    n = int(loo['n_studies'].max()) + 1 if len(loo) else 0
    lines = _pooled_lines(name, pooled, n)
    lines.insert(-1, f'Leave-one-out re-pooled with: {loo_pooling_key(pooled)}')
    return '\n'.join(lines) + \
        loo.to_string(index=False, float_format=lambda v: f'{v:.4f}') + '\n'
