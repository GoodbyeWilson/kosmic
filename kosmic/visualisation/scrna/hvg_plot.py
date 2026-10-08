"""
Highly Variable Genes Plot
==========================
Mean expression against normalised variance (or dispersion) for every
gene, with the genes selected as highly variable highlighted.
"""

import numpy as np

from kosmic.visualisation.style import NS_COLOR, UP_COLOR

# (y column, y-axis label) for each scanpy flavour KOSMIC uses.
_Y_BY_FLAVOR = {
    'seurat_v3': ('variances_norm', 'Normalised variance'),
    'seurat': ('dispersions_norm', 'Normalised dispersion'),
}


def create_hvg_plot(var, flavor='seurat_v3', figsize=None, font_sizes=None):
    """Scatter of each gene's mean against its normalised variability.

    *var* is ``adata.var`` after ``sc.pp.highly_variable_genes``; it needs
    ``means``, ``highly_variable`` and the flavour's variability column.
    The x axis is logarithmic because gene means span several orders of
    magnitude. Returns None when the columns are missing.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ycol, ylabel = _Y_BY_FLAVOR.get(flavor, _Y_BY_FLAVOR['seurat'])
    if not {'means', 'highly_variable', ycol} <= set(var.columns):
        return None
    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    mean = var['means'].to_numpy(dtype=float)
    y = var[ycol].to_numpy(dtype=float)
    hvg = var['highly_variable'].to_numpy(dtype=bool)
    ok = (mean > 0) & np.isfinite(y)

    fig, ax = plt.subplots(figsize=figsize or (6, 4.5))
    ax.scatter(mean[ok & ~hvg], y[ok & ~hvg], s=3, color=NS_COLOR,
               edgecolors='none', label='Other genes')
    ax.scatter(mean[ok & hvg], y[ok & hvg], s=3, color=UP_COLOR,
               edgecolors='none', label='Highly variable')
    ax.set_xscale('log')
    ax.set_xlabel('Mean expression', fontsize=font_sizes['axis_label'])
    ax.set_ylabel(ylabel, fontsize=font_sizes['axis_label'])
    ax.set_title(f'Highly variable genes ({hvg.sum():,} of {len(var):,})',
                 fontsize=font_sizes['title'], fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=font_sizes['tick'])
    ax.legend(frameon=False, fontsize=font_sizes['legend'], markerscale=4)
    fig.tight_layout()
    return fig
