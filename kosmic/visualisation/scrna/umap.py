"""
UMAP/t-SNE Publication Plots
==============================
Publication-quality embedding scatter plots.
"""

import numpy as np


# Distinct, print-friendly default (Tableau-style, saturated hues first).
_DEFAULT_PALETTE = [
    '#4E79A7', '#F28E2B', '#59A14F', '#E15759', '#B07AA1',
    '#76B7B2', '#EDC948', '#FF9DA7', '#9C755F', '#B6992D',
    '#499894', '#D37295', '#8CD17D', '#A0CBE8', '#FFBE7D',
    '#86BCB6', '#F1CE63', '#D4A6C8', '#79706E', '#BAB0AC',
    '#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E',
    '#E6AB02', '#A6761D', '#666666', '#386CB0', '#F0027F',
]
# Selectable categorical schemes: friendly name -> matplotlib qualitative cmap.
_CMAP_SCHEMES = {
    'Tableau 20': 'tab20', 'Bold': 'Set1', 'Dark': 'Dark2',
    'Paired': 'Paired', 'Pastel': 'Set3', 'Classic 10': 'tab10',
}
# Order shown in the UI.
PALETTE_NAMES = ['Default'] + list(_CMAP_SCHEMES.keys())


def _soft_palette(n, scheme='Default'):
    """Return n distinct colours from the named scheme (cycled if n exceeds it)."""
    if scheme in _CMAP_SCHEMES:
        import matplotlib.pyplot as plt
        import matplotlib.colors as mc
        base = [mc.to_hex(c) for c in plt.get_cmap(_CMAP_SCHEMES[scheme]).colors]
    else:
        base = _DEFAULT_PALETTE
    return [base[i % len(base)] for i in range(n)]


def create_embedding_plot(coords, labels, title='', embedding_name='UMAP',
                          figsize=(12, 10), point_size=3, alpha=0.7,
                          dpi=150, rasterized=True, font_sizes=None,
                          scheme='Default'):
    """Create a single embedding scatter plot colored by categorical labels.

    Parameters
    ----------
    coords : numpy.ndarray
        (n_cells, 2) embedding coordinates.
    labels : array-like
        Categorical labels for each cell.
    title : str
        Plot title.
    embedding_name : str
        'UMAP' or 't-SNE' for axis labels.
    figsize : tuple
        Figure size in inches.
    point_size : float
        Scatter point size.
    alpha : float
        Point transparency.
    dpi : int
        Resolution for rasterized elements.
    rasterized : bool
        Rasterize scatter for publication PDFs.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # Keep text/axes as editable vector objects in PDF/SVG (Illustrator);
    # only the dense scatter is rasterized (below).
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    fig, ax = plt.subplots(figsize=figsize)

    labels = np.asarray(labels, dtype=str)
    unique_vals = sorted(set(labels))
    n_unique = len(unique_vals)
    palette = _soft_palette(n_unique, scheme)

    # Honour the requested point size exactly (may be sub-1 for dense
    # embeddings); large cell counts need genuinely tiny dots.
    s = float(point_size)

    for j, val in enumerate(unique_vals):
        mask = labels == val
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            c=[palette[j]], label=f"{val} ({mask.sum():,})",
            s=s, alpha=alpha, rasterized=rasterized, edgecolors='none',
            linewidths=0,
        )

    ax.set_xlabel(f'{embedding_name} 1', fontsize=font_sizes['axis_label'])
    ax.set_ylabel(f'{embedding_name} 2', fontsize=font_sizes['axis_label'])
    if title:
        ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold')
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Adaptive legend — scale from the configured legend size. Legend markers
    # are enlarged independently of the (possibly tiny) plotted point size so
    # they stay visible.
    base_legend = font_sizes['legend']
    ncol = 1 if n_unique <= 20 else 2 if n_unique <= 40 else 3
    legend_fs = (base_legend if n_unique <= 20
                 else base_legend * 0.8 if n_unique <= 40
                 else base_legend * 0.6)
    ax.legend(
        bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False,
        fontsize=legend_fs, ncol=ncol,
        markerscale=max(1.0, 6.0 / max(s, 0.5)),
        handletextpad=0.3, labelspacing=0.35, columnspacing=0.8,
    )

    plt.tight_layout()
    return fig


def create_highlight_plot(coords, labels, highlight_values, title='',
                          embedding_name='UMAP', figsize=(12, 10),
                          point_size=3, alpha=0.8, dpi=150, rasterized=True,
                          font_sizes=None, scheme='Default'):
    """Embedding plot with only the chosen categories coloured, the rest grey.

    Parameters
    ----------
    coords : numpy.ndarray
        (n_cells, 2) embedding coordinates.
    labels : array-like
        Categorical labels.
    highlight_values : list
        Value(s) to colour; every other cell is drawn grey in the background.
    point_size : float
        Size of the highlighted points (background uses a smaller size).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    if font_sizes is None:
        from kosmic.visualisation import default_font_sizes
        font_sizes = default_font_sizes()

    fig, ax = plt.subplots(figsize=figsize)

    labels = np.asarray(labels, dtype=str)
    highlight = [str(v) for v in highlight_values]
    s = float(point_size)
    bg = np.isin(labels, highlight, invert=True)

    # Background (everything not highlighted) in light grey, drawn first.
    ax.scatter(coords[bg, 0], coords[bg, 1], c='#d9d9d9',
               s=max(s * 0.5, 0.3), alpha=0.5, rasterized=rasterized,
               edgecolors='none', linewidths=0)

    # Each highlighted category in its own colour, on top.
    palette = _soft_palette(len(highlight), scheme)
    for j, val in enumerate(highlight):
        m = labels == val
        ax.scatter(coords[m, 0], coords[m, 1], c=[palette[j]], s=s, alpha=alpha,
                   rasterized=rasterized, edgecolors='none', linewidths=0,
                   label=f"{val} ({m.sum():,})")

    ax.set_xlabel(f'{embedding_name} 1', fontsize=font_sizes['axis_label'])
    ax.set_ylabel(f'{embedding_name} 2', fontsize=font_sizes['axis_label'])
    if title:
        ax.set_title(title, fontsize=font_sizes['title'], fontweight='bold')
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False,
              fontsize=font_sizes['legend'],
              markerscale=max(1.0, 6.0 / max(s, 0.5)), handletextpad=0.3)

    plt.tight_layout()
    return fig


def create_overview_figure(coords, obs_df, columns, embedding_name='UMAP',
                           figsize_per_panel=(8, 7), font_sizes=None,
                           scheme='Default'):
    """Create a 2-column overview grid of embedding plots.

    Parameters
    ----------
    coords : numpy.ndarray
        (n_cells, 2) embedding coordinates.
    obs_df : pandas.DataFrame
        Observation metadata.
    columns : list of str
        Column names to plot (max 4).
    embedding_name : str
        'UMAP' or 't-SNE'.
    figsize_per_panel : tuple
        Size per subplot panel.

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

    n_plots = min(len(columns), 4)
    n_rows = (n_plots + 1) // 2
    fig_w = figsize_per_panel[0] * 2
    fig_h = figsize_per_panel[1] * n_rows

    fig, axes = plt.subplots(n_rows, 2, figsize=(fig_w, fig_h))
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    for i, col in enumerate(columns[:4]):
        row, c = i // 2, i % 2
        ax = axes[row, c]
        labels = obs_df[col].astype(str).values
        unique_vals = sorted(set(labels))
        palette = _soft_palette(len(unique_vals), scheme)

        for j, val in enumerate(unique_vals):
            mask = labels == val
            ax.scatter(coords[mask, 0], coords[mask, 1],
                       c=[palette[j]], label=val, s=2, alpha=0.6, rasterized=True)

        ax.set_title(col, fontsize=font_sizes['axis_label'], fontweight='bold')
        ax.set_xlabel(f'{embedding_name}_1', fontsize=font_sizes['tick'])
        ax.set_ylabel(f'{embedding_name}_2', fontsize=font_sizes['tick'])
        ax.set_aspect('equal')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Hide empty axes
    for i in range(n_plots, n_rows * 2):
        row, c = i // 2, i % 2
        axes[row, c].set_visible(False)

    plt.tight_layout()
    return fig


def create_scatter_plots(adata, output_dir, color_by, embedding='umap', dataset_name=''):
    """
    Create and save scatter plots for each color_by column, plus an overview.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have the embedding in obsm (X_umap or X_tsne).
    output_dir : str or Path
        Directory to save PNG files.
    color_by : list of str
        Column names from adata.obs to color by.
    embedding : str
        'umap' or 'tsne'.
    dataset_name : str
        Prefix for saved filenames.

    Returns
    -------
    list of str
        Paths to saved PNG files.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pathlib import Path

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    saved_paths = []

    embed_key = f'X_{embedding}'
    if embed_key not in adata.obsm:
        return []

    coords = adata.obsm[embed_key]

    # Individual plots
    for col in color_by:
        if col not in adata.obs.columns:
            continue

        fig, ax = plt.subplots(figsize=(10, 8))
        values = adata.obs[col].astype(str)
        unique_vals = sorted(values.unique())

        if len(unique_vals) <= 20:
            cmap = plt.cm.get_cmap('tab20' if len(unique_vals) > 10 else 'tab10')
            colors = {v: cmap(j / len(unique_vals)) for j, v in enumerate(unique_vals)}
            for val in unique_vals:
                mask = values == val
                ax.scatter(
                    coords[mask.values, 0], coords[mask.values, 1],
                    c=[colors[val]], label=f"{val} ({mask.sum():,})",
                    s=2, alpha=0.6, rasterized=True,
                )
            ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
                      markerscale=4, fontsize=9)
        else:
            ax.scatter(coords[:, 0], coords[:, 1],
                       c='steelblue', s=1, alpha=0.3, rasterized=True)

        ax.set_xlabel(f'{embedding.upper()}1')
        ax.set_ylabel(f'{embedding.upper()}2')
        ax.set_title(f'{dataset_name} - {col}')
        ax.set_aspect('equal', adjustable='datalim')
        for spine in ax.spines.values():
            spine.set_visible(False)

        plt.tight_layout()
        filename = f"{dataset_name}_{embedding}_{col.replace(' ', '_')}.png"
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()
        saved_paths.append(str(filepath))

    # Overview figure (up to 4 panels)
    if len(color_by) >= 2:
        valid_cols = [c for c in color_by[:4] if c in adata.obs.columns]
        if valid_cols:
            n_cols_grid = min(2, len(valid_cols))
            n_rows_grid = (len(valid_cols) + n_cols_grid - 1) // n_cols_grid

            fig, axes = plt.subplots(n_rows_grid, n_cols_grid,
                                     figsize=(8 * n_cols_grid, 7 * n_rows_grid))
            axes = np.array(axes).flatten() if n_rows_grid * n_cols_grid > 1 else [axes]

            for idx, col in enumerate(valid_cols):
                ax = axes[idx]
                values = adata.obs[col].astype(str)
                unique_vals = sorted(values.unique())

                if len(unique_vals) <= 15:
                    cmap = plt.cm.get_cmap('tab10')
                    for j, val in enumerate(unique_vals):
                        mask = values == val
                        ax.scatter(coords[mask.values, 0], coords[mask.values, 1],
                                   c=[cmap(j % 10)], label=val, s=1, alpha=0.5,
                                   rasterized=True)
                    ax.legend(loc='upper right', fontsize=7, markerscale=3)
                else:
                    ax.scatter(coords[:, 0], coords[:, 1],
                               c='steelblue', s=1, alpha=0.3, rasterized=True)

                ax.set_title(col)
                ax.set_xlabel(f'{embedding.upper()}1')
                ax.set_ylabel(f'{embedding.upper()}2')

            # Hide empty axes
            for idx in range(len(valid_cols), len(axes)):
                axes[idx].set_visible(False)

            plt.tight_layout()
            overview_path = output_dir / f"{dataset_name}_{embedding}_overview.png"
            plt.savefig(overview_path, dpi=150, bbox_inches='tight', facecolor='white')
            plt.close()
            saved_paths.insert(0, str(overview_path))

    return saved_paths
