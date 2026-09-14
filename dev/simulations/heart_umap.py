"""
Generate a heart-shaped UMAP for screenshots.

Usage:
    python simulations/heart_umap.py path/to/your.h5ad output.h5ad

If no input is given, generates a demo dataset with ~5000 cells and 8 clusters.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # repo root: dev/<this dir>/<file>

import numpy as np
import scanpy as sc
import anndata as ad


def heart_curve(t):
    """Parametric heart curve."""
    x = 16 * np.sin(t) ** 3
    y = 13 * np.cos(t) - 5 * np.cos(2 * t) - 2 * np.cos(3 * t) - np.cos(4 * t)
    return x, y


def warp_to_heart(umap_coords, strength=0.85):
    """Warp 2D coordinates onto a heart shape while preserving local structure.

    Parameters
    ----------
    umap_coords : ndarray (n, 2)
        Original UMAP coordinates.
    strength : float
        0.0 = no warping, 1.0 = fully on heart boundary.
        0.85 gives a nice heart shape with visible clusters.
    """
    # Centre and normalise to [-1, 1]
    coords = umap_coords.copy()
    coords -= coords.mean(axis=0)
    scale = np.abs(coords).max()
    coords /= scale

    # Convert to polar
    r = np.sqrt(coords[:, 0] ** 2 + coords[:, 1] ** 2)
    theta = np.arctan2(coords[:, 1], coords[:, 0])

    # Map theta to heart parameter t (shift so top of heart = pi/2)
    t = theta

    # Heart boundary at each angle
    hx, hy = heart_curve(t)
    heart_r = np.sqrt(hx ** 2 + hy ** 2)

    # Normalise heart radius
    heart_r_max = np.sqrt(16 ** 2 + 17 ** 2)  # approximate max
    heart_r_norm = heart_r / heart_r_max

    # Scale original radius to fit within heart boundary
    r_norm = r / r.max()

    # Blend: move points toward heart boundary
    r_new = r_norm * (1 - strength) + r_norm * heart_r_norm * strength

    # Convert back to Cartesian, using heart centre direction
    # Use the heart curve direction instead of original angle for shape
    hx_norm = hx / heart_r_max
    hy_norm = hy / heart_r_max

    # Blend direction: original angle vs heart-curve direction
    dir_x = np.cos(theta) * (1 - strength) + (hx_norm / (heart_r_norm + 1e-10)) * strength
    dir_y = np.sin(theta) * (1 - strength) + (hy_norm / (heart_r_norm + 1e-10)) * strength
    dir_len = np.sqrt(dir_x ** 2 + dir_y ** 2)
    dir_x /= dir_len + 1e-10
    dir_y /= dir_len + 1e-10

    new_coords = np.column_stack([
        dir_x * r_new * scale * 1.2,
        dir_y * r_new * scale * 1.2,
    ])

    return new_coords


def make_demo_dataset(n_cells=5000, n_clusters=8):
    """Generate a demo scRNA-seq dataset with clear clusters."""
    rng = np.random.default_rng(42)

    n_genes = 2000

    # Generate cluster centres in gene space
    cluster_labels = rng.choice(n_clusters, n_cells)

    # Base expression
    X = rng.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(np.float32)

    # Add cluster-specific marker genes
    genes_per_cluster = 50
    for c in range(n_clusters):
        mask = cluster_labels == c
        marker_start = c * genes_per_cluster
        marker_end = marker_start + genes_per_cluster
        X[mask, marker_start:marker_end] += rng.poisson(15, size=(mask.sum(), genes_per_cluster))

    gene_names = [f'Gene_{i}' for i in range(n_genes)]
    cell_names = [f'Cell_{i}' for i in range(n_cells)]

    adata = ad.AnnData(
        X=X,
        obs={'cluster': [f'Cluster {i}' for i in cluster_labels]},
        var=pd.DataFrame(index=gene_names),
    )
    adata.obs.index = cell_names
    adata.obs['cluster'] = adata.obs['cluster'].astype('category')

    return adata


if __name__ == '__main__':
    import pandas as pd

    if len(sys.argv) >= 2 and os.path.exists(sys.argv[1]):
        print(f"Loading {sys.argv[1]}...")
        adata = sc.read_h5ad(sys.argv[1])
        output = sys.argv[2] if len(sys.argv) >= 3 else 'heart_umap.h5ad'
    else:
        print("Generating demo dataset...")
        adata = make_demo_dataset()
        output = sys.argv[1] if len(sys.argv) >= 2 else 'heart_umap.h5ad'

    print(f"  {adata.n_obs} cells, {adata.n_vars} genes")

    # Preprocess if needed
    if 'X_umap' not in adata.obsm:
        print("Running preprocessing pipeline...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=min(2000, adata.n_vars))
        sc.tl.pca(adata, n_comps=min(50, adata.n_obs - 1))
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        if 'cluster' not in adata.obs.columns and 'leiden' not in adata.obs.columns:
            sc.tl.leiden(adata, resolution=0.8)

    print("Warping UMAP to heart shape...")
    original_umap = adata.obsm['X_umap'].copy()
    adata.obsm['X_umap'] = warp_to_heart(original_umap, strength=0.85)

    # Save
    print(f"Saving to {output}...")
    adata.write_h5ad(output)

    print("Done! Load this file in KOSMIC for the screenshot.")
    print("Tip: adjust 'strength' parameter (0.5-0.95) if the shape needs tuning.")
