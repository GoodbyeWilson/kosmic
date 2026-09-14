# Leiden clustering + resolution sweep.
#
# Always uses 'flavor='igraph'' + 'n_iterations=2' + 'directed=False'
# so results are reproducible across scanpy versions.
from __future__ import annotations

from collections import Counter
from typing import List, Sequence



def cluster_leiden(adata, resolution: float = 1.0, random_state: int = 0):
    """Run Leiden clustering once. Writes 'adata.obs['leiden']'.

    Requires a neighbors graph ('compute_neighbors' first).
    """
    import scanpy as sc

    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=random_state,
        flavor='igraph',
        n_iterations=2,
        directed=False,
    )
    return adata


def sweep_leiden_resolutions(
    adata,
    resolutions: Sequence[float],
    random_state: int = 0,
) -> List[dict]:
    """Run Leiden at multiple resolutions; return per-resolution summary.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have a neighbors graph.
    resolutions : sequence of float
    random_state : int
        Random seed for reproducibility.

    Returns
    -------
    list of dict
        Each dict has 'resolution', 'n_clusters', 'min_size',
        'labels' (per-cell label array).
    """
    import scanpy as sc

    results = []
    for r in resolutions:
        adata_copy = adata.copy()
        sc.tl.leiden(
            adata_copy, resolution=r, random_state=random_state,
            flavor='igraph', n_iterations=2, directed=False,
        )
        labels = adata_copy.obs['leiden'].values
        counts = Counter(labels)
        results.append({
            'resolution': r,
            'n_clusters': len(counts),
            'min_size': min(counts.values()),
            'labels': labels,
        })
    return results
