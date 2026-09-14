# Nearest-neighbor graph for clustering.
#
# Wraps 'sc.pp.neighbors' with one piece of domain logic: when no
# 'use_rep' is specified, prefer 'X_pca_harmony' over 'X_pca' so the
# graph is built on the batch-corrected representation if Harmony has
# been run.
from __future__ import annotations

from typing import Optional


def neighbors_are_current(adata, n_neighbors: int, n_pcs: int,
                          use_rep: str) -> bool:
    """True when 'adata' already holds a graph built these parameters.

    scanpy records what it used in ``uns['neighbors']['params']``; a
    graph built from a different representation (say 'X_pca' before
    Harmony ran) or a different neighbour count has to be rebuilt.
    """
    params = (adata.uns.get('neighbors') or {}).get('params')
    if not isinstance(params, dict) or 'connectivities' not in adata.obsp:
        return False
    return (params.get('n_neighbors') == n_neighbors
            and params.get('use_rep') == use_rep
            and params.get('n_pcs') == n_pcs)


def compute_neighbors(
    adata,
    n_neighbors: int = 15,
    n_pcs: int = 30,
    use_rep: Optional[str] = None,
    force: bool = False,
):
    """Compute the kNN graph used by Leiden / UMAP.

    Parameters
    ----------
    adata : anndata.AnnData
    n_neighbors : int
        Passed to 'sc.pp.neighbors'.
    n_pcs : int
        Number of PCs to use. Capped at the size of the chosen
        representation.
    use_rep : str, optional
        'obsm' key to read PCs from. When 'None', picks
        'X_pca_harmony' if present, else 'X_pca'.

    force : bool
        Recompute even when a graph built with the same parameters is
        already present. The graph is the most expensive step in the
        clustering run -- ~35s per 40k cells, minutes on a combined
        master -- and it does not depend on the Leiden resolution, so
        re-running it to try a different resolution is pure waste.

    Returns
    -------
    anndata.AnnData
        Modified in place.
    """
    import scanpy as sc

    if use_rep is None:
        use_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'

    max_pcs = adata.obsm[use_rep].shape[1]
    n_pcs = min(n_pcs, max_pcs)

    if not force and neighbors_are_current(adata, n_neighbors, n_pcs, use_rep):
        return adata

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)
    return adata
