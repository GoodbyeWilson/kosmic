# Clustree-style cluster stability analysis.
#
# Runs Leiden at multiple resolutions and computes a barycentric layout
# (parents below, children above) so transitions between resolutions are
# visualisable as a Sankey-like diagram. The companion renderer lives in
# 'kosmic/visualisation/clustree.py'.
from __future__ import annotations

from collections import Counter
from typing import Sequence

import numpy as np


def compute_clustree_data(adata, resolutions: Sequence[float], random_state: int = 0) -> dict:
    """Compute per-level Leiden + barycentric x-positions + edges.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have neighbors graph computed.
    resolutions : sequence of float
        Resolutions to evaluate (sorted ascending recommended).
    random_state : int

    Returns
    -------
    dict
        'leiden_per_level' -- list of '(resolution, labels_array)'
        'level_positions' -- list of '{cluster_id: x_in_[0,1]}'
        'edges' -- list of dicts with 'source_level', 'target_level',
        'source_cluster', 'target_cluster', 'count', 'proportion',
        'is_dominant'.
        'node_sizes' -- list of '{cluster_id: cell_count}'
    """
    import scanpy as sc

    leiden_per_level = []
    for r in resolutions:
        adata_copy = adata.copy()
        sc.tl.leiden(
            adata_copy, resolution=r, random_state=random_state,
            flavor='igraph', n_iterations=2, directed=False,
        )
        labels = adata_copy.obs['leiden'].values.astype(str)
        leiden_per_level.append((r, labels))

    # Barycentric layout: each cluster's x position is the mean x of its
    # parent cells at the previous level. Re-rank to minimise crossings.
    level_positions = []
    for level_idx, (_r, labels) in enumerate(leiden_per_level):
        unique_clusters = sorted(set(labels), key=lambda x: int(x) if x.isdigit() else x)
        n_clusters = len(unique_clusters)

        if level_idx == 0 or n_clusters <= 1:
            if n_clusters <= 1:
                pos = {c: 0.5 for c in unique_clusters}
            else:
                pos = {c: i / (n_clusters - 1) for i, c in enumerate(unique_clusters)}
        else:
            prev_labels = leiden_per_level[level_idx - 1][1]
            prev_pos = level_positions[level_idx - 1]
            cell_x = np.array([prev_pos.get(pl, 0.5) for pl in prev_labels])
            means = {}
            for c in unique_clusters:
                mask = labels == c
                means[c] = cell_x[mask].mean() if mask.any() else 0.5

            ranked = sorted(unique_clusters, key=lambda c: means[c])
            pos = {c: i / (n_clusters - 1) for i, c in enumerate(ranked)}

        level_positions.append(pos)

    edges = []
    node_sizes = []
    for level_idx, (_r, labels) in enumerate(leiden_per_level):
        counts = Counter(labels)
        node_sizes.append(dict(counts))

        if level_idx > 0:
            prev_labels = leiden_per_level[level_idx - 1][1]
            for prev_cluster in set(prev_labels):
                mask = prev_labels == prev_cluster
                child_labels = labels[mask]
                child_counts = Counter(child_labels)
                total = len(child_labels)
                for child_cluster, count in child_counts.items():
                    proportion = count / total
                    edges.append({
                        'source_level': level_idx - 1,
                        'target_level': level_idx,
                        'source_cluster': prev_cluster,
                        'target_cluster': child_cluster,
                        'count': count,
                        'proportion': proportion,
                        'is_dominant': proportion > 0.5,
                    })

    return {
        'leiden_per_level': leiden_per_level,
        'level_positions': level_positions,
        'edges': edges,
        'node_sizes': node_sizes,
    }
