# Per-cluster differential expression for manual annotation.
#
# Wraps 'sc.tl.rank_genes_groups' with the defaults that make it useful
# for manual cluster naming:
#
# - 'method='wilcoxon'' -- robust for scRNA, the field-standard choice.
# - 'use_raw=True' when 'adata.raw' is available -- so HVG subsetting
#   during clustering doesn't hide non-HVG marker genes (e.g. classic
#   cell-type markers that didn't make the variable-gene cut).
# - 'pts=True' so we get 'pct_in' / 'pct_out' -- crucial for
#   judging marker specificity ("expressed in 95% of T cells, 3% of
#   others" beats "log2FC = 4 from one outlier").
#
# Two helpers:
#
# - 'compute_cluster_marker_genes' -- runs 'rank_genes_groups'
#   once and caches the result in 'adata.uns['rank_genes_groups']'.
# - 'get_top_marker_genes' -- pulls the top-N for a single
#   cluster out of the cached result as a clean DataFrame, ready to
#   render in a Qt table.
#
# A combined master atlas keeps raw counts in X and is never normalised on
# disk. Wilcoxon scores are unaffected by that (ranks within a gene), but
# Scanpy computes the log fold change by exponentiating X as if it were
# log1p data, so on counts it overflows. When X holds counts the rows that
# are ranked are copied and log1p(CP10K) normalised first; the caller's X
# is not modified.
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd


def compute_cluster_marker_genes(
    adata,
    cluster_col: str = 'leiden',
    method: str = 'wilcoxon',
    use_raw: Optional[bool] = None,
    max_cells_per_cluster: Optional[int] = None,
    seed: int = 0,
):
    """Compute per-cluster differential expression. Caches in 'uns'.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have 'cluster_col' in 'obs' (e.g. 'leiden').
    cluster_col : str
        Column to group by. Default 'leiden'.
    method : str
        Passed to 'sc.tl.rank_genes_groups'. Default 'wilcoxon'.
    use_raw : bool, optional
        Passed to scanpy. Default False: X is the full gene list and is
        log-normalised, which is what the Wilcoxon test expects; '.raw'
        in KOSMIC files held raw counts, so ranking on it confounded
        every gene with sequencing depth. When X holds raw counts (a
        combined master atlas), the ranked rows are log1p(CP10K)
        normalised in a copy first, so the log fold changes are valid;
        'adata.X' is not modified.
    max_cells_per_cluster : int, optional
        Rank on at most this many cells from each cluster, drawn at
        random with 'seed'. None (default) uses the configured
        MARKER_MAX_CELLS_PER_CLUSTER; 0 uses every cell. The Wilcoxon
        test ranks every cell for every gene, so its memory grows with
        the cell count -- on a 400,000-nucleus study it needs about 14 GB
        over the study itself, and a million-cell atlas does not fit.
        A few thousand cells per cluster rank the same markers; the
        result records 'max_cells_per_cluster' and 'n_cells_used' in
        its 'params'.

    Returns
    -------
    anndata.AnnData
        Same adata; modified in place. Result lives in
        'adata.uns['rank_genes_groups']'.

    Raises
    ------
    ValueError
        If 'cluster_col' is missing or has only one cluster.
    """
    import scanpy as sc

    if cluster_col not in adata.obs.columns:
        raise ValueError(
            f"Cluster column {cluster_col!r} not in adata.obs. "
            "Run the Cluster tab first.")

    n_clusters = adata.obs[cluster_col].nunique()
    if n_clusters < 2:
        raise ValueError(
            f"Need at least 2 clusters in {cluster_col!r}, found {n_clusters}.")

    # scanpy.tl.rank_genes_groups requires a categorical groupby column.
    # If the column is object/str dtype (common after manual edits or
    # certain h5ad round-trips), coerce in place -- otherwise scanpy
    # raises "can only use .cat accessor with a categorical".
    if not isinstance(adata.obs[cluster_col].dtype, pd.CategoricalDtype):
        adata.obs[cluster_col] = adata.obs[cluster_col].astype('category')

    if use_raw is None:
        use_raw = False
    if max_cells_per_cluster is None:
        from kosmic import MARKER_MAX_CELLS_PER_CLUSTER
        max_cells_per_cluster = MARKER_MAX_CELLS_PER_CLUSTER

    idx = _subsample_per_cluster(adata.obs[cluster_col], max_cells_per_cluster, seed)
    counts = not use_raw and _x_holds_counts(adata)
    if idx is None and counts:
        idx_all = np.arange(adata.n_obs)
        result = _rank_rows(adata, idx_all, cluster_col, method, normalise=True)
        adata.uns['rank_genes_groups'] = result
        return adata
    if idx is None:
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_col,
            method=method,
            use_raw=use_raw,
            pts=True,
        )
        return adata

    result = _rank_rows(adata, idx, cluster_col, method, normalise=counts)
    result['params']['max_cells_per_cluster'] = int(max_cells_per_cluster)
    result['params']['n_cells_used'] = int(len(idx))
    adata.uns['rank_genes_groups'] = result
    return adata


def _x_holds_counts(adata) -> bool:
    """True when X holds raw counts: no 'log1p' record in 'uns' and the
    stored values are integers."""
    from kosmic.scrna.cluster.hvg import _looks_like_counts
    if 'log1p' in getattr(adata, 'uns', {}):
        return False
    return _looks_like_counts(adata)


def _log1p_cp10k(X, target_sum: float = 1e4):
    """log1p(counts per 'target_sum') of a counts matrix, as a new float32
    matrix. Each cell's total is taken over all genes, as in
    'kosmic.scrna.counts.log_normalised'."""
    import scipy.sparse as sp

    X = sp.csr_matrix(X, dtype=np.float32) if sp.issparse(X) else np.array(X, dtype=np.float32)
    totals = np.asarray(X.sum(axis=1), dtype=float).ravel()
    scale = np.where(totals > 0, target_sum / np.where(totals > 0, totals, 1), 0.0)
    if sp.issparse(X):
        X = sp.csr_matrix(sp.diags(scale.astype(np.float32)) @ X)
        X.data = np.log1p(X.data)
        return X
    return np.log1p(X * scale[:, None].astype(np.float32))


def _rank_rows(adata, idx, cluster_col: str, method: str, normalise: bool) -> dict:
    """Run 'rank_genes_groups' on the rows 'idx' of 'adata' and return the
    result. The rows are copied into a small AnnData (log-normalised first
    when 'normalise'); nothing about 'adata' is touched."""
    import anndata as ad
    import scanpy as sc

    X = adata.X[idx]
    if normalise:
        X = _log1p_cp10k(X)
    sub = ad.AnnData(
        X=X,
        obs=adata.obs.iloc[idx][[cluster_col]].copy(),
        var=pd.DataFrame(index=adata.var_names),
    )
    if normalise:
        sub.uns['log1p'] = {'base': None}
    sub.obs[cluster_col] = sub.obs[cluster_col].cat.remove_unused_categories()
    sc.tl.rank_genes_groups(
        sub, groupby=cluster_col, method=method, use_raw=False, pts=True)
    result = sub.uns['rank_genes_groups']
    if normalise:
        result['params']['normalised'] = 'log1p(CP10K) of counts in X'
    return result


def _subsample_per_cluster(groups: pd.Series, cap: Optional[int], seed: int):
    """Row indices with at most 'cap' cells per cluster, or None when no
    cluster exceeds the cap (so the caller ranks on every cell)."""
    if not cap or cap <= 0:
        return None
    counts = groups.value_counts()
    if counts.max() <= cap:
        return None
    rng = np.random.default_rng(seed)
    codes = groups.astype(str).values
    keep = []
    for label, n in counts.items():
        rows = np.flatnonzero(codes == str(label))
        if n > cap:
            rows = rng.choice(rows, size=cap, replace=False)
        keep.append(rows)
    return np.sort(np.concatenate(keep))


def get_top_marker_genes(
    adata,
    cluster_id,
    n: int = 20,
) -> pd.DataFrame:
    """Pull the top-N marker genes for one cluster from cached results.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have 'adata.uns['rank_genes_groups']' populated by
        'compute_cluster_marker_genes'.
    cluster_id : str or int
        The cluster label (matches the keys in the cached result).
    n : int
        How many top genes to return. Default 20.

    Returns
    -------
    pandas.DataFrame
        Columns: 'gene', 'log2fc', 'pval', 'pval_adj',
        'pct_in', 'pct_out', 'score'. Sorted by score descending
        (= most up-regulated in this cluster vs others).
    """
    if 'rank_genes_groups' not in adata.uns:
        raise ValueError(
            "rank_genes_groups not cached. Run "
            "compute_cluster_marker_genes first.")

    result = adata.uns['rank_genes_groups']
    cluster_id = str(cluster_id)
    available_groups = list(result['names'].dtype.names)
    if cluster_id not in available_groups:
        raise KeyError(
            f"Cluster {cluster_id!r} not in rank_genes_groups results "
            f"(have: {available_groups[:5]}...).")

    names = list(result['names'][cluster_id][:n])
    df = pd.DataFrame({
        'gene': names,
        'score': result['scores'][cluster_id][:n],
        'log2fc': result['logfoldchanges'][cluster_id][:n],
        'pval': result['pvals'][cluster_id][:n],
        'pval_adj': result['pvals_adj'][cluster_id][:n],
    })

    # pts / pts_rest are stored as DataFrames keyed by gene name when
    # 'pts=True'. Look up percentages for our top-N.
    if 'pts' in result and 'pts_rest' in result:
        pts = result['pts']
        pts_rest = result['pts_rest']
        # pts/pts_rest are pandas DataFrames in modern scanpy
        if hasattr(pts, 'loc'):
            df['pct_in'] = pts.loc[names, cluster_id].values
            df['pct_out'] = pts_rest.loc[names, cluster_id].values
        else:
            # numpy structured array fallback
            df['pct_in'] = np.array([pts[g][cluster_id] for g in names])
            df['pct_out'] = np.array([pts_rest[g][cluster_id] for g in names])
    else:
        df['pct_in'] = np.nan
        df['pct_out'] = np.nan

    return df
