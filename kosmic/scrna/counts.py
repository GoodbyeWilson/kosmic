# Where a study's raw counts live, resolved in one place.
#
# After QC the original counts are in ``layers['counts']`` and ``X`` is
# log-normalised; that is the data contract (CODEBASE.md, "Raw counts are
# required and preserved"). Earlier versions also snapshotted the counts
# into ``.raw`` at normalisation, and files from other tools (CellxGene)
# carry counts only in ``.raw``. Every reader that needs counts resolves
# them here, in this order:
#
#   1. ``layers[counts_layer]`` when a layer is named (DecontX output);
#   2. ``layers['counts']``;
#   3. ``.raw.X`` -- older KOSMIC files and foreign deposits;
#   4. ``X`` -- a file QC has not touched yet, where X still is the counts.
#
# ``X`` is never subset to highly variable genes (the cluster step keeps a
# mask), so ``adata.var_names`` is always the full gene list and the gene
# names for a count source are ``adata.var_names`` unless the source is
# ``.raw``, which carries its own.
from __future__ import annotations

from typing import Optional


def count_source(adata, counts_layer: Optional[str] = None):
    """The raw-count matrix and its gene names.

    Returns
    -------
    (matrix, var_names, where)
        ``matrix`` is cells x genes, sparse or dense, not copied;
        ``var_names`` a list of str; ``where`` one of
        ``'layers:<name>'``, ``'layers:counts'``, ``'raw'``, ``'X'``.

    Raises
    ------
    ValueError
        If ``counts_layer`` is named but absent.
    """
    if counts_layer is not None:
        if counts_layer not in adata.layers:
            raise ValueError(
                f"Counts layer '{counts_layer}' not found in adata.layers.")
        return adata.layers[counts_layer], list(map(str, adata.var_names)), f"layers:{counts_layer}"
    if 'counts' in adata.layers:
        return adata.layers['counts'], list(map(str, adata.var_names)), 'layers:counts'
    raw = getattr(adata, 'raw', None)
    if raw is not None:
        return raw.X, list(map(str, raw.var_names)), 'raw'
    return adata.X, list(map(str, adata.var_names)), 'X'


def count_var_names(adata) -> list:
    """Gene names of the count source (see :func:`count_source`)."""
    return count_source(adata)[1]


def has_counts_layer(adata) -> bool:
    """True when ``layers['counts']`` or ``.raw`` holds the original counts."""
    return 'counts' in adata.layers or getattr(adata, 'raw', None) is not None


def counts_adata(adata, genes=None, counts_layer: Optional[str] = None,
                 copy: bool = True):
    """An AnnData whose ``X`` is the raw counts, sharing ``obs``.

    Use for code that wants a plain AnnData of counts (scoring, pathway
    scoring, export) rather than reaching into layers itself. ``uns``,
    ``obsm`` and ``layers`` are not carried: this is a count view, not
    the study.

    ``copy=True`` (default) gives the caller its own matrix, so
    normalising or scaling it in place cannot touch the study's counts.
    Pass ``copy=False`` only for read-only use, where sharing the matrix
    saves a full copy on a large study. ``genes`` narrows to those
    present, always as a copy.
    """
    import anndata as ad
    import pandas as pd

    X, names, where = count_source(adata, counts_layer)
    var = adata.raw.var.copy() if where == 'raw' else adata.var.copy()
    if list(var.index) != names:
        var = pd.DataFrame(index=pd.Index(names))
    if genes is not None:
        keep = [g for g in genes if g in var.index]
        idx = [names.index(g) for g in keep]
        sub = X[:, idx]
        sub = sub.copy() if hasattr(sub, 'copy') else sub
        return ad.AnnData(X=sub, obs=adata.obs, var=var.loc[keep].copy())
    if copy:
        X = X.copy()
    return ad.AnnData(X=X, obs=adata.obs, var=var)


def looks_log_normalised(matrix) -> bool:
    """True when the stored values look like log1p-scaled expression
    (maximum below 20). log1p of a 10x count rarely exceeds ~12; raw
    counts run into the thousands."""
    try:
        data = matrix.data if hasattr(matrix, 'data') and hasattr(matrix, 'nnz') else matrix
        if getattr(data, 'size', 0) == 0:
            return True
        return float(data.max()) < 20.0
    except (AttributeError, TypeError, ValueError):
        return False


def log_normalised(adata, target_sum: float = 1e4):
    """Log-normalised expression on the full gene list, with its gene names.

    Returns ``(matrix, var_names)``. Uses ``X`` when it is already
    log-normalised (the state after the QC step); otherwise builds
    log1p(counts per ``target_sum``) from the count source without
    touching the study. Annotation and marker scoring read expression
    through this so they never score raw counts by mistake.
    """
    import numpy as np
    import scipy.sparse as sp

    if 'log1p' in getattr(adata, 'uns', {}) or looks_log_normalised(adata.X):
        return adata.X, list(map(str, adata.var_names))
    X, names, _ = count_source(adata)
    X = sp.csr_matrix(X, dtype=np.float32) if sp.issparse(X) else np.asarray(X, dtype=np.float32).copy()
    totals = np.asarray(X.sum(axis=1)).ravel()
    scale = np.where(totals > 0, target_sum / totals, 0.0)
    X = sp.diags(scale) @ X if sp.issparse(X) else X * scale[:, None]
    if sp.issparse(X):
        X = sp.csr_matrix(X)
        X.data = np.log1p(X.data)
    else:
        X = np.log1p(X)
    return X, names
