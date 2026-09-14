# Highly variable gene selection.
#
# Wraps 'sc.pp.highly_variable_genes' with two non-obvious quality-of-life
# fixes:
#
# - Cap 'n_hvg' at 'adata.n_vars - 1' so the call never raises on
#   small datasets.
# - Detach 'adata.raw' for the duration of the call. Without this scanpy
#   creates a view + full copy of the raw matrix during HVG computation
#   and OOMs on large datasets. The original '.raw' is restored before
#   the function returns.
from __future__ import annotations

from typing import Optional


def _looks_like_counts(adata, n_sample: int = 2000) -> bool:
    """True when 'X' holds raw counts rather than log-normalised values.

    'seurat_v3' models the mean-variance relationship of counts and gives
    nonsense on log data, so the flavour has to follow the matrix. Uses
    the same integer test as the import-time matrix chooser: log1p of a
    UMI count rarely exceeds ~12, and never lands on integers.
    """
    import numpy as np
    import scipy.sparse as sp

    # Not 'hasattr(X, "data")': a dense ndarray has '.data' too, but it is
    # a memoryview of the buffer, not the stored values.
    X = adata.X
    data = X.data if sp.issparse(X) else np.asarray(X).ravel()
    if data.size == 0:
        return False
    sample = np.asarray(data[:n_sample])
    return bool(np.allclose(sample, np.round(sample)))


def find_hvg(adata, n_hvg: int = 2000, flavor: str = 'auto',
             subset: bool = False, batch_key: Optional[str] = None):
    """Annotate 'var['highly_variable']'; optionally subset to HVGs.

    Parameters
    ----------
    adata : anndata.AnnData
        Normalised + log-transformed data.
    n_hvg : int
        Target number of HVGs. Capped at 'adata.n_vars - 1'.
    flavor : str
        'auto' (default) picks 'seurat_v3' when 'X' still holds raw counts
        and scikit-misc is available, else 'seurat'. Anything else is
        passed through to 'sc.pp.highly_variable_genes' unchanged.

        The two rank genes differently. 'seurat' bins log-normalised means
        and takes the most dispersed gene per bin; 'seurat_v3' fits a loess
        mean-variance curve to the *counts* and ranks by standardised
        variance. The latter is what integration pipelines default to,
        because it is not thrown by the depth differences between studies
        -- which is exactly the situation a combined master is in. It
        needs counts, so it has to run before normalisation.
    subset : bool
        If True, return 'adata[:, var['highly_variable']].copy()' so
        downstream PCA / clustering operates on the HVG subset only.
        Default False (annotate-only): scanpy's PCA already auto-uses
        the HVG mask, and keeping the full gene set lets later marker /
        DE lookups find any gene.
    batch_key : str, optional
        Select HVGs within each batch and rank genes by how many batches
        call them variable, rather than pooling every cell first. Without
        it, on a combined master, genes that are variable only because
        the studies differ score highly -- the selection then encodes the
        batch effect it is meant to be blind to, and the study with the
        most cells or the deepest sequencing dominates. Ignored when
        absent from 'obs' or when it has only one level.

    Returns
    -------
    anndata.AnnData
        Modified in-place when 'subset=False'; a new view+copy when
        'subset=True'. Always return a value -- callers should rebind.
    """
    import scanpy as sc

    n_hvg = min(n_hvg, adata.n_vars - 1)

    raw_backup = adata.raw.to_adata() if adata.raw is not None else None
    if raw_backup is not None:
        del adata.raw

    # Only pass batch_key when it would actually do something: scanpy
    # errors on a single-level key, and a per-cell-unique column would
    # make every gene its own batch.
    if flavor == 'auto':
        flavor = 'seurat_v3' if _looks_like_counts(adata) else 'seurat'

    use_batch = None
    if batch_key and batch_key in adata.obs.columns:
        n_levels = int(adata.obs[batch_key].nunique())
        if 2 <= n_levels <= 100:
            use_batch = batch_key

    try:
        sc.pp.highly_variable_genes(
            adata, n_top_genes=n_hvg, flavor=flavor, batch_key=use_batch)
    except ImportError:
        # seurat_v3 delegates its loess fit to scikit-misc. Fall back
        # rather than fail the whole clustering run.
        sc.pp.highly_variable_genes(
            adata, n_top_genes=n_hvg, flavor='seurat', batch_key=use_batch)

    if raw_backup is not None:
        adata.raw = raw_backup

    if subset:
        adata = adata[:, adata.var['highly_variable']].copy()

    return adata
