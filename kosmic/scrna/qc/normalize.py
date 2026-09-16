# Library-size normalisation + log1p for scRNA-seq counts.
#
# Stores raw counts in 'adata.layers['counts']' before normalising so any
# downstream step (DE, MA) can recover the originals. That layer is the
# one copy of the counts a study keeps; '.raw' is no longer written
# (it duplicated the layer -- 8 GB of memory and disk on a large study).
from __future__ import annotations


def normalize_adata(adata, target_sum: int = 10000, log_transform: bool = True):
    """Normalise an AnnData object (library-size + optional log1p).

    Parameters
    ----------
    adata : anndata.AnnData
        Raw count data.
    target_sum : int
        Target sum for library-size normalisation.
    log_transform : bool
        Whether to apply log1p after normalisation.

    Returns
    -------
    anndata.AnnData
        The normalised AnnData (modified in-place and returned).
    """
    import scanpy as sc

    if 'counts' not in adata.layers:
        adata.layers['counts'] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=target_sum)
    if log_transform:
        sc.pp.log1p(adata)

    return adata
