# Harmony batch correction over PCA.
#
# Reads 'adata.obsm['X_pca']' and writes 'adata.obsm['X_pca_harmony']'.
# Downstream neighbors / Leiden / UMAP should auto-pick the harmony output
# via 'kosmic.scrna.cluster.neighbors.compute_neighbors'.
#
# The condition-aware option needs care, and the description it used to
# carry here was wrong. Passing a condition key runs harmony with
# 'vars_use=[batch, condition]' and 'theta=[2, 0]'. In harmonypy 'theta'
# weights only the diversity term of the *clustering* step; the correction
# step ('moe_correct_ridge') regresses out a fitted offset for every level
# of every variable in the design matrix, theta notwithstanding. So
# 'theta=0' does not mean "no correction along the condition axis" -- it
# means "do not also push the clusters to mix conditions". Condition-
# specific shifts are removed from the embedding either way.
#
# That is the right thing when the aim is to label cell types: you want
# one cardiomyocyte cluster, not a diseased one and a healthy one. It is
# the wrong thing when the aim is to find disease-specific cell states,
# because those are precisely what gets regressed away. It never affects
# differential expression, which runs on pseudobulk raw counts and never
# reads 'X_pca_harmony'.
#
# Leave the condition key blank for the standard single-variable run.
from __future__ import annotations

from typing import Optional, Sequence, Union


def run_harmony(
    adata,
    batch_key: str,
    condition_key: Optional[str] = None,
    theta: Optional[Union[float, Sequence[float]]] = None,
    max_iter: int = 20,
):
    """Run Harmony batch correction. Writes 'X_pca_harmony'.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have 'X_pca' in 'obsm' already.
    batch_key : str
        Column in 'adata.obs' identifying batch / sample to correct.
    condition_key : str, optional
        Second variable to include in the model (e.g. disease vs
        control). Its levels are corrected like any other -- 'theta=0'
        only stops the clustering step being pushed to mix them, it does
        not spare them from correction. Include it to merge conditions
        into shared cell-type clusters; leave it out to keep
        disease-specific states apart. Ignored if equal to 'batch_key'
        or absent from 'obs'.
    theta : float or sequence of float, optional
        Diversity weight per variable in the clustering step (not the
        correction step). When 'None' and a condition key is given,
        defaults to '[2, 0]'; otherwise scalar default.
    max_iter : int
        Harmony iterations cap.

    Returns
    -------
    anndata.AnnData
        Modified in place; returned for chaining.

    Raises
    ------
    ValueError
        If 'batch_key' is missing or 'X_pca' hasn't been computed.
    ImportError
        If 'harmonypy' is not installed.
    """
    if 'X_pca' not in adata.obsm:
        raise ValueError("X_pca not found in adata.obsm; run PCA first.")
    if batch_key not in adata.obs.columns:
        raise ValueError(
            f"Batch key {batch_key!r} not found in adata.obs.columns")

    import harmonypy as hm

    use_condition = (
        condition_key
        and condition_key != batch_key
        and condition_key in adata.obs.columns
    )
    if use_condition:
        vars_use = [batch_key, condition_key]
        if theta is None:
            theta = [2, 0]
    else:
        vars_use = [batch_key]

    ho = hm.run_harmony(
        adata.obsm['X_pca'],
        adata.obs,
        vars_use,
        theta=theta,
        max_iter_harmony=max_iter,
    )
    # Orient to (n_cells, n_pcs). Older harmonypy returns Z_corr as
    # (n_pcs, n_cells); 0.2.x already returns (n_cells, n_pcs) -- so key off the
    # cell axis rather than assuming a transpose.
    import numpy as np
    zc = np.asarray(ho.Z_corr)
    if zc.shape[0] != adata.n_obs:
        zc = zc.T
    adata.obsm['X_pca_harmony'] = zc
    return adata
