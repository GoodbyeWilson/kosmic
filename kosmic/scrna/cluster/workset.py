# The object the embedding and clustering steps actually work on.
#
# A study in memory holds two matrices of the same size -- the counts in
# layers['counts'] and the log-normalised X -- plus whatever excluded arm
# it carries. PCA, Harmony, the neighbour graph, Leiden, marker genes and
# UMAP need one matrix, for the included cells, and the obs/var/obsm
# slots. Deep-copying the whole study for them, as the workers once did
# (a full copy, then a second copy of the included rows), tripled the
# memory of a 650,000-nucleus study and put a million-cell atlas out of
# reach on a 64 GB machine.
#
# 'working_subset' builds an AnnData with exactly what a step needs:
# the included rows of one matrix (the matrix itself is shared, not
# copied, when every cell is included and it will only be read), fresh
# obs and var frames, and the requested obsm / uns / obsp entries.
# 'apply_results' writes the step's outputs back onto the full study,
# through 'scatter_results' when cells were excluded and directly
# otherwise, so the caller never touches two copies of the study.
from __future__ import annotations

from typing import Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from kosmic.scrna.inspect.roles import scatter_results


def working_subset(adata, mask: Optional[np.ndarray] = None, *,
                   copy_matrix: bool = False,
                   obsm_keys: Sequence[str] = (),
                   uns_keys: Sequence[str] = (),
                   obsp_keys: Sequence[str] = ()) -> ad.AnnData:
    """One-matrix AnnData of the included cells for an analysis step.

    Parameters
    ----------
    adata
        The full study. Never modified here.
    mask
        Boolean row mask of the cells to work on, or None for all of them.
    copy_matrix
        Give the subset its own matrix even when ``mask`` is None. Set it
        when the step will modify X in place (normalising counts); a
        row-masked matrix is always a new array.
    obsm_keys, uns_keys, obsp_keys
        Entries to carry across. obsm and obsp are row-subset; uns entries
        are shared as they are.

    Returns
    -------
    AnnData
        No layers, no raw. ``obs`` and ``var`` are copies, so columns the
        step adds (``leiden``, ``highly_variable``) do not appear on the
        full study until :func:`apply_results` puts them there.
    """
    if mask is None:
        X = adata.X.copy() if copy_matrix else adata.X
        obs = adata.obs.copy()
    else:
        idx = np.flatnonzero(mask)
        X = adata.X[idx]
        if sp.issparse(X) and not isinstance(X, sp.csr_matrix):
            X = sp.csr_matrix(X)
        obs = adata.obs.iloc[idx].copy()
        # As AnnData's own subsetting does: a sample that is excluded in
        # full must not linger as an empty category, or a per-batch HVG
        # selection meets an empty batch and fails.
        for col in obs.columns:
            if isinstance(obs[col].dtype, pd.CategoricalDtype):
                obs[col] = obs[col].cat.remove_unused_categories()
    work = ad.AnnData(X=X, obs=obs, var=adata.var.copy())
    for key in obsm_keys:
        if key in adata.obsm:
            arr = np.asarray(adata.obsm[key])
            work.obsm[key] = arr if mask is None else arr[np.flatnonzero(mask)]
    for key in uns_keys:
        if key in adata.uns:
            work.uns[key] = adata.uns[key]
    for key in obsp_keys:
        if key in adata.obsp:
            m = adata.obsp[key]
            if mask is None:
                work.obsp[key] = m
            else:
                idx = np.flatnonzero(mask)
                work.obsp[key] = sp.csr_matrix(m)[idx][:, idx]
    return work


def apply_results(full, work, mask: Optional[np.ndarray], *,
                  obsm_keys: Sequence[str] = (),
                  obs_cols: Sequence[str] = (),
                  var_cols: Sequence[str] = (),
                  uns_keys: Sequence[str] = (),
                  obsp_keys: Sequence[str] = (),
                  drop_uns: Sequence[str] = ()):
    """Write a step's outputs from ``work`` onto ``full``.

    With a mask this is :func:`scatter_results` (excluded cells get NaN
    rows in obsm and no category in obs; the graph is dropped from
    ``full`` because it is meaningless off the subset). Without one the
    entries are assigned directly, graph included. ``var_cols`` are
    copied by gene name either way, and ``drop_uns`` removes entries a
    re-run has made stale.
    """
    for col in var_cols:
        if col in work.var.columns:
            full.var[col] = work.var[col].reindex(full.var_names).values
    if mask is not None:
        scatter_results(full, work, mask, obsm_keys=obsm_keys,
                        obs_cols=obs_cols, uns_keys=uns_keys,
                        obsp_keys=obsp_keys)
    else:
        for key in obsm_keys:
            if key in work.obsm:
                full.obsm[key] = work.obsm[key]
        for col in obs_cols:
            if col in work.obs.columns:
                full.obs[col] = work.obs[col].values
        for key in uns_keys:
            if key in work.uns:
                full.uns[key] = work.uns[key]
        for key in obsp_keys:
            if key in work.obsp:
                full.obsp[key] = work.obsp[key]
    for key in drop_uns:
        full.uns.pop(key, None)
    return full


def apply_annotations(full, work) -> list[str]:
    """Copy annotation outputs from ``work`` back onto ``full``.

    An annotation step adds obs columns (labels, scores, per-type scores)
    and uns entries (annotation details, parameters) and touches nothing
    else, so every obs column of ``work`` that is new or differs from
    ``full``'s is written across, and every uns entry of ``work`` is set
    on ``full``. Returns the obs columns written. Both objects must have
    the same cells in the same order.
    """
    if len(work.obs) != len(full.obs) or not work.obs_names.equals(full.obs_names):
        raise ValueError("annotation results do not line up with the study's cells")
    written = []
    for col in work.obs.columns:
        if col not in full.obs.columns or not full.obs[col].equals(work.obs[col]):
            full.obs[col] = work.obs[col].values if not isinstance(
                work.obs[col].dtype, pd.CategoricalDtype) else work.obs[col].copy().values
            written.append(col)
    for key, value in work.uns.items():
        full.uns[key] = value
    return written
