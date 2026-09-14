# Ambient RNA correction (SoupX-style) — vectorised sparse subtraction.
#
# Originally wrapped 'soupx-python''s 'adjustCounts' end-to-end, but
# that path goes through 'expandClusters' which has a Python
# 'for gene_idx in range(n_genes)' loop nested inside a per-cluster
# loop -- empirically ~30 minutes on 150k cells x 30k genes. We replace
# it with a single vectorised sparse subtraction over 'X.data' that
# runs the same dataset in seconds.
#
# The math is the standard SoupX subtraction (R paper eq. 5):
#
#     soup_profile[g]   = sum over cells / total counts
#     expected_soup[c, g] = cell_total[c] * contamination * soup_profile[g]
#     corrected[c, g]   = max(X[c, g] - expected_soup[c, g], 0)
#
# We only operate on non-zero entries of 'X' (since
# 'max(0 - expected, 0) == 0' for zero entries) which keeps memory in
# the tens of MB even on the largest scRNA matrices.
#
# Cluster-aware redistribution (the part that makes 'adjustCounts'
# slow) is dropped: it gates how soup is *redistributed* between cells
# of the same type, not how much is subtracted in total. The
# contamination estimate dominates the result; the redistribution only
# adds noise when the soup profile itself is degraded (e.g. our
# filtered-only fallback). Users who care about cluster-aware behaviour
# can run R SoupX externally.
#
# Auto-estimation via 'soupx.autoEstCont' is still supported when both
# a raw droplet matrix is supplied and a quick Leiden has been run
# beforehand, but we no longer kick off an internal Leiden ourselves --
# that's the user's job. See 'run_soupx' for the contract.
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import pandas as pd
from scipy import sparse as sp


def _get_raw_counts(adata):
    """Return integer count matrix, preferring 'layers['counts']'.

    Raises if neither X nor the counts layer looks like raw integer data.
    """
    if 'counts' in adata.layers:
        X = adata.layers['counts']
    else:
        X = adata.X

    if hasattr(X, 'data'):
        max_val = float(X.data.max()) if X.data.size else 0
    else:
        max_val = float(X.max()) if X.size else 0
    if max_val < 20:
        raise ValueError(
            "SoupX needs raw integer counts but adata.X looks "
            "log-normalised (max < 20). Either run SoupX before "
            "Normalize, or store raw counts in adata.layers['counts'].")
    return X


def _build_soup_profile(X_filtered, raw_adata, filtered_var_names) -> np.ndarray:
    """Return the per-gene soup probability distribution (sums to 1).

    With 'raw_adata' (CellRanger 'raw_feature_bc_matrix') the soup
    is the raw droplet sum -- the canonical SoupX definition. Without
    raw, we fall back to summing over filtered cells (less accurate
    but still useful for bulk subtraction).
    """
    if raw_adata is None:
        gene_total = np.asarray(X_filtered.sum(axis=0)).ravel()
    else:
        raw_X = raw_adata.layers['counts'] if 'counts' in raw_adata.layers else raw_adata.X
        if not sp.issparse(raw_X):
            raw_X = sp.csr_matrix(raw_X)
        # Align raw to the filtered gene set
        if list(raw_adata.var_names) != list(filtered_var_names):
            common = filtered_var_names.intersection(raw_adata.var_names)
            if len(common) < len(filtered_var_names):
                raise ValueError(
                    f"Raw matrix is missing {len(filtered_var_names) - len(common)} "
                    "genes that are in the filtered adata.")
            gene_idx = [list(raw_adata.var_names).index(g) for g in filtered_var_names]
            raw_X = raw_X[:, gene_idx]
        gene_total = np.asarray(raw_X.sum(axis=0)).ravel()

    total = gene_total.sum()
    if total <= 0:
        raise ValueError("Soup profile is empty (zero total counts).")
    return gene_total / total


def _subtract_soup_vectorised(
    X_csr: 'sp.csr_matrix',
    soup_profile: np.ndarray,
    contamination: float,
    round_to_int: bool,
    rng: Optional[np.random.Generator] = None,
) -> 'sp.csr_matrix':
    """Apply SoupX subtraction in a single vectorised pass over X.data.

    Runs in O(nnz) numpy ops. No per-cell or per-gene Python loop.
    """
    n_cells, n_genes = X_csr.shape
    cell_total = np.asarray(X_csr.sum(axis=1)).ravel()
    cell_soup_budget = cell_total * contamination

    indptr = X_csr.indptr
    indices = X_csr.indices  # gene index per nonzero entry
    data = X_csr.data.astype(np.float64, copy=True)

    # Vectorised cell index per nonzero entry: repeat 0..n_cells-1
    # by the number of nonzeros in each row.
    cells_per_row = np.diff(indptr)
    cell_idx_per_entry = np.repeat(np.arange(n_cells), cells_per_row)

    expected = cell_soup_budget[cell_idx_per_entry] * soup_profile[indices]
    data_corrected = np.maximum(data - expected, 0.0)

    if round_to_int:
        rng = rng or np.random.default_rng(0)
        int_part = np.floor(data_corrected)
        frac_part = data_corrected - int_part
        round_up = rng.random(len(data_corrected)) < frac_part
        data_corrected = (int_part + round_up).astype(np.int64)

    return sp.csr_matrix(
        (data_corrected, indices.copy(), indptr.copy()),
        shape=X_csr.shape)


def run_soupx(
    adata,
    raw_adata=None,
    contamination_fraction: Optional[float] = None,
    round_to_int: bool = True,
    progress_callback=None,
    # method / tfidf_min / soup_quantile are no-ops in the vectorised
    # path; only auto-estimation (no fraction supplied + raw_adata
    # given) consults them via the R SoupX wrapper.
    method: str = 'subtraction',
    tfidf_min: float = 1.0,
    soup_quantile: float = 0.90,
) -> Tuple[object, dict]:
    """Run SoupX-style ambient-RNA subtraction on an AnnData.

    Parameters
    ----------
    adata : anndata.AnnData
        Filtered cells. Counts are read from 'layers['counts']' if
        present, else 'X'.
    raw_adata : anndata.AnnData, optional
        Raw droplet matrix (CellRanger 'raw_feature_bc_matrix',
        including empty droplets). When 'None' the soup profile is
        approximated from the filtered cells -- degraded but fast.
    contamination_fraction : float, optional
        Fixed contamination fraction (0..1). Required when
        'raw_adata' is None (auto-estimation has nothing to estimate
        against). Typical heart-tissue values: 0.05 -- 0.15.
    round_to_int : bool
        Round corrected counts back to integers via stochastic rounding
        (the same scheme R SoupX uses). Default True; needed for any
        downstream count-based DE engine.
    progress_callback : callable, optional
        Called with a status string at each pipeline step.
    method, tfidf_min, soup_quantile :
        Auto-estimation parameters (consulted only when
        'contamination_fraction' is None and 'raw_adata' is
        supplied). Ignored by the vectorised correction path.

    Returns
    -------
    tuple of (anndata.AnnData, dict)
        Modified adata (in-place) with corrected 'X' and
        'layers['counts_pre_soupx']' preserving the original counts.
        Info dict has 'contamination_fraction',
        'used_filtered_as_raw', 'auto_estimated', 'message'.
    """
    def _emit(msg: str):
        if progress_callback:
            progress_callback(msg)

    _emit("Reading raw counts...")
    X = _get_raw_counts(adata)
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.tocsr()

    used_filtered_as_raw = raw_adata is None
    if used_filtered_as_raw:
        _emit("No raw droplet matrix supplied; using filtered cells as soup approximation (degraded)")

    auto_estimated = contamination_fraction is None
    if auto_estimated:
        if raw_adata is None:
            raise ValueError(
                "Auto-estimation requires a raw droplet matrix. "
                "Supply raw_adata or pass a manual contamination_fraction "
                "(typical: 0.05 - 0.15).")
        _emit("Estimating contamination fraction (auto)...")
        contamination_fraction = _auto_estimate_via_soupx(
            adata, raw_adata, tfidf_min=tfidf_min, soup_quantile=soup_quantile)

    contamination_fraction = float(contamination_fraction)
    _emit(f"Building soup profile (contamination = {contamination_fraction:.1%})...")
    soup_profile = _build_soup_profile(X, raw_adata, adata.var_names)

    _emit(f"Subtracting soup from {X.nnz:,} non-zero entries...")
    corrected = _subtract_soup_vectorised(
        X, soup_profile, contamination_fraction, round_to_int)

    _emit("Writing corrected counts...")
    if 'counts_pre_soupx' not in adata.layers:
        adata.layers['counts_pre_soupx'] = adata.X.copy()
    adata.X = corrected
    if 'counts' in adata.layers:
        adata.layers['counts'] = corrected.copy()

    info = {
        'contamination_fraction': contamination_fraction,
        'method': 'subtraction',
        'used_filtered_as_raw': used_filtered_as_raw,
        'auto_estimated': auto_estimated,
        'message': (
            f"SoupX complete -- contamination {contamination_fraction:.1%} "
            f"removed{' (filtered-only fallback)' if used_filtered_as_raw else ''}"
            f"{'' if auto_estimated else ' (manual)'}"
        ),
    }
    _emit(info['message'])
    return adata, info


def _auto_estimate_via_soupx(adata, raw_adata, tfidf_min: float, soup_quantile: float) -> float:
    """Use 'soupx.autoEstCont' to estimate contamination from raw + filtered.

    Requires 'adata.obs['leiden']' (or any cluster column) -- the
    caller is responsible for clustering first. We don't run a quick
    Leiden internally any more.
    """
    import soupx

    if 'leiden' not in adata.obs.columns:
        raise RuntimeError(
            "Auto-estimation needs cluster labels in adata.obs['leiden']. "
            "Run the Cluster tab first, or supply a manual contamination_fraction.")

    X = _get_raw_counts(adata)
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    toc = X.T.tocsr().astype(np.float64)

    raw_X = raw_adata.layers['counts'] if 'counts' in raw_adata.layers else raw_adata.X
    if not sp.issparse(raw_X):
        raw_X = sp.csr_matrix(raw_X)
    if list(raw_adata.var_names) != list(adata.var_names):
        common = adata.var_names.intersection(raw_adata.var_names)
        if len(common) < adata.n_vars:
            raise ValueError(
                f"Raw matrix is missing {adata.n_vars - len(common)} genes.")
        gene_idx = [list(raw_adata.var_names).index(g) for g in adata.var_names]
        raw_X = raw_X[:, gene_idx]
    tod = raw_X.T.tocsr().astype(np.float64)

    sc_obj = soupx.SoupChannel(
        tod=tod, toc=toc,
        metaData=pd.DataFrame(index=adata.obs_names))
    sc_obj.setClusters(adata.obs['leiden'].astype(str).values)

    try:
        sc_obj = soupx.autoEstCont(
            sc_obj, tfidfMin=tfidf_min, soupQuantile=soup_quantile)
    except ValueError as exc:
        raise RuntimeError(
            f"SoupX auto-estimation failed: {exc}. "
            "Try lowering tfidf_min or supply a manual contamination_fraction."
        ) from exc

    contam = float(getattr(sc_obj, 'contFrac', np.nan))
    if not np.isfinite(contam):
        raise RuntimeError(
            "SoupX auto-estimation produced no contFrac estimate. "
            "Supply a manual contamination_fraction (typical: 0.05 - 0.15).")
    return contam
