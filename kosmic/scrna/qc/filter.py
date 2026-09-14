# QC filtering for scRNA-seq data.
#
# Single-source-of-truth pipeline for the QC tab: detect MT genes,
# calculate QC metrics, filter cells by gene/count thresholds, filter
# genes by min-cells, drop high-mitochondrial cells. The GUI's
# '_apply_qc_fixed' is a thin wrapper around 'run_qc_pipeline'
# that adds Qt status updates and disk persistence.
from __future__ import annotations

from typing import Optional

import numpy as np


def detect_mitochondrial_genes(adata, species: Optional[str] = None):
    """Add 'adata.var['mt']' boolean column for mitochondrial genes.

    Parameters
    ----------
    adata : anndata.AnnData
    species : {'human', 'mouse', None}, optional
        - 'human' -> match prefix 'MT-'.
        - 'mouse' -> match prefix 'mt-'.
        - 'None' (default) -> auto-detect: try lowercase 'mt-'
          first; if no matches, try uppercase 'MT-'.

    Returns
    -------
    anndata.AnnData
        Modified in-place; returned for chaining.
    """
    if species == 'human':
        mask = adata.var_names.str.startswith('MT-')
    elif species == 'mouse':
        mask = adata.var_names.str.lower().str.startswith('mt-')
    else:
        # Auto-detect: try lowercase first (mouse convention), fall
        # back to uppercase (human/standard convention).
        mask = adata.var_names.str.lower().str.startswith('mt-')
        if mask.sum() == 0:
            mask = adata.var_names.str.startswith('MT-')

    # Also flag mitochondrial genes carrying legacy symbols (ND1, CYTB, CO1 ...)
    # so mito-% QC works even when gene names were never harmonised to MT-*.
    from kosmic.scrna.inspect.gene_names import _MITO_ALIAS_MAP
    legacy = {k.upper() for k in _MITO_ALIAS_MAP}
    legacy_mask = adata.var_names.str.upper().isin(legacy)
    adata.var['mt'] = mask | legacy_mask
    return adata


def fix_nan_in_sparse(adata):
    """Replace NaN / Inf entries in 'adata.X' with 0.

    Handles sparse and dense matrices alike. Returns the modified
    adata for chaining.
    """
    from scipy import sparse as sp

    if sp.issparse(adata.X):
        data = adata.X.data
        data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
        adata.X.data = data
    else:
        adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
    return adata


def run_qc_pipeline(adata, params: Optional[dict] = None):
    """Run the QC filter pipeline on a scRNA-seq AnnData.

    Mutates 'adata' in place. 'sc.pp.filter_cells' /
    'sc.pp.filter_genes' are already in-place; the 'max_counts' and
    'max_mt' slicing steps create a new view-then-copy, so the
    *returned* object may be a different instance than the input.
    Callers should always rebind: 'adata, stats = run_qc_pipeline(adata, ...)'
    and treat the original reference as invalid.

    Filter steps run in this order so the user's thresholds compose
    predictably:

      1. Annotate 'var['mt']' (auto-detect species).
      2. Replace NaN / Inf in 'X'.
      3. Compute QC metrics via 'sc.pp.calculate_qc_metrics' (skipped
         when 'n_genes_by_counts' / 'pct_counts_mt' are already in
         'obs' -- the GUI's histogram render path always populates
         them first, so the apply button doesn't pay for it twice).
      4. (optional) preserve original counts in 'adata.layers['counts']'.
      5. Filter cells by 'min_genes' / 'max_genes'.
      6. Filter cells by 'min_counts' / 'max_counts'.
      7. Filter genes by 'min_cells'.
      8. Recompute QC metrics on the filtered set.
      9. Filter cells by 'max_mt' (mitochondrial percent).

    Parameters
    ----------
    adata : anndata.AnnData
    params : dict, optional
        Threshold knobs (any missing key falls back to its default):

          'min_genes' (int, default 200) -- min genes/cell
          'max_genes' (int, default 0)   -- max genes/cell; 0 = disabled
          'min_counts' (int, default 0)  -- min total counts/cell
          'max_counts' (int, default 0)  -- max total counts/cell
          'min_cells' (int, default 3)   -- min cells/gene
          'max_mt' (float, default 0)    -- max %MT; 0 = disabled
          'preserve_counts_layer' (bool, default True)
                                           -- keep raw counts in
                                              'layers['counts']' so
                                              normalisation can be
                                              undone later
          'signature_panels' (dict, default None)
                                           -- {key: [genes]} ambient-
                                              contamination panels; each
                                              is scored into
                                              'obs['pct_counts_<key>']'
                                              on the pre-filter counts so
                                              the metric survives cell/gene
                                              filtering and later subsetting

    Returns
    -------
    tuple of (filtered_adata, stats_dict)
        stats_dict has before/after cell + gene counts, percentages,
        and the auto-detected species.
    """
    import scanpy as sc
    from kosmic.scrna.inspect.detection import detect_species

    p = dict(params or {})
    min_genes = p.get('min_genes', 200)
    max_genes = p.get('max_genes', 0)
    min_counts = p.get('min_counts', 0)
    max_counts = p.get('max_counts', 0)
    min_cells = p.get('min_cells', 3)
    max_mt = p.get('max_mt', 0)
    preserve_counts_layer = p.get('preserve_counts_layer', True)
    signature_panels = p.get('signature_panels', None)

    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    species = detect_species(list(adata.var_names))
    detect_mitochondrial_genes(adata, species)
    fix_nan_in_sparse(adata)

    # Score ambient-contamination panels on the pre-filter counts (all
    # genes still present, incl. dominant-cell-type genes) and stash the
    # result in obs. Stored per-cell, it rides through cell/gene filtering
    # and any later subset -- the percent.mito lifecycle.
    sig_matched = {}
    if signature_panels:
        from kosmic.scrna.qc.signatures import add_signature_metrics
        sig_matched = add_signature_metrics(adata, signature_panels, use_raw=False)

    metrics_already_present = (
        'n_genes_by_counts' in adata.obs.columns
        and 'pct_counts_mt' in adata.obs.columns
    )
    if not metrics_already_present:
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Preserve original counts before any per-cell filtering or later
    # normalisation. Idempotent: only writes the layer if not present.
    if preserve_counts_layer and 'counts' not in adata.layers:
        adata.layers['counts'] = adata.X.copy()

    if min_genes > 0:
        sc.pp.filter_cells(adata, min_genes=min_genes)
    if max_genes > 0:
        sc.pp.filter_cells(adata, max_genes=max_genes)
    if min_counts > 0:
        sc.pp.filter_cells(adata, min_counts=min_counts)
    if max_counts > 0:
        adata = adata[adata.obs['total_counts'] <= max_counts, :].copy()
    if min_cells > 0:
        sc.pp.filter_genes(adata, min_cells=min_cells)

    if 'mt' in adata.var.columns:
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    if max_mt > 0 and 'pct_counts_mt' in adata.obs.columns:
        adata = adata[adata.obs['pct_counts_mt'] < max_mt, :].copy()

    stats = {
        'n_cells_before': n_cells_before,
        'n_genes_before': n_genes_before,
        'n_cells_after': int(adata.n_obs),
        'n_genes_after': int(adata.n_vars),
        'cells_removed': int(n_cells_before - adata.n_obs),
        'genes_removed': int(n_genes_before - adata.n_vars),
        'pct_cells_removed': (n_cells_before - adata.n_obs)
                             / max(n_cells_before, 1) * 100,
        'pct_genes_removed': (n_genes_before - adata.n_vars)
                             / max(n_genes_before, 1) * 100,
        'species': species,
        'signature_matched': sig_matched,
    }

    return adata, stats
