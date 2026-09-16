# Ambient-contamination signature scoring for scRNA/snRNA-seq.
#
# A `percent.mito`-style per-cell metric: the fraction of a cell's counts
# coming from a panel of genes specific to a dominant, RNA-rich cell type
# (e.g. cardiomyocytes). Scored *inside a different* cell type it estimates
# how much ambient RNA from the dominant type has leaked in -- the quantity
# you regress against to tell real DE from soup-driven DE.
#
# Pure logic (no Qt). The QC pipeline calls `compute_signature_pct` to
# stash `obs['pct_counts_<key>']`; because it lives in `obs` it survives
# subsetting (e.g. to endothelial cells) and flows through to DE, exactly
# like `pct_counts_mt`.
from __future__ import annotations


import numpy as np


def _format_panel(adata, genes):
    """Format panel gene symbols to the dataset's species convention."""
    from kosmic.scrna.inspect.detection import (
        detect_species, format_gene_for_species)

    from kosmic.scrna.counts import count_var_names
    species = detect_species(count_var_names(adata))
    return [format_gene_for_species(g, species) for g in genes]


def _counts_source(adata, use_raw=True):
    """Pick the best raw-count source and its var_names.

    The counts (layers['counts'], else .raw, else X) when `use_raw`,
    otherwise `adata.X`. Returns `(matrix, var_names_list)`.
    """
    from kosmic.scrna.counts import count_source
    if use_raw:
        X, names, _ = count_source(adata)
        return X, names
    return adata.X, list(map(str, adata.var_names))


def flag_signature_genes(adata, genes, key):
    """Add a boolean `adata.var[key]` flagging panel genes present in `X`.

    Mirrors `detect_mitochondrial_genes`; usable as a `qc_vars` entry for
    `sc.pp.calculate_qc_metrics`. Returns the list of matched gene names.
    """
    formatted = set(_format_panel(adata, genes))
    matched = [g for g in adata.var_names if g in formatted]
    adata.var[key] = adata.var_names.isin(matched)
    return matched


def compute_signature_pct(adata, genes, key=None, use_raw=True, store=True):
    """Per-cell percent of counts from a contamination panel.

    ``pct_i = 100 * sum(counts in panel genes) / total counts`` for cell
    *i*, computed on raw counts so it is independent of normalisation and
    HVG subsetting. This is the on-demand path that works on any object,
    reading the panel genes from the counts.

    Parameters
    ----------
    adata : anndata.AnnData
    genes : list of str
        Panel gene symbols (species-formatted internally).
    key : str, optional
        If given and `store`, result is written to
        ``obs['pct_counts_<key>']``.
    use_raw : bool
        Read the counts rather than X.
    store : bool
        Write the metric into `adata.obs`.

    Returns
    -------
    numpy.ndarray
        Per-cell contamination percent, shape ``(n_obs,)``.
    """
    formatted = set(_format_panel(adata, genes))
    X, var_names = _counts_source(adata, use_raw=use_raw)

    idx = [i for i, g in enumerate(var_names) if g in formatted]

    total = np.asarray(X.sum(axis=1)).ravel().astype(float)
    total = np.clip(total, 1.0, None)

    if idx:
        panel_sum = np.asarray(X[:, idx].sum(axis=1)).ravel().astype(float)
    else:
        panel_sum = np.zeros(X.shape[0], dtype=float)

    pct = panel_sum / total * 100.0

    if store and key:
        adata.obs[f'pct_counts_{key}'] = pct
    return pct


def add_signature_metrics(adata, panels, use_raw=True):
    """Compute `pct_counts_<key>` for each panel in `panels`.

    Parameters
    ----------
    adata : anndata.AnnData
    panels : dict
        ``{key: [gene_list]}``.
    use_raw : bool
        Passed to `compute_signature_pct`.

    Returns
    -------
    dict
        ``{key: n_matched_genes}`` -- how many panel genes were found in
        the count source (0 means the metric is all-zero / uninformative).
    """
    matched_counts = {}
    for key, genes in (panels or {}).items():
        formatted = set(_format_panel(adata, genes))
        _, var_names = _counts_source(adata, use_raw=use_raw)
        matched_counts[key] = sum(1 for g in var_names if g in formatted)
        compute_signature_pct(adata, genes, key=key, use_raw=use_raw, store=True)
    return matched_counts
