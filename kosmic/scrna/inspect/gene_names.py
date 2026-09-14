# Gene name harmonisation using the HGNC approved symbol table.
#
# Maps previous symbols and aliases to current HGNC-approved symbols,
# ensuring consistency across datasets from different genome annotation builds.

import os
from dataclasses import dataclass, field
from importlib.resources import files
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

_DEFAULT_HGNC_PATH = str(files("kosmic.reference.hgnc") / "hgnc_symbols.tsv")
_DEFAULT_LOCUS_PATH = str(
    files("kosmic.reference.hgnc") / "hgnc_locus_groups.tsv.gz")

# The 13 mtDNA protein-coding genes are routinely labelled with legacy symbols
# (ND1, CYTB, ATP6, COI, ...) whose HGNC aliases collide with unrelated nuclear
# genes -- e.g. 'ND1' is an alias of both MT-ND1 and IVNS1ABP. In single-cell
# data these legacy names are always the mitochondrial gene, so map them
# explicitly to MT-*. This both fixes the label and hard-blocks the collision
# (an old 'ND1' can never be silently relabelled to a nuclear gene). The
# cyclooxygenase-colliding forms (COX1/COX2 -> PTGS1/PTGS2) are deliberately
# excluded and left to the ambiguity guard.
_MITO_ALIAS_MAP = {
    'ND1': 'MT-ND1', 'ND2': 'MT-ND2', 'ND3': 'MT-ND3', 'ND4': 'MT-ND4',
    'ND4L': 'MT-ND4L', 'ND5': 'MT-ND5', 'ND6': 'MT-ND6',
    'MTND1': 'MT-ND1', 'MTND2': 'MT-ND2', 'MTND3': 'MT-ND3', 'MTND4': 'MT-ND4',
    'MTND4L': 'MT-ND4L', 'MTND5': 'MT-ND5', 'MTND6': 'MT-ND6',
    'CYTB': 'MT-CYB', 'CYB': 'MT-CYB', 'MTCYB': 'MT-CYB',
    'ATP6': 'MT-ATP6', 'ATP8': 'MT-ATP8', 'MTATP6': 'MT-ATP6', 'MTATP8': 'MT-ATP8',
    'COI': 'MT-CO1', 'COII': 'MT-CO2', 'COIII': 'MT-CO3',
    'CO1': 'MT-CO1', 'CO2': 'MT-CO2', 'CO3': 'MT-CO3',
    'MTCO1': 'MT-CO1', 'MTCO2': 'MT-CO2', 'MTCO3': 'MT-CO3',
}


@dataclass
class HarmonisationReport:
    """Report of gene name harmonisation applied to a dataset."""
    total_genes: int = 0
    renamed: int = 0
    duplicates_merged: int = 0
    unchanged: int = 0
    mappings: List[Tuple[str, str, str]] = field(default_factory=list)
    """List of (old_name, new_name, reason) tuples."""
    merged_genes: List[Tuple[str, List[str]]] = field(default_factory=list)
    """List of (approved_symbol, [old_names_that_merged]) tuples."""
    locus_before: Dict[str, int] = field(default_factory=dict)
    """HGNC locus-group counts of the names as they stand, from
    :func:`count_locus_groups`. Empty if not computed."""
    locus_after: Dict[str, int] = field(default_factory=dict)
    """The same counts as they would be after harmonisation."""
    unrecognised: List[str] = field(default_factory=list)
    """Names HGNC does not know, after harmonisation. Retired clone-based
    names and bare Ensembl IDs -- they match nothing in another study, so
    they are lost to any cross-study join."""


def load_hgnc_lookup(
    path: str = None,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Load HGNC symbol table and build old-name → approved-symbol lookup.

    Parameters
    ----------
    path : str, optional
        Path to HGNC TSV file. Defaults to bundled kosmic/reference/hgnc/hgnc_symbols.tsv.

    Returns
    -------
    tuple of (lookup, reasons)
        'lookup': old symbol -> current approved symbol.
        'reasons': old symbol -> 'previous_symbol' or 'alias'.
    """
    if path is None:
        path = _DEFAULT_HGNC_PATH

    hgnc = pd.read_csv(path, sep='\t')

    # Build set of all approved symbols for conflict detection
    approved_set = set(hgnc['Approved symbol'].dropna())

    lookup = {}
    reasons = {}

    # Process previous symbols first (they take priority), then aliases. Within
    # each level, only map an old name when a SINGLE approved symbol claims it --
    # an old name claimed by two genes (e.g. 'ND1', which is an alias of both
    # IVNS1ABP and MT-ND1) is ambiguous, so remapping it would silently corrupt
    # one of them. Ambiguous names are left unchanged rather than guessed.
    for col, reason in [('Previous symbols', 'previous_symbol'),
                        ('Alias symbols', 'alias')]:
        # Explode comma-separated values into individual rows
        mask = hgnc[col].notna()
        if not mask.any():
            continue
        subset = hgnc.loc[mask, ['Approved symbol', col]].copy()
        # Split comma-separated strings into lists, then explode
        subset['_old'] = subset[col].str.split(', ')
        exploded = subset.explode('_old')
        exploded['_old'] = exploded['_old'].str.strip()

        # Collect all approved symbols claiming each old name at this level.
        claimants = {}
        for old, approved in zip(exploded['_old'], exploded['Approved symbol']):
            if not old or old == approved or old in approved_set:
                continue
            claimants.setdefault(old, set()).add(approved)

        for old, apprs in claimants.items():
            if old in lookup:
                continue  # already mapped by a higher-priority column
            if len(apprs) == 1:
                lookup[old] = next(iter(apprs))
                reasons[old] = reason
            # len > 1: ambiguous -> do not remap (avoids ND1->IVNS1ABP-style
            # corruption); the name stays as-is in the data.

    # Mitochondrial legacy names win outright over any ambiguous alias handling.
    for old, mt in _MITO_ALIAS_MAP.items():
        if mt in approved_set and old not in approved_set:
            lookup[old] = mt
            reasons[old] = 'mitochondrial'

    return lookup, reasons


def load_locus_groups(path: str = None) -> Dict[str, str]:
    """Load the HGNC approved symbol -> locus group table.

    Locus groups are HGNC's own coarse classification: ``'protein-coding
    gene'``, ``'non-coding RNA'``, ``'pseudogene'`` or ``'other'``.

    Parameters
    ----------
    path : str, optional
        Path to the locus-group TSV. Defaults to the bundled
        kosmic/reference/hgnc/hgnc_locus_groups.tsv.gz.

    Returns
    -------
    dict
        Approved symbol -> locus group.
    """
    if path is None:
        path = _DEFAULT_LOCUS_PATH
    table = pd.read_csv(path, sep='	')
    return dict(zip(table['Approved symbol'], table['Locus group']))


def count_locus_groups(var_names, locus_groups: Dict[str, str] = None,
                       lookup: Dict[str, str] = None) -> Dict[str, int]:
    """Count how many of ``var_names`` fall into each HGNC locus group.

    A gene counts towards a group only if its name is an HGNC approved
    symbol; everything else lands in ``'unrecognised'``. That bucket is
    the interesting one -- it holds retired clone-based names
    (``RP11-34P13.3``), bare Ensembl IDs left on the index where the
    depositor had no symbol, and genuine typos, all of which fail to
    match the same gene in another study.

    Parameters
    ----------
    var_names : iterable of str
        Gene names, e.g. ``adata.var_names``.
    locus_groups : dict, optional
        From :func:`load_locus_groups`. Loaded if omitted.
    lookup : dict, optional
        Old-symbol -> approved-symbol map from :func:`load_hgnc_lookup`.
        When given, names are counted as they would be *after*
        harmonisation, so the caller can show what harmonising recovers.

    Returns
    -------
    dict
        Keys ``'protein-coding gene'``, ``'non-coding RNA'``,
        ``'pseudogene'``, ``'other'``, ``'unrecognised'`` and ``'total'``.
        Counts are of distinct names, since harmonisation merges
        duplicates.
    """
    if locus_groups is None:
        locus_groups = load_locus_groups()

    names = {lookup.get(n, n) for n in var_names} if lookup else set(var_names)
    counts = {'protein-coding gene': 0, 'non-coding RNA': 0,
              'pseudogene': 0, 'other': 0, 'unrecognised': 0}
    for name in names:
        counts[locus_groups.get(name, 'unrecognised')] += 1
    counts['total'] = len(names)
    return counts


def harmonise_adata(adata, lookup: Dict[str, str] = None,
                    reasons: Dict[str, str] = None,
                    hgnc_path: str = None,
                    progress_callback=None) -> Tuple:
    """Harmonise gene names in an AnnData object to current HGNC symbols.

    Parameters
    ----------
    adata : AnnData
        The dataset to harmonise. Modified in place.
    lookup : dict, optional
        Pre-built old→new mapping. If None, loaded from HGNC table.
    reasons : dict, optional
        Pre-built old→reason mapping.
    hgnc_path : str, optional
        Path to HGNC TSV (only used if lookup is None).

    Returns
    -------
    (adata, HarmonisationReport)
        The modified AnnData and a report of what changed.
    """
    if lookup is None:
        lookup, reasons = load_hgnc_lookup(hgnc_path)
    if reasons is None:
        reasons = {}

    def _progress(msg):
        if progress_callback is not None:
            progress_callback(msg)

    report = HarmonisationReport(total_genes=adata.n_vars)

    _progress(f"Mapping {adata.n_vars:,} gene names...")
    var_names = list(adata.var_names)
    new_names = []
    mapping_list = []

    for name in var_names:
        if name in lookup:
            new_name = lookup[name]
            reason = reasons.get(name, 'alias')
            mapping_list.append((name, new_name, reason))
            new_names.append(new_name)
        else:
            new_names.append(name)

    report.renamed = len(mapping_list)
    report.mappings = mapping_list
    _progress(f"Found {report.renamed:,} genes to rename")

    # Check for duplicates (multiple old names → same approved symbol)
    from collections import Counter
    name_counts = Counter(new_names)
    duplicates = {name: count for name, count in name_counts.items() if count > 1}

    if duplicates:
        import scipy.sparse as sp

        name_to_positions = {}
        for i, name in enumerate(new_names):
            name_to_positions.setdefault(name, []).append(i)

        drop_indices = set()
        merge_info = []
        n_to_merge = sum(1 for p in name_to_positions.values() if len(p) > 1)

        # Convert to CSC for fast column operations if sparse
        X_was_csr = False
        if sp.issparse(adata.X) and sp.isspmatrix_csr(adata.X):
            _progress("Converting matrix for column operations...")
            adata.X = adata.X.tocsc()
            X_was_csr = True

        merged_count = 0
        for name, positions in name_to_positions.items():
            if len(positions) <= 1:
                continue
            merged_count += 1
            if merged_count % 10 == 0 or merged_count == n_to_merge:
                _progress(f"Merging duplicates {merged_count}/{n_to_merge}")
            primary = positions[0]
            merge_info.append((name, [var_names[p] for p in positions]))
            for secondary in positions[1:]:
                if sp.issparse(adata.X):
                    adata.X[:, primary] = adata.X[:, primary] + adata.X[:, secondary]
                else:
                    adata.X[:, primary] += adata.X[:, secondary]
                drop_indices.add(secondary)

        # Convert back to CSR
        if X_was_csr:
            _progress("Converting matrix back...")
            adata.X = adata.X.tocsr()

        report.merged_genes = merge_info
        report.duplicates_merged = len(merge_info)

        # Drop secondary columns
        _progress("Dropping duplicate columns...")
        keep_mask = np.ones(adata.n_vars, dtype=bool)
        keep_mask[list(drop_indices)] = False
        adata = adata[:, keep_mask]
        new_names = [n for i, n in enumerate(new_names) if keep_mask[i]]
    else:
        report.merged_genes = []
        report.duplicates_merged = 0

    # Apply the new names
    adata.var_names = pd.Index(new_names)
    adata.var_names_make_unique()

    # Sync any companion gene-symbol columns so converters.py doesn't
    # revert the index on reload (e.g. cellxgene 'feature_name').
    for col in ('feature_name', 'gene_symbol', 'gene_name', 'gene_symbols', 'Symbol', 'symbol'):
        if col in adata.var.columns:
            adata.var[col] = list(adata.var_names)

    # Also sync in adata.raw if present (CellTypist reads from there).
    # AnnData's raw.var is sometimes read-only depending on the backing
    # store; tolerate that, but let other failures propagate.
    if getattr(adata, 'raw', None) is not None:
        try:
            raw_var = adata.raw.var
            for col in ('feature_name', 'gene_symbol', 'gene_name', 'gene_symbols', 'Symbol', 'symbol'):
                if col in raw_var.columns:
                    raw_names = list(raw_var.index)
                    new_raw = [lookup.get(n, n) for n in raw_names]
                    raw_var[col] = new_raw
        except (AttributeError, ValueError):
            pass

    report.unchanged = report.total_genes - report.renamed
    return adata, report


def download_hgnc_table(output_path: str = None) -> str:
    """Download the latest HGNC symbol table.

    Parameters
    ----------
    output_path : str, optional
        Where to save. Defaults to kosmic/reference/hgnc/hgnc_symbols.tsv.

    Returns
    -------
    str : path to the saved file.
    """
    import urllib.request

    if output_path is None:
        output_path = _DEFAULT_HGNC_PATH

    url = (
        'https://www.genenames.org/cgi-bin/download/custom'
        '?col=gd_app_sym&col=gd_prev_sym&col=gd_aliases'
        '&status=Approved&hgnc_datea=&hgnc_dateb='
        '&order_by=gd_app_sym_sort&format=text&submit=submit'
    )

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    urllib.request.urlretrieve(url, output_path)
    return output_path
