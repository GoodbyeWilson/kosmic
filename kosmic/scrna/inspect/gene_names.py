# Gene name harmonisation against the current official symbol table of the
# dataset's species: HGNC for human, MGI for mouse.
#
# Maps earlier symbols and aliases to current symbols, ensuring consistency
# across datasets from different genome annotation builds. The species is
# detected from the gene names (kosmic.scrna.inspect.detection): mouse
# symbols are mixed case. Mouse data harmonised against the human table
# would have mouse symbols renamed to unrelated human genes wherever the
# upper-case form happens to be an HGNC alias (Hc -> CLTC), so the two
# tables are never mixed.

import os
from dataclasses import dataclass, field
from importlib.resources import files
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

_DEFAULT_HGNC_PATH = str(files("kosmic.reference.hgnc") / "hgnc_symbols.tsv")
_DEFAULT_LOCUS_PATH = str(
    files("kosmic.reference.hgnc") / "hgnc_locus_groups.tsv.gz")
_DEFAULT_MGI_PATH = str(files("kosmic.reference.mgi") / "mgi_symbols.tsv.gz")
_MGI_VERSION_PATH = str(files("kosmic.reference.mgi") / "mgi_symbols.json")

#: Name of the symbol table used for each species, for display and records.
REFERENCE_NAMES = {'human': 'HGNC', 'mouse': 'MGI'}

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
    """Names the symbol table does not know, after harmonisation. Retired
    clone-based names and bare Ensembl IDs -- they match nothing in another
    study, so they are lost to any cross-study join."""
    species: str = 'human'
    """'human' or 'mouse': which symbol table was used."""
    reference: str = 'HGNC'
    """The table's name, from :data:`REFERENCE_NAMES`."""


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


def load_mgi_lookup(
    path: str = None,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Load the MGI mouse symbol table and build old-name -> symbol lookup.

    MGI lists each gene's synonyms, which include its earlier symbols,
    without saying which are which, so every mapping has the reason
    ``'synonym'``. As for HGNC, a name is mapped only when exactly one
    current gene claims it and it is not itself a current symbol; an
    ambiguous name is left unchanged rather than guessed.

    Parameters
    ----------
    path : str, optional
        Path to the MGI TSV. Defaults to the bundled
        kosmic/reference/mgi/mgi_symbols.tsv.gz.

    Returns
    -------
    tuple of (lookup, reasons)
        'lookup': old symbol -> current MGI symbol.
        'reasons': old symbol -> 'synonym'.
    """
    if path is None:
        path = _DEFAULT_MGI_PATH
    mgi = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    current = set(mgi['Symbol'])

    claimants = {}
    for symbol, synonyms in zip(mgi['Symbol'], mgi['Synonyms']):
        for old in synonyms.split('|'):
            old = old.strip()
            if not old or old == symbol or old in current:
                continue
            claimants.setdefault(old, set()).add(symbol)

    lookup = {old: next(iter(s)) for old, s in claimants.items() if len(s) == 1}
    reasons = dict.fromkeys(lookup, 'synonym')
    return lookup, reasons


def load_locus_groups(path: str = None, species: str = 'human') -> Dict[str, str]:
    """Load the current symbol -> locus group table for a species.

    Locus groups are HGNC's own coarse classification: ``'protein-coding
    gene'``, ``'non-coding RNA'``, ``'pseudogene'`` or ``'other'``. For
    mouse, MGI's feature types are reduced to the same four groups, so
    both species are counted the same way.

    Parameters
    ----------
    path : str, optional
        Path to the table. Defaults to the bundled HGNC locus-group table
        (human) or MGI symbol table (mouse).
    species : str
        'human' or 'mouse'.

    Returns
    -------
    dict
        Current symbol -> locus group.
    """
    if species == 'mouse':
        table = pd.read_csv(path or _DEFAULT_MGI_PATH, sep='\t', dtype=str,
                            keep_default_na=False)
        return dict(zip(table['Symbol'], table['Locus group']))
    table = pd.read_csv(path or _DEFAULT_LOCUS_PATH, sep='\t')
    return dict(zip(table['Approved symbol'], table['Locus group']))


def load_symbol_lookup(species: str = 'human'
                       ) -> Tuple[Dict[str, str], Dict[str, str]]:
    """The bundled old-name -> current-symbol lookup for a species:
    :func:`load_mgi_lookup` for ``'mouse'``, :func:`load_hgnc_lookup`
    otherwise. Returns ``(lookup, reasons)``."""
    return load_mgi_lookup() if species == 'mouse' else load_hgnc_lookup()


def reference_version(species: str = 'human') -> Dict[str, str]:
    """Which release of the bundled symbol table is in use, for the record.

    Mouse: the MGI report's download and build dates from
    ``mgi_symbols.json`` (written by ``dev/scripts/build_mgi_symbols.py``).
    Both species: a short SHA-1 of the table file, which identifies the
    exact table even where no date was recorded (the HGNC table) or the
    table has been updated in place.
    """
    import hashlib
    import json

    path = _DEFAULT_MGI_PATH if species == 'mouse' else _DEFAULT_HGNC_PATH
    info = {'table': REFERENCE_NAMES.get(species, 'HGNC')}
    if species == 'mouse':
        try:
            with open(_MGI_VERSION_PATH, encoding='utf-8') as f:
                meta = json.load(f)
            info.update({k: str(meta[k]) for k in ('source', 'downloaded', 'built')
                         if k in meta})
        except (OSError, ValueError):
            pass
    try:
        with open(path, 'rb') as f:
            info['sha1'] = hashlib.sha1(f.read()).hexdigest()[:12]
    except OSError:
        pass
    return info


def detect_gene_species(var_names) -> str:
    """'human' or 'mouse', from the case of the gene names.

    See :func:`kosmic.scrna.inspect.detection.detect_species`.
    """
    from kosmic.scrna.inspect.detection import detect_species
    return detect_species([str(n) for n in var_names])


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


def _remap_csr_columns_inplace(X, target, n_cols: int, chunk: int = 50_000_000) -> None:
    """Map every column index of CSR ``X`` through ``target`` and sum
    entries that land on the same (row, column). ``X`` ends up with
    ``n_cols`` columns. Works on the existing index and data arrays;
    the only transient allocation is scipy's trimmed copy after the sum.
    """
    idx = X.indices
    for start in range(0, idx.size, chunk):
        stop = min(start + chunk, idx.size)
        idx[start:stop] = target[idx[start:stop]]
    X._shape = (X.shape[0], int(n_cols))
    X.has_sorted_indices = False
    X.has_canonical_format = False
    X.sum_duplicates()


def harmonise_adata(adata, lookup: Dict[str, str] = None,
                    reasons: Dict[str, str] = None,
                    hgnc_path: str = None,
                    progress_callback=None,
                    species: str = None) -> Tuple:
    """Harmonise gene names in an AnnData object to current symbols.

    Parameters
    ----------
    adata : AnnData
        The dataset to harmonise. Modified in place.
    lookup : dict, optional
        Pre-built old→new mapping. If None, loaded for the species: the
        HGNC table for human, the MGI table for mouse.
    reasons : dict, optional
        Pre-built old→reason mapping.
    hgnc_path : str, optional
        Path to HGNC TSV (only used if lookup is None and the data is human).
    species : str, optional
        'human' or 'mouse'. Detected from the gene names when omitted.

    Returns
    -------
    (adata, HarmonisationReport)
        The modified AnnData and a report of what changed.
    """
    if species is None:
        species = detect_gene_species(adata.var_names[:100])
    if lookup is None:
        if species == 'mouse':
            lookup, reasons = load_mgi_lookup()
        else:
            lookup, reasons = load_hgnc_lookup(hgnc_path)
    if reasons is None:
        reasons = {}

    def _progress(msg):
        if progress_callback is not None:
            progress_callback(msg)

    report = HarmonisationReport(total_genes=adata.n_vars, species=species,
                                 reference=REFERENCE_NAMES.get(species, 'HGNC'))

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

    # Check for duplicates (multiple old names -> same approved symbol).
    # Merging is one sparse matmul with a column-merge matrix G (old genes
    # x kept genes, one 1 per row): every column of X lands in its kept
    # column and duplicates sum. The previous per-duplicate column
    # assignment rebuilt the whole sparse matrix once per merge, which on
    # a 220k x 34k, 630M-non-zero dataset took the better part of an hour.
    from collections import Counter
    name_counts = Counter(new_names)
    duplicates = {name: count for name, count in name_counts.items() if count > 1}

    if duplicates:
        import anndata as ad
        import scipy.sparse as sp

        first_position = {}
        kept_positions = []          # old column index of each kept column, in order
        target = np.empty(adata.n_vars, dtype=np.int64)
        for i, name in enumerate(new_names):
            if name not in first_position:
                first_position[name] = len(kept_positions)
                kept_positions.append(i)
            target[i] = first_position[name]
        merge_info = [(name, [var_names[i] for i, n in enumerate(new_names) if n == name])
                      for name in duplicates]

        n_kept = len(kept_positions)
        _progress(f"Merging {len(duplicates):,} duplicate symbols "
                  f"({adata.n_vars - n_kept:,} columns)...")

        def _merge(M):
            """Relabel columns to their kept index and sum what lands together.

            In place on a sparse matrix: the merge is only a change of
            column index (secondary -> primary) followed by summing the
            entries that now share a cell, which scipy does on the
            existing arrays. Nothing the size of the matrix is
            allocated -- on Reichart (881,081 cells, 1.35 billion
            non-zeros) a product against a merge matrix needed the
            original and the result live at once and exhausted 64 GB.
            """
            if M is None:
                return None
            if sp.issparse(M):
                M = M.tocsr() if not sp.isspmatrix_csr(M) else M
                _remap_csr_columns_inplace(M, target, n_kept)
                return M
            dense = np.asarray(M)
            out = np.zeros((dense.shape[0], n_kept), dtype=dense.dtype)
            np.add.at(out, (slice(None), target), dense)
            return out

        kept_mask = np.zeros(adata.n_vars, dtype=bool)
        kept_mask[kept_positions] = True
        merged = ad.AnnData(
            X=_merge(adata.X),
            obs=adata.obs,
            var=adata.var.iloc[kept_positions].copy(),
            obsm=dict(adata.obsm) if adata.obsm is not None else None,
            obsp=dict(adata.obsp) if adata.obsp is not None else None,
            uns=dict(adata.uns),
            layers={k: _merge(v) for k, v in adata.layers.items() if k is not None},
        )
        for k in list(adata.varm.keys()):
            merged.varm[k] = np.asarray(adata.varm[k])[kept_mask]
        if getattr(adata, 'raw', None) is not None:
            merged.raw = adata.raw
        adata = merged
        new_names = [new_names[i] for i in kept_positions]

        report.merged_genes = merge_info
        report.duplicates_merged = len(merge_info)
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
