# Streamed import of large 10x MatrixMarket deposits.
#
# Some deposits contain every barcode of the 10x whitelist rather than the
# cells Cell Ranger called: an unfiltered 'raw_feature_bc_matrix', or a
# 'cellranger aggr' run on raw matrices. A v3 whitelist is 6,794,880
# barcodes per sample (v2: 737,280), so GSE234052 (six samples aggregated)
# has 40.8 million barcode columns, of which about 22,000 are cells. Reading
# such a file with mmread and one obs row per barcode exhausts memory.
#
# This module reads the matrix in chunks instead. A first pass totals the
# UMIs and detected genes of every barcode. When most barcodes hold fewer
# than EMPTY_UMI counts the file is taken to be unfiltered, and only the
# barcodes that pass a UMI and gene threshold are kept; otherwise every
# barcode is kept. A second pass builds the matrix of the kept barcodes, so
# memory is bounded by the cells and not by the whitelist. The rule applied
# and the numbers kept and dropped are written to uns['kosmic_import'],
# which the import records in the study's provenance.
#
# Barcodes of an aggregated run carry the sample's position as a suffix
# ('-1' .. '-6'); split_by_barcode_suffix turns those into samples.
# See docs/ADR-013-unfiltered-10x-import.md.
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import pandas as pd

#: Matrices with more barcode columns than this are read in chunks.
LARGE_BARCODES = 200_000
#: A barcode below this many UMIs is an empty droplet for the unfiltered test.
EMPTY_UMI = 100
#: Default thresholds for a barcode to be kept from an unfiltered matrix.
MIN_UMI = 500
MIN_GENES = 200

_CHUNK = 5_000_000


def _open_text(path):
    path = Path(path)
    return gzip.open(path, 'rt') if path.suffix == '.gz' else open(path)


def mtx_header(path) -> tuple[int, int, int, int]:
    """(n_rows, n_cols, n_entries, header_lines) of a MatrixMarket file.

    Reads only the banner, comments and size line. 10x files have genes
    as rows and barcodes as columns.
    """
    with _open_text(path) as f:
        banner = f.readline()
        if not banner.startswith('%%MatrixMarket'):
            raise ValueError(f"{Path(path).name} is not a MatrixMarket file")
        if 'coordinate' not in banner or 'pattern' in banner:
            raise ValueError(f"{Path(path).name}: only coordinate matrices with "
                             "values can be read in chunks")
        n_lines = 1
        for line in f:
            n_lines += 1
            if not line.startswith('%'):
                n_rows, n_cols, n_entries = (int(x) for x in line.split()[:3])
                return n_rows, n_cols, n_entries, n_lines
    raise ValueError(f"{Path(path).name} has no size line")


def is_large(path) -> bool:
    """Whether the matrix has more barcode columns than LARGE_BARCODES."""
    return mtx_header(path)[1] > LARGE_BARCODES


def _entries(path, header_lines: int):
    """The (row, col, value) entries in chunks, as 0-based indices."""
    reader = pd.read_csv(
        path, sep=r'\s+', header=None, skiprows=header_lines,
        names=['r', 'c', 'v'], dtype={'r': np.int64, 'c': np.int64, 'v': np.float64},
        chunksize=_CHUNK)
    for chunk in reader:
        yield (chunk['r'].to_numpy() - 1, chunk['c'].to_numpy() - 1,
               chunk['v'].to_numpy())


def barcode_totals(path, progress: Optional[Callable[[str], None]] = None):
    """UMIs and detected genes of every barcode column, in one pass.

    Returns two int32 arrays of length n_cols.
    """
    _, n_cols, n_entries, header = mtx_header(path)
    umi = np.zeros(n_cols, dtype=np.int32)
    genes = np.zeros(n_cols, dtype=np.int32)
    done = 0
    for _, c, v in _entries(path, header):
        cols, inverse = np.unique(c, return_inverse=True)
        umi[cols] += np.rint(np.bincount(inverse, weights=v)).astype(np.int32)
        genes[cols] += np.bincount(inverse).astype(np.int32)
        done += len(c)
        if progress:
            progress(f"  Totalling counts per barcode: {done / n_entries:.0%}")
    return umi, genes


def _read_features(path) -> tuple[list, list]:
    """(gene names, gene ids) from a 10x features or genes file."""
    df = pd.read_csv(path, sep='\t', header=None, dtype=str)
    ids = df[0].tolist()
    if df.shape[1] >= 2:
        names = [n if isinstance(n, str) and n else i for n, i in zip(df[1], ids)]
    else:
        names = ids
    return names, ids


def _read_barcodes(path, keep_idx: np.ndarray) -> list:
    """The barcodes at the (sorted) indices 'keep_idx', read in chunks."""
    out, start = [], 0
    for chunk in pd.read_csv(path, sep='\t', header=None, usecols=[0], dtype=str,
                             chunksize=_CHUNK):
        stop = start + len(chunk)
        lo, hi = np.searchsorted(keep_idx, [start, stop])
        out.extend(chunk[0].to_numpy()[keep_idx[lo:hi] - start])
        start = stop
    return out


def load_large_10x(matrix, features, barcodes, *, min_umi: int = MIN_UMI,
                   min_genes: int = MIN_GENES,
                   progress: Optional[Callable[[str], None]] = None):
    """Read a 10x MatrixMarket set in chunks, keeping cells only if unfiltered.

    Parameters
    ----------
    matrix, features, barcodes : str or Path
        The '.mtx[.gz]', features/genes and barcodes files.
    min_umi, min_genes : int
        Thresholds a barcode must reach to be kept from an unfiltered matrix.
    progress : callable(str), optional

    Returns
    -------
    anndata.AnnData
        Cells x genes, float32 counts, with uns['kosmic_import'] describing
        the barcodes kept.
    """
    from anndata import AnnData
    from scipy.sparse import coo_matrix

    matrix = Path(matrix)
    n_genes, n_cols, _, header = mtx_header(matrix)
    say = progress or (lambda _msg: None)
    say(f"  {matrix.name}: {n_cols:,} barcodes x {n_genes:,} genes; "
        "reading in chunks")

    umi, genes = barcode_totals(matrix, progress)
    empty_fraction = float(np.mean(umi < EMPTY_UMI))
    unfiltered = empty_fraction > 0.5
    if unfiltered:
        keep = (umi >= min_umi) & (genes >= min_genes)
    else:
        keep = np.ones(n_cols, dtype=bool)
    keep_idx = np.flatnonzero(keep)
    if keep_idx.size == 0:
        raise RuntimeError(
            f"{matrix.name} holds {n_cols:,} barcodes and none has at least "
            f"{min_umi} UMIs and {min_genes} genes.")
    del genes
    if unfiltered:
        say(f"  {empty_fraction:.1%} of barcodes have fewer than {EMPTY_UMI} UMIs: "
            f"the file holds every barcode, not only cells. Keeping "
            f"{keep_idx.size:,} barcodes with >= {min_umi} UMIs and "
            f">= {min_genes} genes.")

    new_col = np.full(n_cols, -1, dtype=np.int64)
    new_col[keep_idx] = np.arange(keep_idx.size)
    rows, cells, vals = [], [], []
    for r, c, v in _entries(matrix, header):
        nc = new_col[c]
        m = nc >= 0
        rows.append(r[m].astype(np.int32))
        cells.append(nc[m])
        vals.append(v[m].astype(np.float32))
    del new_col
    X = coo_matrix((np.concatenate(vals), (np.concatenate(cells), np.concatenate(rows))),
                   shape=(keep_idx.size, n_genes)).tocsr()

    names, ids = _read_features(features)
    adata = AnnData(X=X)
    adata.obs_names = _read_barcodes(barcodes, keep_idx)
    adata.var_names = names
    adata.var['gene_ids'] = ids
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    record = {'source_file': matrix.name, 'barcodes_in_file': n_cols}
    if unfiltered:
        record.update({
            'cell_filter': (f"barcodes with >= {min_umi} UMIs and >= {min_genes} "
                            f"genes; the file holds every barcode ({empty_fraction:.1%} "
                            f"with < {EMPTY_UMI} UMIs)"),
            'barcodes_kept': int(keep_idx.size),
            'barcodes_dropped': int(n_cols - keep_idx.size),
        })
    adata.uns['kosmic_import'] = record
    return adata


def split_by_barcode_suffix(adata, name: str) -> None:
    """Set obs['sample'], one sample per barcode suffix when there are several.

    'cellranger aggr' appends each input sample's position to its barcodes
    ('-1', '-2', ...). When the barcodes carry more than one numeric suffix,
    each becomes a sample named '{name}_{N}'; otherwise every cell is
    sample 'name'. An existing obs['sample'] is kept.
    """
    if 'sample' in adata.obs.columns:
        return
    suffix = pd.Series(adata.obs_names, index=adata.obs_names).str.rsplit('-', n=1).str[-1]
    has_suffix = pd.Series(adata.obs_names).str.contains('-').all()
    if has_suffix and suffix.str.isdigit().all() and suffix.nunique() > 1:
        adata.obs['sample'] = (name + '_' + suffix).to_numpy()
        record = adata.uns.setdefault('kosmic_import', {})
        record['samples_from'] = 'barcode suffix (cellranger aggr)'
    else:
        adata.obs['sample'] = name
