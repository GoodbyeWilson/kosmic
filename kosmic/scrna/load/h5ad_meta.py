# Reading a study's metadata without loading its matrices.
#
# anndata's backed mode ('read_h5ad(path, backed="r")') keeps only X on
# disk: obsm, obsp, uns and every layer are read into memory. KOSMIC keeps
# the counts in layers['counts'], so a "cheap" backed read of a study
# loads the whole counts matrix -- 15 GB to fetch 189 MB of obs on the DCM
# atlas. Anything that only needs obs, the shape or the gene names goes
# through these readers instead; 'read_var_names' in gene_overlap.py is
# the same idea for var.
#
# Figure Export draws from a study's file when no workspace holds the
# study (ADR-014): 'read_figure_meta' reads everything but the matrices,
# 'read_genes' a few genes for every cell and 'read_cells' every gene for
# some cells. Matrices are read a chunk at a time, so the memory a read
# needs is set by 'chunk_bytes', not by the size of the study.
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd

#: Bytes of matrix (values plus indices) read at a time. On the DCM atlas
#: (1.17 million cells, 15 GB file) a 69-gene marker panel reads in 9.6 s
#: with a 1.0 GB peak at 256 MB, and 14 s with a 2.6 GB peak at 1 GB.
DEFAULT_CHUNK_BYTES = 256 << 20

#: Embeddings the figures draw. 'X_pca' and the graphs stay in the file.
FIGURE_EMBEDDINGS = ('X_umap', 'X_tsne')


def read_obs(h5ad_path) -> pd.DataFrame:
    """The obs frame of an h5ad, read on its own through h5py."""
    import h5py
    from anndata.io import read_elem

    with h5py.File(str(h5ad_path), 'r') as f:
        if 'obs' not in f:
            return pd.DataFrame()
        return read_elem(f['obs'])


def read_obs_columns(h5ad_path) -> list[str]:
    """The obs column names alone (no values read)."""
    import h5py

    with h5py.File(str(h5ad_path), 'r') as f:
        if 'obs' not in f:
            return []
        g = f['obs']
        idx = g.attrs.get('_index', '_index')
        if isinstance(idx, bytes):
            idx = idx.decode('utf-8')
        return [k for k in g.keys() if k != idx]


def read_shape(h5ad_path) -> Tuple[int, int]:
    """(n_obs, n_vars) from the file's X without reading it."""
    import h5py

    with h5py.File(str(h5ad_path), 'r') as f:
        X = f['X']
        if isinstance(X, h5py.Group):
            shape = X.attrs['shape']
            return int(shape[0]), int(shape[1])
        return int(X.shape[0]), int(X.shape[1])


def read_layer_names(h5ad_path) -> list[str]:
    """The names of the file's layers (none of them is read)."""
    import h5py

    with h5py.File(str(h5ad_path), 'r') as f:
        return list(f['layers'].keys()) if 'layers' in f else []


def write_obs(h5ad_path, obs: pd.DataFrame) -> None:
    """Replace the obs group of an h5ad in place; nothing else is touched.

    Object columns are written as categoricals (with missing values kept
    missing), the encoding anndata uses for strings, so the file reads
    back as a normal AnnData.
    """
    import h5py
    from anndata.io import write_elem

    obs = obs.copy()
    for col in obs.columns:
        if obs[col].dtype == object:
            obs[col] = pd.Categorical(obs[col].where(obs[col].notna(), None))
    with h5py.File(str(h5ad_path), 'r+') as f:
        write_elem(f, 'obs', obs)


def read_counts_adata(h5ad_path, *, counts_layer: Optional[str] = None,
                      uns: bool = True):
    """An AnnData holding only what pseudobulk analysis needs: obs, var and
    the counts as X.

    The counts come from ``layers['counts']`` when the file has it (a
    study after QC), otherwise from ``X`` (a file QC has not touched) --
    the same precedence as ``kosmic.scrna.counts.count_source``. A named
    ``counts_layer`` (``'decontX_counts'``) is read instead, and must be
    in the file. The normalised matrix, embeddings, graphs and other
    layers are never read, so a study costs half of what a full load does
    (Reichart-LV 16.7 GB -> 8 GB; the DCM atlas 27 GB -> 14 GB). ``uns`` is
    small and carries the role map, so it is read unless told otherwise.

    The result is a dataset in the "counts in X, no layers" state every
    KOSMIC reader already understands; ``uns['counts_loaded_from']`` says
    which counts X holds.
    """
    import anndata as ad
    import h5py
    from anndata.io import read_elem

    with h5py.File(str(h5ad_path), 'r') as f:
        obs = read_elem(f['obs'])
        var = read_elem(f['var'])
        layers = f['layers'] if 'layers' in f else {}
        if counts_layer is not None:
            if counts_layer not in layers:
                raise ValueError(
                    f"{Path(str(h5ad_path)).name} has no layer '{counts_layer}'.")
            X = read_elem(layers[counts_layer])
            source = f'layers:{counts_layer}'
        elif 'counts' in layers:
            X = read_elem(layers['counts'])
            source = 'layers:counts'
        else:
            X = read_elem(f['X'])
            source = 'X'
        extra = read_elem(f['uns']) if uns and 'uns' in f else {}
    a = ad.AnnData(X=X, obs=obs, var=var, uns=extra)
    a.uns['counts_loaded_from'] = source
    return a


# -- A study's file for Figure Export (ADR-014) ---------------------------------

def registered_h5ad(study_dir) -> Optional[Path]:
    """The study's working h5ad: the manifest's file if it still exists,
    else the newest h5ad in processed_data/. None if there is none.

    The manifest names the file relative to processed_data/, so a project
    folder copied to another computer resolves the same file.
    """
    from kosmic.manifest import study_overview
    from kosmic.paths import processed_data_dir

    study = Path(study_dir)
    processed = processed_data_dir(study)
    dataset = study_overview(study)["dataset"] or {}
    if dataset.get('file'):
        named = processed / dataset['file']
        if named.is_file():
            return named
    h5ads = sorted(processed.glob('*.h5ad'), key=lambda p: p.stat().st_mtime)
    return h5ads[-1] if h5ads else None


def read_figure_meta(h5ad_path):
    """Everything in a study's file except its matrices: obs, var, uns and
    the embeddings figures draw ('X_umap', 'X_tsne').

    ``X`` is None rather than an empty matrix, so code that reads
    expression from this object fails instead of drawing zeros. Expression
    comes from :func:`read_genes` and :func:`read_cells`. On the DCM atlas
    (1.17 million cells) this holds about 1.3 GB, most of it obs.
    """
    import anndata as ad
    import h5py
    from anndata.io import read_elem

    with h5py.File(str(h5ad_path), 'r') as f:
        obs = read_elem(f['obs'])
        var = read_elem(f['var'])
        uns = read_elem(f['uns']) if 'uns' in f else {}
        obsm = {}
        if 'obsm' in f:
            for key in FIGURE_EMBEDDINGS:
                if key in f['obsm']:
                    obsm[key] = np.asarray(read_elem(f['obsm'][key]))
    return ad.AnnData(X=None, obs=obs, var=var, uns=uns, obsm=obsm)


def read_genes(h5ad_path, genes: Sequence[str], *, source: str = 'log',
               counts_layer: Optional[str] = None,
               chunk_bytes: int = DEFAULT_CHUNK_BYTES, target_sum: float = 1e4):
    """A few genes for every cell of a study, read from its file.

    Parameters
    ----------
    h5ad_path : path
        The study's h5ad.
    genes : sequence of str
        Genes to read, in order; genes not in the file are left out.
    source : {'log', 'counts'}
        ``'log'``: log-normalised expression, the values
        :func:`kosmic.scrna.counts.log_normalised_genes` gives on the loaded
        study. ``X`` is read as it is when it is log-normalised; otherwise
        log1p(counts per ``target_sum``) of the count source, each cell's
        total taken over all genes in the same pass. ``uns['log1p']`` is
        set, so the result is read as it is. ``'counts'``: the counts, from
        ``counts_layer`` when given, else in the order of
        :func:`kosmic.scrna.counts.count_source`; ``uns['counts_loaded_from']``
        names where from.
    chunk_bytes : int
        Bytes of matrix read at a time.

    Returns
    -------
    AnnData
        All cells, the genes found, obs and the genes' var rows.
    """
    import anndata as ad
    import h5py
    import scipy.sparse as sp
    from anndata.io import read_elem

    if source not in ('log', 'counts'):
        raise ValueError(f"source must be 'log' or 'counts', not {source!r}")
    with h5py.File(str(h5ad_path), 'r') as f:
        obs = read_elem(f['obs'])
        if source == 'log' and _x_is_log(f, h5ad_path, chunk_bytes):
            where, normalise = 'X', False
        else:
            where = _count_where(f, counts_layer if source == 'counts' else None)
            normalise = source == 'log'
        node, var = _matrix(f, where)
        pos = {n: i for i, n in enumerate(map(str, var.index))}
        found = [str(g) for g in genes if str(g) in pos]
        cols = np.array([pos[g] for g in found], dtype=np.int64)
        sub, totals = _read_columns(node, cols, chunk_bytes, totals=normalise)
    if normalise:
        scale = np.where(totals > 0, target_sum / np.where(totals > 0, totals, 1), 0.0)
        if sp.issparse(sub):
            sub = sp.csr_matrix(sp.diags(scale) @ sub.astype(np.float32))
            sub.data = np.log1p(sub.data)
        else:
            sub = np.log1p(np.asarray(sub, dtype=np.float32) * scale[:, None])
    a = ad.AnnData(X=sub, obs=obs, var=var.iloc[cols].copy())
    if source == 'log':
        a.uns['log1p'] = {'base': None}
    else:
        a.uns['counts_loaded_from'] = 'X' if where == 'X' else where.replace(
            'layer:', 'layers:')
    return a


def read_cells(h5ad_path, rows, *, source: str = 'log',
               chunk_bytes: int = DEFAULT_CHUNK_BYTES, target_sum: float = 1e4):
    """Every gene for some cells of a study, read from its file.

    For figures that need the whole gene list, such as module scores whose
    control genes are drawn from all genes. ``rows`` are cell positions;
    they are read sorted and without repeats. ``source='log'`` gives the
    values :func:`kosmic.scrna.counts.log_normalised` gives on the loaded
    cells; ``'counts'`` the count source as it is.

    Returns
    -------
    AnnData
        The cells read (in ascending position), all genes of the source.
    """
    import anndata as ad
    import h5py
    import scipy.sparse as sp
    from anndata.io import read_elem

    if source not in ('log', 'counts'):
        raise ValueError(f"source must be 'log' or 'counts', not {source!r}")
    rows = np.unique(np.asarray(rows, dtype=np.int64))
    with h5py.File(str(h5ad_path), 'r') as f:
        obs = read_elem(f['obs']).iloc[rows]
        if source == 'log' and _x_is_log(f, h5ad_path, chunk_bytes):
            where, normalise = 'X', False
        else:
            where, normalise = _count_where(f, None), source == 'log'
        node, var = _matrix(f, where)
        X = _read_rows(node, rows, chunk_bytes)
    if normalise:
        X = sp.csr_matrix(X, dtype=np.float32) if sp.issparse(X) else np.asarray(
            X, dtype=np.float32).copy()
        totals = np.asarray(X.sum(axis=1)).ravel()
        with np.errstate(divide='ignore'):
            scale = np.where(totals > 0, target_sum / totals, 0.0)
        X = sp.diags(scale) @ X if sp.issparse(X) else X * scale[:, None]
        if sp.issparse(X):
            X = sp.csr_matrix(X)
            X.data = np.log1p(X.data)
        else:
            X = np.log1p(X)
    a = ad.AnnData(X=X, obs=obs, var=var.copy())
    if source == 'log':
        a.uns['log1p'] = {'base': None}
    return a


# Whether a file's X is log-normalised, by (path, mtime, size): deciding
# can mean reading all of X's values once.
_X_IS_LOG: dict = {}


def _x_is_log(f, h5ad_path, chunk_bytes: int) -> bool:
    """The test :func:`kosmic.scrna.counts.log_normalised_genes` applies to a
    loaded study -- ``uns['log1p']`` present, or every stored value below 20
    -- applied to the file. Values are read a chunk at a time and the scan
    stops at the first value of 20 or more, so a count matrix is decided
    from its first chunk."""
    if 'uns' in f and 'log1p' in f['uns']:
        return True
    st = Path(str(h5ad_path)).stat()
    key = (str(Path(str(h5ad_path)).resolve()), st.st_mtime_ns, st.st_size)
    if key in _X_IS_LOG:
        return _X_IS_LOG[key]
    import h5py
    node = f['X']
    values = node['data'] if isinstance(node, h5py.Group) else node
    is_log = True
    if values.size:
        step = max(1, chunk_bytes // max(values.dtype.itemsize, 1))
        if isinstance(node, h5py.Group):
            for i in range(0, values.shape[0], step):
                if float(np.max(values[i:i + step])) >= 20.0:
                    is_log = False
                    break
        else:
            n_rows = max(1, step // max(values.shape[1], 1))
            for i in range(0, values.shape[0], n_rows):
                block = values[i:i + n_rows]
                if block.size and float(np.max(block)) >= 20.0:
                    is_log = False
                    break
    _X_IS_LOG[key] = is_log
    return is_log


def _count_where(f, counts_layer: Optional[str]) -> str:
    """Where the counts are, in the order of ``count_source``: the named
    layer, ``layers['counts']``, ``raw``, then ``X``."""
    layers = f['layers'] if 'layers' in f else {}
    if counts_layer is not None:
        if counts_layer not in layers:
            raise ValueError(f"Counts layer '{counts_layer}' not found in the file.")
        return f'layer:{counts_layer}'
    if 'counts' in layers:
        return 'layer:counts'
    if 'raw' in f and 'X' in f['raw']:
        return 'raw'
    return 'X'


def _matrix(f, where: str):
    """(matrix node, var frame) for 'X', 'raw' or 'layer:<name>'."""
    from anndata.io import read_elem

    if where == 'raw':
        return f['raw']['X'], read_elem(f['raw']['var'])
    var = read_elem(f['var'])
    if where.startswith('layer:'):
        return f['layers'][where.split(':', 1)[1]], var
    return f['X'], var


def _shape_kind(node) -> tuple:
    import h5py

    if isinstance(node, h5py.Group):
        enc = node.attrs.get('encoding-type', 'csr_matrix')
        enc = enc.decode() if isinstance(enc, bytes) else str(enc)
        shape = tuple(int(s) for s in node.attrs['shape'])
        return shape, ('csc' if enc.startswith('csc') else 'csr')
    return tuple(node.shape), 'dense'


def _row_spans(indptr: np.ndarray, per_entry: int, chunk_bytes: int):
    """(first, stop) row ranges whose stored entries fit in 'chunk_bytes';
    a single row larger than that is a range of its own."""
    n = len(indptr) - 1
    budget = max(1, chunk_bytes // max(per_entry, 1))
    r0 = 0
    while r0 < n:
        r1 = int(np.searchsorted(indptr, indptr[r0] + budget, side='right')) - 1
        r1 = min(max(r1, r0 + 1), n)
        yield r0, r1
        r0 = r1


def _csr_block(node, indptr, r0: int, r1: int, n_cols: int):
    import scipy.sparse as sp

    a, b = int(indptr[r0]), int(indptr[r1])
    return sp.csr_matrix(
        (node['data'][a:b], node['indices'][a:b], indptr[r0:r1 + 1] - a),
        shape=(r1 - r0, n_cols))


def _read_columns(node, cols: np.ndarray, chunk_bytes: int, *, totals: bool):
    """(the columns 'cols' of every row, each row's total over all columns
    or None). Sparse in, CSR out; dense in, dense out."""
    import scipy.sparse as sp

    (n_rows, n_cols), kind = _shape_kind(node)
    row_totals = np.zeros(n_rows, dtype=float) if totals else None
    if kind == 'csr':
        indptr = node['indptr'][:].astype(np.int64)
        per = node['data'].dtype.itemsize + node['indices'].dtype.itemsize
        parts = []
        for r0, r1 in _row_spans(indptr, per, chunk_bytes):
            block = _csr_block(node, indptr, r0, r1, n_cols)
            parts.append(block[:, cols])
            if totals:
                row_totals[r0:r1] = np.asarray(block.sum(axis=1), dtype=float).ravel()
        sub = sp.vstack(parts, format='csr') if parts else sp.csr_matrix(
            (n_rows, len(cols)), dtype=node['data'].dtype)
        return sub, row_totals
    if kind == 'csc':
        indptr = node['indptr'][:].astype(np.int64)
        data, indices = node['data'], node['indices']
        pieces = [sp.csc_matrix(
            (data[indptr[c]:indptr[c + 1]], indices[indptr[c]:indptr[c + 1]],
             np.array([0, indptr[c + 1] - indptr[c]])), shape=(n_rows, 1))
            for c in cols]
        sub = (sp.hstack(pieces, format='csr') if pieces
               else sp.csr_matrix((n_rows, 0), dtype=data.dtype))
        if totals:
            step = max(1, chunk_bytes // (data.dtype.itemsize + indices.dtype.itemsize))
            for i in range(0, data.shape[0], step):
                row_totals += np.bincount(indices[i:i + step], weights=data[i:i + step],
                                          minlength=n_rows)
        return sub, row_totals
    step = max(1, chunk_bytes // max(n_cols * node.dtype.itemsize, 1))
    parts = []
    for r0 in range(0, n_rows, step):
        block = node[r0:r0 + step]
        parts.append(block[:, cols])
        if totals:
            row_totals[r0:r0 + len(block)] = np.asarray(block.sum(axis=1), dtype=float)
    sub = np.vstack(parts) if parts else np.zeros((n_rows, len(cols)), dtype=node.dtype)
    return sub, row_totals


def _read_rows(node, rows: np.ndarray, chunk_bytes: int):
    """The rows 'rows' (sorted, unique) of a matrix node, all columns."""
    import scipy.sparse as sp

    (n_rows, n_cols), kind = _shape_kind(node)
    if kind == 'csr':
        indptr = node['indptr'][:].astype(np.int64)
        per = node['data'].dtype.itemsize + node['indices'].dtype.itemsize
        parts = []
        for r0, r1 in _row_spans(indptr, per, chunk_bytes):
            want = rows[(rows >= r0) & (rows < r1)]
            if want.size:
                parts.append(_csr_block(node, indptr, r0, r1, n_cols)[want - r0])
        return sp.vstack(parts, format='csr') if parts else sp.csr_matrix(
            (0, n_cols), dtype=node['data'].dtype)
    if kind == 'csc':
        from anndata.io import read_elem
        return sp.csr_matrix(read_elem(node))[rows]
    step = max(1, chunk_bytes // max(n_cols * node.dtype.itemsize, 1))
    parts = []
    for r0 in range(0, n_rows, step):
        want = rows[(rows >= r0) & (rows < r0 + step)]
        if want.size:
            parts.append(node[r0:r0 + step][want - r0])
    return np.vstack(parts) if parts else np.zeros((0, n_cols), dtype=node.dtype)
