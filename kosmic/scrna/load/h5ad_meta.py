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

import os
import threading
from contextlib import contextmanager
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

    with file_lock(h5ad_path).reading(), h5py.File(str(h5ad_path), 'r') as f:
        if 'obs' not in f:
            return pd.DataFrame()
        return read_elem(f['obs'])


def read_obs_columns(h5ad_path) -> list[str]:
    """The obs column names alone (no values read)."""
    import h5py

    with file_lock(h5ad_path).reading(), h5py.File(str(h5ad_path), 'r') as f:
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

    with file_lock(h5ad_path).reading(), h5py.File(str(h5ad_path), 'r') as f:
        X = f['X']
        if isinstance(X, h5py.Group):
            shape = X.attrs['shape']
            return int(shape[0]), int(shape[1])
        return int(X.shape[0]), int(X.shape[1])


def read_layer_names(h5ad_path) -> list[str]:
    """The names of the file's layers (none of them is read)."""
    import h5py

    with file_lock(h5ad_path).reading(), h5py.File(str(h5ad_path), 'r') as f:
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
    with file_lock(h5ad_path).writing(), h5py.File(str(h5ad_path), 'r+') as f:
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

    with file_lock(h5ad_path).reading(), h5py.File(str(h5ad_path), 'r') as f:
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


# -- Reading and saving one file from two threads ---------------------------------
#
# On Windows a file open for reading can neither be overwritten nor
# replaced, and a reader that opens a file while a save writes it reads a
# half-written file. Figure Export reads a study's file in a worker thread
# while a step can be saving it, so both go through the per-file lock
# below: a save waits for the chunk being read (about a second), and a
# reader that waited for a save stops with 'StudyFileChanged' rather than
# read a mix of two versions. The lock is for this process only.

class FileLock:
    """Readers share the file; a save has it alone. A waiting save goes
    before readers that arrive after it, so a stream of chunk reads cannot
    hold a save off. Not reentrant."""

    def __init__(self):
        self._cond = threading.Condition()
        self._readers = 0
        self._writing = False
        self._writers_waiting = 0

    @contextmanager
    def reading(self):
        with self._cond:
            self._cond.wait_for(lambda: not self._writing and not self._writers_waiting)
            self._readers += 1
        try:
            yield
        finally:
            with self._cond:
                self._readers -= 1
                self._cond.notify_all()

    @contextmanager
    def writing(self):
        with self._cond:
            self._writers_waiting += 1
            self._cond.wait_for(lambda: not self._writing and self._readers == 0)
            self._writers_waiting -= 1
            self._writing = True
        try:
            yield
        finally:
            with self._cond:
                self._writing = False
                self._cond.notify_all()


_LOCKS: dict = {}
_LOCKS_GUARD = threading.Lock()


class StudyFileChanged(RuntimeError):
    """The file was saved again while it was being read."""


def file_lock(path) -> FileLock:
    """The lock shared by everything in this process that reads or saves
    the file at 'path'."""
    key = os.path.normcase(os.path.abspath(str(path)))
    with _LOCKS_GUARD:
        return _LOCKS.setdefault(key, FileLock())


def save_h5ad(adata, path, *, atomic: bool = False, **write_kwargs) -> None:
    """Write 'adata' to 'path' while no reader has the file open.

    ``atomic=True`` writes a temporary file beside it and renames it into
    place, so a failed write leaves the previous file intact; this needs
    room for both files at once. ``atomic=False`` writes the file in place,
    as 'AnnData.write_h5ad' does. Extra keyword arguments (``compression``)
    go to 'write_h5ad'.
    """
    target = Path(str(path))
    with file_lock(target).writing():
        if not atomic:
            adata.write_h5ad(str(target), **write_kwargs)
            return
        tmp = target.with_name(target.name + '.tmp')
        adata.write_h5ad(str(tmp), **write_kwargs)
        os.replace(tmp, target)


class _StudyFile:
    """A study's file, opened afresh for each read under its lock.

    The file's modification time and size are taken when the first read
    starts; a later open that finds them changed raises 'StudyFileChanged'.
    """

    def __init__(self, path):
        self.path = Path(str(path))
        self.lock = file_lock(self.path)
        st = self.path.stat()
        self.signature = (st.st_mtime_ns, st.st_size)

    @contextmanager
    def open(self):
        import h5py

        with self.lock.reading():
            st = self.path.stat()
            if (st.st_mtime_ns, st.st_size) != self.signature:
                raise StudyFileChanged(f"{self.path.name} was saved while it was being read.")
            with h5py.File(str(self.path), 'r') as f:
                yield f


def read_figure_meta(h5ad_path):
    """Everything in a study's file except its matrices: obs, var, uns and
    the embeddings figures draw ('X_umap', 'X_tsne').

    ``X`` is None rather than an empty matrix, so code that reads
    expression from this object fails instead of drawing zeros. Expression
    comes from :func:`read_genes` and :func:`read_cells`. On the DCM atlas
    (1.17 million cells) this holds about 1.3 GB, most of it obs.
    """
    import anndata as ad
    from anndata.io import read_elem

    with _StudyFile(h5ad_path).open() as f:
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
        Bytes of matrix read at a time; the file is opened once per chunk.

    Returns
    -------
    AnnData
        All cells, the genes found, obs and the genes' var rows.

    Raises
    ------
    StudyFileChanged
        The file was saved again during the read.
    """
    import anndata as ad
    import scipy.sparse as sp
    from anndata.io import read_elem

    if source not in ('log', 'counts'):
        raise ValueError(f"source must be 'log' or 'counts', not {source!r}")
    sf = _StudyFile(h5ad_path)
    with sf.open() as f:
        obs = read_elem(f['obs'])
    if source == 'log' and _x_is_log(sf, chunk_bytes):
        where, normalise = 'X', False
    else:
        with sf.open() as f:
            where = _count_where(f, counts_layer if source == 'counts' else None)
        normalise = source == 'log'
    m = _Matrix(sf, where)
    pos = {n: i for i, n in enumerate(map(str, m.var.index))}
    found = [str(g) for g in genes if str(g) in pos]
    cols = np.array([pos[g] for g in found], dtype=np.int64)
    sub, totals = m.columns(cols, chunk_bytes, totals=normalise)
    if normalise:
        scale = np.where(totals > 0, target_sum / np.where(totals > 0, totals, 1), 0.0)
        if sp.issparse(sub):
            sub = sp.csr_matrix(sp.diags(scale) @ sub.astype(np.float32))
            sub.data = np.log1p(sub.data)
        else:
            sub = np.log1p(np.asarray(sub, dtype=np.float32) * scale[:, None])
    a = ad.AnnData(X=sub, obs=obs, var=m.var.iloc[cols].copy())
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

    Raises
    ------
    StudyFileChanged
        The file was saved again during the read.
    """
    import anndata as ad
    import scipy.sparse as sp
    from anndata.io import read_elem

    if source not in ('log', 'counts'):
        raise ValueError(f"source must be 'log' or 'counts', not {source!r}")
    rows = np.unique(np.asarray(rows, dtype=np.int64))
    sf = _StudyFile(h5ad_path)
    with sf.open() as f:
        obs = read_elem(f['obs']).iloc[rows]
    if source == 'log' and _x_is_log(sf, chunk_bytes):
        where, normalise = 'X', False
    else:
        with sf.open() as f:
            where = _count_where(f, None)
        normalise = source == 'log'
    m = _Matrix(sf, where)
    X = m.rows(rows, chunk_bytes)
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
    a = ad.AnnData(X=X, obs=obs, var=m.var.copy())
    if source == 'log':
        a.uns['log1p'] = {'base': None}
    return a


# Whether a file's X is log-normalised, by (path, mtime, size): deciding
# can mean reading all of X's values once.
_X_IS_LOG: dict = {}


def _x_is_log(sf: _StudyFile, chunk_bytes: int) -> bool:
    """The test :func:`kosmic.scrna.counts.log_normalised_genes` applies to a
    loaded study -- ``uns['log1p']`` present, or every stored value below 20
    -- applied to the file. Values are read a chunk at a time and the scan
    stops at the first value of 20 or more, so a count matrix is decided
    from its first chunk."""
    import h5py

    key = (os.path.normcase(str(sf.path.resolve())),) + sf.signature
    if key in _X_IS_LOG:
        return _X_IS_LOG[key]
    with sf.open() as f:
        if 'uns' in f and 'log1p' in f['uns']:
            _X_IS_LOG[key] = True
            return True
        node = f['X']
        sparse = isinstance(node, h5py.Group)
        values = node['data'] if sparse else node
        n, itemsize = values.shape[0], values.dtype.itemsize
        width = 1 if sparse else max(values.shape[1], 1)
    step = max(1, chunk_bytes // (itemsize * width))
    is_log = True
    for i in range(0, n, step):
        with sf.open() as f:
            block = f['X']['data'][i:i + step] if sparse else f['X'][i:i + step]
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


def _node(f, where: str):
    if where == 'raw':
        return f['raw']['X']
    if where.startswith('layer:'):
        return f['layers'][where.split(':', 1)[1]]
    return f['X']


class _Matrix:
    """One matrix of a study's file ('X', 'raw' or 'layer:<name>'), read a
    chunk at a time with the file opened afresh for each chunk."""

    def __init__(self, sf: _StudyFile, where: str):
        import h5py
        from anndata.io import read_elem

        self.sf, self.where = sf, where
        with sf.open() as f:
            node = _node(f, where)
            self.var = read_elem(f['raw']['var'] if where == 'raw' else f['var'])
            if isinstance(node, h5py.Group):
                enc = node.attrs.get('encoding-type', 'csr_matrix')
                enc = enc.decode() if isinstance(enc, bytes) else str(enc)
                self.kind = 'csc' if enc.startswith('csc') else 'csr'
                self.shape = tuple(int(s) for s in node.attrs['shape'])
                self.indptr = node['indptr'][:].astype(np.int64)
                self.dtype = node['data'].dtype
                self.per_entry = node['data'].dtype.itemsize + node['indices'].dtype.itemsize
            else:
                self.kind, self.shape, self.dtype = 'dense', tuple(node.shape), node.dtype

    def _spans(self, chunk_bytes: int):
        """(first, stop) ranges along the stored axis -- rows of a CSR or
        dense matrix, entries of a CSC one -- that fit in 'chunk_bytes'."""
        if self.kind == 'dense':
            step = max(1, chunk_bytes // max(self.shape[1] * self.dtype.itemsize, 1))
            return [(r, min(r + step, self.shape[0])) for r in range(0, self.shape[0], step)]
        if self.kind == 'csc':
            n, step = int(self.indptr[-1]), max(1, chunk_bytes // self.per_entry)
            return [(i, min(i + step, n)) for i in range(0, n, step)]
        return list(_row_spans(self.indptr, self.per_entry, chunk_bytes))

    def _csr_block(self, r0: int, r1: int):
        import scipy.sparse as sp

        a, b = int(self.indptr[r0]), int(self.indptr[r1])
        with self.sf.open() as f:
            node = _node(f, self.where)
            data, indices = node['data'][a:b], node['indices'][a:b]
        return sp.csr_matrix((data, indices, self.indptr[r0:r1 + 1] - a),
                             shape=(r1 - r0, self.shape[1]))

    def _dense_block(self, r0: int, r1: int):
        with self.sf.open() as f:
            return _node(f, self.where)[r0:r1]

    def columns(self, cols: np.ndarray, chunk_bytes: int, *, totals: bool):
        """(the columns 'cols' of every row, each row's total over all
        columns or None). Sparse in, CSR out; dense in, dense out."""
        import scipy.sparse as sp

        n_rows = self.shape[0]
        row_totals = np.zeros(n_rows, dtype=float) if totals else None
        if self.kind == 'csr':
            parts = []
            for r0, r1 in self._spans(chunk_bytes):
                block = self._csr_block(r0, r1)
                parts.append(block[:, cols])
                if totals:
                    row_totals[r0:r1] = np.asarray(block.sum(axis=1), dtype=float).ravel()
            sub = sp.vstack(parts, format='csr') if parts else sp.csr_matrix(
                (n_rows, len(cols)), dtype=self.dtype)
            return sub, row_totals
        if self.kind == 'csc':
            pieces = []
            with self.sf.open() as f:
                node = _node(f, self.where)
                for c in cols:
                    a, b = int(self.indptr[c]), int(self.indptr[c + 1])
                    pieces.append(sp.csc_matrix(
                        (node['data'][a:b], node['indices'][a:b], np.array([0, b - a])),
                        shape=(n_rows, 1)))
            sub = (sp.hstack(pieces, format='csr') if pieces
                   else sp.csr_matrix((n_rows, 0), dtype=self.dtype))
            if totals:
                for i, j in self._spans(chunk_bytes):
                    with self.sf.open() as f:
                        node = _node(f, self.where)
                        row_totals += np.bincount(node['indices'][i:j],
                                                  weights=node['data'][i:j],
                                                  minlength=n_rows)
            return sub, row_totals
        parts = []
        for r0, r1 in self._spans(chunk_bytes):
            block = self._dense_block(r0, r1)
            parts.append(block[:, cols])
            if totals:
                row_totals[r0:r1] = np.asarray(block.sum(axis=1), dtype=float)
        sub = np.vstack(parts) if parts else np.zeros((n_rows, len(cols)), dtype=self.dtype)
        return sub, row_totals

    def rows(self, rows: np.ndarray, chunk_bytes: int):
        """The rows 'rows' (sorted, unique), all columns. A CSC matrix --
        KOSMIC writes none -- is read whole."""
        import scipy.sparse as sp

        n_cols = self.shape[1]
        if self.kind == 'csc':
            from anndata.io import read_elem
            with self.sf.open() as f:
                return sp.csr_matrix(read_elem(_node(f, self.where)))[rows]
        parts = []
        for r0, r1 in self._spans(chunk_bytes):
            want = rows[(rows >= r0) & (rows < r1)]
            if not want.size:
                continue
            if self.kind == 'csr':
                parts.append(self._csr_block(r0, r1)[want - r0])
            else:
                parts.append(self._dense_block(r0, r1)[want - r0])
        if self.kind == 'csr':
            return sp.vstack(parts, format='csr') if parts else sp.csr_matrix(
                (0, n_cols), dtype=self.dtype)
        return np.vstack(parts) if parts else np.zeros((0, n_cols), dtype=self.dtype)


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
