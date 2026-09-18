# Reading a study's metadata without loading its matrices.
#
# anndata's backed mode ('read_h5ad(path, backed="r")') keeps only X on
# disk: obsm, obsp, uns and every layer are read into memory. KOSMIC keeps
# the counts in layers['counts'], so a "cheap" backed read of a study
# loads the whole counts matrix -- 15 GB to fetch 189 MB of obs on the DCM
# atlas. Anything that only needs obs, the shape or the gene names goes
# through these readers instead; 'read_var_names' in gene_overlap.py is
# the same idea for var.
from __future__ import annotations

from typing import Tuple

import pandas as pd


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


def read_counts_adata(h5ad_path, *, uns: bool = True):
    """An AnnData holding only what pseudobulk analysis needs: obs, var and
    the counts as X.

    The counts come from ``layers['counts']`` when the file has it (a
    study after QC), otherwise from ``X`` (a file QC has not touched) --
    the same precedence as ``kosmic.scrna.counts.count_source``. The
    normalised matrix, embeddings, graphs and other layers are never
    read, so a study costs half of what a full load does (Reichart-LV
    16.7 GB -> 8 GB; the DCM atlas 27 GB -> 14 GB). ``uns`` is small and
    carries the role map, so it is read unless told otherwise.

    The result is a dataset in the "counts in X, no layers" state every
    KOSMIC reader already understands.
    """
    import anndata as ad
    import h5py
    from anndata.io import read_elem

    with h5py.File(str(h5ad_path), 'r') as f:
        obs = read_elem(f['obs'])
        var = read_elem(f['var'])
        if 'layers' in f and 'counts' in f['layers']:
            X = read_elem(f['layers']['counts'])
            source = 'layers:counts'
        else:
            X = read_elem(f['X'])
            source = 'X'
        extra = read_elem(f['uns']) if uns and 'uns' in f else {}
    a = ad.AnnData(X=X, obs=obs, var=var, uns=extra)
    a.uns['counts_loaded_from'] = source
    return a
