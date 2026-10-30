"""Stage 1 of ADR-014: reading what figures need from a study's file.

Each reader is checked, bit for bit, against what the same code gives on
the study loaded whole, in both forms a study's file takes
(log-normalised X with the counts in a layer; counts in X). The chunk
size is set far below the matrix size so every read crosses chunk
boundaries.
"""
from __future__ import annotations

import os

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp

from kosmic.scrna.counts import (
    counts_adata, log_normalised, log_normalised_genes, looks_log_normalised,
)
from kosmic.scrna.load import h5ad_meta
from kosmic.scrna.load.h5ad_meta import (
    read_cells, read_figure_meta, read_genes, registered_h5ad,
)
from tests.study_folder import MARKER_GENES, make_study_folder

SMALL_CHUNK = 1500          # bytes: a handful of rows per chunk


def _dense(m):
    return m.toarray() if sp.issparse(m) else np.asarray(m)


@pytest.fixture(params=['log', 'counts'])
def study(request, tmp_path):
    folder = make_study_folder(tmp_path, form=request.param)
    path = registered_h5ad(folder)
    return path, ad.read_h5ad(path)


def test_figure_meta_holds_obs_var_uns_and_embeddings_but_no_matrix(study):
    path, full = study
    meta = read_figure_meta(path)
    assert meta.X is None
    assert meta.shape == full.shape
    assert meta.obs.equals(full.obs)
    assert meta.var.equals(full.var)
    np.testing.assert_array_equal(meta.uns['pca']['variance_ratio'],
                                  full.uns['pca']['variance_ratio'])
    assert meta.uns['hvg'] == full.uns['hvg']
    np.testing.assert_array_equal(meta.obsm['X_umap'], full.obsm['X_umap'])
    assert 'X_pca' not in meta.obsm             # not drawn: left in the file
    assert len(meta.obsp) == 0 and len(meta.layers) == 0


def test_x_is_judged_log_normalised_as_the_loaded_study_is(study):
    path, full = study
    h5ad_meta._X_IS_LOG.clear()
    sf = h5ad_meta._StudyFile(path)
    assert h5ad_meta._x_is_log(sf, SMALL_CHUNK) == looks_log_normalised(full.X)


def test_log_genes_equal_the_loaded_studys(study):
    path, full = study
    genes = ['VWF', 'NOT_A_GENE', 'PTPRC', 'G3', 'VWF']
    want, found = log_normalised_genes(full, genes)
    got = read_genes(path, genes, source='log', chunk_bytes=SMALL_CHUNK)
    assert list(got.var_names) == found == ['VWF', 'PTPRC', 'G3', 'VWF']
    np.testing.assert_array_equal(_dense(got.X), _dense(want))
    assert 'log1p' in got.uns                   # read as it is downstream
    assert got.obs.equals(full.obs)
    # What the marker summary reads from it is the same values again.
    again, _ = log_normalised_genes(got, found)
    np.testing.assert_array_equal(_dense(again), _dense(got.X))


def test_count_genes_equal_the_loaded_studys(study):
    path, full = study
    genes = MARKER_GENES[:4] + ['NOT_A_GENE', 'G7']
    want = counts_adata(full, genes=genes)
    got = read_genes(path, genes, source='counts', chunk_bytes=SMALL_CHUNK)
    assert list(got.var_names) == list(want.var_names)
    np.testing.assert_array_equal(_dense(got.X), _dense(want.X))
    assert got.uns['counts_loaded_from'] in ('layers:counts', 'X')


def test_decontx_layer_is_read_when_named(tmp_path):
    path = registered_h5ad(make_study_folder(tmp_path, form='log'))
    full = ad.read_h5ad(path)
    got = read_genes(path, ['CD14', 'G1'], source='counts', counts_layer='decontX_counts',
                     chunk_bytes=SMALL_CHUNK)
    want = counts_adata(full, genes=['CD14', 'G1'], counts_layer='decontX_counts')
    np.testing.assert_array_equal(_dense(got.X), _dense(want.X))
    assert got.uns['counts_loaded_from'] == 'layers:decontX_counts'
    with pytest.raises(ValueError):
        read_genes(path, ['CD14'], source='counts', counts_layer='nope')


def test_cells_equal_the_loaded_studys(study):
    path, full = study
    rows = np.array([200, 3, 5, 3, 117, 0])
    got = read_cells(path, rows, source='log', chunk_bytes=SMALL_CHUNK)
    keep = np.unique(rows)
    want, names = log_normalised(full[keep])
    assert list(got.obs_names) == list(full.obs_names[keep])
    assert list(got.var_names) == names
    np.testing.assert_array_equal(_dense(got.X), _dense(want))


def test_a_dense_matrix_is_read_in_row_chunks(tmp_path):
    from tests.study_folder import make_adata
    a = make_adata('counts')
    a.X = a.X.toarray()
    path = tmp_path / 'dense.h5ad'
    a.write_h5ad(path)
    got = read_genes(path, ['LYZ', 'G2'], source='log', chunk_bytes=SMALL_CHUNK)
    want, _ = log_normalised_genes(a, ['LYZ', 'G2'])
    np.testing.assert_array_equal(_dense(got.X), _dense(want))
    cells = read_cells(path, [4, 9], source='counts', chunk_bytes=SMALL_CHUNK)
    np.testing.assert_array_equal(_dense(cells.X), a.X[[4, 9]])


def test_registered_file_is_the_manifests_else_the_newest(tmp_path):
    study = make_study_folder(tmp_path, form='counts')
    named = registered_h5ad(study)
    assert named.name == 'StudyA.h5ad'
    newer = named.with_name('backup.h5ad')
    newer.write_bytes(named.read_bytes())
    os.utime(newer, (named.stat().st_atime, named.stat().st_mtime + 60))
    assert registered_h5ad(study) == named          # the manifest wins
    named.unlink()
    assert registered_h5ad(study) == newer          # else the newest
    newer.unlink()
    assert registered_h5ad(study) is None
