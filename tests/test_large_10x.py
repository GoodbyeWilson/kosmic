"""Import of 10x matrices that hold every whitelist barcode.

An unfiltered or 'cellranger aggr' deposit (GSE234052: 40.8 million
barcodes, about 22,000 cells) is read in chunks, keeps only the barcodes
that look like cells, and is split into samples by barcode suffix. The
synthetic matrices here stand in for the whitelist with a lowered
LARGE_BARCODES threshold.
"""
from __future__ import annotations

import gzip
import shutil

import numpy as np
import pytest
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

from kosmic.scrna.load import large_10x
from kosmic.scrna.load.converters import load_10x_folder, load_prefixed_10x

N_GENES = 300
N_BARCODES = 50          # per sample
CELLS = [3, 4, 5]        # per sample: barcode positions that are cells


def _write_10x(folder, *, empty_counts=True, prefix=''):
    """Two aggregated samples ('-1', '-2'); cells have 300 genes x 2 UMIs."""
    rng = np.random.default_rng(0)
    X = np.zeros((N_GENES, 2 * N_BARCODES), dtype=np.int64)
    for s in range(2):
        for b in range(N_BARCODES):
            col = s * N_BARCODES + b
            if b in CELLS:
                X[:, col] = 2
            elif empty_counts and b % 2:
                X[rng.integers(0, N_GENES, 3), col] = 1   # ambient: 3 UMIs
            elif not empty_counts:
                X[:, col] = 1                             # every barcode a cell
    folder.mkdir(parents=True, exist_ok=True)
    mmwrite(str(folder / 'tmp.mtx'), csr_matrix(X))
    with open(folder / 'tmp.mtx', 'rb') as src, \
            gzip.open(folder / f'{prefix}matrix.mtx.gz', 'wb') as dst:
        shutil.copyfileobj(src, dst)
    (folder / 'tmp.mtx').unlink()
    with gzip.open(folder / f'{prefix}features.tsv.gz', 'wt') as f:
        f.writelines(f"ENSMUSG{i:05d}\tGene{i}\tGene Expression\n" for i in range(N_GENES))
    with gzip.open(folder / f'{prefix}barcodes.tsv.gz', 'wt') as f:
        f.writelines(f"BC{b:05d}-{s + 1}\n" for s in range(2) for b in range(N_BARCODES))
    return X


@pytest.fixture
def small_whitelist(monkeypatch):
    monkeypatch.setattr(large_10x, 'LARGE_BARCODES', 20)
    monkeypatch.setattr(large_10x, '_CHUNK', 1000)       # several chunks


def test_unfiltered_matrix_keeps_the_cells(tmp_path, small_whitelist):
    X = _write_10x(tmp_path)
    adata = load_10x_folder(tmp_path, cells_only=True)
    assert adata.n_obs == 2 * len(CELLS) and adata.n_vars == N_GENES
    assert list(adata.obs_names[:3]) == ['BC00003-1', 'BC00004-1', 'BC00005-1']
    np.testing.assert_array_equal(adata.X.toarray(), X[:, [3, 4, 5, 53, 54, 55]].T)
    record = adata.uns['kosmic_import']
    assert record['barcodes_in_file'] == 100 and record['barcodes_kept'] == 6
    assert record['barcodes_dropped'] == 94 and '>= 500 UMIs' in record['cell_filter']


def test_aggregated_barcodes_become_samples(tmp_path, small_whitelist):
    _write_10x(tmp_path, prefix='GSE1_')
    info = {'matrix': tmp_path / 'GSE1_matrix.mtx.gz',
            'features': tmp_path / 'GSE1_features.tsv.gz',
            'barcodes': tmp_path / 'GSE1_barcodes.tsv.gz'}
    adata = load_prefixed_10x(info, cells_only=True)
    large_10x.split_by_barcode_suffix(adata, 'GSE1')
    assert adata.obs['sample'].value_counts().to_dict() == {'GSE1_1': 3, 'GSE1_2': 3}
    assert adata.uns['kosmic_import']['samples_from'].startswith('barcode suffix')


def test_a_large_filtered_matrix_keeps_every_barcode(tmp_path, small_whitelist):
    _write_10x(tmp_path, empty_counts=False)
    adata = load_10x_folder(tmp_path, cells_only=True)
    assert adata.n_obs == 100
    assert 'cell_filter' not in adata.uns['kosmic_import']


def test_without_cells_only_every_barcode_is_read(tmp_path, small_whitelist):
    """The SoupX raw-matrix read needs the empty droplets."""
    _write_10x(tmp_path)
    assert load_10x_folder(tmp_path).n_obs == 100


def test_one_suffix_is_one_sample():
    from anndata import AnnData
    adata = AnnData(np.zeros((2, 1)))
    adata.obs_names = ['AAAC-1', 'AAAG-1']
    large_10x.split_by_barcode_suffix(adata, 'GSM1')
    assert set(adata.obs['sample']) == {'GSM1'}
    adata.obs['sample'] = 'kept'
    large_10x.split_by_barcode_suffix(adata, 'other')
    assert set(adata.obs['sample']) == {'kept'}
