"""A study file written under pandas 3 can be saved again under pandas 2.

pandas 3 makes text columns a 'string' dtype, and anndata writes them as
nullable-string arrays. Read back under pandas 2 they stay 'string', and
anndata refuses to write them unless 'allow_write_nullable_strings' is
on. Every study saved during the pandas-3 period was in that format, so
under pandas 2 every save failed -- and a save that overwrote the study
file in place left it truncated.
"""
from __future__ import annotations

import os

import h5py
import numpy as np
import pandas as pd
import pytest

import kosmic  # noqa: F401  -- importing the package applies the setting

ad = pytest.importorskip("anndata")


@pytest.fixture
def pandas3_file(tmp_path):
    """A small h5ad whose names are stored as nullable-string arrays."""
    a = ad.AnnData(X=np.array([[1, 0, 2], [0, 3, 0], [4, 0, 0]], dtype=np.float32))
    a.obs_names = pd.Index(['c1', 'c2', 'c3'], dtype='string')
    a.var_names = pd.Index(['GENE1', 'GENE2', 'GENE3'], dtype='string')
    a.obs['sample'] = pd.array(['s1', 's2', None], dtype='string')
    path = tmp_path / 'pandas3.h5ad'
    with ad.settings.override(allow_write_nullable_strings=True):
        a.write_h5ad(path)
    with h5py.File(path, 'r') as f:
        assert f['obs/_index'].attrs['encoding-type'] == 'nullable-string-array'
    return path


def test_importing_kosmic_enables_writing_nullable_strings():
    assert os.environ.get('ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS') == '1'
    assert ad.settings.allow_write_nullable_strings is True


def test_a_pandas3_era_file_saves_again(pandas3_file, tmp_path):
    a = ad.read_h5ad(pandas3_file)
    out = tmp_path / 'resaved.h5ad'
    a.write_h5ad(out)

    b = ad.read_h5ad(out)
    assert list(b.obs_names) == ['c1', 'c2', 'c3']
    assert list(b.var_names) == ['GENE1', 'GENE2', 'GENE3']
    assert pd.isna(b.obs['sample'].iloc[2])
    np.testing.assert_array_equal(b.X, a.X)

