"""The Combine dialog's study scan counts the cells the roles exclude."""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from kosmic.gui.combine.workspace import _read_obs_metadata


def _write(path: Path, roles):
    n = len(roles)
    a = ad.AnnData(X=sp.csr_matrix(np.ones((n, 3), dtype=np.float32)),
                   obs=pd.DataFrame({'_role': pd.Categorical(roles),
                                     'condition': pd.Categorical(['x'] * n)},
                                    index=[f'c{i}' for i in range(n)]))
    a.write_h5ad(path)
    return path


def test_scan_reports_excluded_cells(tmp_path):
    m = _read_obs_metadata(_write(tmp_path / 's.h5ad', ['disease'] * 3 + ['exclude'] * 2 + ['control']))
    assert m['n_obs'] == 6 and m['n_excluded'] == 2


def test_scan_without_roles_reports_none_excluded(tmp_path):
    n = 4
    a = ad.AnnData(X=sp.csr_matrix(np.ones((n, 3), dtype=np.float32)),
                   obs=pd.DataFrame({'condition': pd.Categorical(['x'] * n)},
                                    index=[f'c{i}' for i in range(n)]))
    a.write_h5ad(tmp_path / 's.h5ad')
    m = _read_obs_metadata(tmp_path / 's.h5ad')
    assert m['n_obs'] == 4 and m['n_excluded'] == 0
