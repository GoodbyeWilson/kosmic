"""Tests for kosmic.scrna.inspect.roles.de_readiness: the scRNA -> DE handoff gate."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from kosmic.scrna.inspect.roles import de_readiness


def _adata(X, obs=None, uns=None):
    n_obs = X.shape[0]
    obs = obs if obs is not None else pd.DataFrame(index=[f'c{i}' for i in range(n_obs)])
    a = ad.AnnData(X=X, obs=obs)
    if uns:
        a.uns.update(uns)
    return a


def test_none_adata_not_ready():
    ready, reason = de_readiness(None)
    assert not ready
    assert reason


def test_raw_counts_not_normalised_blocks_de():
    """Un-normalised data (large integer counts) must not be handed to DE
    even if a role_map is present -- the '_role' handoff we care about
    doesn't matter if pseudobulk would be computed on the wrong scale."""
    X = np.random.default_rng(0).integers(0, 500, size=(10, 5)).astype(np.float32)
    obs = pd.DataFrame({'condition': ['A'] * 5 + ['B'] * 5},
                        index=[f'c{i}' for i in range(10)])
    adata = _adata(X, obs, uns={'role_map': {'A': 'control', 'B': 'disease'}})
    ready, reason = de_readiness(adata)
    assert not ready
    assert 'normalis' in reason.lower()


def test_normalised_no_role_assignment_blocks_de():
    X = np.random.default_rng(0).random((10, 5)).astype(np.float32) * 5
    adata = _adata(X)
    ready, reason = de_readiness(adata)
    assert not ready
    assert 'role' in reason.lower()


def test_normalised_with_role_map_is_ready():
    X = np.random.default_rng(0).random((10, 5)).astype(np.float32) * 5
    obs = pd.DataFrame({'condition': ['NF'] * 5 + ['DCM'] * 5},
                        index=[f'c{i}' for i in range(10)])
    adata = _adata(X, obs, uns={'role_map': {'NF': 'control', 'DCM': 'disease'}})
    ready, reason = de_readiness(adata)
    assert ready
    assert reason == ""


def test_normalised_with_resolved_role_column_is_ready():
    X = np.random.default_rng(0).random((10, 5)).astype(np.float32) * 5
    obs = pd.DataFrame({'_role': ['control'] * 5 + ['disease'] * 5},
                        index=[f'c{i}' for i in range(10)])
    adata = _adata(X, obs)
    ready, reason = de_readiness(adata)
    assert ready
    assert reason == ""


def test_role_map_missing_one_side_not_ready():
    """A role_map with only control (or only disease) values assigned
    isn't a usable two-group contrast yet."""
    X = np.random.default_rng(0).random((10, 5)).astype(np.float32) * 5
    obs = pd.DataFrame({'condition': ['NF'] * 10}, index=[f'c{i}' for i in range(10)])
    adata = _adata(X, obs, uns={'role_map': {'NF': 'control'}})
    ready, reason = de_readiness(adata)
    assert not ready
