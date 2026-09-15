"""Tests for the per-donor transcript-depth floor ('min_counts').

'min_cells' guards a pseudobulk profile by cell number; 'min_counts'
guards it by summed transcripts. These tests pin down that the floor is
decided in one place ('profile_table'), that 'create_pseudobulk', the
detection mask and the per-cell-type planner all apply it identically,
that depth is summed over every gene rather than the tested subset, and
that the 'total_counts' column round-trips through the pseudobulk CSV
without being mistaken for a gene.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")

from kosmic.de import batch as B  # noqa: E402
from kosmic.de.de_analysis import (  # noqa: E402
    cell_depths, create_pseudobulk, donors_meeting_min_cells, profile_table,
)
from kosmic.meta_analysis.cc_permutation import load_pseudobulk  # noqa: E402


def _adata(depth_per_cell, cells_per_donor=20, n_genes=30, seed=0):
    """Three disease + three control donors; *depth_per_cell* maps donor
    to the Poisson mean of every count, so a donor's depth is controllable."""
    rng = np.random.default_rng(seed)
    rows, blocks = [], []
    for role, donors in (('disease', ['D1', 'D2', 'D3']),
                         ('control', ['C1', 'C2', 'C3'])):
        for donor in donors:
            rows += [(donor, role)] * cells_per_donor
            blocks.append(rng.poisson(depth_per_cell[donor],
                                      size=(cells_per_donor, n_genes)))
    obs = pd.DataFrame(rows, columns=['donor', '_role'])
    obs['condition'] = obs['_role'].map({'disease': 'DCM', 'control': 'Donor'})
    obs['cell_type'] = 'Endothelial'
    obs.index = [f'cell{i}' for i in range(len(obs))]
    X = np.vstack(blocks).astype(np.float32)
    var = pd.DataFrame(index=[f'GENE{i}' for i in range(n_genes)])
    return anndata.AnnData(X=X, obs=obs, var=var)


UNIFORM = {d: 5.0 for d in ('D1', 'D2', 'D3', 'C1', 'C2', 'C3')}
ONE_THIN = {**UNIFORM, 'C2': 0.2}   # C2 sums to ~120 transcripts, others ~3,000


# --- profile_table ----------------------------------------------------------

def test_profile_table_reports_every_donor_with_depth():
    adata = _adata(UNIFORM)
    table = profile_table(adata, 'donor', min_cells=10, min_counts=0)
    assert table['sample'].tolist() == ['D1', 'D2', 'D3', 'C1', 'C2', 'C3']
    assert (table['n_cells'] == 20).all()
    expected = np.asarray(adata.X.sum(axis=1)).ravel()
    per_donor = pd.Series(expected).groupby(adata.obs['donor'].values).sum()
    assert np.allclose(table.set_index('sample')['total_counts'],
                       per_donor.loc[table['sample']].round())
    assert table['kept'].all()
    assert (table['reason'] == '').all()


def test_profile_table_floor_off_keeps_thin_donor():
    table = profile_table(_adata(ONE_THIN), 'donor', min_cells=10, min_counts=0)
    assert table['kept'].all()


def test_profile_table_floor_drops_thin_donor_with_reason():
    table = profile_table(_adata(ONE_THIN), 'donor', min_cells=10, min_counts=1000)
    dropped = table[~table['kept']]
    assert dropped['sample'].tolist() == ['C2']
    assert 'transcripts < 1,000' in dropped['reason'].iloc[0]


def test_profile_table_cell_floor_reported_before_depth():
    """A donor failing both floors is reported for the cell count."""
    table = profile_table(_adata(ONE_THIN), 'donor', min_cells=21, min_counts=1000)
    assert not table['kept'].any()
    assert all('cells < 21' in r for r in table['reason'])


def test_cell_depths_prefers_the_counts_layer():
    adata = _adata(UNIFORM)
    adata.layers['counts'] = adata.X.copy()
    adata.X = np.log1p(adata.X)          # a normalised X must not be summed
    counts_total = np.asarray(adata.layers['counts'].sum(axis=1)).ravel()
    assert np.allclose(cell_depths(adata), counts_total)
    assert not np.allclose(cell_depths(adata), np.asarray(adata.X.sum(axis=1)).ravel())


def test_cell_depths_unknown_layer_raises():
    with pytest.raises(ValueError, match="not found"):
        cell_depths(_adata(UNIFORM), counts_layer='nope')


# --- create_pseudobulk ------------------------------------------------------

def test_create_pseudobulk_default_keeps_thin_donor_and_reports_depth():
    adata = _adata(ONE_THIN)
    matrix, sample_df, genes = create_pseudobulk(
        adata, list(adata.var_names), 'donor', 'condition',
        min_cells=10, aggregate='sum')
    assert sample_df['sample'].tolist() == ['D1', 'D2', 'D3', 'C1', 'C2', 'C3']
    assert 'total_counts' in sample_df.columns
    # With every gene requested, the profile sum is the reported depth.
    assert np.allclose(matrix.sum(axis=1), sample_df['total_counts'])


def test_create_pseudobulk_floor_drops_thin_donor():
    adata = _adata(ONE_THIN)
    matrix, sample_df, _ = create_pseudobulk(
        adata, list(adata.var_names), 'donor', 'condition',
        min_cells=10, min_counts=1000, aggregate='sum')
    assert 'C2' not in sample_df['sample'].tolist()
    assert matrix.shape[0] == 5


def test_depth_is_over_all_genes_not_the_requested_subset():
    """Asking for two genes must not make every donor look thin."""
    adata = _adata(UNIFORM)
    matrix, sample_df, genes = create_pseudobulk(
        adata, ['GENE0', 'GENE1'], 'donor', 'condition',
        min_cells=10, min_counts=1000, aggregate='sum')
    assert genes == ['GENE0', 'GENE1']
    assert matrix.shape == (6, 2)
    assert (sample_df['total_counts'] > matrix.sum(axis=1)).all()


# --- detection mask and planner agree with create_pseudobulk ---------------

def test_detection_mask_excludes_the_same_donor():
    adata = _adata(ONE_THIN)
    assert donors_meeting_min_cells(adata, 'donor', 10, min_counts=0) is None
    mask = donors_meeting_min_cells(adata, 'donor', 10, min_counts=1000)
    assert mask is not None
    assert set(adata.obs.loc[~mask, 'donor']) == {'C2'}


def test_plan_cell_types_applies_the_depth_floor():
    adata = _adata(ONE_THIN)
    plan_off, = B.plan_cell_types(adata, 'cell_type', 'donor',
                                  min_cells=10, min_counts=0)
    plan_on, = B.plan_cell_types(adata, 'cell_type', 'donor',
                                 min_cells=10, min_counts=1000)
    assert plan_off.n_control_samples == 3
    assert plan_on.n_control_samples == 2
    assert plan_on.eligible


def test_plan_reason_names_the_depth_floor():
    adata = _adata(ONE_THIN)
    plan, = B.plan_cell_types(adata, 'cell_type', 'donor',
                              min_cells=10, min_counts=1_000_000)
    assert plan.eligible is False
    assert '1,000,000 transcripts' in plan.reason


# --- CSV round trip ---------------------------------------------------------

def test_pseudobulk_frame_places_total_counts_before_role():
    sample_df = pd.DataFrame({
        'sample': ['D1', 'C1'], 'condition': ['DCM', 'Donor'],
        'n_cells': [10, 12], 'total_counts': [50_000, 60_000],
        'role': ['disease', 'control'],
    })
    df = B.pseudobulk_frame(np.array([[1.0, 2.0], [3.0, 4.0]]),
                            sample_df, ['A', 'B'])
    assert list(df.columns) == [
        'condition', 'n_cells', 'total_counts', 'role', 'A', 'B']


def test_load_pseudobulk_does_not_read_total_counts_as_a_gene(tmp_path):
    sample_df = pd.DataFrame({
        'sample': ['D1', 'C1'], 'condition': ['DCM', 'Donor'],
        'n_cells': [10, 12], 'total_counts': [50_000, 60_000],
        'role': ['disease', 'control'],
    })
    path = tmp_path / 'X_pseudobulk.csv'
    B.pseudobulk_frame(np.array([[1.0, 2.0], [3.0, 4.0]]),
                       sample_df, ['A', 'B']).to_csv(path)
    expr, conds, genes, samples, n_cells, roles = load_pseudobulk(path)
    assert genes == ['A', 'B']
    assert expr.shape == (2, 2)
    assert n_cells.tolist() == [10, 12]


def test_load_pseudobulk_still_reads_files_without_total_counts(tmp_path):
    sample_df = pd.DataFrame({
        'sample': ['D1'], 'condition': ['DCM'], 'n_cells': [10]})
    path = tmp_path / 'X_pseudobulk.csv'
    B.pseudobulk_frame(np.array([[1.0]]), sample_df, ['A']).to_csv(path)
    expr, _, genes, _, _, _ = load_pseudobulk(path)
    assert genes == ['A'] and expr.shape == (1, 1)
