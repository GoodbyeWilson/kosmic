"""Tests for kosmic.scrna.inspect.gene_names: HGNC symbol harmonisation."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from kosmic.scrna.inspect.gene_names import (
    harmonise_adata,
    load_hgnc_lookup,
)


def test_hgnc_lookup_excludes_approved_symbols():
    """An approved symbol like GAPDH must not appear as a key -- otherwise
    a current symbol that happens to be an alias for another would get
    silently re-mapped, corrupting the count matrix."""
    lookup, _ = load_hgnc_lookup()
    assert 'GAPDH' not in lookup
    assert lookup.get('FAM206A') == 'ABITRAM'   # known previous-symbol mapping


def test_mito_legacy_names_resolve_to_mt_not_nuclear():
    """A legacy mitochondrial label ('ND1') must resolve to MT-ND1, never to the
    colliding nuclear gene (IVNS1ABP, which also lists 'ND1' as an alias)."""
    lookup, reasons = load_hgnc_lookup()
    assert lookup.get('ND1') == 'MT-ND1'
    assert reasons.get('ND1') == 'mitochondrial'
    assert lookup.get('CYTB') == 'MT-CYB'

    rng = np.random.default_rng(0)
    X = rng.integers(0, 50, size=(6, 4)).astype(np.float32)
    adata = ad.AnnData(
        X=X, obs=pd.DataFrame(index=[f'c{i}' for i in range(6)]),
        var=pd.DataFrame(index=['ND1', 'CYTB', 'GAPDH', 'ACTB']))
    result, _ = harmonise_adata(adata)
    assert 'MT-ND1' in result.var_names and 'MT-CYB' in result.var_names
    assert 'IVNS1ABP' not in result.var_names


def test_ambiguous_nonmito_alias_is_not_remapped():
    """A non-mitochondrial alias claimed by two approved symbols is left alone
    (never guessed), so it can't silently corrupt either gene."""
    lookup, _ = load_hgnc_lookup()
    assert 'MRP2' not in lookup  # claimed by multiple approved symbols


def test_qc_flags_legacy_mito_names_without_harmonisation():
    """mito-% QC must catch mitochondrial genes under legacy names (ND1, CYTB),
    not just MT-* -- otherwise skipping harmonisation silently breaks mito
    filtering."""
    from kosmic.scrna.qc.filter import detect_mitochondrial_genes
    X = np.ones((3, 5), dtype=np.float32)
    adata = ad.AnnData(
        X=X, obs=pd.DataFrame(index=[f'c{i}' for i in range(3)]),
        var=pd.DataFrame(index=['ND1', 'CYTB', 'MT-CO1', 'GAPDH', 'ACTB']))
    detect_mitochondrial_genes(adata)
    flagged = set(adata.var_names[adata.var['mt']])
    assert flagged == {'ND1', 'CYTB', 'MT-CO1'}


def test_harmonise_adata_renames_aliases_and_preserves_data():
    """Aliased gene names get rewritten; non-aliased names pass through;
    the count matrix is unchanged."""
    rng = np.random.default_rng(42)
    X = rng.integers(0, 100, size=(10, 5)).astype(np.float32)
    var = pd.DataFrame(index=['GAPDH', 'FAM206A', 'ACTB', 'FAM175A', 'TP53'])
    obs = pd.DataFrame(index=[f'cell_{i}' for i in range(10)])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    result, report = harmonise_adata(adata)

    assert {'GAPDH', 'ACTB', 'TP53', 'ABITRAM', 'ABRAXAS1'} == set(result.var_names)
    assert report.renamed == 2
    np.testing.assert_array_equal(result.X, X)


def test_harmonise_merges_duplicates_when_two_aliases_share_canonical():
    """Two aliases pointing at the same approved symbol must merge into one
    column with summed counts -- otherwise downstream DE sees a phantom
    duplicate gene."""
    lookup = {'ALIAS_A': 'MERGED', 'ALIAS_B': 'MERGED'}
    reasons = {'ALIAS_A': 'alias', 'ALIAS_B': 'alias'}

    X = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    var = pd.DataFrame(index=['ALIAS_A', 'ALIAS_B', 'NORMAL'])
    obs = pd.DataFrame(index=['c1', 'c2'])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    result, report = harmonise_adata(adata, lookup=lookup, reasons=reasons)

    assert result.n_vars == 2
    assert 'MERGED' in result.var_names
    assert report.duplicates_merged == 1
    merged_idx = list(result.var_names).index('MERGED')
    np.testing.assert_array_equal(result.X[:, merged_idx], [3, 9])
