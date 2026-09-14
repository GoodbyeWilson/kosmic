"""Tests for kosmic.scrna.load.converters: format ingestion + orientation
detection.

Focus: the genuinely tricky bits — orientation (cells×genes vs genes×cells)
and the per-format round-trip into h5ad. Format-extension detection is a
simple if-else and not exercised here.
"""
from __future__ import annotations


import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy.io import mmwrite
from scipy.sparse import csr_matrix


_BARCODES = [
    'AAACCTGAGAAACCAT-1', 'AAACCTGCATTTGCTT-1', 'AAACGGGAGTACGATA-1',
    'AAACGGGTCAGCAACT-1', 'AACCATGAGCCAACAG-1',
]


# ----------------------------------------------------------------------
# Orientation detection — the canonical scRNA loader gotcha
# ----------------------------------------------------------------------
class TestOrientationDetection:
    def test_cells_on_rows_recognised_via_barcode_index(self):
        from kosmic.scrna.load.converters import detect_matrix_orientation
        df = pd.DataFrame(
            np.random.rand(5, 3),
            index=_BARCODES, columns=['GAPDH', 'HK1', 'BRCA1'],
        )
        assert detect_matrix_orientation(df) == 'cells_x_genes'

    def test_genes_on_rows_recognised_via_gene_index(self):
        from kosmic.scrna.load.converters import detect_matrix_orientation
        df = pd.DataFrame(
            np.random.rand(3, 5),
            index=['GAPDH', 'HK1', 'BRCA1'], columns=_BARCODES,
        )
        assert detect_matrix_orientation(df) == 'genes_x_cells'

    def test_ambiguous_defaults_to_genes_on_rows(self):
        """When neither axis looks like barcodes, default to the more
        common scientific convention (genes×cells)."""
        from kosmic.scrna.load.converters import detect_matrix_orientation
        df = pd.DataFrame(
            np.random.rand(3, 3),
            index=['A', 'B', 'C'], columns=['X', 'Y', 'Z'],
        )
        assert detect_matrix_orientation(df) == 'genes_x_cells'


def test_looks_like_barcode_distinguishes_barcodes_from_gene_names():
    """The barcode heuristic powers orientation detection."""
    from kosmic.scrna.load.converters import looks_like_barcode
    assert looks_like_barcode(['AAACCTGAGAAACCAT-1', 'AAACCTGCATTTGCTT-1'])
    assert not looks_like_barcode(['GAPDH', 'HK1', 'BRCA1'])
    assert not looks_like_barcode([])


# ----------------------------------------------------------------------
# CSV / TSV ingest: both orientations + featureCounts edge case
# ----------------------------------------------------------------------
class TestCSVTSVRoundTrip:
    def test_csv_genes_on_rows(self, tmp_path):
        from kosmic.scrna.load.converters import csv_to_h5ad
        genes = ['GAPDH', 'HK1', 'MT-ND1']
        df = pd.DataFrame(np.random.rand(3, 5), index=genes, columns=_BARCODES)
        path = tmp_path / 'test.csv'
        df.to_csv(path)
        adata = csv_to_h5ad(str(path), str(tmp_path / 'out.h5ad'))
        assert (adata.n_obs, adata.n_vars) == (5, 3)

    def test_csv_cells_on_rows(self, tmp_path):
        from kosmic.scrna.load.converters import csv_to_h5ad
        df = pd.DataFrame(
            np.random.rand(3, 5),
            index=_BARCODES[:3],
            columns=['GAPDH', 'HK1', 'MT-ND1', 'BRCA1', 'TP53'],
        )
        path = tmp_path / 'test.csv'
        df.to_csv(path)
        adata = csv_to_h5ad(str(path), str(tmp_path / 'out.h5ad'))
        assert (adata.n_obs, adata.n_vars) == (3, 5)

    def test_tsv_featurecounts_format(self, tmp_path):
        """featureCounts output has Geneid/Chr/Start/End/Strand/Length
        annotation columns followed by per-sample count columns. The
        loader must drop the annotation columns and treat the rest as
        cells."""
        from kosmic.scrna.load.converters import tsv_to_h5ad
        df = pd.DataFrame({
            'Geneid': ['Gene1', 'Gene2', 'Gene3'],
            'Chr': ['chr1', 'chr2', 'chr3'],
            'Start': [100, 200, 300],
            'End': [200, 300, 400],
            'Strand': ['+', '-', '+'],
            'Length': [100, 100, 100],
            'Cell_A': [10, 20, 30],
            'Cell_B': [40, 50, 60],
        })
        path = tmp_path / 'fc.tsv'
        df.to_csv(path, sep='\t', index=False)
        adata = tsv_to_h5ad(str(path))
        assert (adata.n_obs, adata.n_vars) == (2, 3)


# ----------------------------------------------------------------------
# MTX + 10x folder ingest
# ----------------------------------------------------------------------
class TestMTXAnd10X:
    def test_mtx_round_trip(self, tmp_path):
        from kosmic.scrna.load.converters import mtx_to_h5ad
        n_genes, n_cells = 10, 5
        X = csr_matrix(np.random.rand(n_genes, n_cells).astype(np.float32))
        mmwrite(str(tmp_path / 'matrix.mtx'), X)
        with open(tmp_path / 'genes.tsv', 'w') as f:
            for i in range(n_genes):
                f.write(f"ENSG{i:04d}\tGene{i}\n")
        with open(tmp_path / 'barcodes.tsv', 'w') as f:
            for i in range(n_cells):
                f.write(f"AAACCTGA{i:08d}-1\n")
        adata = mtx_to_h5ad(str(tmp_path / 'matrix.mtx'),
                             str(tmp_path / 'out.h5ad'))
        assert (adata.n_obs, adata.n_vars) == (n_cells, n_genes)

    def test_10x_folder_loads_and_orients_cells_on_rows(self, tmp_path):
        from kosmic.scrna.load.converters import load_10x_folder
        n_genes, n_cells = 5, 3
        X = csr_matrix(np.random.rand(n_genes, n_cells).astype(np.float32))
        mmwrite(str(tmp_path / 'matrix.mtx'), X)
        with open(tmp_path / 'genes.tsv', 'w') as f:
            for i in range(n_genes):
                f.write(f"ENSG{i}\tGene{i}\n")
        with open(tmp_path / 'barcodes.tsv', 'w') as f:
            for i in range(n_cells):
                f.write(f"CELL{i}-1\n")
        adata = load_10x_folder(tmp_path)
        assert (adata.n_obs, adata.n_vars) == (n_cells, n_genes)

    def test_10x_folder_missing_matrix_raises(self, tmp_path):
        from kosmic.scrna.load.converters import load_10x_folder
        with pytest.raises(FileNotFoundError):
            load_10x_folder(tmp_path)


# ----------------------------------------------------------------------
# R / Seurat export reassembly
# ----------------------------------------------------------------------
def test_assemble_from_r_export_normalises_orig_ident_to_sample(tmp_path):
    """Seurat's ``orig.ident`` column must be standardised to ``sample``
    so the rest of the pipeline can find it."""
    from kosmic.scrna.load.converters import assemble_h5ad_from_r_export
    prefix, n_genes, n_cells = 'sample1', 4, 3
    X = csr_matrix(np.random.rand(n_genes, n_cells).astype(np.float32))
    mmwrite(str(tmp_path / f'{prefix}_expression.mtx'), X)
    cells = [f'Cell{i}' for i in range(n_cells)]
    (tmp_path / f'{prefix}_genes.txt').write_text(
        '\n'.join(f'Gene{i}' for i in range(n_genes)))
    (tmp_path / f'{prefix}_cells.txt').write_text('\n'.join(cells))
    pd.DataFrame(
        {'orig.ident': ['s1'] * n_cells, 'nCount_RNA': [100, 200, 300]},
        index=cells,
    ).to_csv(tmp_path / f'{prefix}_metadata.csv')
    adata = assemble_h5ad_from_r_export(
        tmp_path, prefix, tmp_path / f'{prefix}.h5ad')
    assert (adata.n_obs, adata.n_vars) == (n_cells, n_genes)
    assert 'sample' in adata.obs.columns


# ----------------------------------------------------------------------
# Sample-name extraction: a real-world bug source for GSM-prefixed files
# ----------------------------------------------------------------------
def test_extract_sample_name_from_gsm_prefix():
    from kosmic.scrna.load.converters import _extract_sample_name
    assert _extract_sample_name('GSM8761493_matrix_gene_ND178.tsv.gz') == 'ND178'


def test_extract_sample_name_falls_back_to_parent_folder(tmp_path):
    """When the filename itself is generic ('matrix.tsv'), the parent
    folder name supplies the sample identifier."""
    from kosmic.scrna.load.converters import _extract_sample_name
    sub = tmp_path / 'MyProject'
    sub.mkdir()
    path = sub / 'matrix.tsv'
    path.touch()
    assert _extract_sample_name(str(path)) == 'MyProject'


# ----------------------------------------------------------------------
# AnnData summary: distinguishes raw counts from normalised data
# ----------------------------------------------------------------------
class TestSummaryDetectsNormalisation:
    def test_normalised_data_flagged_as_normalised(self):
        from kosmic.scrna.load.converters import build_adata_summary
        # log-normalised values are typically O(1)-O(5)
        adata = AnnData(X=np.random.rand(10, 5).astype(np.float32) * 5)
        assert bool(build_adata_summary(adata)['is_normalized']) is True

    def test_raw_count_data_flagged_as_not_normalised(self):
        from kosmic.scrna.load.converters import build_adata_summary
        X = csr_matrix(np.random.poisson(100, (10, 5)).astype(np.float32))
        adata = AnnData(X=X)
        assert bool(build_adata_summary(adata)['is_normalized']) is False


# ----------------------------------------------------------------------
# h5ad load: Ensembl-indexed files swap to symbols but keep the stable ID
# ----------------------------------------------------------------------
def test_ensembl_indexed_h5ad_swaps_to_symbols_and_preserves_id(tmp_path):
    from kosmic.scrna.load.converters import load_h5ad_with_summary

    # cellxgene-style: var indexed by Ensembl ID, symbol in 'feature_name'
    var = pd.DataFrame(
        {'feature_name': ['HMGB2', 'ACTA2', 'MYH11']},
        index=['ENSG00000164104', 'ENSG00000107796', 'ENSG00000133392'],
    )
    adata = AnnData(X=np.random.rand(6, 3).astype(np.float32), var=var)
    adata.raw = adata.copy()
    path = tmp_path / 'ensembl.h5ad'
    adata.write_h5ad(path)

    loaded, _ = load_h5ad_with_summary(path)

    # Index now holds symbols; the gene picker can find HMGB2
    assert 'HMGB2' in set(loaded.var_names)
    assert not any(str(v).startswith('ENSG') for v in loaded.var_names)
    # Stable Ensembl IDs are not lost — stashed for provenance
    assert 'ensembl_id' in loaded.var.columns
    assert 'ENSG00000164104' in set(loaded.var['ensembl_id'])
    # raw slot mirrors the swap and retains its Ensembl IDs too
    assert 'HMGB2' in set(loaded.raw.var_names)
    assert 'ENSG00000164104' in set(loaded.raw.var['ensembl_id'])
