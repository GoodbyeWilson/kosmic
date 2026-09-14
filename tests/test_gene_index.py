"""Gene-identifier normalisation on import.

Studies are combined by joining on ``var_names``. Depositors disagree about
which identifier goes on the index -- CellxGene mandates Ensembl IDs, most
other sources use symbols -- and the mismatch is silent and total: an
Ensembl-indexed study shares no genes at all with a symbol-indexed one, so an
inner join across them returns nothing. These sentinels cover the three
layouts actually seen in the DCM cohort.
"""
import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from kosmic.scrna.load.converters import normalise_gene_index

ENSG = ["ENSG00000155657", "ENSG00000175206", "ENSG00000092054"]
SYMBOLS = ["TTN", "NPPA", "MYH7"]


def _adata(var):
    return ad.AnnData(X=csr_matrix(np.eye(len(var), dtype=np.float32)), var=var)


def test_cellxgene_layout_swaps_index_to_symbols():
    """Ensembl index + feature_name column -- reichart's layout."""
    a = _adata(pd.DataFrame({"feature_name": SYMBOLS}, index=ENSG))
    normalise_gene_index(a)
    assert list(a.var_names) == SYMBOLS
    assert list(a.var["ensembl_id"]) == ENSG


def test_cellranger_layout_surfaces_ensembl_column():
    """Symbol index + gene_ids column -- chaffin's layout."""
    a = _adata(pd.DataFrame({"gene_ids": ENSG}, index=SYMBOLS))
    normalise_gene_index(a)
    assert list(a.var_names) == SYMBOLS      # index untouched
    assert list(a.var["ensembl_id"]) == ENSG


def test_seurat_layout_invents_nothing():
    """Symbols only -- Seurat objects do not retain Ensembl IDs."""
    a = _adata(pd.DataFrame(index=SYMBOLS))
    normalise_gene_index(a)
    assert list(a.var_names) == SYMBOLS
    assert "ensembl_id" not in a.var.columns


def test_ensembl_index_without_symbols_is_left_alone():
    """No symbol column to swap to: better an Ensembl index than a broken one."""
    a = _adata(pd.DataFrame(index=ENSG))
    normalise_gene_index(a)
    assert list(a.var_names) == ENSG


def test_existing_ensembl_id_column_is_not_overwritten():
    a = _adata(pd.DataFrame(
        {"gene_ids": ENSG, "ensembl_id": ["keep", "these", "values"]},
        index=SYMBOLS))
    normalise_gene_index(a)
    assert list(a.var["ensembl_id"]) == ["keep", "these", "values"]


def test_symbol_index_is_not_mistaken_for_ensembl():
    """A stray ENSG-looking symbol must not flip the whole index."""
    var = pd.DataFrame(index=["TTN", "NPPA", "ENSG00000999999"])
    a = _adata(var)
    normalise_gene_index(a)
    assert list(a.var_names) == ["TTN", "NPPA", "ENSG00000999999"]


def test_duplicate_symbols_are_made_unique():
    """Two Ensembl IDs can map to one symbol; the index must stay unique."""
    a = _adata(pd.DataFrame(
        {"feature_name": ["TTN", "TTN", "MYH7"]}, index=ENSG))
    normalise_gene_index(a)
    assert len(set(a.var_names)) == 3
    assert list(a.var["ensembl_id"]) == ENSG
