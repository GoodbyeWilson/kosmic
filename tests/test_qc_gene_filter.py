"""The QC filter keeps every gene unless told otherwise, and reports the
gene count either way.

A gene seen in a handful of cells is a measured near-zero; dropping it
made each study's gene list depend on itself (Reichart-LV lost 3,126
genes to the old default of 3). The default is now 0 and the stats
always carry genes before and after.
"""
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from kosmic.scrna.qc.filter import run_qc_pipeline


def _adata():
    rng = np.random.default_rng(0)
    X = rng.poisson(2, size=(100, 30)).astype(np.float32)
    X[:, 25:] = 0                     # five genes seen in no cell
    X[0, 25] = 1                      # one of them seen in exactly one cell
    genes = [f'G{i}' for i in range(30)]
    a = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(index=[f'c{i}' for i in range(100)]),
                   var=pd.DataFrame(index=genes))
    return a


def test_default_keeps_every_gene():
    a, stats = run_qc_pipeline(_adata(), {'min_genes': 0, 'max_mt': 0})
    assert a.n_vars == 30
    assert stats['n_genes_before'] == 30 and stats['n_genes_after'] == 30


def test_explicit_min_cells_drops_and_reports():
    a, stats = run_qc_pipeline(_adata(), {'min_genes': 0, 'max_mt': 0, 'min_cells': 3})
    assert a.n_vars == 25
    assert stats['n_genes_before'] == 30 and stats['n_genes_after'] == 25
    assert stats['genes_removed'] == 5
