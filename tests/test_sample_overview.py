# Sentinel for the Inspect Samples overview: per-sample counts, roles,
# and QC metrics derived from X when scanpy QC columns are absent.
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from kosmic.scrna.inspect.batch import sample_overview


def _toy_adata():
    # 6 cells, 3 genes (one mitochondrial), 2 samples across 2 roles
    X = sp.csr_matrix(np.array([
        [5, 0, 5],   # s1
        [4, 1, 5],   # s1
        [0, 2, 0],   # s1
        [3, 3, 0],   # s2
        [2, 2, 2],   # s2
        [1, 0, 1],   # s2
    ], dtype=float))
    obs = pd.DataFrame({
        'sample': ['s1', 's1', 's1', 's2', 's2', 's2'],
        '_role': ['disease'] * 3 + ['control'] * 3,
    })
    var = pd.DataFrame(index=['GENE1', 'MT-CO1', 'GENE2'])
    return ad.AnnData(X=X, obs=obs, var=var)


def test_sample_overview_counts_and_metrics():
    result = sample_overview(_toy_adata())
    assert result['n_cells'] == 6
    assert result['n_genes'] == 3
    assert result['n_samples'] == 2
    assert result['n_conditions'] == 2

    table = result['table'].set_index('sample')
    assert table.loc['s1', 'n_cells'] == 3
    assert table.loc['s1', 'condition'] == 'disease'
    assert table.loc['s2', 'condition'] == 'control'
    # s2 genes detected per cell: 2, 3, 2 -> median 2
    assert table.loc['s2', 'median_genes'] == 2
    # Mito derived from the MT- gene without scanpy QC columns
    assert table.loc['s1', 'pct_mito'] > 0


def test_sample_overview_degrades_without_sample_column():
    adata = _toy_adata()
    adata.obs = adata.obs.drop(columns=['sample'])
    result = sample_overview(adata)
    assert result['n_cells'] == 6
    assert result['n_samples'] == 0
    assert result['table'].empty
