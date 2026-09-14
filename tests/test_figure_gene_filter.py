"""Sentinel: per-gene figures (in-app and exported) restrict to the DE-tested
gene set via the shared 'tested_gene_set' helper, so a gene the DE table filtered
out (e.g. expressed in one donor) never appears in a heatmap. Applying the filter
in shared code keeps the in-app and exported plots identical."""
from __future__ import annotations

import numpy as np
import pandas as pd

# Import the module (not the symbol) so pytest doesn't collect the
# 'tested_gene_set' function as a test case.
from kosmic.de import de_analysis


def test_no_de_results_returns_none():
    assert de_analysis.tested_gene_set(None) is None


def test_only_tested_genes_pass_the_filter():
    dr = pd.DataFrame({
        'names': ['ENO1', 'PKM', 'NDUFA9', 'SDHB'],
        'filter_status': ['tested', 'tested', 'independent_filter', 'cooks_outlier'],
    })
    tested = de_analysis.tested_gene_set(dr)
    assert tested == {'ENO1', 'PKM'}
    assert 'NDUFA9' not in tested  # one-donor gene DESeq2 independent-filtered


def test_no_filter_status_column_falls_back_to_all_names():
    dr = pd.DataFrame({'names': ['ENO1', 'PKM']})
    assert de_analysis.tested_gene_set(dr) == {'ENO1', 'PKM'}


def test_prepare_heatmap_data_honours_allowed_genes():
    """prepare_heatmap_data drops genes outside allowed_genes, so both the in-app
    and exported heatmaps show only DE-tested genes."""
    import anndata as ad
    from kosmic.visualisation.de.gene_heatmap import prepare_heatmap_data

    rng = np.random.default_rng(0)
    genes = ['ENO1', 'PKM', 'LDHA', 'NDUFA9']
    n = 40
    X = rng.integers(1, 50, size=(n, len(genes))).astype(np.float32)
    obs = pd.DataFrame({
        'sample': ['d%d' % (i % 4) for i in range(n)],
        'condition': (['control'] * (n // 2)) + (['disease'] * (n // 2)),
    }, index=[f'c{i}' for i in range(n)])
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))

    gs = {'Glycolysis': ['ENO1', 'PKM', 'LDHA', 'NDUFA9']}
    _, avail = prepare_heatmap_data(
        adata, gs, 'sample', 'condition', min_cells=1,
        allowed_genes={'ENO1', 'PKM', 'LDHA'})
    assert set(avail['Glycolysis']) == {'ENO1', 'PKM', 'LDHA'}
    assert 'NDUFA9' not in avail['Glycolysis']
