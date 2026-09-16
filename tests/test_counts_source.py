"""Where the raw counts are read from (kosmic.scrna.counts).

After QC a study holds its original counts in layers['counts'] and a
log-normalised X; .raw is no longer written. Older KOSMIC files and
foreign deposits (CellxGene) carry counts in .raw; a file QC has not
touched has them in X. Every reader resolves the counts through
count_source in that order, and counts_adata hands back a copy unless
told otherwise, so normalising the copy cannot touch the study.
"""
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from kosmic.scrna.counts import (
    count_source, count_var_names, counts_adata, has_counts_layer,
    log_normalised, looks_log_normalised,
)


def _study(n=20, g=6, seed=0):
    rng = np.random.default_rng(seed)
    X = sp.csr_matrix(rng.poisson(4, size=(n, g)).astype(np.float32))
    return ad.AnnData(X=X, obs=pd.DataFrame(index=[f'c{i}' for i in range(n)]),
                      var=pd.DataFrame(index=[f'G{j}' for j in range(g)]))


def test_fresh_import_reads_x():
    a = _study()
    X, names, where = count_source(a)
    assert where == 'X' and names == list(a.var_names)
    assert not has_counts_layer(a)


def test_after_qc_reads_counts_layer_not_x():
    a = _study()
    a.layers['counts'] = a.X.copy()
    a.X = sp.csr_matrix(np.log1p(a.X.toarray()))
    X, names, where = count_source(a)
    assert where == 'layers:counts'
    assert (X != a.layers['counts']).nnz == 0
    assert has_counts_layer(a)


def test_named_layer_wins_and_missing_layer_raises():
    import pytest
    a = _study()
    a.layers['counts'] = a.X.copy()
    a.layers['decontX_counts'] = a.X.copy() * 0.9
    assert count_source(a, 'decontX_counts')[2] == 'layers:decontX_counts'
    with pytest.raises(ValueError, match='decontX_missing'):
        count_source(a, 'decontX_missing')


def test_legacy_raw_is_read_when_no_layer():
    a = _study()
    a.raw = a.copy()
    a.X = sp.csr_matrix(np.log1p(a.X.toarray()))
    X, names, where = count_source(a)
    assert where == 'raw'
    assert count_var_names(a) == list(a.raw.var_names)


def test_counts_adata_copies_so_normalising_it_leaves_the_study_alone():
    a = _study()
    a.layers['counts'] = a.X.copy()
    before = a.layers['counts'].copy()
    work = counts_adata(a)
    work.X = work.X * 100.0                      # in-place style mutation
    assert (a.layers['counts'] != before).nnz == 0
    shared = counts_adata(a, copy=False)
    assert shared.X is a.layers['counts']


def test_counts_adata_gene_subset():
    a = _study()
    a.layers['counts'] = a.X.copy()
    sub = counts_adata(a, genes=['G3', 'G1', 'nope'])
    assert list(sub.var_names) == ['G3', 'G1']
    assert (sub.X != a.layers['counts'][:, [3, 1]]).nnz == 0


def test_log_normalised_uses_x_when_already_log_and_rebuilds_otherwise():
    a = _study()
    a.X = a.X * 30                                 # real counts run into the hundreds
    X, names = log_normalised(a)                  # X holds counts -> rebuilt
    assert looks_log_normalised(X) and not looks_log_normalised(a.X)
    assert names == list(a.var_names)
    a.layers['counts'] = a.X.copy()
    a.X = X
    a.uns['log1p'] = {}
    X2, _ = log_normalised(a)
    assert X2 is a.X
