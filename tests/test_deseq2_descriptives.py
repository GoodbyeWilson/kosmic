"""The DESeq2 table's descriptive columns (group means, variances, Cohen's d)
are on normalised counts. On raw pseudobulk counts they scaled with each
sample's depth: Youness's TPM3 showed a disease variance of 1.3 million
where the normalised values were tightly grouped."""
from __future__ import annotations

import numpy as np
import pandas as pd


def test_depth_differences_do_not_inflate_the_spread():
    from kosmic.de.de_analysis import _run_deseq2
    rng = np.random.default_rng(0)
    n_genes = 200
    base = rng.uniform(50, 500, n_genes)
    depth = np.array([1, 10, 1, 10, 1, 10, 1, 10], dtype=float)
    counts = rng.poisson(base * depth[:, None] * 5).astype(float)
    genes = [f'g{i}' for i in range(n_genes)]
    sample_df = pd.DataFrame({'sample': [f's{i}' for i in range(8)],
                              'role': ['disease', 'control'] * 4})

    res = _run_deseq2(counts, sample_df, genes).set_index('names')

    # Same profile in every sample: the normalised spread is Poisson-sized,
    # far below the spread of the raw counts.
    raw = pd.DataFrame(counts, columns=genes)
    raw_cv = raw.std(ddof=1) / raw.mean()
    cv = np.sqrt(res['disease_var']) / res['disease_mean']
    assert cv.median() < 0.05 < raw_cv.median()
    assert res['cohens_d'].abs().median() < 1
    # The means are on the normalised scale, between the two depths.
    assert (res['disease_mean'] / (base * 5)).between(1, 10).all()
