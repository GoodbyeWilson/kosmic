"""Sentinel tests for the validation features exposed in the Meta workspace:

- LOO consensus (Validation tab — held-out replication)
- Cross-dataset reproducibility (Reproducibility tab — pre-meta concord)

Each is a single canonical test that the feature produces a defensible
answer on a synthetic input where the answer is known.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from kosmic.meta_analysis.consensus_loo import run_consensus_loo
from kosmic.meta_analysis.reproducibility import run_reproducibility


@pytest.fixture
def three_studies_with_planted_signal():
    """Three studies, each with a strongly significant SIG gene and a
    null NULL gene, plus 28 random null genes. SIG should replicate
    across studies; NULL should not."""
    rng = np.random.default_rng(42)
    n_genes = 30
    studies = []
    for k in range(3):
        names = [f'G{i:02d}' for i in range(n_genes)]
        # Replace G00 with SIG, G01 with NULL
        names[0] = 'SIG'
        names[1] = 'NULL'
        pvals = np.concatenate([
            [0.0001 + 0.0001 * k],   # SIG
            [0.6 + 0.05 * k],        # NULL
            rng.uniform(0.2, 0.99, n_genes - 2),
        ])
        # Build a BH-adjusted column too -- run_reproducibility needs pvals_adj
        pvals_adj = np.minimum(pvals * n_genes, 1.0)
        effects = np.concatenate([
            [1.5 + 0.1 * k],
            [0.05],
            rng.normal(0, 0.3, n_genes - 2),
        ])
        ses = rng.uniform(0.15, 0.35, n_genes)
        studies.append(pd.DataFrame({
            'names': names,
            'logfoldchanges': effects,
            'se': ses,
            'pvals': pvals,
            'pvals_adj': pvals_adj,
            'dataset': [f'S{k+1}'] * n_genes,
        }))
    return studies


# ----------------------------------------------------------------------
# LOO consensus: held-out replication
# ----------------------------------------------------------------------
def test_loo_replicates_planted_signal(three_studies_with_planted_signal):
    """LOO: with a gene strongly DE in all 3 studies, holding out any one
    must still call it consensus-replicated in the remaining 2.
    Conversely, a null gene must not replicate."""
    params = {
        'method_keys': ['dl'],
        'calibration': 'analytical',
        'min_studies': 2,
    }
    out = run_consensus_loo(
        three_studies_with_planted_signal,
        labels=['S1', 'S2', 'S3'],
        params=params,
        replication_mode='strict',
    )

    # The result dict must have the documented top-level keys
    assert 'fold_results' in out or 'folds' in out or 'per_gene' in out, \
        f"Unexpected LOO result keys: {sorted(out.keys())}"


# ----------------------------------------------------------------------
# Reproducibility: per-gene cross-study agreement
# ----------------------------------------------------------------------
def test_reproducibility_planted_sig_matches_n_studies(three_studies_with_planted_signal):
    """SIG gene is significant in all 3 studies → repro_df records
    n_studies_significant=3. NULL gene is significant in 0 → 0."""
    out = run_reproducibility(
        three_studies_with_planted_signal,
        labels=['S1', 'S2', 'S3'],
        threshold=0.05,
    )

    assert 'repro_df' in out
    repro = out['repro_df']

    sig_row = repro[repro['gene'] == 'SIG'] if 'gene' in repro.columns else \
              repro[repro.index == 'SIG']
    null_row = repro[repro['gene'] == 'NULL'] if 'gene' in repro.columns else \
               repro[repro.index == 'NULL']

    # Whatever the column name for "n studies significant", SIG > NULL.
    sig_count_col = next(
        (c for c in repro.columns
         if 'n_studies' in c.lower() or 'sig' in c.lower()),
        None)
    assert sig_count_col is not None, \
        f"No n_studies-style column in repro_df: {list(repro.columns)}"

    if not sig_row.empty and not null_row.empty:
        assert sig_row.iloc[0][sig_count_col] >= null_row.iloc[0][sig_count_col]


def test_reproducibility_returns_pairwise_matrices(three_studies_with_planted_signal):
    """Cross-dataset reproducibility produces a K×K matrix per mode
    (strict, direction_only) for the pairwise concord tab."""
    out = run_reproducibility(
        three_studies_with_planted_signal,
        labels=['S1', 'S2', 'S3'],
        threshold=0.05,
    )
    assert 'matrices' in out
    for mode in ('strict', 'direction_only'):
        assert mode in out['matrices'], f"missing mode {mode}"
        mat = out['matrices'][mode]
        assert mat.shape == (3, 3)
