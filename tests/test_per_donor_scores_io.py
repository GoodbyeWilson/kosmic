"""Sentinel: per-donor pathway scores round-trip through disk unchanged, so the
figure-export patient dotplot can reload the exact values a scoring run produced
instead of recomputing them."""
from __future__ import annotations

import numpy as np
import pandas as pd

from kosmic.de.pathway_scoring import (
    save_per_donor_scores, load_per_donor_scores,
)


def test_per_donor_scores_round_trip(tmp_path):
    scores = {
        'names': ['Glycolysis', 'TCA_Cycle', 'Oxidative_Phosphorylation'],
        'matrix': np.array([[0.1, -0.2, 0.3],
                            [0.4, 0.5, -0.6],
                            [-0.7, 0.8, 0.9]], dtype=float),
        'sample_df': pd.DataFrame({
            'sample': ['d1', 'd2', 'd3'],
            'condition': ['control', 'disease', 'disease'],
            'n_cells': [100, 120, 90],
        }),
    }
    path = save_per_donor_scores(scores, tmp_path / 'pathway_scores_per_donor.csv')
    assert path.exists()

    back = load_per_donor_scores(path)
    assert back['names'] == scores['names']
    np.testing.assert_allclose(back['matrix'], scores['matrix'])
    assert list(back['sample_df']['condition']) == ['control', 'disease', 'disease']
    assert list(back['sample_df']['sample']) == ['d1', 'd2', 'd3']


def test_load_missing_file_returns_none(tmp_path):
    assert load_per_donor_scores(tmp_path / 'nope.csv') is None
