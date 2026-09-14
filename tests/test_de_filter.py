"""Sentinel for the sample-level expression-prevalence gene filter.

The point of this filter is donor-awareness: a gene carried by a single
dominant donor must not survive, because it can't be reproducible across the
cohort. Auto mode requires detection in at least the smaller group's size.
"""
from __future__ import annotations

import numpy as np

from kosmic.de import filter_by_expression


def test_filter_by_expression_is_donor_aware():
    # 4 disease + 4 control donors, three genes.
    roles = np.array(['disease'] * 4 + ['control'] * 4)
    counts = np.array([
        [500, 800, 0],   # g0 broad, g1 one dominant donor, g2 absent
        [500, 0, 0],
        [500, 0, 0],
        [500, 0, 0],
        [500, 0, 0],
        [500, 0, 0],
        [500, 0, 0],
        [500, 0, 0.0],
    ])

    # Auto: smaller group = 4, so g1 (one donor) is dropped; g0 kept; g2 dropped.
    mask = filter_by_expression(counts, roles, min_count=10)
    assert mask.tolist() == [True, False, False]

    # Overriding min_samples down to 1 recovers the single-donor gene.
    mask1 = filter_by_expression(counts, roles, min_count=10, min_samples=1)
    assert mask1.tolist() == [True, True, False]

    # A near-zero total-count gene is dropped by the total-count floor even if a
    # couple of donors have a lone read.
    low = np.zeros((8, 1))
    low[0, 0] = 3
    low[4, 0] = 3
    assert filter_by_expression(low, roles, min_count=10,
                                min_total_count=15, min_samples=1).tolist() == [False]
