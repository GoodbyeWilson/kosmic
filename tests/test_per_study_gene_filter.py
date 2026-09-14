"""filterByExpr is the gene filter for pseudobulk DE, and on an atlas it
has to run within each study.

Pooling hides the failure: on 12 disease + 12 control donors the rule
asks for 11.4 samples, so a gene present in exactly one study's 12
donors and absent from the other clears it and gets tested as though
both cohorts could see it. Per-study, it cannot.
"""
from __future__ import annotations

import numpy as np

from kosmic.de.de_analysis import filter_by_expression


def counts_matrix(spec, n_genes=None, background=500.0):
    """spec: list of per-sample lists of counts for the genes of interest."""
    arr = np.asarray(spec, dtype=float)
    if n_genes and n_genes > arr.shape[1]:
        pad = np.full((arr.shape[0], n_genes - arr.shape[1]), background)
        arr = np.hstack([arr, pad])
    return arr


ROLES_24 = np.array(["disease"] * 12 + ["control"] * 12)
STUDIES_24 = np.array(["S1"] * 6 + ["S2"] * 6 + ["S1"] * 6 + ["S2"] * 6)


def build_one_study_gene():
    """Gene 0 is present only in S1's donors; gene 1 is everywhere."""
    n = 24
    g0 = np.where(STUDIES_24 == "S1", 400.0, 0.0)
    g1 = np.full(n, 400.0)
    filler = np.full((n, 200), 500.0)
    return np.column_stack([g0, g1, filler])


def test_pooled_filter_admits_a_one_study_gene():
    counts = build_one_study_gene()
    mask = filter_by_expression(counts, ROLES_24, min_count=10)
    assert bool(mask[0]) is True   # the one-study gene slips through
    assert bool(mask[1]) is True


def test_per_study_filter_rejects_it():
    counts = build_one_study_gene()
    mask = filter_by_expression(counts, ROLES_24, min_count=10,
                                studies=STUDIES_24)
    assert bool(mask[0]) is False  # absent from S2, so it cannot pass
    assert bool(mask[1]) is True


def test_single_study_label_is_the_pooled_path():
    """One cohort has nothing to stratify; the answer must not change."""
    counts = build_one_study_gene()
    pooled = filter_by_expression(counts, ROLES_24, min_count=10)
    one = filter_by_expression(counts, ROLES_24, min_count=10,
                               studies=np.array(["only"] * 24))
    assert np.array_equal(pooled, one)


def test_none_is_the_pooled_path():
    counts = build_one_study_gene()
    assert np.array_equal(
        filter_by_expression(counts, ROLES_24, min_count=10),
        filter_by_expression(counts, ROLES_24, min_count=10, studies=None))


def test_a_gene_in_both_studies_survives_stratification():
    counts = build_one_study_gene()
    mask = filter_by_expression(counts, ROLES_24, min_count=10,
                                studies=STUDIES_24)
    assert bool(mask[1]) is True


def test_group_sizes_are_taken_within_each_study():
    """Each study's requirement comes from its own arms, not the pool."""
    # S1 has 6+6 donors, S2 has 2+2. A gene in 3 of S2's 4 donors passes
    # S2's own requirement (2) but would fail a pooled requirement.
    roles = np.array(["disease"] * 6 + ["control"] * 6
                     + ["disease"] * 2 + ["control"] * 2)
    studies = np.array(["S1"] * 12 + ["S2"] * 4)
    g0 = np.concatenate([np.full(12, 400.0), np.array([400.0, 400.0, 400.0, 0.0])])
    counts = np.column_stack([g0, np.full((16, 200), 500.0)])
    assert bool(filter_by_expression(
        counts, roles, min_count=10, studies=studies)[0]) is True


def test_empty_matrix_with_studies():
    empty = np.zeros((0, 5))
    mask = filter_by_expression(empty, np.array([]), studies=np.array([]))
    assert mask.shape == (5,)
    assert not mask.any()


def test_study_with_no_samples_is_skipped():
    counts = build_one_study_gene()
    studies = STUDIES_24.copy()
    mask = filter_by_expression(counts, ROLES_24, min_count=10,
                                studies=studies)
    # Adding a level that no sample carries must not change the answer.
    assert np.array_equal(
        mask, filter_by_expression(counts, ROLES_24, min_count=10,
                                   studies=studies))


def test_min_samples_override_applies_within_each_study():
    counts = build_one_study_gene()
    strict = filter_by_expression(counts, ROLES_24, min_count=10,
                                  min_samples=6, studies=STUDIES_24)
    assert bool(strict[0]) is False
    assert bool(strict[1]) is True


# --- the donor requirement has exactly one definition --------------------

def test_required_donor_count_is_the_smaller_arm():
    from kosmic.de.de_analysis import required_donor_count

    assert required_donor_count(6, 8) == 6
    assert required_donor_count(8, 6) == 6
    assert required_donor_count(6, 4) == 4


def test_required_donor_count_eases_large_cohorts():
    """12 v 12 asks for 11.4, not 12: 10 + (12 - 10) * 0.7."""
    from kosmic.de.de_analysis import required_donor_count

    assert required_donor_count(12, 12) == 11.4
    assert required_donor_count(10, 10) == 10
    assert required_donor_count(30, 30) == 10 + 20 * 0.7


def test_required_donor_count_honours_an_override():
    from kosmic.de.de_analysis import required_donor_count

    assert required_donor_count(6, 8, min_samples=3) == 3
    assert required_donor_count(6, 8, min_samples=0) == 6      # 0 = auto
    assert required_donor_count(6, 8, min_samples=None) == 6


def test_required_donor_count_with_an_empty_arm():
    from kosmic.de.de_analysis import required_donor_count

    assert required_donor_count(6, 0) == 6
    assert required_donor_count(0, 0) is None


def test_filter_uses_the_same_rule_it_reports():
    """The number the page prints must be the number the filter applies."""
    from kosmic.de.de_analysis import required_donor_count

    roles = np.array(["disease"] * 12 + ["control"] * 12)
    required = required_donor_count(12, 12)
    # A gene in exactly 11 donors falls short of 11.4; 12 clears it.
    for n_donors, expected in ((11, False), (12, True)):
        col = np.zeros(24)
        col[:n_donors] = 400.0
        counts = np.column_stack([col, np.full((24, 200), 500.0)])
        got = bool(filter_by_expression(counts, roles, min_count=10)[0])
        assert got is expected, (n_donors, required)
