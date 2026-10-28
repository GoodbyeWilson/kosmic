"""The methods text states each pooling method in words (#61), and the
HKSJ t-test's smallest reachable p-value is known for a study count (#63).
"""
from __future__ import annotations

import pytest

from kosmic.meta_analysis.pooling_dispatch import (
    describe_pooling_method, hksj_smallest_p,
)


def test_sumrank_top50_is_the_top_r_sum_not_a_gene_filter():
    text = describe_pooling_method('sumrank_top50', 5)
    assert text == ("SumRank over the best ceil(k/2) study ranks per gene "
                    "(top-r sum, r = 3 of k = 5)")


def test_effect_size_top50_and_hksj_are_both_stated():
    text = describe_pooling_method('dl_top50_hksj', 5)
    assert text.startswith('DerSimonian-Laird random effects')
    assert 'most precise studies' in text and 'Hartung-Knapp' in text


@pytest.mark.parametrize('key', ['sumrank', 'gwop', 'fisher', 'stouffer_weighted', 'reml'])
def test_plain_methods_have_a_name(key):
    text = describe_pooling_method(key)
    assert text != key and 'ceil' not in text and 'Hartung' not in text


def test_hksj_with_five_studies_cannot_reach_genome_wide_fdr():
    """BH over 11,267 genes needs the strongest p below 0.05 / 11,267."""
    assert hksj_smallest_p(5) > 0.05 / 11_267
    assert hksj_smallest_p(20) < 0.05 / 11_267
    assert hksj_smallest_p(1) == 1.0
