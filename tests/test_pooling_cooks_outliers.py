"""A gene a study's DE gave no raw p-value (a DESeq2 Cook's outlier: fold
change reported, test withheld) is missing from that study for every pooling
method, so all methods pool the same genes, as the permutation null does."""
import numpy as np
import pandas as pd
import pytest

from kosmic.meta_analysis.pooling_dispatch import get_analytical_pool_fn


def _study(i, cooks=False):
    rng = np.random.default_rng(i)
    names = [f"G{j}" for j in range(30)]
    p = rng.uniform(0.001, 0.9, 30)
    if cooks:
        p[0] = np.nan                       # G0: Cook's outlier in this study
    return pd.DataFrame({"names": names, "logfoldchanges": rng.normal(0, 1, 30),
                         "se": rng.uniform(0.1, 0.5, 30), "pvals": p,
                         "pvals_adj": np.minimum(p * 3, 1), "dataset": f"S{i}"})


STUDIES = [_study(i, cooks=(i == 2)) for i in range(4)]


@pytest.mark.parametrize("method", ["reml", "dl", "stouffer", "sumrank", "fisher"])
def test_a_cooks_outlier_counts_as_missing_for_every_method(method):
    out = get_analytical_pool_fn(method, 4)(STUDIES, min_studies=4).set_index("names")
    assert "G0" not in out.index                 # in only 3 of the 4 studies
    assert len(out) == 29
    three = get_analytical_pool_fn(method, 3)(STUDIES, min_studies=3).set_index("names")
    assert three.at["G0", "n_studies"] == 3


def test_a_missing_adjusted_p_value_alone_does_not_drop_the_gene():
    """DESeq2's independent filter leaves padj NaN with a valid raw test."""
    studies = [d.assign(pvals_adj=np.where(d.names == "G1", np.nan, d.pvals_adj))
               for d in STUDIES]
    out = get_analytical_pool_fn("reml", 4)(studies, min_studies=4).set_index("names")
    assert "G1" in out.index
