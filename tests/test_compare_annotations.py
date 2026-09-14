"""Independent per-study annotation against atlas-propagated labels.

Both survive propagation on purpose: agreement between them is a result,
and disagreement is something to understand before building on either.
"""
import anndata as ad
import numpy as np
import pandas as pd

from kosmic.scrna.inspect.gene_group import compare_annotations


def _adata(own, atlas):
    n = len(own)
    return ad.AnnData(
        X=np.zeros((n, 2), dtype=np.float32),
        obs=pd.DataFrame({"cell_type": pd.Categorical(own),
                          "cell_type_atlas": pd.Categorical(atlas)},
                         index=[f"c{i}" for i in range(n)]))


def test_identical_labels_agree_completely():
    r = compare_annotations(_adata(["CM", "Fib"], ["CM", "Fib"]))
    assert r["agreement"] == 1.0
    assert r["n_compared"] == 2


def test_a_disagreement_is_counted_and_named():
    r = compare_annotations(_adata(["CM", "Fib"], ["CM", "Pericyte"]))
    assert r["agreement"] == 0.5
    assert r["disagreements"][0] == ("Fib", "Pericyte", 1)


def test_disagreements_come_worst_first():
    r = compare_annotations(
        _adata(["A", "B", "B", "B"], ["X", "Y", "Y", "Y"]))
    assert r["disagreements"][0][2] == 3


def test_cells_the_atlas_never_saw_are_not_disagreements():
    """An excluded cell has no atlas label. Counting it as a mismatch
    would understate agreement in proportion to how much was excluded."""
    r = compare_annotations(_adata(["CM", "CM", "CM"], ["CM", None, None]))
    assert r["n_compared"] == 1
    assert r["agreement"] == 1.0
    assert r["n_unlabelled"] == 2


def test_a_missing_column_is_empty_not_an_error():
    a = _adata(["CM"], ["CM"])
    del a.obs["cell_type_atlas"]
    r = compare_annotations(a)
    assert r["agreement"] is None
    assert r["n_compared"] == 0


def test_nothing_comparable_reports_none_not_zero():
    """Zero agreement and 'no overlap to judge' are different findings."""
    r = compare_annotations(_adata(["CM", "Fib"], [None, None]))
    assert r["agreement"] is None
    assert r["n_unlabelled"] == 2


def test_the_crosstab_is_left_by_right():
    r = compare_annotations(_adata(["CM", "Fib"], ["CM", "CM"]))
    assert r["crosstab"].loc["Fib", "CM"] == 1


def test_relabelling_shows_as_a_clean_block():
    """A pure rename ('Fib' -> 'Fibroblast') is 0% agreement but one
    disagreement, not scattered noise -- worth telling apart."""
    r = compare_annotations(
        _adata(["Fib"] * 5, ["Fibroblast"] * 5))
    assert r["agreement"] == 0.0
    assert len(r["disagreements"]) == 1
