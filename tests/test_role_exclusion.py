"""Honouring '_role == exclude' upstream of differential expression.

'exclude' began as a DE-only idea, so an excluded arm still drove HVG
selection, PCA, the neighbour graph, the clusters and the UMAP. On
GSE292067 that was 44% of the cells -- a doxorubicin cohort shaping the
embedding behind a DCM-vs-donor question. These cover the mask that opts
those cells out and the scatter that writes a subset analysis back without
deleting them.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from kosmic.scrna.inspect.roles import included_mask, scatter_results


def _adata(roles):
    n = len(roles)
    return ad.AnnData(
        X=np.zeros((n, 3), dtype=np.float32),
        var=pd.DataFrame(index=list("ABC")),
        obs=pd.DataFrame({"_role": roles},
                         index=[f"c{i}" for i in range(n)]))


# --- the mask ---------------------------------------------------------------

def test_no_role_column_means_nothing_to_honour():
    a = ad.AnnData(X=np.zeros((4, 3), dtype=np.float32))
    assert included_mask(a) is None


def test_no_excluded_cells_means_nothing_to_honour():
    """Callers can then keep the fast path instead of subsetting."""
    assert included_mask(_adata(["disease", "control"])) is None


def test_excluded_cells_are_masked_out():
    mask = included_mask(_adata(["disease", "exclude", "control", "exclude"]))
    assert list(mask) == [True, False, True, False]


def test_case_is_ignored():
    assert list(included_mask(_adata(["Disease", "EXCLUDE"]))) == [True, False]


def test_everything_excluded_is_refused():
    """Analysing nothing is never what was meant."""
    assert included_mask(_adata(["exclude", "exclude"])) is None


# --- scattering results back ------------------------------------------------

def _analysed_subset(full, mask):
    sub = full[mask].copy()
    sub.obsm["X_pca"] = np.arange(sub.n_obs * 2, dtype=np.float32).reshape(-1, 2)
    sub.obs["leiden"] = pd.Categorical(["0"] * sub.n_obs)
    sub.uns["neighbors"] = {"params": {"n_neighbors": 15}}
    return sub


def test_cell_count_is_unchanged():
    """Excluded cells stay in the file -- deleting them is Filter Dataset."""
    full = _adata(["disease", "exclude", "control"])
    mask = included_mask(full)
    scatter_results(full, _analysed_subset(full, mask), mask,
                    obsm_keys=("X_pca",), obs_cols=("leiden",))
    assert full.n_obs == 3


def test_included_rows_get_their_values():
    full = _adata(["disease", "exclude", "control"])
    mask = included_mask(full)
    sub = _analysed_subset(full, mask)
    scatter_results(full, sub, mask, obsm_keys=("X_pca",))
    assert np.array_equal(full.obsm["X_pca"][mask], sub.obsm["X_pca"])


def test_excluded_rows_get_nan_not_zero():
    """Zero is a coordinate; NaN says the analysis never saw this cell."""
    full = _adata(["disease", "exclude", "control"])
    mask = included_mask(full)
    scatter_results(full, _analysed_subset(full, mask), mask,
                    obsm_keys=("X_pca",))
    assert np.isnan(full.obsm["X_pca"][1]).all()


def test_excluded_rows_get_no_cluster_label():
    full = _adata(["disease", "exclude", "control"])
    mask = included_mask(full)
    scatter_results(full, _analysed_subset(full, mask), mask,
                    obs_cols=("leiden",))
    assert pd.isna(full.obs["leiden"].iloc[1])
    assert full.obs["leiden"].iloc[0] == "0"


def test_a_stale_graph_is_dropped_not_left_behind():
    """obsp is cell-by-cell; off the subset it is meaningless."""
    full = _adata(["disease", "exclude", "control"])
    full.obsp["connectivities"] = np.eye(3, dtype=np.float32)
    mask = included_mask(full)
    scatter_results(full, _analysed_subset(full, mask), mask,
                    obsp_keys=("connectivities",))
    assert "connectivities" not in full.obsp


def test_uns_is_carried_across():
    full = _adata(["disease", "exclude", "control"])
    mask = included_mask(full)
    scatter_results(full, _analysed_subset(full, mask), mask,
                    uns_keys=("neighbors",))
    assert full.uns["neighbors"]["params"]["n_neighbors"] == 15


@pytest.mark.parametrize("key", ["X_umap", "not_present"])
def test_a_missing_key_is_skipped_quietly(key):
    full = _adata(["disease", "exclude", "control"])
    mask = included_mask(full)
    scatter_results(full, _analysed_subset(full, mask), mask,
                    obsm_keys=(key,))
    assert key not in full.obsm


# --- cluster labels after exclusion -----------------------------------------

from kosmic.scrna.inspect.roles import cluster_labels     # noqa: E402


def _clustered(labels):
    n = len(labels)
    return ad.AnnData(
        X=np.zeros((n, 2), dtype=np.float32),
        obs=pd.DataFrame({"leiden": pd.Categorical(labels)},
                         index=[f"c{i}" for i in range(n)]))


def test_nan_clusters_are_not_returned():
    """The crash: NaN is in unique(), but `series == NaN` matches nothing,
    so the caller got an empty group and value_counts().index[0] raised."""
    a = _clustered(["0", "1", None, None])
    assert cluster_labels(a) == ["0", "1"]


def test_every_returned_label_has_cells():
    a = _clustered(["0", "1", None])
    for label in cluster_labels(a):
        assert int((a.obs["leiden"] == label).sum()) > 0


def test_an_emptied_category_is_dropped():
    """A categorical keeps categories after their last cell is filtered
    away; those fail exactly like NaN does."""
    a = _clustered(["0", "1"])
    a.obs["leiden"] = a.obs["leiden"].cat.add_categories(["9"])
    assert "9" not in cluster_labels(a)


def test_numeric_labels_sort_numerically():
    a = _clustered(["10", "2", "1"])
    assert cluster_labels(a) == ["1", "2", "10"]


def test_mixed_labels_put_numbers_first():
    a = _clustered(["2", "Fibroblast", "1"])
    assert cluster_labels(a) == ["1", "2", "Fibroblast"]


def test_a_missing_column_is_empty_not_an_error():
    assert cluster_labels(_clustered(["0"]), "not_a_column") == []


def test_all_cells_unclustered_gives_nothing():
    a = _clustered([None, None])
    assert cluster_labels(a) == []


# --- labels are not always strings ------------------------------------------

from kosmic.scrna.annotate.score import collapse_celltypist_label  # noqa: E402


def test_collapse_passes_a_real_label_through_the_map():
    assert collapse_celltypist_label("Tcm/Naive helper T cells").startswith(
        "T") is not None


def test_collapse_survives_nan():
    """Cells left out of the embedding carry no label, and
    Categorical.map calls the mapper with np.nan to decide what NA
    becomes -- so this is reached with a float even when every real
    label is a string."""
    assert np.isnan(collapse_celltypist_label(np.nan))


def test_collapse_survives_none():
    assert collapse_celltypist_label(None) is None


def test_mapping_a_categorical_with_nulls_does_not_raise():
    labels = pd.Categorical(["Fibroblasts", "Endothelial cells", None])
    out = pd.Series(labels).map(collapse_celltypist_label)
    assert out.isna().sum() == 1
    assert out.notna().sum() == 2


# --- plots must show the same donors the contrast used -------------------

def _adata_with_excluded_arm():
    """3 disease, 3 control, 2 excluded donors."""
    import anndata
    rows = []
    for role, donors, cond in (("disease", ["D1", "D2", "D3"], "DCM"),
                               ("control", ["C1", "C2", "C3"], "Donor"),
                               ("exclude", ["X1", "X2"], "DoxCM")):
        for d in donors:
            rows += [(d, role, cond)] * 40
    obs = pd.DataFrame(rows, columns=["sample", "_role", "condition"])
    obs.index = [f"c{i}" for i in range(len(obs))]
    rng = np.random.default_rng(0)
    X = rng.poisson(8, size=(len(obs), 30)).astype(np.float32)
    return anndata.AnnData(
        X=X, obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(30)]))


def test_expression_matrix_drops_the_excluded_arm():
    """Plots sit beside DE results; they must not add donors the DE never saw."""
    from kosmic.de.de_analysis import pseudobulk_expression_matrix

    adata = _adata_with_excluded_arm()
    df = pseudobulk_expression_matrix(
        adata, "sample", "condition", normalization="cpm", min_cells=10)

    assert len(df) == 6
    assert set(df["sample"]) == {"D1", "D2", "D3", "C1", "C2", "C3"}
    assert set(df["condition"].astype(str)) == {"DCM", "Donor"}
    assert "DoxCM" not in set(df["condition"].astype(str))


def test_expression_matrix_without_roles_keeps_everyone():
    """No _role column means nothing has been excluded."""
    from kosmic.de.de_analysis import pseudobulk_expression_matrix

    adata = _adata_with_excluded_arm()
    del adata.obs["_role"]
    df = pseudobulk_expression_matrix(
        adata, "sample", "condition", normalization="cpm", min_cells=10)
    assert len(df) == 8
