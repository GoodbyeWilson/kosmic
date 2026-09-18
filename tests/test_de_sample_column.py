"""The sample column must identify donors, and a column with hundreds of
donors must be offered as one (kosmic.de.de_analysis.check_sample_column,
DE Setup candidates)."""
import numpy as np
import pandas as pd
import pytest

from kosmic.de.de_analysis import check_sample_column


def _obs(n_donors, per=10):
    donors = [f"d{i}" for i in range(n_donors) for _ in range(per)]
    cond = ["DCM" if int(d[1:]) % 2 else "NF" for d in donors]
    return pd.DataFrame({"sample": donors, "condition": cond,
                         "_role": ["disease" if c == "DCM" else "control" for c in cond]})


def test_condition_like_column_used_as_sample_is_refused():
    obs = _obs(20)
    obs["group"] = obs["condition"]          # a second column with one value per condition
    with pytest.raises(ValueError, match="identifies conditions, not donors"):
        check_sample_column(obs, "group", "condition")


def test_same_column_for_both_is_refused():
    obs = _obs(20)
    with pytest.raises(ValueError, match="both 'condition'"):
        check_sample_column(obs, "condition", "condition")


def test_six_condition_strings_as_samples_is_refused():
    """The DCM atlas case: 141 donors, six condition strings, and the
    condition column chosen as the sample column."""
    obs = _obs(141)
    obs["condition"] = np.random.default_rng(0).choice(
        ["DCM", "Donor", "NICM", "NF", "dilated cardiomyopathy", "normal"], size=len(obs))
    obs["group"] = obs["condition"]
    with pytest.raises(ValueError, match="6"):
        check_sample_column(obs, "group", "condition")


def test_real_sample_column_passes():
    check_sample_column(_obs(141), "sample", "condition")
    check_sample_column(_obs(4), "sample", "condition")


def test_setup_offers_a_sample_column_with_hundreds_of_donors():
    import anndata as ad
    from kosmic.gui.de_analysis.pages.setup_page import SetupPage as DESetupPage
    obs = _obs(300)
    a = ad.AnnData(X=np.ones((len(obs), 3), dtype=np.float32), obs=obs)
    a.obs["sample"] = a.obs["sample"].astype("category")
    a.obs["condition"] = a.obs["condition"].astype("category")
    cands = DESetupPage._sample_column_candidates(a)
    assert "sample" in cands
    assert "condition" in cands           # one arm per value, so it is not excluded by the role rule alone


def test_pseudobulk_sparse_aggregation_matches_the_dense_loop():
    """One indicator-matrix product replaces a per-donor .toarray() loop;
    sums and means must be identical, donors in the same order."""
    import anndata as ad
    import scipy.sparse as sp
    from kosmic.de.de_analysis import create_pseudobulk
    rng = np.random.default_rng(0)
    n, g = 400, 50
    obs = _obs(20, per=20)
    X = sp.csr_matrix(rng.poisson(1.5, size=(n, g)).astype(np.float32))
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"G{j}" for j in range(g)]))
    genes = [f"G{j}" for j in range(0, g, 3)]
    pb, sdf, used = create_pseudobulk(a, genes, "sample", "condition", min_cells=1, aggregate="sum")
    dense = X[:, [int(x[1:]) for x in used]].toarray()
    expected = np.vstack([dense[(obs["sample"] == s).to_numpy()].sum(axis=0) for s in sdf["sample"]])
    assert np.allclose(pb, expected)
    assert (sdf["n_cells"] == 20).all()
    pb_mean, _, _ = create_pseudobulk(a, genes, "sample", "condition", min_cells=1, aggregate="mean")
    assert np.allclose(pb_mean, expected / 20)
