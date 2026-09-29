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


def test_setup_page_combo_offers_a_sample_column_with_hundreds_of_donors(tmp_path):
    """Both the candidate rule and the combo loop must admit it."""
    import os
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    import anndata as ad
    from PyQt6.QtCore import QEvent
    from PyQt6.QtWidgets import QApplication
    from kosmic.gui.de_analysis.workspace import DEWorkspace
    app = QApplication.instance() or QApplication([])
    obs = _obs(141)
    a = ad.AnnData(X=np.ones((len(obs), 3), dtype=np.float32), obs=obs)
    for c in ("sample", "condition", "_role"):
        a.obs[c] = a.obs[c].astype("category")
    ws = DEWorkspace()
    ws.current_adata = a
    ws.setup_page._refresh_from_adata()
    items = [ws.setup_page.sample_col_combo.itemText(i) for i in range(ws.setup_page.sample_col_combo.count())]
    assert "sample" in items
    assert ws.sample_col == "sample"
    ws.close()
    ws.deleteLater()
    app.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    app.processEvents()


def test_pseudobulk_takes_genes_in_the_requested_order_without_copying_the_matrix():
    """The requested genes are picked out of the per-sample sums, so the
    columns follow the request (not the matrix), genes absent from the
    matrix are dropped, and the count matrix itself is untouched."""
    import anndata as ad
    import scipy.sparse as sp
    from kosmic.de.de_analysis import create_pseudobulk
    rng = np.random.default_rng(1)
    n, g = 200, 30
    obs = _obs(10, per=20)
    X = sp.csr_matrix(rng.poisson(1.0, size=(n, g)).astype(np.float32))
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"G{j}" for j in range(g)]))
    genes = ["G7", "G2", "NOT_A_GENE", "G29", "G0"]
    pb, sdf, used = create_pseudobulk(a, genes, "sample", "condition", min_cells=1, aggregate="sum")
    assert used == ["G7", "G2", "G29", "G0"]
    dense = X.toarray()
    for j, gname in enumerate(used):
        col = int(gname[1:])
        expected = [dense[(obs["sample"] == s).to_numpy(), col].sum() for s in sdf["sample"]]
        assert np.allclose(pb[:, j], expected)
    assert a.X is X                                     # no copy was made or swapped in


def test_pct_expressing_per_var_is_chunked_and_exact():
    """Row-chunked non-zero counts per gene equal the dense computation,
    across chunk boundaries and with explicit zeros stored in the matrix."""
    import anndata as ad
    import scipy.sparse as sp
    from kosmic.de.de_analysis import annotate_pct_expressing, pct_expressing_per_var
    rng = np.random.default_rng(2)
    n, g = 333, 12
    obs = _obs(37, per=9)
    dense = rng.poisson(0.7, size=(n, g)).astype(np.float32)
    X = sp.csr_matrix(dense)
    X.data[:5] = 0.0                                    # stored zeros must not count
    dense = X.toarray()
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"G{j}" for j in range(g)]))
    pd_, pc_ = pct_expressing_per_var(a, chunk_cells=50)
    is_d = (obs["_role"] == "disease").to_numpy()
    assert np.allclose(pd_, (dense[is_d] != 0).mean(axis=0))
    assert np.allclose(pc_, (dense[~is_d] != 0).mean(axis=0))
    de = pd.DataFrame({"names": ["G3", "G0", "MISSING"]})
    out = annotate_pct_expressing(de, a, (pd_, pc_))
    assert np.isclose(out.pct_disease[0], pd_[3]) and np.isclose(out.pct_control[1], pc_[0])
    assert np.isnan(out.pct_disease[2])
    dense_a = ad.AnnData(X=dense.copy(), obs=obs, var=a.var.copy())
    pdd, pcd = pct_expressing_per_var(dense_a, chunk_cells=100)
    assert np.allclose(pdd, pd_) and np.allclose(pcd, pc_)
