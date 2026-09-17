"""Writing master labels back onto each study.

A study can be short of the master for two reasons that want opposite
handling. A cell cap leaves an arbitrary subsample, so the absent cells
should be projected by kNN. Cells left out because their role is
'exclude' were deliberately omitted -- projecting them would invent
labels, from a space they were never part of, for cells already declared
outside the analysis.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from kosmic.combine import concat_studies, propagate_labels

GENES = ["TTN", "NPPA", "MYH7"]


def _study(path, roles):
    path.parent.mkdir(parents=True, exist_ok=True)
    n = len(roles)
    a = ad.AnnData(
        X=csr_matrix(np.ones((n, len(GENES)), dtype=np.float32)),
        var=pd.DataFrame(index=GENES),
        obs=pd.DataFrame({"_role": roles},
                         index=[f"cell{i}" for i in range(n)]))
    a.write_h5ad(path)
    return path


def _master_from(tmp_path, study_paths, **kw):
    # concat_on_disk raises on a single input file, so every fixture
    # combines at least two studies -- which is the real use case anyway.
    out = concat_studies(study_paths, tmp_path / "master.h5ad", **kw)
    m = ad.read_h5ad(out)
    # Stand in for clustering + annotation.
    m.obs["leiden"] = pd.Categorical(["0"] * m.n_obs)
    m.obs["cell_type"] = pd.Categorical(["Cardiomyocyte"] * m.n_obs)
    m.write_h5ad(out)
    return out


def test_every_cell_present_labels_all_of_them(tmp_path):
    s1 = _study(tmp_path / "S1" / "processed_data" / "S1.h5ad",
                ["disease", "control", "disease"])
    s2 = _study(tmp_path / "S2" / "processed_data" / "S2.h5ad", ["control"])
    master = _master_from(tmp_path, [s1, s2])
    propagate_labels(master, s1)
    out = ad.read_h5ad(s1)
    assert out.obs["cell_type"].notna().all()


def test_excluded_cells_stay_unlabelled(tmp_path):
    """The master omitted them on purpose; inventing labels for them
    would be worse than leaving them blank."""
    s1 = _study(tmp_path / "S1" / "processed_data" / "S1.h5ad",
                ["disease", "control", "exclude", "exclude"])
    s2 = _study(tmp_path / "S2" / "processed_data" / "S2.h5ad", ["control"])
    master = _master_from(tmp_path, [s1, s2], drop_excluded=True)
    propagate_labels(master, s1)
    out = ad.read_h5ad(s1)
    labelled = out.obs["cell_type"].notna()
    assert int(labelled.sum()) == 2
    assert not labelled[out.obs["_role"].astype(str) == "exclude"].any()


def test_the_study_keeps_all_its_cells(tmp_path):
    s1 = _study(tmp_path / "S1" / "processed_data" / "S1.h5ad",
                ["disease", "exclude"])
    s2 = _study(tmp_path / "S2" / "processed_data" / "S2.h5ad", ["control"])
    master = _master_from(tmp_path, [s1, s2], drop_excluded=True)
    propagate_labels(master, s1)
    assert ad.read_h5ad(s1).n_obs == 2


def test_two_studies_get_their_own_labels(tmp_path):
    """Barcode names repeat across studies; the join must not cross them."""
    s1 = _study(tmp_path / "S1" / "processed_data" / "S1.h5ad",
                ["disease", "exclude"])
    s2 = _study(tmp_path / "S2" / "processed_data" / "S2.h5ad",
                ["control", "control"])
    out = concat_studies([s1, s2], tmp_path / "master.h5ad",
                         drop_excluded=True)
    m = ad.read_h5ad(out)
    m.obs["leiden"] = pd.Categorical(["0"] * m.n_obs)
    m.obs["cell_type"] = pd.Categorical(
        ["TypeA"] + ["TypeB"] * (m.n_obs - 1))
    m.write_h5ad(out)

    propagate_labels(out, s1)
    propagate_labels(out, s2)
    assert ad.read_h5ad(s1).obs["cell_type"].iloc[0] == "TypeA"
    assert ad.read_h5ad(s2).obs["cell_type"].notna().all()


def test_a_missing_label_column_is_refused(tmp_path):
    s1 = _study(tmp_path / "S1" / "processed_data" / "S1.h5ad", ["disease"])
    s2 = _study(tmp_path / "S2" / "processed_data" / "S2.h5ad", ["control"])
    master = concat_studies([s1, s2], tmp_path / "master.h5ad")
    with pytest.raises(ValueError, match="missing"):
        propagate_labels(master, s1)


# --- not clobbering the study's own annotation ------------------------------

def _annotated_study(path, roles, own_type):
    """A study that was already clustered and annotated on its own."""
    path.parent.mkdir(parents=True, exist_ok=True)
    n = len(roles)
    a = ad.AnnData(
        X=csr_matrix(np.ones((n, len(GENES)), dtype=np.float32)),
        var=pd.DataFrame(index=GENES),
        obs=pd.DataFrame({
            "_role": roles,
            "leiden": pd.Categorical(["7"] * n),
            "cell_type": pd.Categorical([own_type] * n),
        }, index=[f"cell{i}" for i in range(n)]))
    a.write_h5ad(path)
    return path


def test_the_studys_own_labels_survive(tmp_path):
    """Comparing per-study DE under independent labels against the same
    DE under shared atlas labels needs both to exist."""
    s1 = _annotated_study(tmp_path / "S1" / "processed_data" / "S1.h5ad",
                          ["disease", "control"], "MyOwnCall")
    s2 = _annotated_study(tmp_path / "S2" / "processed_data" / "S2.h5ad",
                          ["control"], "MyOwnCall")
    master = _master_from(tmp_path, [s1, s2])
    propagate_labels(master, s1)
    out = ad.read_h5ad(s1)
    assert (out.obs["cell_type"] == "MyOwnCall").all()
    assert (out.obs["leiden"].astype(str) == "7").all()


def test_atlas_labels_land_in_suffixed_columns(tmp_path):
    s1 = _annotated_study(tmp_path / "S1" / "processed_data" / "S1.h5ad",
                          ["disease", "control"], "MyOwnCall")
    s2 = _annotated_study(tmp_path / "S2" / "processed_data" / "S2.h5ad",
                          ["control"], "MyOwnCall")
    master = _master_from(tmp_path, [s1, s2])
    propagate_labels(master, s1)
    out = ad.read_h5ad(s1)
    assert (out.obs["cell_type_atlas"] == "Cardiomyocyte").all()
    assert "leiden_atlas" in out.obs.columns


def test_an_unannotated_study_still_gets_a_usable_cell_type(tmp_path):
    """Nothing to preserve, so the plain column is filled too."""
    s1 = _study(tmp_path / "S1" / "processed_data" / "S1.h5ad",
                ["disease", "control"])
    s2 = _study(tmp_path / "S2" / "processed_data" / "S2.h5ad", ["control"])
    master = _master_from(tmp_path, [s1, s2])
    propagate_labels(master, s1)
    out = ad.read_h5ad(s1)
    assert (out.obs["cell_type"] == "Cardiomyocyte").all()
    assert (out.obs["cell_type_atlas"] == "Cardiomyocyte").all()


def test_a_custom_suffix_is_honoured(tmp_path):
    s1 = _annotated_study(tmp_path / "S1" / "processed_data" / "S1.h5ad",
                          ["disease"], "MyOwnCall")
    s2 = _annotated_study(tmp_path / "S2" / "processed_data" / "S2.h5ad",
                          ["control"], "MyOwnCall")
    master = _master_from(tmp_path, [s1, s2])
    propagate_labels(master, s1, suffix="_v2")
    assert "cell_type_v2" in ad.read_h5ad(s1).obs.columns


# --- obs columns that survive combining ---------------------------------------

def test_atlas_labels_sex_and_age_are_carried_onto_a_rebuilt_master(tmp_path):
    """After propagation a study holds cell_type_atlas; re-creating the
    master (e.g. uncapped, for the mega-analysis DE) must keep it, and the
    sex the Inspect step designated must arrive as one female/male column."""
    s1 = _study(tmp_path / "S1" / "processed_data" / "S1.h5ad",
                ["disease", "control", "disease"])
    s2 = _study(tmp_path / "S2" / "processed_data" / "S2.h5ad", ["control", "control"])
    a1 = ad.read_h5ad(s1)
    a1.obs["cell_type_atlas"] = pd.Categorical(["CM", "Fib", "CM"])
    a1.obs["leiden_atlas"] = pd.Categorical(["0", "1", "0"])
    a1.obs["Sex"] = ["F", "M", "F"]              # recorded under a different name
    a1.obs["age"] = [61.0, 45.0, 70.0]
    a1.uns["sex_column"] = "Sex"
    a1.write_h5ad(s1)
    a2 = ad.read_h5ad(s2)
    a2.obs["sex_inferred"] = pd.Categorical(["male", "female"])
    a2.uns["sex_column"] = "sex_inferred"
    a2.write_h5ad(s2)

    m = ad.read_h5ad(concat_studies([s1, s2], tmp_path / "master.h5ad"))
    assert list(m.obs["cell_type_atlas"].astype(str)[:3]) == ["CM", "Fib", "CM"]
    assert list(m.obs["leiden_atlas"].astype(str)[:3]) == ["0", "1", "0"]
    assert list(m.obs["sex"].astype(str)) == ["female", "male", "female", "male", "female"]
    assert list(m.obs["age"][:3]) == [61.0, 45.0, 70.0]
    assert "Sex" not in m.obs.columns and "sex_inferred" not in m.obs.columns


def test_no_sex_designation_writes_no_sex_column(tmp_path):
    s1 = _study(tmp_path / "S1" / "processed_data" / "S1.h5ad", ["disease"])
    s2 = _study(tmp_path / "S2" / "processed_data" / "S2.h5ad", ["control"])
    m = ad.read_h5ad(concat_studies([s1, s2], tmp_path / "master.h5ad"))
    assert "sex" not in m.obs.columns


def test_text_age_in_one_study_and_float_in_another_concatenate(tmp_path):
    """Koenig recorded age as a categorical of strings ('68', ''), Chaffin
    as float; mixed, the on-disk concat could not write the column."""
    s1 = _study(tmp_path / "S1" / "processed_data" / "S1.h5ad", ["disease", "control", "disease"])
    s2 = _study(tmp_path / "S2" / "processed_data" / "S2.h5ad", ["control", "control"])
    a1 = ad.read_h5ad(s1)
    a1.obs["age"] = pd.Categorical(["68", "", "53"])
    a1.write_h5ad(s1)
    a2 = ad.read_h5ad(s2)
    a2.obs["age"] = [56.0, np.nan]
    a2.write_h5ad(s2)
    m = ad.read_h5ad(concat_studies([s1, s2], tmp_path / "master.h5ad"))
    assert m.obs["age"].dtype.kind == "f"
    vals = m.obs["age"].tolist()
    assert vals[0] == 68.0 and np.isnan(vals[1]) and vals[2] == 53.0 and vals[3] == 56.0 and np.isnan(vals[4])
