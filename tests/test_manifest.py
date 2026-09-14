# Sentinels for the study manifest (ADR-001): write/read round-trip,
# tolerant absence/corruption, and disk-derived overview merging.
import json

from kosmic.manifest import (
    manifest_path, read_manifest, study_overview, write_manifest,
)
from kosmic.paths import scaffold_study


def test_roundtrip(tmp_path):
    study = tmp_path / "GSE1_EC"
    scaffold_study(study)
    write_manifest(
        study,
        dataset={"file": "GSE1_EC.h5ad", "n_cells": 100, "n_genes": 2000},
        semantics={
            "sample_column": "sample",
            "condition_column": "condition",
            "cell_type_column": None,  # dropped
            "role_map": {"Donor": "control", "DCM": "disease"},
        },
    )
    m = read_manifest(study)
    assert m["version"] == 1
    assert m["dataset"]["n_cells"] == 100
    assert m["semantics"]["role_map"]["DCM"] == "disease"
    assert "cell_type_column" not in m["semantics"]


def test_absent_and_corrupt_tolerated(tmp_path):
    assert read_manifest(tmp_path) is None
    manifest_path(tmp_path).write_text("{not json", encoding="utf-8")
    assert read_manifest(tmp_path) is None
    manifest_path(tmp_path).write_text(json.dumps([1, 2]), encoding="utf-8")
    assert read_manifest(tmp_path) is None  # wrong shape -> unconfigured


def test_overview_peeks_h5ad_when_manifest_absent(tmp_path):
    import anndata as ad
    import numpy as np
    import pandas as pd

    study = tmp_path / "GSE3_FB"
    scaffold_study(study)
    adata = ad.AnnData(
        X=np.zeros((4, 3)),
        obs=pd.DataFrame({
            'sample': ['s1', 's1', 's2', 's2'],
            '_role': pd.Categorical(
                ['control', 'control', 'disease', 'disease'],
                categories=['control', 'disease', 'exclude']),
        }, index=[f"c{i}" for i in range(4)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(3)]))
    adata.uns['role_map'] = {'Donor': 'control', 'DCM': 'disease'}
    adata.uns['role_condition_col'] = 'condition'
    adata.write_h5ad(study / "processed_data" / "GSE3_FB.h5ad")

    overview = study_overview(study)
    assert overview["dataset"]["n_cells"] == 4
    assert overview["dataset"]["n_genes"] == 3
    assert overview["dataset"]["file"] == "GSE3_FB.h5ad"
    # Semantics peeked from the h5ad header without a manifest
    sem = overview["semantics"]
    assert sem["sample_column"] == "sample"
    assert sem["condition_column"] == "condition"
    assert sem["role_map"] == {'Donor': 'control', 'DCM': 'disease'}


def test_overview_merges_disk_and_manifest(tmp_path):
    study = tmp_path / "GSE2_CM"
    scaffold_study(study)
    # A processed h5ad on disk flips the derived flag with no manifest
    (study / "processed_data" / "GSE2_CM.h5ad").write_bytes(b"x")

    before = study_overview(study)
    assert before["accession"] == "GSE2_CM"
    assert before["status"]["processed"] is True
    assert before["semantics"] is None  # unconfigured, not an error

    write_manifest(study, semantics={"sample_column": "donor"})
    after = study_overview(study)
    assert after["semantics"]["sample_column"] == "donor"
    # Asset flags stay disk-derived, never cached in the manifest
    assert after["status"]["raw_data"] is False
