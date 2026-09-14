"""Sentinel tests for the project / study scaffolding contract.

`scaffold_project` and `scaffold_study` are how KOSMIC creates the
canonical on-disk layout when a user picks a folder or adds a study.
If they break, the whole add-flow breaks silently.
"""

from kosmic.paths import (
    de_pathway_scoring_dir,
    de_stats_dir,
    de_figures_dir,
    de_result_path,
    list_studies,
    meta_output_dir,
    processed_data_dir,
    pseudobulk_path,
    raw_data_dir,
    scaffold_project,
    scaffold_study,
    scrna_figures_dir,
    study_status,
)


def test_scaffold_project_creates_meta_analysis(tmp_path):
    """A project root scaffolds with meta_analysis/ ready for outputs."""
    project = tmp_path / "MyProject"
    scaffold_project(project)

    assert project.is_dir()
    assert meta_output_dir(project).is_dir()


def test_scaffold_project_is_idempotent(tmp_path):
    """Calling twice on an existing project does not raise."""
    project = tmp_path / "MyProject"
    scaffold_project(project)
    scaffold_project(project)


def test_scaffold_study_creates_all_canonical_subdirs(tmp_path):
    """A study folder scaffolds with every directory the workspaces need."""
    study = tmp_path / "GSE12345"
    scaffold_study(study)

    assert study.is_dir()
    assert raw_data_dir(study).is_dir()
    assert processed_data_dir(study).is_dir()
    assert scrna_figures_dir(study).is_dir()
    assert de_figures_dir(study).is_dir()
    assert de_stats_dir(study).is_dir()
    assert de_pathway_scoring_dir(study).is_dir()


def test_scaffold_study_preserves_existing_files(tmp_path):
    """Scaffolding an already-populated study folder leaves data untouched."""
    study = tmp_path / "GSE12345"
    scaffold_study(study)
    sentinel = raw_data_dir(study) / "existing.txt"
    sentinel.write_text("preserved")

    scaffold_study(study)

    assert sentinel.read_text() == "preserved"


def test_list_studies_excludes_meta_analysis_and_dotfiles(tmp_path):
    """Study enumeration sees real studies, hides reserved + hidden folders."""
    project = tmp_path / "Project"
    scaffold_project(project)
    scaffold_study(project / "GSE12345")
    scaffold_study(project / "GSE99999")
    (project / ".cache").mkdir()
    (project / "notes.txt").write_text("scratch")

    assert list_studies(project) == ["GSE12345", "GSE99999"]


def test_list_studies_returns_empty_for_missing_or_empty_project(tmp_path):
    """Missing / empty project directories report no studies, no error."""
    assert list_studies(tmp_path / "does-not-exist") == []

    project = tmp_path / "Empty"
    scaffold_project(project)
    assert list_studies(project) == []


def test_list_studies_rejects_study_folder_opened_as_project(tmp_path):
    """If user picks a study folder by mistake, its substructure isn't listed.

    Without this guard, opening ``GSE12345/`` (a study) as a project would
    report ``raw_data``, ``processed_data``, ``results`` as three studies.
    """
    study = tmp_path / "GSE12345"
    scaffold_study(study)

    assert list_studies(study) == []


def test_study_status_pseudobulk_handles_celltype_suffixed_accession(tmp_path):
    """Pseudobulk badge picks up files where the accession includes a suffix.

    Real KOSMIC projects sometimes use ``{folder}_{celltype}`` as the
    accession in filenames (e.g. ``GSE12345_DCMEC_pseudobulk.csv``)
    even though the folder name is just ``GSE12345``.
    """
    study = tmp_path / "GSE12345"
    scaffold_study(study)
    pseudobulk_path(study, "GSE12345_DCMEC").write_text("gene\n")

    status = study_status(study)
    assert status["pseudobulk"] is True


def test_study_status_reflects_actual_files(tmp_path):
    """Status flags flip true only when the underlying files exist."""
    study = tmp_path / "GSE12345"
    scaffold_study(study)

    initial = study_status(study)
    assert initial == {
        "raw_data": False,
        "processed": False,
        "de_results": False,
        "pseudobulk": False,
    }

    (raw_data_dir(study) / "GSE12345_RAW.tar").write_bytes(b"x")
    (study / "processed_data" / "clustered.h5ad").write_bytes(b"x")
    de_result_path(study, "GSE12345", "deseq2_fast").write_text("gene\n")
    pseudobulk_path(study, "GSE12345").write_text("gene\n")

    assert study_status(study) == {
        "raw_data": True,
        "processed": True,
        "de_results": True,
        "pseudobulk": True,
    }
