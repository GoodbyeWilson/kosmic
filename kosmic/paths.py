"""Canonical project/study/meta path helpers for KOSMIC.

KOSMIC uses a consistent on-disk layout that's spread across three
workspaces (scRNA / DE / Meta).  The helpers here centralise it so
every caller builds paths the same way.

Layout
------

::

    project_dir/                     <- Meta workspace root
        meta_analysis/               <- consensus MA outputs land here
        {accession_1}/               <- one study; a "study_dir"
            raw_data/                <- downloaded source files
            processed_data/          <- pipeline intermediate h5ad files
            figures/                 <- scRNA figure outputs
            results/
                figures/             <- DE workspace figure exports
                metabolic_DE_results_all.csv
                metabolic_DE_results_significant.csv
                de_analysis/
                    statistics/
                        {accession_1}_pseudobulk.csv
                        {accession_1}_DE_deseq2_fast.csv
                        {accession_1}_DE_welch_cpm.csv
                        ...
                    pathway_scoring/
                        pathway_de_results.csv
        {accession_2}/
            ...

The DE workspace operates on a single study at a time -- its
'project_dir' is one of the '{accession}' folders (a "study_dir").
The Meta workspace's ``project_dir`` is the *parent* of many study_dirs.
These two concepts have the same variable name in callers today; the
helpers disambiguate via typed function names (``study_dir`` argument
for per-study paths; ``project_dir`` for the meta-level root).

Filename conventions
--------------------

DE-result CSVs follow ``{accession}_DE_{method}.csv`` where ``method``
is one of ``deseq2``, ``deseq2_fast``, ``welch_cpm``, ``welch_cpm_eb``,
``welch_raw``, ``welch_ttest``, ``welch_ttest_eb``.  The regex that
parses this lives in :mod:`kosmic.meta_analysis.io`; helpers here
produce paths that match it.
"""
from __future__ import annotations

from pathlib import Path
from typing import Union

PathLike = Union[str, Path]


# ---------------------------------------------------------------------------
# Per-study (DE workspace / scRNA workspace)
# ---------------------------------------------------------------------------

def raw_data_dir(study_dir: PathLike) -> Path:
    """``{study_dir}/raw_data`` -- downloaded source files."""
    return Path(study_dir) / "raw_data"


def processed_data_dir(study_dir: PathLike) -> Path:
    """``{study_dir}/processed_data`` -- pipeline intermediate h5ad files."""
    return Path(study_dir) / "processed_data"


def scrna_figures_dir(study_dir: PathLike) -> Path:
    """``{study_dir}/figures`` -- scRNA figure outputs.

    Distinct from :func:`de_figures_dir` (the DE workspace exports
    under ``results/figures/``).
    """
    return Path(study_dir) / "figures"


def results_dir(study_dir: PathLike) -> Path:
    """``{study_dir}/results`` -- parent of DE + figures + metabolic CSVs."""
    return Path(study_dir) / "results"


def de_analysis_dir(study_dir: PathLike) -> Path:
    """``{study_dir}/results/de_analysis`` -- DE analysis root."""
    return results_dir(study_dir) / "de_analysis"


def de_stats_dir(study_dir: PathLike) -> Path:
    """``{study_dir}/results/de_analysis/statistics`` -- per-method CSVs."""
    return de_analysis_dir(study_dir) / "statistics"


def de_pathway_scoring_dir(study_dir: PathLike) -> Path:
    """``{study_dir}/results/de_analysis/pathway_scoring`` -- pathway DE CSVs."""
    return de_analysis_dir(study_dir) / "pathway_scoring"


def de_result_path(study_dir: PathLike, accession: str, method: str) -> Path:
    """Path to one DE-method result CSV for a study.

    Produces ``{study_dir}/results/de_analysis/statistics/{accession}_DE_{method}.csv``.
    ``method`` must match the regex in
    :data:`kosmic.meta_analysis.io._DE_RESULTS_PATTERN` for the
    Meta workspace to discover it.
    """
    return de_stats_dir(study_dir) / f"{accession}_DE_{method}.csv"


def significant_path(study_dir, accession):
    """``{study_dir}/results/de_analysis/statistics/{accession}_significant_genes.csv``.

    Deliberately not named ``*_DE_*``: the meta-analysis discovers study
    inputs by that pattern, and a significant-only table is a summary for
    reading, not a pooling input.
    """
    return de_stats_dir(study_dir) / f"{accession}_significant_genes.csv"


def pseudobulk_path(study_dir: PathLike, accession: str) -> Path:
    """Path to the per-study pseudobulk counts CSV.

    Produces ``{study_dir}/results/de_analysis/statistics/{accession}_pseudobulk.csv``.
    """
    return de_stats_dir(study_dir) / f"{accession}_pseudobulk.csv"


def de_figures_dir(study_dir: PathLike) -> Path:
    """``{study_dir}/results/figures`` -- DE workspace figure exports.

    Distinct from :func:`scrna_figures_dir` (``{study_dir}/figures``).
    """
    return results_dir(study_dir) / "figures"


def metabolic_de_all_path(study_dir: PathLike) -> Path:
    """``{study_dir}/results/metabolic_DE_results_all.csv`` -- default export path."""
    return results_dir(study_dir) / "metabolic_DE_results_all.csv"


def metabolic_de_significant_path(study_dir: PathLike) -> Path:
    """``{study_dir}/results/metabolic_DE_results_significant.csv`` -- default export path."""
    return results_dir(study_dir) / "metabolic_DE_results_significant.csv"


def processed_h5ad_path(study_dir: PathLike, name: str) -> Path:
    """``{study_dir}/processed_data/{name}`` -- canonical h5ad slot.

    ``name`` is typically ``"annotated.h5ad"`` / ``"clustered.h5ad"`` /
    ``"pca.h5ad"`` for the three scRNA pipeline checkpoints.
    """
    return processed_data_dir(study_dir) / name


# ---------------------------------------------------------------------------
# Project-level (Meta workspace)
# ---------------------------------------------------------------------------

def study_dir(project_dir: PathLike, accession: str) -> Path:
    """``{project_dir}/{accession}`` -- one study inside a meta project."""
    return Path(project_dir) / accession


def meta_output_dir(project_dir: PathLike) -> Path:
    """``{project_dir}/meta_analysis`` -- consensus MA output lands here."""
    return Path(project_dir) / "meta_analysis"


MASTER_ACCESSION = "_master"


def master_study_dir(project_dir: PathLike) -> Path:
    """``{project_dir}/_master`` -- synthetic study holding the Combine output."""
    return Path(project_dir) / MASTER_ACCESSION


def master_h5ad_path(project_dir: PathLike) -> Path:
    """``{project_dir}/_master/processed_data/master.h5ad``."""
    return processed_data_dir(master_study_dir(project_dir)) / "master.h5ad"


# ---------------------------------------------------------------------------
# Helpers for the Meta workspace that crosses study boundaries
# ---------------------------------------------------------------------------

def meta_de_result_path(project_dir: PathLike, accession: str, method: str) -> Path:
    """DE-result path for ``accession`` inside a meta project.

    Equivalent to ``de_result_path(study_dir(project_dir, accession), accession, method)``.
    """
    return de_result_path(study_dir(project_dir, accession), accession, method)


def meta_pseudobulk_path(project_dir: PathLike, accession: str) -> Path:
    """Pseudobulk path for ``accession`` inside a meta project."""
    return pseudobulk_path(study_dir(project_dir, accession), accession)


# ---------------------------------------------------------------------------
# Scaffolding -- create the canonical layout on disk
# ---------------------------------------------------------------------------

def scaffold_project(project_dir: PathLike) -> None:
    """Create the project root and its ``meta_analysis/`` subdir if missing.

    Idempotent. Safe to call on existing projects.
    """
    Path(project_dir).mkdir(parents=True, exist_ok=True)
    meta_output_dir(project_dir).mkdir(parents=True, exist_ok=True)


def scaffold_study(study_dir: PathLike) -> None:
    """Create a study folder and all its canonical subdirectories.

    Creates ``raw_data/``, ``processed_data/``, ``figures/``,
    ``results/figures/``, ``results/de_analysis/statistics/``, and
    ``results/de_analysis/pathway_scoring/``. Idempotent.
    """
    study_dir = Path(study_dir)
    study_dir.mkdir(parents=True, exist_ok=True)
    raw_data_dir(study_dir).mkdir(parents=True, exist_ok=True)
    processed_data_dir(study_dir).mkdir(parents=True, exist_ok=True)
    scrna_figures_dir(study_dir).mkdir(parents=True, exist_ok=True)
    de_figures_dir(study_dir).mkdir(parents=True, exist_ok=True)
    de_stats_dir(study_dir).mkdir(parents=True, exist_ok=True)
    de_pathway_scoring_dir(study_dir).mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Discovery -- enumerate studies and report their status
# ---------------------------------------------------------------------------

# Reserved subfolder names inside a project_dir that are not studies.
_PROJECT_RESERVED_DIRS = frozenset({"meta_analysis"})


def _looks_like_study(path: Path) -> bool:
    """Heuristic: does ``path`` contain any of a study's canonical subdirs?

    Required because users may pick a study folder by mistake instead of
    its parent project. Studies that have been scaffolded (or populated)
    will have at least one of raw_data/, processed_data/, or results/.
    """
    return (
        (path / "raw_data").is_dir()
        or (path / "processed_data").is_dir()
        or (path / "results").is_dir()
    )


def list_studies(project_dir: PathLike) -> list[str]:
    """Return accession names of study subfolders inside ``project_dir``.

    A direct subdirectory counts as a study only if it (a) is not
    ``meta_analysis``, (b) does not start with a dot, and (c) contains
    at least one of ``raw_data/``, ``processed_data/``, or ``results/``
    (so a study folder accidentally opened as a project shows zero
    studies, not three pseudo-studies named after its subdirs).
    Result is sorted lexicographically.
    """
    root = Path(project_dir)
    if not root.is_dir():
        return []
    return sorted(
        child.name for child in root.iterdir()
        if child.is_dir()
        and child.name not in _PROJECT_RESERVED_DIRS
        and not child.name.startswith(".")
        and _looks_like_study(child)
    )


def study_status(study_dir: PathLike) -> dict[str, bool]:
    """Return per-stage completion flags for a study folder.

    Keys: ``raw_data`` (any file in raw_data/), ``processed`` (any
    .h5ad in processed_data/), ``de_results`` (any *_DE_*.csv in
    de_analysis/statistics/), ``pseudobulk`` (any *_pseudobulk.csv in
    de_analysis/statistics/). Filename-based checks tolerate accession
    prefixes that include a cell-type suffix (e.g. ``GSE123_DCMEC``)
    that doesn't match the folder name.
    """
    study = Path(study_dir)
    raw = raw_data_dir(study)
    processed = processed_data_dir(study)
    stats = de_stats_dir(study)
    return {
        "raw_data": raw.is_dir() and any(raw.iterdir()),
        "processed": processed.is_dir() and any(processed.glob("*.h5ad")),
        "de_results": stats.is_dir() and any(stats.glob("*_DE_*.csv")),
        "pseudobulk": stats.is_dir() and any(stats.glob("*_pseudobulk.csv")),
    }


__all__ = [
    # study-level
    "raw_data_dir",
    "processed_data_dir",
    "scrna_figures_dir",
    "results_dir",
    "de_analysis_dir",
    "de_stats_dir",
    "de_pathway_scoring_dir",
    "de_result_path",
    "pseudobulk_path",
    "de_figures_dir",
    "metabolic_de_all_path",
    "metabolic_de_significant_path",
    "processed_h5ad_path",
    # meta-level
    "study_dir",
    "meta_output_dir",
    "meta_de_result_path",
    "meta_pseudobulk_path",
    # combine-level
    "MASTER_ACCESSION",
    "master_study_dir",
    "master_h5ad_path",
    # scaffolding
    "scaffold_project",
    "scaffold_study",
    # discovery
    "list_studies",
    "study_status",
]
