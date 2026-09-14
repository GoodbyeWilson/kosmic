"""Project workspace workflow: the five steps and when each is complete.

Qt-free. The Project page renders these rules; this module owns them so
they can be tested and so the documentation can state them exactly
(ADR-002).

A project folder holds three kinds of subfolder that look alike on
disk and in the study table:

- a **study** -- one published dataset, named by accession
- the **atlas** -- ``_master``, the Combine output
- a **subset** -- a cell-type slice of a study or of the atlas, written
  by the scRNA Subset step as a new sibling folder

The kind is read from the provenance record (``source`` points at the
parent h5ad for a subset), never guessed from the name.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from kosmic import provenance
from kosmic.paths import MASTER_ACCESSION, processed_data_dir

STEP_PROJECT, STEP_STUDIES, STEP_DATA, STEP_ATLAS, STEP_REVIEW = range(5)

# (name, hint) pairs in the shape 'ExplorerPane.set_workflow_steps' takes.
STEPS = (
    ("Project", "Choose the folder that holds this analysis."),
    ("Studies", "Add each dataset you want to compare as a study."),
    ("Data", "Import or download data for every study."),
    ("Shared Atlas (optional)",
     "Combine studies so cell types are annotated consistently."),
    ("Review", "Check every study is ready for analysis."),
)

# Two studies is the floor for both pooling and the shared-gene count.
MIN_STUDIES = 2


def study_kind(study_dir) -> tuple[str, str | None]:
    """Classify a study folder as ``('study'|'atlas'|'subset', parent)``.

    ``parent`` is the accession a subset was cut from, else None.
    """
    study = Path(study_dir)
    if study.name == MASTER_ACCESSION:
        return "atlas", None
    rec = provenance.load(processed_data_dir(study)) or {}
    source = rec.get("source")
    if source:
        # source is '<project>/<parent>/processed_data/<file>.h5ad'
        return "subset", Path(source).parent.parent.name
    return "study", None


@dataclass
class ProjectProgress:
    """What the five steps see, computed from per-folder status flags."""

    project_open: bool = False
    n_studies: int = 0            # real datasets only (kind == 'study')
    n_rows: int = 0               # every folder in the table
    n_with_data: int = 0          # studies with a processed dataset
    missing_data: list[str] = field(default_factory=list)
    has_atlas: bool = False
    n_ready: int = 0              # rows with DE results + pseudobulk
    complete: tuple[bool, bool, bool, bool, bool] = (False,) * 5

    @property
    def status(self) -> tuple[str, str]:
        """One line for the explorer STATUS panel: ``(text, state)``."""
        if not self.project_open:
            return "Open or create a project to begin", "warning"
        if self.n_studies < MIN_STUDIES:
            n = self.n_studies
            return (f"{n} {'study' if n == 1 else 'studies'} — add at least "
                    f"{MIN_STUDIES} to compare", "warning")
        if self.missing_data:
            k = len(self.missing_data)
            return (f"{k} {'study' if k == 1 else 'studies'} without data",
                    "warning")
        if self.n_ready < MIN_STUDIES:
            return (f"{self.n_ready} of {self.n_rows} ready for "
                    f"meta-analysis", "info")
        return (f"{self.n_ready} of {self.n_rows} ready for meta-analysis",
                "success")


def project_progress(rows: dict[str, dict] | None,
                     project_open: bool = True) -> ProjectProgress:
    """Compute step completion.

    ``rows`` maps accession -> ``{'kind': str, 'status': study_status
    dict}``; pass None (or ``project_open=False``) for no project.
    """
    if rows is None or not project_open:
        return ProjectProgress()

    studies = [a for a, r in rows.items() if r.get("kind") == "study"]
    missing = [a for a in studies
               if not (rows[a].get("status") or {}).get("processed")]
    ready = sum(1 for a, r in rows.items()
                if a != MASTER_ACCESSION
                and (r.get("status") or {}).get("de_results")
                and (r.get("status") or {}).get("pseudobulk"))
    has_atlas = bool((rows.get(MASTER_ACCESSION, {}).get("status") or {})
                     .get("processed"))

    n_studies = len(studies)
    complete = (
        True,
        n_studies >= MIN_STUDIES,
        n_studies >= 1 and not missing,
        has_atlas,
        ready >= MIN_STUDIES,
    )
    return ProjectProgress(
        project_open=True,
        n_studies=n_studies,
        n_rows=len(rows),
        n_with_data=n_studies - len(missing),
        missing_data=missing,
        has_atlas=has_atlas,
        n_ready=ready,
        complete=complete,
    )
