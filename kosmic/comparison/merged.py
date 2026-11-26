"""Merged pseudobulk: one DESeq2 model on the shared atlas (ADR-017 §3).

The merged result is the atlas's (``_master``) saved gene DE run for the
selection's counts and cell type when that run's design matches the
per-study designs (:func:`kosmic.comparison.design.differences`). When
there is no saved run, or the user chooses to refit, :func:`run_merged`
fits it with the DE workspace's own batch function
(``kosmic.de.batch.run_de_by_cell_type``) and saves it where the DE
workspace saves the atlas's runs, so the DE workspace lists it too.

The refit uses the minimum cells and counts per donor that the studies'
runs recorded. The studies' runs do not record their gene filters
(detection thresholds, DESeq2 independent and Cook's filtering, gene
groups excluded), so the refit uses KOSMIC's defaults for those, as a
new DE run does, and filters genes within each study (``study``), as the
per-study runs filtered within their own donors.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

import pandas as pd

from kosmic.comparison.design import ROLE_COLUMN, STUDY_COLUMN, Design, differences

#: What the refit cannot take from the per-study records, shown with it.
DEFAULT_FILTERS_NOTE = (
    "The per-study runs do not record their gene filters, so merged pseudobulk "
    "is fitted with KOSMIC's default gene filters, applied within each study.")


@dataclass(frozen=True)
class SavedMerged:
    """The atlas's saved gene DE run for a selection, and how it differs
    from the design (empty ``differences`` when it matches)."""

    path: Path
    record: Optional[dict]
    counts: Optional[str]
    differences: tuple

    @property
    def matches(self) -> bool:
        return not self.differences

    def read(self) -> pd.DataFrame:
        from kosmic.meta_analysis.io import load_de_results
        return load_de_results(str(self.path))


def find_saved_merged(project, design: Design, method: str = 'deseq2') -> Optional[SavedMerged]:
    """The atlas's saved run on the design's counts and cell type, or None."""
    from kosmic.de.saved import list_gene_de_runs, read_run_record
    from kosmic.paths import master_study_dir

    folder = master_study_dir(project)
    if not folder.is_dir():
        return None
    want = (design.cell_type or '').replace('_', ' ').lower()
    for run in list_gene_de_runs(folder, counts=design.counts):
        if run.method != method or (run.cell_type or '').replace('_', ' ').lower() != want:
            continue
        if design.counts and run.counts not in (None, design.counts):
            continue
        record = read_run_record(run.path)
        return SavedMerged(path=run.path, record=record, counts=run.counts,
                           differences=tuple(differences(record, design, run.counts)))
    return None


def atlas_counts_layer(atlas_path, counts: Optional[str]) -> tuple:
    """``(layer, problem)``: the layer to read the design's counts from the
    atlas (None for the atlas's own counts), or a problem when the atlas
    does not hold those counts."""
    from kosmic.paths import counts_folder_name
    from kosmic.scrna.load.h5ad_meta import combined_from, read_layer_names

    own = counts_folder_name(combined_from(atlas_path) or 'raw')
    if counts in (None, own):
        return None, None
    layer = counts if counts != counts_folder_name('raw') else 'counts'
    if layer in read_layer_names(atlas_path):
        return layer, None
    return None, (f"The atlas holds {own.replace('_', ' ')}, not "
                  f"{counts.replace('_', ' ')}, so merged pseudobulk cannot be fitted "
                  "on the studies' counts.")


def run_merged(project, design: Design, *, progress: Optional[Callable[[str], None]] = None,
               cancelled: Optional[Callable[[], bool]] = None) -> SavedMerged:
    """Fit merged pseudobulk on the atlas with ``design`` and save it.

    Raises ValueError when the design cannot be fitted (a covariate with no
    atlas column, an atlas without the studies' counts) or the fit fails.
    """
    from kosmic import DE_MIN_CELLS, DE_MIN_COUNTS
    from kosmic.de.batch import run_de_by_cell_type
    from kosmic.de.saved import build_run_record
    from kosmic.paths import MASTER_ACCESSION, master_h5ad_path, master_study_dir
    from kosmic.scrna.load.h5ad_meta import read_counts_adata

    if not design.can_fit:
        raise ValueError("Merged pseudobulk cannot be fitted: the atlas has no column for "
                         + ', '.join(design.unresolved) + '.')
    atlas = master_h5ad_path(project)
    layer, problem = atlas_counts_layer(atlas, design.counts)
    if problem:
        raise ValueError(problem)
    say = progress or (lambda _m: None)
    say(f"Reading the atlas's counts ({atlas.name})...")
    adata = read_counts_adata(atlas, counts_layer=layer)
    min_cells = DE_MIN_CELLS if design.min_cells is None else int(design.min_cells)
    min_counts = DE_MIN_COUNTS if design.min_counts is None else int(design.min_counts)
    covariates = [*design.covariates, STUDY_COLUMN]
    record = build_run_record(
        adata, method='deseq2', sample_col=design.sample_column, condition_col=ROLE_COLUMN,
        disease_label='disease', control_label='control', min_cells=min_cells,
        min_counts=min_counts, covariates=covariates, h5ad_path=atlas)
    say(f"Fitting merged pseudobulk: {design.formula}")
    runs = run_de_by_cell_type(
        adata, design.cell_type_column, design.sample_column, ROLE_COLUMN,
        master_study_dir(project), MASTER_ACCESSION, cell_types=[design.cell_type],
        de_method='deseq2', covariates=covariates, min_cells=min_cells,
        min_counts=min_counts, counts_layer=None, detection_study_col=STUDY_COLUMN,
        progress_callback=say, cancelled=cancelled, run_record=record)
    run = runs[0] if runs else None
    if run is None or getattr(run, 'de_path', None) is None:
        reason = (run.message if run is not None else '') or 'no result'
        raise ValueError(f"Merged pseudobulk was not fitted: {reason}")
    saved = find_saved_merged(project, design)
    if saved is None:
        raise ValueError("Merged pseudobulk was fitted but its saved table was not found.")
    return saved
