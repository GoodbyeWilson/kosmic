"""The merged pseudobulk design, matched to the per-study designs (ADR-017 §4).

Merged pseudobulk and a meta-analysis can be compared only when the
merged model adjusts for the same covariates as the per-study models,
plus study. The studies and the atlas name some variables differently:
a study's model may use ``sex_inferred`` while the atlas holds ``sex``,
which Combine writes from each study's designated sex column
(``uns['sex_column']``); a study may identify its donors by ``donor_id``
while the atlas uses ``sample``; the per-study runs may select the cell
type from ``cell_type_atlas`` (labels propagated from the atlas) while
the atlas's own column is ``cell_type``.

This module reads the per-study DE run records, maps each covariate to
an atlas column, and checks the values behind each mapping rather than
the names: a sex covariate is accepted when, for every donor, it agrees
with the column the atlas's ``sex`` was written from; a different sample
column is accepted when it identifies the same donors as ``sample``.
Every problem found is returned as a warning; no covariate is dropped.

:func:`resolve_design` works on data frames and is what the tests call;
:func:`load_design` reads the files of a project and calls it.
:func:`differences` compares a saved merged run's record with a design.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

#: The atlas columns Combine always writes (``kosmic.combine.concat_studies``).
STUDY_COLUMN = 'study'
SAMPLE_COLUMN = 'sample'
SEX_COLUMN = 'sex'
ROLE_COLUMN = '_role'

#: Sex covariates written by KOSMIC itself, besides recorded sex columns.
_INFERRED_SEX_COLUMNS = ('sex_inferred',)


@dataclass(frozen=True)
class StudyDesign:
    """One study's gene DE run, as its record describes it.

    ``covariates`` are the study's own column names; ``counts`` is the
    counts folder ('decontX_counts', 'original_counts').
    """

    study: str
    de_path: Path
    covariates: tuple
    sample_column: str
    counts: Optional[str]
    cell_type: Optional[str]
    cell_type_column: Optional[str]
    min_cells: Optional[int]
    min_counts: Optional[int]
    method: str = 'deseq2'
    filters: Optional[dict] = None
    fdr_scope: Optional[str] = None
    gene_sets: Optional[dict] = None


@dataclass
class Design:
    """The merged model's design and what was found while building it.

    ``covariates`` are atlas column names, in the order of the per-study
    covariates, without ``study``; ``formula`` adds ``study`` and
    ``condition``. ``mapping`` gives each per-study covariate's atlas
    column (None when the atlas has none). ``unresolved`` lists the
    per-study covariates with no atlas column; while it is not empty the
    merged model cannot be fitted to match the studies.
    """

    covariates: list
    mapping: dict
    sample_column: str
    cell_type_column: Optional[str]
    cell_type: Optional[str]
    counts: Optional[str]
    min_cells: Optional[int]
    min_counts: Optional[int]
    warnings: list = field(default_factory=list)
    unresolved: list = field(default_factory=list)
    filters: Optional[dict] = None
    fdr_scope: Optional[str] = None
    gene_sets: Optional[dict] = None

    @property
    def formula(self) -> str:
        return '~' + ' + '.join([*self.covariates, STUDY_COLUMN, 'condition'])

    @property
    def can_fit(self) -> bool:
        return not self.unresolved


def study_design(study: str, de_path, record: dict) -> StudyDesign:
    """A :class:`StudyDesign` from a run's record (``kosmic.de.saved.read_run_record``)."""
    from kosmic.de.saved import run_counts
    return StudyDesign(
        study=study, de_path=Path(de_path),
        covariates=tuple(str(c) for c in record.get('covariates') or ()),
        sample_column=str(record.get('sample_column') or SAMPLE_COLUMN),
        counts=run_counts(de_path),
        cell_type=record.get('cell_type'),
        cell_type_column=record.get('cell_type_column'),
        min_cells=record.get('min_cells'), min_counts=record.get('min_counts'),
        method=str(record.get('method') or 'deseq2'),
        filters=record.get('filters'), fdr_scope=record.get('fdr_scope'),
        gene_sets=record.get('pathway_gene_sets') or None)


#: Filters that differ between the runs by design: a study filters genes over
#: its own donors, the merged model within each study.
_FILTERS_NOT_COMPARED = ('detection_study_col',)


def comparable_filters(filters: Optional[dict]) -> Optional[dict]:
    """``filters`` without the keys that differ between the runs by design."""
    if filters is None:
        return None
    return {k: (sorted(v) if isinstance(v, (list, tuple)) else v)
            for k, v in filters.items() if k not in _FILTERS_NOT_COMPARED}


def uses_gene_sets(fdr_scope: Optional[str]) -> bool:
    """Whether a run's BH correction ran over the committed gene sets only."""
    return bool(fdr_scope) and 'gene_sets' in fdr_scope


def _is_sex_covariate(name: str, designated: Optional[str]) -> bool:
    from kosmic.scrna.inspect.sex import SEX_COLUMN_CANDIDATES
    return name in SEX_COLUMN_CANDIDATES or name in _INFERRED_SEX_COLUMNS or name == designated


def _donor_values(obs: pd.DataFrame, donor_col: str, col: str) -> pd.Series:
    """Each donor's value of ``col`` ('' where missing); several values per
    donor are joined with '|', so a donor-level mismatch stays visible."""
    values = obs[col].astype(object).where(obs[col].notna(), '').astype(str)
    return values.groupby(obs[donor_col].astype(str).values).agg(
        lambda v: '|'.join(sorted(set(v))))


def _sex_disagreements(obs: pd.DataFrame, donor_col: str, used: str, source: str) -> list:
    """Donors whose sex in ``used`` differs from ``source`` (both normalised)."""
    from kosmic.scrna.inspect.sex import normalise_sex_value
    a = _donor_values(obs, donor_col, used).map(lambda v: normalise_sex_value(v) or v)
    b = _donor_values(obs, donor_col, source).map(lambda v: normalise_sex_value(v) or v)
    return sorted(a.index[a != b.reindex(a.index)])


def _same_donors(obs: pd.DataFrame, a: str, b: str) -> bool:
    """Whether columns ``a`` and ``b`` split the cells into the same donors."""
    pairs = obs[[a, b]].astype(str).drop_duplicates()
    return pairs[a].is_unique and pairs[b].is_unique


def _atlas_cell_type_column(column: Optional[str], atlas_columns) -> Optional[str]:
    """The atlas column holding the labels a study selected its cell type by:
    the same column, or for a propagated ``<name>_atlas`` column the
    atlas's own ``<name>``."""
    if not column or column in atlas_columns and not column.endswith('_atlas'):
        return column
    stem = column[:-len('_atlas')] if column.endswith('_atlas') else None
    if stem and stem in atlas_columns:
        return stem
    return column if column in atlas_columns else None


def _donor_table(atlas_obs: pd.DataFrame, columns, cell_type_column=None, cell_type=None):
    """One row per atlas donor (study and sample) with the given columns."""
    obs = atlas_obs
    if cell_type_column and cell_type and cell_type_column in obs:
        obs = obs[obs[cell_type_column].astype(str) == str(cell_type)]
    keep = [c for c in dict.fromkeys([STUDY_COLUMN, SAMPLE_COLUMN, ROLE_COLUMN, *columns])
            if c in obs]
    table = obs[keep].astype(object).astype(str).drop_duplicates(
        [c for c in (STUDY_COLUMN, SAMPLE_COLUMN) if c in keep])
    return table.reset_index(drop=True)


def _confounding_warnings(donors: pd.DataFrame, covariate: str) -> list:
    """Why ``covariate`` cannot be estimated beside study and condition, if it cannot."""
    if covariate not in donors:
        return []
    values = donors[covariate]
    out = []
    missing = (values == '') | (values.str.lower() == 'nan')
    if missing.any():
        out.append(f"'{covariate}' is missing for {int(missing.sum())} of "
                   f"{len(donors)} atlas donors.")
    present = donors[~missing]
    if present[covariate].nunique() < 2:
        return out + [f"'{covariate}' has one value across the atlas donors, so it "
                      "cannot be estimated."]
    varying = sorted(s for s, v in present.groupby(STUDY_COLUMN)[covariate]
                     if v.nunique() > 1) if STUDY_COLUMN in present else []
    if STUDY_COLUMN in present and not varying:
        out.append(f"'{covariate}' does not vary within any study, so it is "
                   "confounded with study.")
    elif len(varying) == 1:
        out.append(f"'{covariate}' varies within one study only ({varying[0]}).")
    if ROLE_COLUMN in present:
        by_level = present.groupby(covariate)[ROLE_COLUMN].nunique()
        if (by_level < 2).all():
            out.append(f"'{covariate}' is confounded with condition: each of its "
                       "values occurs in one condition only.")
    return out


def resolve_design(studies: list, study_obs: dict, sex_columns: dict,
                   atlas_obs: pd.DataFrame) -> Design:
    """The merged design that matches ``studies``.

    Parameters
    ----------
    studies : list of StudyDesign
        The pooled studies' gene DE runs.
    study_obs : dict
        ``{study: obs}`` of each study's dataset. A study without one is
        not checked by value, and a warning says so.
    sex_columns : dict
        ``{study: designated sex column or None}``, the ``uns['sex_column']``
        Combine wrote the atlas's ``sex`` from.
    atlas_obs : DataFrame
        The atlas's obs (``study``, ``sample``, ``_role`` and the covariates).
    """
    from kosmic.scrna.inspect.sex import recorded_sex_column

    warnings: list = []
    atlas_cols = set(atlas_obs.columns)
    used = {s.study: s.covariates for s in studies}
    covariates = list(dict.fromkeys(c for s in studies for c in s.covariates))
    if len({tuple(sorted(c)) for c in used.values()}) > 1:
        detail = '; '.join(f"{s}: {', '.join(c) or 'none'}" for s, c in used.items())
        warnings.append(f"The studies used different covariates ({detail}).")

    mapping: dict = {}
    for cov in covariates:
        if cov in atlas_cols:
            mapping[cov] = cov
        elif SEX_COLUMN in atlas_cols and any(
                _is_sex_covariate(cov, sex_columns.get(s.study)) for s in studies):
            mapping[cov] = SEX_COLUMN
        else:
            mapping[cov] = None
    unresolved = [c for c, m in mapping.items() if m is None]
    for cov in unresolved:
        warnings.append(f"The atlas has no column for the per-study covariate '{cov}'.")

    for s in studies:
        obs = study_obs.get(s.study)
        if obs is None:
            warnings.append(f"{s.study}: its dataset could not be read, so its covariates "
                            "and donors were not checked against the atlas.")
            continue
        donor_col = s.sample_column if s.sample_column in obs else None
        if s.sample_column != SAMPLE_COLUMN:
            if SAMPLE_COLUMN not in obs or donor_col is None:
                warnings.append(f"{s.study}: donors are identified by '{s.sample_column}' "
                                f"and the atlas uses '{SAMPLE_COLUMN}', which the study "
                                "does not have.")
            elif not _same_donors(obs, s.sample_column, SAMPLE_COLUMN):
                warnings.append(f"{s.study}: '{s.sample_column}' and '{SAMPLE_COLUMN}' "
                                "identify different donors, so the atlas's donors are "
                                "not the study's.")
        for cov in s.covariates:
            if mapping.get(cov) != SEX_COLUMN:
                continue
            source = sex_columns.get(s.study) or recorded_sex_column(obs)
            if source == cov or donor_col is None or cov not in obs:
                continue
            if source is None or source not in obs:
                warnings.append(f"{s.study}: the model used '{cov}', and the study has no "
                                "designated sex column for the atlas's 'sex'.")
                continue
            bad = _sex_disagreements(obs, donor_col, cov, source)
            if bad:
                shown = ', '.join(bad[:5]) + (', ...' if len(bad) > 5 else '')
                warnings.append(f"{s.study}: the model used '{cov}', but the atlas's 'sex' "
                                f"comes from '{source}', which differs for {len(bad)} "
                                f"donor(s) ({shown}).")

    ref = studies[0] if studies else None
    for attr, label in (('counts', 'counts'), ('min_cells', 'minimum cells per donor'),
                        ('min_counts', 'minimum counts per donor'),
                        ('cell_type', 'cell type')):
        values = {getattr(s, attr) for s in studies}
        if len(values) > 1:
            warnings.append(f"The studies' runs differ in {label}: "
                            + ', '.join(sorted(map(str, values))) + '.')

    filters, scope, gene_sets = None, None, None
    recorded = [comparable_filters(s.filters) for s in studies]
    if studies and all(f is not None for f in recorded):
        if all(f == recorded[0] for f in recorded):
            filters = dict(studies[0].filters)
        else:
            detail = '; '.join(f"{s.study}: {f}" for s, f in zip(studies, recorded))
            warnings.append(f"The studies' runs used different gene filters ({detail}).")
    scopes = {uses_gene_sets(s.fdr_scope) for s in studies if s.fdr_scope}
    if len(scopes) > 1:
        warnings.append("Some studies corrected for multiple testing over the committed "
                        "gene sets and others over every gene.")
    elif scopes and studies and all(s.fdr_scope for s in studies):
        scope = studies[0].fdr_scope
        gene_sets = studies[0].gene_sets if uses_gene_sets(scope) else None

    atlas_covs = list(dict.fromkeys(m for m in mapping.values() if m))
    ct_col = _atlas_cell_type_column(ref.cell_type_column if ref else None, atlas_cols)
    if ref and ref.cell_type and ct_col is None:
        warnings.append(f"The atlas has no column for the cell type labels "
                        f"'{ref.cell_type_column}'.")
    donors = _donor_table(atlas_obs, atlas_covs, ct_col, ref.cell_type if ref else None)
    for cov in atlas_covs:
        warnings.extend(_confounding_warnings(donors, cov))

    return Design(covariates=atlas_covs, mapping=mapping, sample_column=SAMPLE_COLUMN,
                  cell_type_column=ct_col, cell_type=ref.cell_type if ref else None,
                  counts=ref.counts if ref else None,
                  min_cells=ref.min_cells if ref else None,
                  min_counts=ref.min_counts if ref else None,
                  warnings=warnings, unresolved=unresolved,
                  filters=filters, fdr_scope=scope, gene_sets=gene_sets)


def differences(record: Optional[dict], design: Design, counts: Optional[str] = None) -> list:
    """How a saved merged run differs from ``design``; empty when it matches.

    ``record`` is the run's record; ``counts`` its counts folder
    (``kosmic.de.saved.run_counts``). Covariates are compared as sets of
    atlas columns, so the order of the terms does not matter.
    """
    if not record:
        return ["The saved merged run has no record of its design."]
    out = []
    saved = set(map(str, record.get('covariates') or ()))
    wanted = {*design.covariates, STUDY_COLUMN}
    if STUDY_COLUMN not in saved:
        out.append("The saved merged run does not include study as a covariate.")
    if saved - wanted - {STUDY_COLUMN}:
        out.append("The saved merged run adjusts for covariates the studies did not use: "
                   + ', '.join(sorted(saved - wanted)) + '.')
    if wanted - saved - {STUDY_COLUMN}:
        out.append("The saved merged run does not adjust for: "
                   + ', '.join(sorted(wanted - saved - {STUDY_COLUMN})) + '.')
    if design.counts and counts and counts != design.counts:
        out.append(f"The saved merged run used {counts}; the studies used {design.counts}.")
    for key, label in (('min_cells', 'minimum cells per donor'),
                       ('min_counts', 'minimum counts per donor')):
        want = getattr(design, key)
        if want is not None and record.get(key) is not None and record.get(key) != want:
            out.append(f"The saved merged run's {label} is {record.get(key)}; "
                       f"the studies used {want}.")
    want_filters = comparable_filters(design.filters)
    if want_filters is not None:
        saved_filters = comparable_filters(record.get('filters'))
        if saved_filters is None:
            out.append("The saved merged run did not record its gene filters, so they "
                       "cannot be checked against the studies'.")
        elif saved_filters != want_filters:
            changed = sorted(k for k in set(saved_filters) | set(want_filters)
                             if saved_filters.get(k) != want_filters.get(k))
            out.append("The saved merged run's gene filters differ from the studies': "
                       + ', '.join(changed) + '.')
    saved_scope = record.get('fdr_scope')
    if design.fdr_scope and saved_scope and \
            uses_gene_sets(saved_scope) != uses_gene_sets(design.fdr_scope):
        over = {True: 'the gene sets', False: 'every gene'}
        out.append("The saved merged run corrected for multiple testing over "
                   f"{over[uses_gene_sets(saved_scope)]}; the studies over "
                   f"{over[uses_gene_sets(design.fdr_scope)]}.")
    sample = record.get('sample_column')
    if sample and sample != design.sample_column:
        out.append(f"The saved merged run identifies donors by '{sample}', "
                   f"not '{design.sample_column}'.")
    return out


def _designated_sex_column(h5ad_path) -> Optional[str]:
    """A study file's ``uns['sex_column']``, read without loading the file."""
    import h5py
    from kosmic.scrna.load.h5ad_meta import file_lock
    with file_lock(h5ad_path).reading(), h5py.File(str(h5ad_path), 'r') as f:
        if 'uns' not in f or 'sex_column' not in f['uns']:
            return None
        value = f['uns']['sex_column'][()]
    return value.decode() if isinstance(value, bytes) else str(value)


def load_design(project, de_paths: dict) -> Design:
    """:func:`resolve_design` for a project.

    ``de_paths`` maps each pooled study to its gene DE table. The studies'
    obs and designated sex columns and the atlas's obs are read from disk.
    """
    from kosmic.de.saved import read_run_record
    from kosmic.paths import master_h5ad_path, processed_h5ad_path, study_dir
    from kosmic.scrna.load.h5ad_meta import read_obs, registered_h5ad

    studies, study_obs, sex_columns = [], {}, {}
    for study, path in de_paths.items():
        record = read_run_record(path) or {}
        studies.append(study_design(study, path, record))
        folder = study_dir(project, study)
        h5ad = processed_h5ad_path(folder, record['dataset_file']) \
            if record.get('dataset_file') else None
        if h5ad is None or not h5ad.is_file():
            h5ad = registered_h5ad(folder)
        if h5ad is not None:
            try:
                study_obs[study] = read_obs(h5ad)
                sex_columns[study] = _designated_sex_column(h5ad)
            except OSError:
                pass
    atlas = master_h5ad_path(project)
    atlas_obs = read_obs(atlas) if atlas.is_file() else pd.DataFrame()
    design = resolve_design(studies, study_obs, sex_columns, atlas_obs)
    if not atlas.is_file():
        design.warnings.insert(0, "The project has no combined atlas (_master), so merged "
                                  "pseudobulk cannot be fitted.")
        design.unresolved = design.unresolved or ['atlas']
    return design
