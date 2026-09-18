# Cross-study dataset assembly for the Combine workspace.
#
# Two operations:
#   - concat_studies(): build a master h5ad from per-study h5ads using
#     disk-backed concat so peak RAM stays bounded to one study at a
#     time. Uses each study's raw integer counts (from .raw.to_adata()
#     when available, else .X). Optional per-study subsample.
#   - propagate_labels(): after the user clusters + annotates the master
#     in the scRNA pipeline, write the cell-type labels back to each
#     per-study h5ad. Uses obs-only barcode join when every cell is
#     present in the master, falls back to sc.tl.ingest otherwise.
from __future__ import annotations

import gc
import shutil
import tempfile
from pathlib import Path
from typing import Callable, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd


# obs columns carried from each study onto the master. The atlas labels
# ('cell_type_atlas', 'leiden_atlas') are there so an atlas re-created
# after propagation -- for example uncapped, as the container for the
# mega-analysis DE -- still carries the shared labels; 'sex' and 'age' so
# they are available as DE covariates.
# Included cells absent from the master, as a fraction of the study's
# included cells, above which they are projected by kNN (a cap was used)
# rather than left unlabelled (the master's own QC removed a few).
PROJECT_MIN_FRACTION = 0.01

_OBS_KEEP_CANDIDATES = ('sample', 'condition', 'cell_type', '_role',
                        'cell_type_atlas', 'leiden_atlas', 'sex', 'age')


def _kept_columns_across(study_h5ad_paths) -> list[str]:
    """The kept obs columns any of the studies has (obs read only, no X).

    'sex' counts as present when the study designates a sex column in
    uns['sex_column'] or carries a recognisable one, since _unify_sex_column
    will write it.
    """
    import h5py
    from kosmic.scrna.inspect.sex import SEX_COLUMN_CANDIDATES
    found: set[str] = set()
    for path in study_h5ad_paths:
        with h5py.File(path, 'r') as h:
            cols = set(h['obs'].keys()) if 'obs' in h else set()
            found |= {c for c in _OBS_KEEP_CANDIDATES if c in cols}
            has_designation = 'uns' in h and 'sex_column' in h['uns']
            if has_designation or any(c in cols for c in SEX_COLUMN_CANDIDATES):
                found.add('sex')
    return [c for c in _OBS_KEEP_CANDIDATES if c in found]


def _unify_sex_column(a, designated: Optional[str]) -> None:
    """Write obs['sex'] as 'female' / 'male' (or '') from the designated column.

    'designated' is the study's uns['sex_column'] -- the column the Inspect
    step settled on ('sex_inferred' when the sex was called from XIST and
    Y genes). With no designation, an existing recorded sex column is used,
    normalised the same way; with neither, no column is written.
    """
    from kosmic.scrna.inspect.sex import normalise_sex_value, recorded_sex_column
    col = designated
    if not col or col not in a.obs.columns:
        col = recorded_sex_column(a.obs)
    if not col:
        return
    values = a.obs[col].astype(object).map(normalise_sex_value)
    a.obs['sex'] = values.fillna('').astype(str).values


def concat_studies(
    study_h5ad_paths: Sequence[Path],
    output_path: Path,
    cap_per_study: Optional[int] = None,
    seed: int = 0,
    role_map: Optional[dict] = None,
    drop_genes: Optional[Sequence[str]] = None,
    drop_excluded: bool = False,
    progress_callback: Optional[Callable[[str], None]] = None,
) -> Path:
    """Build a master h5ad from per-study h5ads using disk-backed concat.

    For each input h5ad:
      1. Load it (one at a time -- peak RAM = single study).
      2. Put the raw counts in '.X' (from layers['counts'], else '.raw',
         else X itself) so the temp holds integer counts.
      3. Optional subsample to 'cap_per_study'.
      4. Strip obsm / varm / uns / layers; keep a tight obs subset.
      5. Write to a temp h5ad.

    Then 'anndata.experimental.concat_on_disk' joins the temps into the
    master without loading them all at once. Inner join on genes.

    Parameters
    ----------
    study_h5ad_paths : sequence of Path
        Source h5ad files (typically each study's raw_data/*.h5ad).
    output_path : Path
        Where to write the master h5ad.
    cap_per_study : int, optional
        Random subsample cap. None keeps every cell.
    seed : int
        Subsample RNG seed.
    role_map : dict, optional
        Per-study (column, {value: role}) mapping that becomes the
        unified '_role' obs column on the master. Shape:
        '{accession: (condition_col_name, {raw_value: role})}'. Roles
        are 'control' / 'disease' / 'exclude'. When omitted, '_role' is
        carried through from each study's existing obs only if present
        there; otherwise the master has no '_role' column.
    drop_excluded : bool
        Leave out cells whose ``obs['_role']`` is ``'exclude'``. Off by
        default so the behaviour is unchanged unless asked for. A study's
        off-target arm is often a large fraction of it -- across the two
        DCM studies it is ~49% of the cells -- and those cells never enter
        a contrast, so carrying them into the master means clustering,
        annotating and storing them for nothing.
    drop_genes : sequence of str, optional
        Gene names to exclude from every study before the join. Intended
        for genes that harmonisation merged inconsistently across studies
        (see 'kosmic.combine.gene_merge'), where the same column name is a
        sum of two features in one cohort and a single feature in another.
        Names absent from a study are ignored.
    progress_callback : callable, optional
        Status string callback.

    Returns
    -------
    Path
        The output path (now written to disk).
    """
    rng = np.random.default_rng(seed)
    drop = set(drop_genes or ())
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    temp_dir = Path(tempfile.mkdtemp(prefix='kosmic_combine_'))
    temp_files: list[Path] = []
    study_names: list[str] = []

    # concat_on_disk joins obs columns by intersection, so a column that
    # only some studies carry (age, or atlas labels on a study added
    # later) would vanish. Every study gets the union of the kept
    # columns, blank where it has none.
    union_cols = _kept_columns_across(study_h5ad_paths)

    try:
        for path in study_h5ad_paths:
            path = Path(path)
            accession = path.parent.parent.name
            if progress_callback:
                progress_callback(
                    f"Extracting raw counts: {accession} (from {path.parent.name}/)...")

            a = ad.read_h5ad(path)
            sex_col = a.uns.get('sex_column') if isinstance(a.uns, dict) else None
            from kosmic.scrna.counts import counts_adata
            a = counts_adata(a, copy=False)

            if cap_per_study is not None and a.n_obs > cap_per_study:
                idx = np.sort(
                    rng.choice(a.n_obs, size=cap_per_study, replace=False))
                a = a[idx, :].copy()

            if drop_excluded:
                from kosmic.scrna.inspect.roles import included_mask
                mask = included_mask(a)
                if mask is not None:
                    n_out = int((~mask).sum())
                    a = a[mask].copy()
                    if progress_callback:
                        progress_callback(
                            f"  {accession}: leaving out {n_out:,} excluded "
                            f"cell(s), keeping {a.n_obs:,}")

            # Drop excluded genes here rather than after the join: the
            # temp files are what concat_on_disk reads, and trimming now
            # keeps them smaller too.
            if drop:
                keep = [g for g in a.var_names if g not in drop]
                if len(keep) < a.n_vars:
                    a = a[:, keep].copy()

            # One 'sex' column in female/male form, from whichever column
            # the study's Inspect step designated (recorded or inferred).
            _unify_sex_column(a, sex_col)

            # Strip per-study cruft; PCA/UMAP/etc. recomputed on master.
            a.obsm = None
            a.varm = None
            a.uns = {}
            a.layers = None

            # Apply the user-supplied role map (if any) BEFORE the obs
            # column allowlist trims, so we can read the source condition
            # column even if it's not in _OBS_KEEP_CANDIDATES.
            if role_map and accession in role_map:
                cond_col, value_map = role_map[accession]
                if cond_col in a.obs.columns:
                    src = a.obs[cond_col].astype(str)
                    a.obs['_role'] = src.map(value_map).fillna('').astype(str)

            keep_cols = [c for c in _OBS_KEEP_CANDIDATES if c in a.obs.columns]
            a.obs = a.obs[keep_cols].copy()
            for c in union_cols:
                if c not in a.obs.columns:
                    a.obs[c] = np.nan if c == 'age' else ''
            a.obs = a.obs[[c for c in _OBS_KEEP_CANDIDATES if c in a.obs.columns]]
            # One dtype per column across studies, or the on-disk concat
            # cannot write it: age is numeric (a study that recorded it as
            # text, '68' or '', becomes 68.0 / NaN), the rest are strings.
            if 'age' in a.obs.columns:
                a.obs['age'] = pd.to_numeric(
                    a.obs['age'].astype(object).replace('', np.nan), errors='coerce').astype(float)
            for c in a.obs.columns:
                if c != 'age' and (isinstance(a.obs[c].dtype, pd.CategoricalDtype) or a.obs[c].dtype == object):
                    a.obs[c] = a.obs[c].astype(object).where(a.obs[c].notna(), '').astype(str)

            temp_path = temp_dir / f"{accession}.h5ad"
            a.write_h5ad(temp_path)
            temp_files.append(temp_path)
            study_names.append(accession)
            del a
            gc.collect()

        if progress_callback:
            progress_callback(
                f"Concatenating {len(temp_files)} studies on disk...")
        # concat_on_disk refuses to overwrite; remove an existing master.
        if output_path.exists():
            output_path.unlink()
        ad.experimental.concat_on_disk(
            in_files=[str(p) for p in temp_files],
            out_file=str(output_path),
            label='study',
            keys=study_names,
            join='inner',
            index_unique='-',
        )
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    _ensure_int64_sparse(output_path, progress_callback)

    if progress_callback:
        from kosmic.scrna.load.h5ad_meta import read_shape
        n_obs, n_vars = read_shape(output_path)
        progress_callback(f"Master: {n_obs:,} cells x {n_vars:,} genes")

    return output_path


def _ensure_int64_sparse(
    path: Path,
    progress_callback: Optional[Callable[[str], None]] = None,
) -> None:
    """Cast a sparse master's indices / indptr to int64 if needed.

    concat_on_disk writes int32 index arrays. When total nnz exceeds
    2**31, scipy's C routines (eliminate_zeros, axis sums) raise
    'Output dtype not compatible with inputs'. Detection is cheap --
    peek at the X/indices HDF5 dataset's dtype and length via h5py
    without loading anything. The fix loads + rewrites the master once.
    """
    import h5py
    import scipy.sparse as sp

    with h5py.File(path, 'r') as h:
        if 'X' not in h or 'indices' not in h.get('X', {}):
            return  # dense X or non-sparse layout; nothing to do
        indices_dtype = h['X/indices'].dtype
        nnz = h['X/indices'].shape[0]

    if indices_dtype == np.int64:
        return
    if nnz <= 2**31 - 1:
        return  # int32 indices are valid for this nnz

    if progress_callback:
        progress_callback(
            f"Total nnz = {nnz:,} exceeds int32; rewriting sparse "
            f"indices as int64...")
    a = ad.read_h5ad(path)
    X = a.X
    a.X = sp.csr_matrix(
        (X.data, X.indices.astype(np.int64), X.indptr.astype(np.int64)),
        shape=X.shape,
    )
    a.write_h5ad(path)


def _write_label(obs, col: str, values, suffix: str) -> None:
    """Store an atlas label in an obs frame without clobbering the study's own.

    Always writes '<col><suffix>'. Also fills the plain column when the
    study does not already have one, so a study that was never annotated
    individually is still usable downstream.
    """
    obs[f"{col}{suffix}"] = values
    if col not in obs.columns:
        obs[col] = values


def propagate_labels(
    master_h5ad_path: Path,
    study_h5ad_path: Path,
    label_cols: Sequence[str] = ('leiden', 'cell_type'),
    suffix: str = '_atlas',
    progress_callback: Optional[Callable[[str], None]] = None,
) -> int:
    """Write master-derived labels into a per-study h5ad, alongside its own.

    Labels land in suffixed columns -- ``cell_type_atlas``,
    ``leiden_atlas`` -- rather than over ``cell_type`` and ``leiden``.
    A study's own annotation is a result in its own right: comparing
    per-study DE under independent labels against the same DE under
    shared atlas labels is a question you can only ask while both
    survive, and overwriting made it unanswerable.

    The unsuffixed columns are written only when the study has none, so
    a study that was never annotated individually still comes out with a
    usable ``cell_type`` for the steps that assume one.

    Fast path (every per-study cell is in the master):
      Read master.obs only, join by obs_names suffix, write labels back.
      No matrix loaded for the master; per-study h5ad opened in backed
      mode and rewritten only after the obs columns are set.

    Deliberate-omission path (every per-study cell absent from the
    master has ``_role == 'exclude'``):
      Barcode-join the cells that are present; the excluded ones are left
      unlabelled rather than projected, since they were left out on
      purpose.

    Fallback path (per-study has cells not present in the master for some
    other reason, e.g. the master was built with a cap):
      sc.tl.ingest to project labels via kNN in master's PCA space.

    Parameters
    ----------
    master_h5ad_path : Path
        Master h5ad. Must have the requested label columns in obs.
    study_h5ad_path : Path
        Per-study h5ad. Rewritten with new obs columns.
    label_cols : sequence of str
        obs columns to transfer.

    Returns
    -------
    int
        Number of cells in the study after label transfer.
    """
    accession = study_h5ad_path.parent.parent.name

    # Read master.obs only (no X) -- a few MB max.
    if progress_callback:
        progress_callback(f"Reading master labels for {accession}...")
    # obs only: a backed read would load the master's counts layer.
    from kosmic.scrna.load.h5ad_meta import read_obs, write_obs
    _master_obs_full = read_obs(master_h5ad_path)
    missing = [c for c in label_cols if c not in _master_obs_full.columns]
    if missing:
        raise ValueError(
            f"Master obs is missing {missing}. Annotate before propagating.")
    master_obs = _master_obs_full[list(label_cols)].copy()
    master_obs_names = _master_obs_full.index.to_numpy()
    del _master_obs_full

    # Cells in the master are indexed as '<original>-<accession>'
    # (concat_studies uses index_unique='-'). Recover the original
    # per-study cell names by stripping that. Named for what it is, not
    # 'suffix' -- that is the caller's column suffix, and shadowing it
    # here silently produced 'cell_type-S1' instead of 'cell_type_atlas'.
    barcode_suffix = f"-{accession}"
    master_for_study = master_obs.copy()
    master_for_study.index = [
        n[:-len(barcode_suffix)] if n.endswith(barcode_suffix) else n
        for n in master_obs_names
    ]

    # Only obs is needed to join labels by barcode; the matrix stays on
    # disk. A 590,000-cell study is a few hundred MB this way against
    # 13 GB loaded whole -- which it was, until this was written.
    if progress_callback:
        progress_callback(f"Reading cells of {study_h5ad_path.name}...")
    study_obs = read_obs(study_h5ad_path)
    study_n_obs = len(study_obs)
    study_obs_names = study_obs.index

    study_cells = set(study_obs_names)
    master_cells_for_study = set(master_for_study.index)
    coverage = len(master_cells_for_study & study_cells) / max(1, len(study_cells))

    # A study can be short of the master for two very different reasons,
    # and they want opposite handling. A cap leaves an arbitrary
    # subsample, so the absent cells should be projected. Cells left out
    # because their role is 'exclude' were deliberately omitted, and
    # projecting them would invent labels -- from a space they were never
    # part of -- for cells you have already said are not part of the
    # analysis. Barcode-join what is there and leave the rest unlabelled.
    uncovered = study_obs_names[~study_obs_names.isin(master_cells_for_study)]
    if len(uncovered) and '_role' in study_obs.columns:
        roles = study_obs.loc[uncovered, '_role'].astype(str).str.lower()
        n_excluded = int((roles == 'exclude').sum())
    else:
        n_excluded = 0
    # Included cells that are not in the master: either a cap left them
    # out (an arbitrary subsample, worth projecting) or the master's own
    # QC removed a handful (not worth projecting -- and on the DCM atlas
    # 390 such cells once sent a 590,000-cell study through sc.tl.ingest
    # for two hours). Below the threshold they stay unlabelled.
    leftover = uncovered[~uncovered.isin(
        study_obs_names[(study_obs['_role'].astype(str).str.lower() == 'exclude')]
        if '_role' in study_obs.columns else [])]
    n_included = study_n_obs - (
        int((study_obs['_role'].astype(str).str.lower() == 'exclude').sum())
        if '_role' in study_obs.columns else 0)
    project = len(leftover) > PROJECT_MIN_FRACTION * max(1, n_included)

    # Every cell that is in the master gets its label by barcode.
    if progress_callback:
        msg = f"Joining labels by barcode ({coverage*100:.1f}% of cells are in the master"
        if n_excluded:
            msg += f"; {n_excluded:,} excluded cell(s) stay unlabelled"
        if len(leftover) and not project:
            msg += (f"; {len(leftover):,} included cell(s) are not in the master "
                    f"-- removed by its QC -- and stay unlabelled")
        progress_callback(msg + ")...")
    joined = master_for_study.reindex(study_obs_names)
    for col in label_cols:
        _write_label(study_obs, col, joined[col].values, suffix)

    if project:
        # A cap was applied: project the absent included cells, and only
        # those, into the master's PCA space by kNN (sc.tl.ingest).
        if progress_callback:
            progress_callback(
                f"{len(leftover):,} included cell(s) are not in the master "
                f"(a cap was used); projecting them by kNN in the master's "
                f"PCA space...")
        import scanpy as sc
        study = ad.read_h5ad(study_h5ad_path)
        master = ad.read_h5ad(master_h5ad_path)
        if 'X_pca' not in master.obsm or 'PCs' not in master.varm:
            raise ValueError(
                "Master h5ad has no PCA (obsm['X_pca'] and varm['PCs']). "
                "Run PCA on the master before propagating with a subsample cap.")
        common = master.var_names.intersection(study.var_names)
        if len(common) < 100:
            raise ValueError(
                f"Only {len(common)} shared genes between master and study; "
                f"refusing to project.")
        common_list = list(common)
        master_sub = master[:, common_list].copy()
        del master
        study_sub = study[leftover, common_list].copy()
        del study
        sc.tl.ingest(study_sub, master_sub, obs=list(label_cols))
        del master_sub
        for col in label_cols:
            if col in study_sub.obs.columns:
                projected = joined[col].astype(object).copy()
                projected.loc[study_sub.obs_names] = study_sub.obs[col].astype(object).values
                _write_label(study_obs, col, projected.values, suffix)

    # Write obs back in place: the matrix and layers are untouched, so the
    # file is not rewritten.
    if progress_callback:
        progress_callback(f"Writing labels into {study_h5ad_path.name}...")
    write_obs(study_h5ad_path, study_obs)
    return int(study_n_obs)


__all__ = ['concat_studies', 'propagate_labels']
