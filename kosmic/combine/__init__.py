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


_OBS_KEEP_CANDIDATES = ('sample', 'condition', 'cell_type', '_role')


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

    try:
        for path in study_h5ad_paths:
            path = Path(path)
            accession = path.parent.parent.name
            if progress_callback:
                progress_callback(
                    f"Extracting raw counts: {accession} (from {path.parent.name}/)...")

            a = ad.read_h5ad(path)
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
        a = ad.read_h5ad(output_path, backed='r')
        try:
            progress_callback(
                f"Master: {a.n_obs:,} cells x {a.n_vars:,} genes")
        finally:
            a.file.close()

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


def _write_label(study, col: str, values, suffix: str) -> None:
    """Store an atlas label without clobbering the study's own.

    Always writes '<col><suffix>'. Also fills the plain column when the
    study does not already have one, so a study that was never annotated
    individually is still usable downstream.
    """
    study.obs[f"{col}{suffix}"] = values
    if col not in study.obs.columns:
        study.obs[col] = values


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
    master_backed = ad.read_h5ad(master_h5ad_path, backed='r')
    try:
        missing = [c for c in label_cols if c not in master_backed.obs.columns]
        if missing:
            raise ValueError(
                f"Master obs is missing {missing}. Annotate before propagating.")
        master_obs = master_backed.obs[list(label_cols)].copy()
        master_obs_names = master_backed.obs_names.to_numpy()
    finally:
        master_backed.file.close()

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

    if progress_callback:
        progress_callback(f"Loading {study_h5ad_path.name}...")
    study = ad.read_h5ad(study_h5ad_path)

    study_cells = set(study.obs_names)
    master_cells_for_study = set(master_for_study.index)
    coverage = len(master_cells_for_study & study_cells) / max(1, len(study_cells))

    # A study can be short of the master for two very different reasons,
    # and they want opposite handling. A cap leaves an arbitrary
    # subsample, so the absent cells should be projected. Cells left out
    # because their role is 'exclude' were deliberately omitted, and
    # projecting them would invent labels -- from a space they were never
    # part of -- for cells you have already said are not part of the
    # analysis. Barcode-join what is there and leave the rest unlabelled.
    uncovered = study.obs_names[~study.obs_names.isin(master_cells_for_study)]
    excluded_only = False
    if len(uncovered) and '_role' in study.obs.columns:
        roles = study.obs.loc[uncovered, '_role'].astype(str).str.lower()
        excluded_only = bool((roles == 'exclude').all())

    if coverage >= 0.999 or excluded_only:
        # Fast path: barcode join.
        if progress_callback:
            if excluded_only and coverage < 0.999:
                progress_callback(
                    f"Joining labels by barcode; {len(uncovered):,} excluded "
                    f"cell(s) are not in the master and stay unlabelled...")
            else:
                progress_callback(
                    f"Joining labels by barcode ({coverage*100:.1f}% match)...")
        joined = master_for_study.reindex(study.obs_names)
        for col in label_cols:
            _write_label(study, col, joined[col].values, suffix)
    else:
        # Fallback: sc.tl.ingest. Only triggered when a cap was applied.
        if progress_callback:
            progress_callback(
                f"Master covers {coverage*100:.0f}% of cells; using sc.tl.ingest...")
        import scanpy as sc
        master = ad.read_h5ad(master_h5ad_path)
        if 'X_pca' not in master.obsm:
            raise ValueError(
                "Master h5ad has no 'X_pca' in obsm. Run PCA on the master "
                "before propagating with a subsample cap.")
        common = master.var_names.intersection(study.var_names)
        if len(common) < 100:
            raise ValueError(
                f"Only {len(common)} shared genes between master and study; "
                f"refusing to project.")
        common_list = list(common)
        master_sub = master[:, common_list].copy()
        study_sub = study[:, common_list].copy()
        sc.tl.ingest(study_sub, master_sub, obs=list(label_cols))
        for col in label_cols:
            if col in study_sub.obs.columns:
                _write_label(study, col, study_sub.obs[col].values, suffix)

    if progress_callback:
        progress_callback(f"Writing {study_h5ad_path.name}...")
    study.write_h5ad(study_h5ad_path)
    return int(study.n_obs)


__all__ = ['concat_studies', 'propagate_labels']
