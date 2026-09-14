# Supervised ambient-RNA decontamination (DecontX).
#
# DecontX attributes each cell's 'off-profile' transcripts to the other
# cell types present and subtracts them. The correction is only valid on
# the FULL multi-cell-type dataset: the abundant source cells (e.g.
# cardiomyocytes) must be present for their ambient signal to be
# recognised and removed from the cells they contaminate (e.g. endothelium).
# Run it on a single-cell-type subset and there is no contamination source
# to attribute to, so it removes nothing useful.
#
# Placement is therefore post-annotation, pre-subset. Three prerequisites
# enforce this (see 'check_decontx_prerequisites'):
#   1. cell-type labels exist (annotation done),
#   2. at least two cell types present (not already subset),
#   3. a sample/donor column is set (ambient is per-sample).
#
# Supervised mode passes the annotation as 'cluster_key'. Each sample is
# decontaminated on its own and its dense result is sparsified before the
# next sample, so peak memory stays at one sample rather than a dense
# '(n_cells, n_genes)' matrix (which the decontx package allocates whole --
# tens of GB on a full nuclei dataset). Decontaminated counts are written
# to a new sparse layer; raw counts in 'X' / 'layers['counts']' are never
# overwritten, so the step is reversible and DE can compare raw vs decont.
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from scipy import sparse as sp

# Native output names written by the decontx package.
DECONTX_LAYER = 'decontX_counts'
DECONTX_CONTAMINATION = 'decontX_contamination'

CELL_TYPE_COL = 'cell_type'
MIN_CELL_TYPES = 2


def check_decontx_prerequisites(
    adata,
    sample_col: Optional[str],
    cell_type_col: str = CELL_TYPE_COL,
    min_cell_types: int = MIN_CELL_TYPES,
) -> Tuple[bool, str]:
    """Return (ok, reason) for whether DecontX may run on 'adata'.

    'ok' is True only when all three gates pass. When False, 'reason'
    is a one-line, user-facing explanation the GUI shows next to the
    disabled action. When True, 'reason' is an empty string.
    """
    if adata is None:
        return False, "Load and process a dataset first."

    if cell_type_col not in adata.obs.columns:
        return False, "Annotate cell types first (no 'cell_type' column)."
    labels = adata.obs[cell_type_col]
    if labels.isna().all():
        return False, "Annotate cell types first (cell_type is empty)."

    n_types = labels.dropna().astype(str).nunique()
    if n_types < min_cell_types:
        return False, (
            "DecontX needs the full multi-cell-type dataset -- run it "
            "before Subset (only one cell type present).")

    if not sample_col:
        return False, "Select the sample/donor column -- DecontX decontaminates per sample."
    if sample_col not in adata.obs.columns:
        return False, f"Sample column '{sample_col}' is not in the data."
    if adata.obs[sample_col].dropna().nunique() < 1:
        return False, f"Sample column '{sample_col}' has no values."

    return True, ""


def _raw_counts_matrix(adata):
    """Return a raw integer count matrix, preferring 'layers['counts']'.

    DecontX reads counts from '.X', but KOSMIC normalises 'X' in place at
    the QC step and stashes raw counts in 'layers['counts']'. We resolve
    the count source the same way SoupX does and hand DecontX a copy with
    raw counts in 'X', so we never depend on whether Normalize has run.
    """
    if 'counts' in adata.layers:
        X = adata.layers['counts']
    else:
        X = adata.X
    max_val = float(X.data.max()) if sp.issparse(X) else float(np.asarray(X).max())
    if max_val < 20:
        raise ValueError(
            "DecontX needs raw integer counts but the data looks "
            "log-normalised (max < 20) and no raw 'counts' layer was "
            "found. Run DecontX before Normalize, or keep raw counts in "
            "layers['counts'].")
    return X


def run_decontx(
    adata,
    sample_col: str,
    cell_type_col: str = CELL_TYPE_COL,
    *,
    seed: int = 12345,
    max_iter: int = 500,
    round_to_int: bool = True,
    progress_callback=None,
    progress_pct_callback=None,
) -> Tuple[object, dict]:
    """Run supervised, per-sample DecontX on the full annotated dataset.

    Writes 'layers['decontX_counts']' (sparse), 'obs['decontX_contamination']'
    and 'uns['decontX']' onto 'adata' in place. Raw counts in 'X' and
    'layers['counts']' are left untouched. Returns '(adata, info)'.

    Each sample is decontaminated on its own (statistically identical to the
    package's 'batch_key' path, since batches are run independently) and its
    dense result is sparsified immediately. This bounds peak memory to a
    single sample instead of a dense '(n_cells, n_genes)' matrix, which the
    decontx package would otherwise allocate whole -- tens of GB, and a hard
    out-of-memory crash on a full-size dataset. Samples with fewer than two
    cell types cannot be decontaminated (no contamination source) and keep
    their raw counts.

    'round_to_int' rounds corrected counts back to integers (stochastic
    rounding, seeded) so the layer feeds count-based DE directly.

    Raises ValueError if the prerequisites are not met.
    """
    def _emit(msg: str):
        if progress_callback:
            progress_callback(msg)

    def _emit_pct(pct: float):
        if progress_pct_callback:
            progress_pct_callback(int(round(pct)))

    ok, reason = check_decontx_prerequisites(adata, sample_col, cell_type_col)
    if not ok:
        raise ValueError(reason)

    import anndata as ad
    import pandas as pd
    try:
        import decontx
    except ImportError as e:
        raise RuntimeError(
            "The 'decontx' package is not installed. Install it to use "
            "DecontX decontamination (pip install decontx).") from e

    _emit("Preparing raw counts for DecontX...")
    counts = _raw_counts_matrix(adata)
    counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)
    n_cells, n_genes = counts.shape

    cell_types = adata.obs[cell_type_col].astype(str).values
    samples = adata.obs[sample_col].astype(str).values
    unique_samples = list(pd.unique(samples))
    rng = np.random.default_rng(seed)

    contamination = np.zeros(n_cells, dtype=float)
    row_blocks = []          # (original_indices, sparse block) in processing order
    skipped = []

    n_unique = len(unique_samples)
    for i, s in enumerate(unique_samples):
        idx = np.where(samples == s)[0]
        z = cell_types[idx]
        _emit(f"DecontX: sample {i + 1}/{n_unique} "
              f"'{s}' ({idx.size:,} cells)...")
        # Reserve the last 5% for the assembly/write step.
        _emit_pct(95.0 * i / n_unique)

        if np.unique(z).size < 2:
            # No contamination source within this sample -- keep raw counts.
            row_blocks.append((idx, counts[idx].copy()))
            skipped.append(s)
            continue

        sub = ad.AnnData(
            X=counts[idx].copy(),
            obs=pd.DataFrame({cell_type_col: z}, index=idx.astype(str)),
        )
        decontx.decontx(
            sub, cluster_key=cell_type_col, max_iter=max_iter,
            seed=seed, copy=False, verbose=False)

        dec = sub.layers[DECONTX_LAYER]
        dec = dec.toarray() if sp.issparse(dec) else np.asarray(dec)
        if round_to_int:
            dec = _stochastic_round(dec, rng)
        row_blocks.append((idx, sp.csr_matrix(dec)))
        contamination[idx] = np.asarray(sub.obs[DECONTX_CONTAMINATION].values, dtype=float)
        del sub, dec

    _emit("Assembling decontaminated layer...")
    _emit_pct(96.0)
    proc_idx = np.concatenate([idx for idx, _ in row_blocks])
    stacked = sp.vstack([blk for _, blk in row_blocks]).tocsr()
    inverse = np.empty(n_cells, dtype=np.int64)
    inverse[proc_idx] = np.arange(proc_idx.size)
    adata.layers[DECONTX_LAYER] = stacked[inverse]
    adata.obs[DECONTX_CONTAMINATION] = contamination

    n_samples = len(unique_samples)
    n_types = int(np.unique(cell_types).size)
    adata.uns['decontX'] = {
        'sample_col': sample_col,
        'cell_type_col': cell_type_col,
        'n_samples': n_samples,
        'n_cell_types': n_types,
        'n_samples_skipped': len(skipped),
        'rounded': bool(round_to_int),
    }

    mean_contam = float(np.nanmean(contamination)) if contamination.size else float('nan')
    skip_note = f"; {len(skipped)} single-cell-type sample(s) kept raw" if skipped else ""
    info = {
        'mean_contamination': mean_contam,
        'n_samples': n_samples,
        'n_cell_types': n_types,
        'n_samples_skipped': len(skipped),
        'layer': DECONTX_LAYER,
        'message': (
            f"DecontX complete -- mean estimated contamination "
            f"{mean_contam:.1%} across {n_samples} sample(s){skip_note}; "
            f"decontaminated counts in layers['{DECONTX_LAYER}']"),
    }
    _emit_pct(100.0)
    _emit(info['message'])
    return adata, info


def _stochastic_round(arr: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Round to integers stochastically so column/row totals are preserved."""
    floor = np.floor(arr)
    return (floor + (rng.random(arr.shape) < (arr - floor))).astype(np.float32)
