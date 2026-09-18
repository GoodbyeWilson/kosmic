# ADR-003: One dataset in memory, and only the counts

- Status: Proposed (2026-09-18). Part 1 accepted and done (#46, #49).
  Part 2 awaiting Calum and Kirk.

## Problem

A study in KOSMIC is a large object, and the application held more of
it, and more copies of it, than any step needed. Measured on the DCM
project (five studies, 1.09 million included nuclei) in September 2026:

- After QC a study carries two full matrices: the counts in
  `layers['counts']` and the log-normalised copy in `X`. They are the
  same size. Reichart-LV is 16.7 GB loaded; the shared atlas 27 GB.
  The normalised matrix is entirely derivable from the counts.
- The DE workspace loaded a study whole -- both matrices, embeddings,
  graph -- for an analysis that reads obs and the counts (fixed, #49:
  it now loads counts only; the atlas costs 14 GB there).
- Switching study did not free the previous one: the workspace dropped
  its pointer but every tab kept a reference (fixed, #49: `release_dataset`
  clears every tab; the rule is applied across workspaces).
- `read_h5ad(backed='r')` keeps `X` on disk and loads every layer, so a
  "cheap" backed read of a study loaded its counts (fixed, #46: obs,
  shape and gene names are read through `h5ad_meta.py`).
- The cluster-step workers deep-copied the study, twice (fixed, #33);
  the annotation workers copied it once (fixed, #45); the QC step copied
  it for each threshold (fixed, #41).

With those fixes the DCM atlas clusters in about 35 GB and its DE in
about 18 GB. Kirk's laptop has 16 GB. A million-cell atlas, or the two
larger studies, do not run there, and the reason is the second matrix.

## Decision

**A study holds its counts, and only its counts. The log-normalised
matrix is a cache that no longer exists; each step derives the
normalised values it needs, for the genes it needs, from the counts.**

Concretely:

1. `layers['counts']` goes away as a separate copy. After QC, `X` holds
   the counts (integer-valued, as on import). The processed file holds
   one matrix. Files from earlier KOSMIC versions, with both matrices,
   are read as they are and written back in the new form on the next
   save.
2. Normalisation is a recorded choice, not a stored matrix. The QC tab's
   Normalise step writes the parameters to `uns['normalisation']`
   (`target_sum`, `log1p`, and the date) and a per-cell size factor to
   `obs['size_factor']`, computed from the full gene set at that moment.
   Nothing is materialised.
3. One function produces normalised values: `kosmic.scrna.counts.
   log_normalised(adata, genes=None, cells=None)`, reading the recorded
   parameters and the stored size factor. No caller passes its own
   target sum. The result is a matrix for the requested genes and cells
   only -- 2,000 HVG columns for PCA, the sampled cells for marker
   ranking, one gene for the UMAP colouring, per-cluster means for
   reference annotation.
4. The size factor is stored, not recomputed, so a later gene filter or
   subset does not change the normalised values of the cells that
   remain. This is what the cached matrix guaranteed by accident.
5. HVG selection uses the counts (seurat_v3), which it should have done
   all along; the current 'seurat' fallback on log data goes.

What a derived matrix equals is testable: a value produced by
`log_normalised` for a gene subset must equal the value in a fully
materialised matrix, bit for bit. That test is part of the change.

## Why this is safe

Normalisation is a deterministic per-element transform of the counts
given the size factor and the parameters. The three ways it could go
wrong -- different parameters in different places, size factors from a
gene subset, drift between materialised and derived values -- are each
closed by one of the rules above (one recorded parameter set; a stored
per-cell factor from the full gene set; one function with an equality
test).

Every reader of expression already goes through `kosmic.scrna.counts`
(`count_source`, `counts_adata`, `log_normalised`), so the change is
inside that module and in the few places that read `X` directly
(the UMAP gene colouring, reference annotation means, the marker
ranking, the QC histograms).

## Consequences

- Reichart-LV: 16.7 GB -> 8 GB loaded; the DCM atlas 27 GB -> 14 GB; the
  cluster step on the atlas about 20 GB instead of 35. Files and saves
  halve. Reichart-LV and Chaffin fit a 16 GB machine; the atlas is
  close (see ADR-005 for what makes it comfortable).
- The counts contract in CODEBASE.md changes: "X holds the counts;
  normalised values are derived", replacing "counts in layers['counts'],
  log-normalised X".
- `looks_log_normalised` and the "already normalised" guards become
  checks on `uns['normalisation']`, not on the values in `X`.
- Every step that normalises a private copy today (the PCA worker when
  X held counts) stops doing so.
- Files written by earlier versions keep working; a migration on save
  drops the cached `X` and keeps the counts. `raw_data/` copies are
  untouched.
- Estimated effort: two to three days including tests and the help
  pages. It touches what a study file is, so it needs both maintainers'
  agreement before it starts.
