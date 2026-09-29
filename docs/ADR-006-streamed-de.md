# ADR-006: The DE workspace streams counts from the processed file

- Status: Proposed (2026-09-18). Awaiting Calum and Kirk. Builds on the
  block-wise pseudobulk of #55; independent of ADR-003 and ADR-005,
  and the block reader it adds is the one ADR-005 needs. Issue #56.

## Problem

The DE workspace loads a study's counts and keeps them resident for
the session: obs, var and the count matrix as `X` (ADR-003 part 1
removed the normalised matrix, the embeddings and the graph from this
load). On the DCM project that is 14 GB for the atlas and 8 GB for
Reichart-LV. Per-cell-type DE then copies each cell type out of the
resident matrix: 5.3 GB for the atlas's cardiomyocytes. Measured after
#55, the atlas cardiomyocyte run peaks at 21.2 GB. Kirk's laptop has
16 GB.

Nothing in the pseudobulk workflow needs the matrix resident. What
the workspace reads, by feature:

| feature | reads |
|---|---|
| Setup page: columns, donor counts, condition roles | `obs` only |
| pseudobulk DE, whole study or per cell type | one pass over the counts |
| detection pre-filter, expression rates | one pass over the counts |
| gene-expression and pathway-gene plots | one gene's column |
| cell-level Wilcoxon; pathway scoring (scoring mode) | the matrix, per cell |

Since #55 the three passes are already block-wise: `_indicator_sums`
multiplies the indicator matrix by blocks of 25,000 cells;
`pct_expressing_per_var` counts non-zeros in blocks of 10,000. They
take their blocks from memory only because that is where the matrix
is.

## Decision

**The DE workspace holds a study's `obs`, `var` and `uns`, and reads
the counts from the processed file in row blocks when a step needs
them. The count matrix is never resident in the DE workspace unless a
per-cell step asks for it.**

1. Loading a study into DE reads `obs`, `var` and `uns` through
   `kosmic/scrna/load/h5ad_meta.py` (megabytes). `current_adata` is
   that AnnData with no `X`; `h5ad_path` names the file, as now. The
   Setup page, the covariate list and the gene universe work from it
   unchanged.
2. `h5ad_meta.py` gains a row-block reader: it yields CSR blocks of
   the count source (`layers['counts']`, else `X`) for a requested
   range of cells, optionally masked to a cell selection, reading each
   block's `data` and `indices` once from the file. A block is the only
   count data in memory at a time.
3. `run_de_pipeline` and `run_de_by_cell_type` take a *count source*:
   either an in-memory matrix (a study loaded by another path, the
   tests) or the block reader. Pseudobulk sums, the detection
   pre-filter and the expression rates consume blocks from either.
   Per cell type, the cell mask is applied to each block as it is
   read; the subset copy goes.
4. The gene-expression and pathway-gene plots read the one column they
   need by a pass over the file, and cache the last few genes. On the
   atlas that is a scan of 14 GB, tens of seconds, for an occasional
   plot; on a study it is seconds.
5. Cell-level Wilcoxon and pathway scoring need every cell's counts.
   In the streamed workspace they load the matrix when the step is
   run, after telling the user its size, and release it when the step
   finishes. They are not part of the pseudobulk workflow and stay
   opt-in.
6. The one-dataset rule (ADR-003 part 1) is unchanged: loading a study
   into DE still releases what scRNA holds, and a DE workspace holding
   obs and var is what "loaded" means there.

A run on a block reader and a run on the same matrix in memory must
give identical output; that equality is a test, on a study small
enough to hold both.

## Consequences

- Atlas DE on the DCM project: about 1-2 GB in the workspace, from
  21.2 GB today (30.3 before #55). Per-study DE the same. The DE step
  of the paper's worked example runs on a 16 GB laptop, and on the
  atlas the "Run DE Analysis" over all cells costs the same as one
  cell type.
- Loading a study into DE becomes instant, since only obs and var are
  read.
- The block reader is the chunked pass over the file that ADR-005
  (sketch mode) specifies for pseudobulk and annotation means; it is
  written once, here.
- The DE data-contract paragraph of CODEBASE.md changes: "the DE
  workspace holds obs and var; counts are streamed from the processed
  file", replacing "DE loads obs, var and the counts".
- The two per-cell steps get a load prompt they did not have; a user
  of scoring mode on a large study sees where the memory goes.
- Effort: about two days. `kosmic/de/de_analysis.py` (the count-source
  abstraction over what #55 made block-wise), `kosmic/de/batch.py`,
  `h5ad_meta.py` (the reader), the Setup page loader, the two plot call
  sites, tests, CODEBASE.md and the DE help pages. It changes what
  `DEWorkspace.current_adata` is, so it needs both maintainers'
  agreement before it starts.
