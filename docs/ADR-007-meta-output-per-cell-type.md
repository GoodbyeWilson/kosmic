# ADR-007: Meta-analysis outputs are kept per cell type

- Status: Proposed (2026-09-24). Folder naming, mixed selections and
  external imports agreed by Kirk; awaiting Calum. Issue #60.

## Problem

Every meta-analysis run writes its outputs to the project's
`meta_analysis/` folder: `consensus_<methods>_<calibration>.csv/.json`,
`pathway_ma_<method>.csv/.json`, and the meta-analysis `provenance.json`
and `methods.md` that the Methods step renders. The file names depend on
the pooling method, not on what was pooled.

Per-cell-type DE (`kosmic.de.batch`) gives one result per study per cell
type, and the Select Studies step now pools one cell type at a time. On
the DCM project each run of the same method on another cell type
replaced the previous one: the cardiomyocyte REML table was overwritten
by the endothelial one, then by the fibroblast one, within an hour. The
provenance record, which keeps only the latest run of each stage, then
described only the last cell type.

## Decision

**Each meta-analysis selection writes to its own folder under
`meta_analysis/`, named after what was pooled.**

| Selected results | Folder |
|---|---|
| One cell type across studies | `meta_analysis/<cell type>/`, e.g. `meta_analysis/Endothelial_Cell/` |
| Whole-study results only | `meta_analysis/all_cells/` |
| More than one cell type | `meta_analysis/mixed/`, with a warning in the output panel |

1. The folder name is the cell type as the per-cell-type DE file names
   and the Select Studies box spell it, so the three match.
2. Externally imported results have no cell type and do not affect the
   name: adding one to an Endothelial_Cell selection still writes to
   `Endothelial_Cell/`. A selection of external results only writes to
   `all_cells/`.
3. Everything a selection produces goes to its folder: the consensus and
   pathway tables with their settings sidecars, and the meta-analysis
   `provenance.json` and `methods.md`. The Methods step renders the
   record of the current selection's folder.
4. Mixing cell types is allowed, because pooling different cell types
   can be deliberate, but it is named `mixed/` and reported.
5. Files already at the top of `meta_analysis/` from earlier versions
   are left in place. The Methods step falls back to the top-level
   record when the selection's folder has none, so existing projects
   still open.

The folder is derived from the selection in one place
(`kosmic.meta_analysis.io`) and the path is built by
`kosmic.paths.meta_output_dir`, which gains an optional selection
argument; every writer and the Methods step go through them.

## Consequences

- Running the same method on each cell type in turn keeps every result
  and its methods record. On the DCM project that is twelve folders
  after a full pass.
- Two runs on the same cell type with the same method still replace
  each other, as before: the folder separates cell types, not attempts.
- The provenance record of a folder describes that cell type's run
  only, which is what a methods section for that cell type needs.
- Code reading meta-analysis outputs from disk must take the selection
  into account. Today only the Methods step does; Validation and Figure
  Export use the results of the current run held in memory.
- CODEBASE.md's description of the project layout gains the per-selection
  folders.
