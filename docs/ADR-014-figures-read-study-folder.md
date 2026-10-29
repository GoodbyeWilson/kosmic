# ADR-014: Figures read a study's saved data when no workspace holds it

- Status: Proposed (2026-10-29). Agreed by Kirk; awaiting Calum.
  Extends ADR-009 (Figure Export) and ADR-010 (figures read the
  meta-analysis); keeps ADR-001 and ADR-003.

## Problem

Figure Export drew every figure from what a workspace held in memory:
the scRNA figures from the scRNA workspace's dataset, the DE figures
from the DE workspace's attributes, the meta-analysis figures from
`figure_inputs()` (ADR-010). One dataset is kept in memory at a time
(ADR-003), so:

1. After a study is analysed from QC to DE, the scRNA figures were
   unavailable: DE read its own copy of the counts and the scRNA
   workspace released the study. Opening scRNA again reloaded it and
   released DE's copy, which made the Gene Heatmap unavailable; going
   back to DE cleared its results from memory.
2. Figures for any study other than the one the workspaces held could not
   be made at all, which includes the per-study figures once an atlas
   has been built.
3. A user who analyses on a workstation and makes figures on a laptop had
   to run the analysis again on the laptop.

Each step of the analysis already saves its output to the study's
folder: the processed h5ad is the source of truth (ADR-001) and the DE
steps write their tables under `results/de_analysis/`; meta-analysis
runs are saved under `meta_analysis/`.

## Decision

**A figure reads the data it draws from the workspace that holds the
study when one does, and otherwise from the study's folder.** Nothing is
recomputed and no matrix is loaded whole.

1. The study's processed file is read in parts
   (`kosmic/scrna/load/h5ad_meta.py`): `read_figure_meta` reads obs, var,
   uns and the drawn embeddings, and no matrix; `read_genes` reads a few
   genes for every cell and `read_cells` every gene for some cells, a
   chunk of the matrix at a time. Log-normalised values are produced by
   the same rule as `kosmic.scrna.counts`, and the results equal the
   values from the loaded study bit for bit (tested).
2. The DE and meta-analysis outputs are read from the files the steps
   write (`kosmic/de/saved.py`, `kosmic/meta_analysis/saved.py`). Paths
   are rebuilt from the study folder and the file names, never taken from
   inside a saved file, so a project folder copied to another computer
   works.
3. What a figure needs but no step saved is saved from now on, next to
   the step's other output: each gene DE run's settings (count source,
   gene sets, columns) as a `.json` record, over-representation results,
   and the running-enrichment curves of a library GSEA run. Runs made
   before this change cannot be recovered; their figures say which step
   to run again.
4. Figure Export has a study of its own: the application's active study
   by default, or the study chosen with **Open in Figures** on the Project
   page, which loads and releases nothing.

The figures read through one source object in the Figure Export
workspace; the figure pages do not know which source served them.

## Why this keeps ADR-001 and ADR-003

The file and the saved tables are the record the analysis already
treats as authoritative (ADR-001); figures now read the same record.
Figure Export holds obs, var, the embeddings and small tables, never a
matrix, so the one-dataset rule (ADR-003) still holds. Measured on the
DCM atlas (1.17 million cells, 15 GB file) with nothing else loaded:

| Read | Time | Peak memory |
|---|---|---|
| obs, var, uns, embeddings | 1.6 s | 0.64 GB |
| 69 marker genes, every cell (256 MB chunks) | 9.6 s | 1.0 GB |
| every gene, 150,000 sampled cells | 8.6 s | 3.8 GB |

The last is the size of its result, the same object the Marker Score
UMAPs page builds from a loaded study.

## Consequences

- Every figure an analysis produced can be exported without opening the
  workspaces in a particular order, for any study in the project, and on
  another computer from a copy of the project folder (`raw_data/` is not
  needed).
- A study held in a workspace takes precedence over its file, so unsaved
  state such as an open subclustering is what the figure shows.
- A figure from saved results can be older than the study's file (cell
  types relabelled after DE was run); the page says so.
- Reading a file can collide with a save of it on Windows, where a file
  open for reading cannot be replaced. Readers and the save path share a
  lock.
- Follow-ups not decided here: exporting figures from the command line
  (ADR-004), and two computers writing to one project folder at once.
