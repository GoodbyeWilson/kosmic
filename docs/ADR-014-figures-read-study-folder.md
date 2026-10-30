# ADR-014: Figures read a study's saved data when no workspace holds it

- Status: Approved (2026-10-30). Proposed 2026-10-29.
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
   page, which loads and releases nothing and does not change the active
   study.
5. When the study has several saved DE runs (one per cell type and
   method), the user chooses one with a Run selector; the DE workspace's
   run, when it holds one for the study, is listed first. The Meta
   figures have their own selector over the runs saved under
   `meta_analysis/`, after the Meta-Analysis workspace's run. A saved
   meta-analysis run's forest plots load each pooled study's table with
   the workspace's own loader (`load_study_de`), so the numbers are the
   same as in the workspace.
6. A value that a run made before this change did not record is taken
   from the next record that holds it: the disease and control labels
   from the run's record, else the manifest's role map, else the
   dataset's `uns['role_map']`, else the role names when the dataset has
   a `_role` column (an atlas has no manifest).

The figures read through one source object in the Figure Export
workspace (`StudySource`); the figure pages do not know which source
served them. Whether a figure is available is decided from which files
exist, so filling the figure list reads no table.

## Why this keeps ADR-001 and ADR-003

The file and the saved tables are the record the analysis already
treats as authoritative (ADR-001); figures now read the same record.
Figure Export holds obs, var, the embeddings and small tables, never a
whole matrix, so the one-dataset rule (ADR-003) still holds. It keeps
the last four reads of a few genes for every cell, so that changing a
figure's controls does not read the file again; on the atlas the Gene
Heatmap's read of 707 genes is the largest of these (about 1.2 GB). The
sample of every gene that the Marker Score UMAPs read is not kept, only
the scores computed from it.

Measured on the DCM atlas (`_master`: 1.17 million cells, 15 GB file)
with nothing loaded in any workspace, every figure drawn as the
application draws it, in one process limited to 16 GB of committed
memory by a Windows job object (2026-10-30). Peak memory is the
process's private memory, sampled every 20 ms; 1.6 GB of it is Python,
Qt and the libraries before anything is read.

| Step | Time | Peak memory |
|---|---|---|
| obs, var, uns, embeddings | 1.7 s | 2.1 GB |
| 69 marker genes, every cell (256 MB chunks) | 9.5 s | 2.7 GB |
| Gene Heatmap counts: 707 genes of 5 Hallmark sets, every cell | 7.9 s | 3.9 GB |
| Marker Score UMAPs: every gene of 150,000 sampled cells, scored and drawn | 28.1 s | 6.6 GB |
| each other figure with data (UMAP, dot plot, composition, QC violins, variable genes, PCA elbow, volcano, library GSEA, leading edge, the gene meta-analysis figures) | 0.1-7 s | at most 3.8 GB |

Every figure for which the atlas has results drew within the limit. On
Chaffin (339,609 cells, 12 GB file) the metadata read took 2.1 s, the
marker genes 4.8 s, and the peak was 2.3 GB.

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
- Runs made before this change have no record of their gene set
  selection, over-representation results or GSEA curves; the figures
  that need them say which step to run again.
- A saved meta-analysis run's forest plots read each study's DE table as
  it is now. A study whose DE was run again after the meta-analysis
  shows its newer table, and nothing notes this.
- Only the pathway level of pathway scoring is read from a study's
  folder, not the pathway-group level.
- Follow-ups not decided here: exporting figures from the command line
  (ADR-004), and two computers writing to one project folder at once.
