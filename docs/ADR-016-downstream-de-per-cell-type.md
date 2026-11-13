# ADR-016: Enrichment and pathway analyses per cell type, chosen in one bar for the DE workspace

- Status: Accepted (Kirk, 2026-11-12)
- Extends: ADR-015 (its follow-up on `pathway_scoring/`)

## Problem

A study can be analysed per cell type without making a subset per cell type:
the per-cell-type batch run writes one gene DE run per cell type (ADR-015).
The analyses that follow gene DE cannot be run per cell type in the same way.

1. **The cell type is chosen in one place only, and the choice is not
   shown elsewhere.** The cell types of a batch run are offered in a small
   *Cell type* box under the Gene DE volcano. Choosing one replaces the
   workspace's `de_results`, which the Enrichment page then reads. Enrichment
   therefore runs on the cell type last chosen on Gene DE, and the Enrichment
   page does not show which one it is.
2. **Every downstream output is written to one folder.** Over-representation,
   library GSEA, fgsea of the gene-set panel, pathway DE and the per-donor
   pathway scores are all saved in `results/de_analysis/pathway_scoring/`
   under fixed file names. A run for a second cell type, or on the other
   counts, replaces the first. Figure Export and the pathway-level
   meta-analysis read that folder and cannot tell the runs apart.
3. **Pathway DE scores the whole loaded dataset.** In Scoring mode, pathway
   scoring uses every cell of the study, whichever cell type's gene DE is
   shown, so a per-cell-type pathway analysis is not possible.
4. **Provenance does not record the cell type.** An enrichment stage replaces
   the previous entry with the same method and library, so the record of one
   cell type's run is lost when another cell type is run.

## Decision

**The DE workspace holds one selection of the gene DE run being analysed,
made in a bar across the top of the workspace. Every step after gene DE uses
that run, and saves its outputs in a folder for that run's counts and cell
type.**

### 1. The run bar

A bar above the DE workspace's pages, shown on Gene Differential
Expression, Pathway Differential Expression, Enrichment and Methods:

```
Results for  [ Ventricular Cardiomyocyte                ▾ ]
DESeq2 · DecontX counts · 18 control, 23 disease donors · 3,412 significant (FDR < 0.05)
```

- **What it lists.** The study's saved gene DE runs on the loaded counts and
  the chosen DE method: the whole-study run first, if there is one, then
  each cell type by name with its number of significant genes and donors.
  A cell type the batch run skipped, or whose run failed, is listed but
  cannot be chosen, with the reason from the batch summary.
- **When it is shown.** As soon as the study has a gene DE run, saved or
  just made. A study with only a whole-study run shows one entry, so every
  page shows which run its results belong to. The bar is hidden on Load
  Data, Choose Analysis and Select Gene Sets.
- **One selection for the workspace.** The selection is held by the
  workspace (`DEWorkspace.cell_type_runs` and the run shown), not by a page.
  Changing it changes the run on every page. The *Cell type* box under the
  Gene DE volcano is removed.
- **No results shown under another run's name.** When the selection
  changes, each page's results are cleared, then reloaded from the new run's
  folder (§2) where it holds saved results for that page.
- **The second line** gives the run's method, counts, donors per arm and
  number of significant genes under the run's recorded rule (ADR-015 §5).
  The Explorer pane's status line adds the cell type.
- **A new run.** A whole-study Gene DE run or a batch run refreshes the list
  and selects the run just made (the first cell type of a batch run).

The bar is a widget of the DE workspace (`kosmic/gui/de_analysis/`), placed
in the workspace's layout above its page stack. It is not a new page
archetype; the pages keep their archetypes.

### 2. Layout of the downstream outputs

```
results/de_analysis/pathway_scoring/
    <counts>/
        enrichment_results.csv, gsea_library_results.csv (+ _details.npz),
        fgsea_results.csv, pathway_de_results.csv, parent_pathway_de_results.csv,
        pathway_scores_per_donor.csv, pathway_coverage.csv, ...
                                            a whole-study run
        <cell type>/
            the same files                  a run on one cell type
```

- `<counts>` and `<cell type>` are the folder names ADR-015 uses under
  `statistics/` (`original_counts`, `decontX_counts`; the cell type made
  file-safe by `slugify_cell_type`).
- **File names do not change.** The folder identifies the run.
- One helper in `kosmic/paths.py`, `de_downstream_dir(study_dir, counts,
  cell_type=None)`, builds the folder. The workers that write these files
  and every reader (`kosmic/de/saved.py`, Figure Export's `StudySource`,
  the pathway-level discovery in `kosmic/meta_analysis/io.py`) use it.

### 3. Pathway DE per cell type

When the run shown is a cell type's, pathway scoring uses that cell type's
cells only. The run's `.json` record names the column (`cell_type_column`)
and the value (`cell_type`). Donors are kept by the floors the gene DE run
recorded (`min_cells`, `min_counts`). The gene sets are those of the Select
Gene Sets step, as now. fgsea of the gene-set panel and library GSEA rank
the run's own genes, as now.

The per-gene expression views (spillover, the pathway gene plots, the Gene
DE heatmap and expression tabs) also use only the cell type's cells, so
they show the same cells as the DE results beside them.

### 4. Provenance and Methods

The enrichment and pathway stages record `cell_type` (None for a whole-study
run) and `counts`. A rerun replaces only the entry with the same method,
library, cell type and counts, so the entries of each cell type are kept.
The Methods page lists the downstream analyses per cell type.

### 5. Figure Export

The Enrichment, library GSEA, fgsea and pathway figures read the folder of
the run selected in Figure Export's existing Run list, which names cell
type, method and counts (ADR-015 §2). When the DE workspace holds the
study, Figure Export reads the run selected in the bar (ADR-014).

### 6. The pathway-level meta-analysis

Discovery lists a study's pathway DE results once per folder, with the
accession `{study}_{cell type}` and the counts taken from the folder, as for
gene DE runs. A pathway-level selection is then made per cell type and
counts source, by ADR-007 and ADR-015 §3.

## Compatibility

- **Files in the flat `pathway_scoring/` folder** of earlier versions are
  read as the whole-study run's, with counts and cell type marked *not
  recorded*, when the whole-study run's folder has none. Nothing is moved.
  The current DCM project (`New_QC_DecontX`) has no files there.
- The Meta-Analysis workspace's Enrichment step (the same page class, on the
  pooled gene list) has no run bar and saves where it does now, in the
  selection's `meta_analysis/` folder.

## Consequences

- Enrichment and pathway analyses can be made for each cell type of one
  study and kept side by side, on each counts source.
- Every downstream page shows, in the same place, which run its results
  belong to.
- Pathway DE of a cell type subsets the loaded matrix when the selection
  changes; the full matrix stays loaded (ADR-003).
- About 12 modules change: `kosmic/paths.py`, `kosmic/de/saved.py`,
  `kosmic/meta_analysis/io.py`, the DE workspace and its Gene DE,
  Pathway DE, Enrichment and Methods pages and workers, Figure Export's
  `StudySource` and its DE pages, and the help pages for those screens.

## Stages

Each stage ends with `pytest tests/` and `ruff check kosmic/ main.py
tests/` passing.

1. **Paths and readers.** `de_downstream_dir`; `saved.py`, `StudySource`
   and pathway discovery read the new folders and the flat one. Tests: a
   fixture with the flat layout and one with the new layout; two cell types
   with an `enrichment_results.csv` each are read separately; pathway
   discovery gives `{study}_{cell type}` accessions with their counts.
2. **The run bar and the shared selection.** The workspace holds the runs
   and the selection; the bar replaces the volcano box; Gene DE shows the
   selected run. Tests: changing the bar changes `de_results` and the
   significance rule; a skipped cell type is listed and cannot be chosen; a
   new batch run selects its first cell type.
3. **Enrichment and Pathway DE per run.** Writers save to the run's folder;
   pages clear and reload when the selection changes; pathway scoring and
   the expression views use the cell type's cells; provenance records cell
   type and counts. Tests: ORA for two cell types writes two files and
   neither replaces the other; pathway scoring of a cell type uses only its
   cells and the run's donors; enrichment stages for two cell types both
   stay in the provenance record.
4. **Figure Export, Methods and documentation.** Figure Export reads the
   selected run's folder; Methods lists runs per cell type; CODEBASE.md, the
   DE help pages and their screenshots.
