# ADR-015: DE outputs by counts source and cell type, with the significance rule recorded per run

- Status: Accepted (Kirk, 2026-11-09)

## Problem

A study's DE outputs are written to one flat folder,
`results/de_analysis/statistics/`. A per-cell-type run of one study on the
DCM project writes 49 files there: for each of 12 cell types the DE table,
its `.json` record, the pseudobulk counts and the significant genes, plus
one `{accession}_by_cell_type_summary.csv`. Three problems follow.

1. **The folder cannot hold two count sources.** A study can be analysed
   on its original counts or on DecontX-corrected counts (ADR-011). The
   file names do not say which, so a second run on the other counts
   replaces the first. Comparing the two meant keeping two copies of the
   whole project. The same holds for the meta-analysis outputs in
   `meta_analysis/<selection>/` (ADR-007).
2. **The batch summary is incomplete and not portable.** It records the
   number of donors but not per condition, does not say which counts or
   significance rule were used, and stores an absolute `de_path`, which is
   wrong once the project is copied to another computer.
3. **The significance rule is a hidden project setting.** `config.toml`
   sets `lfc_threshold = 0.25`. The batch summary's `n_significant`, every
   `*_significant_genes.csv`, the Gene DE page's counts and top-DE list,
   the input gene list of over-representation analysis, and the Methods
   text all apply it, while Figure Export and the meta-analysis use
   FDR < 0.05 alone. The same run therefore reports different numbers in
   different places (Chaffin cardiomyocytes, October run: 4,135 genes with
   the gate, 4,317 without), and the Methods text can describe a threshold
   the analysis did not use.

## Decision

**A DE run's outputs are kept in a folder named after the counts it used
and the cell type it tested, and each run records the significance rule
it was made with; every place that counts or lists significant genes uses
that recorded rule unless the user changes it on purpose.**

### 1. Layout

```
results/de_analysis/statistics/
    <counts>/
        {study}_by_cell_type_summary.csv
        {study}_DE_{method}.csv                       (+ .json, pseudobulk,
                                                      significant genes)
                                                      a whole-study run
        <cell type>/
            {study}_{cell type}_DE_{method}.csv
            {study}_{cell type}_DE_{method}.json
            {study}_{cell type}_pseudobulk.csv
            {study}_{cell type}_significant_genes.csv
```

- `<counts>` is `original_counts` for the study's original counts and the
  layer name for a corrected layer (`decontX_counts`), the value
  `kosmic.scrna.counts.loaded_count_source` reports (`raw` is written as
  `original_counts`).
- `<cell type>` is the cell type with spaces as underscores, as in the
  file names today.
- **File names do not change.** A file copied out of its folder still
  names its study and cell type, and the name parsing in
  `kosmic/de/saved.py` and `kosmic/meta_analysis/io.py` keeps working.
- One helper in `kosmic/paths.py` builds these paths; the batch writer and
  every reader use it, so the layout is defined in one place.

### 2. What identifies a run

A run is identified by **cell type, method and counts source**. The DE run
picker and Figure Export's Run list name all three
("Ventricular Cardiomyocyte · DESeq2 · DecontX counts").

### 3. Meta-analysis outputs

- The meta-analysis pools runs made on **one counts source** in every
  selected study, and states which. A selection whose studies' runs used
  different counts is refused, naming the studies. In the application a
  **Counts** box on Select Studies chooses the counts, and only runs on
  them are listed, so a mixed selection cannot be made there.
- Its outputs go to `meta_analysis/<counts>/<selection>/`, extending
  ADR-007's per-selection folders, so a DecontX and an original-counts
  meta-analysis of the same cell type are kept side by side.
- Imported external tables have no known counts source; a selection that
  includes them is written under `meta_analysis/unspecified_counts/`.

### 4. The batch summary

One row per cell type, in this order:

| Columns | Meaning |
|---|---|
| `cell_type`, `status`, `message` | what ran; why a run was skipped or failed |
| `method`, `counts`, `design`, `min_cells` | how it was run |
| `n_donors_control`, `n_donors_disease`, `n_donors_excluded`, `n_donors_without_cells` | donors per arm, donors left out for having fewer than `min_cells` cells, and donors in the contrast with none of the cell type's cells (added 2026-11-10), so the three add up to the study's donors |
| `n_cells`, `median_cells_per_donor` | cells tested |
| `n_genes_tested`, `n_significant`, `n_up`, `n_down` | results under the run's significance rule |
| `fdr`, `lfc_threshold` | that rule (`lfc_threshold` 0 when no fold-change gate was set) |
| `saved` | when the run was written; matches its `.json` record |

The summary stores no path (the absolute `de_path` it replaces was wrong
on another computer): each run is in the folder beside it named after its
cell type. A `folder` column naming that folder was dropped on 2026-11-10
as it repeated the cell type.

### 5. Significance is a setting of each run

- The Gene DE page and the batch DE run have a **Min |log2FC|** control
  beside the FDR, **Off by default**: significant means FDR < 0.05 alone.
- The run's `.json` record and its provenance stage store `fdr` and
  `lfc_threshold`.
- The run's `significant_genes.csv`, the summary counts, the Gene DE
  page's counts and top-DE list and the Methods text use the recorded
  rule. With no gate, the Methods text states FDR < 0.05 and no
  fold-change threshold.
- The Enrichment page's fold-change box and Figure Export's
  **Min |log2FC|** start from the run's recorded value. Changing them
  affects only that analysis or figure, which then says so.
- `lfc_threshold` is removed from `config.toml`; `fdr` stays there as the
  starting value of the FDR control.

## Compatibility

- **Readers accept both layouts.** A study written before this change has
  its files directly in `statistics/`; they are listed as runs with the
  counts source taken from the run's `.json` record, or shown as
  *counts not recorded* where there is none. Nothing is moved; existing
  projects (the DCM `Clean_Data_Run`) open and pool as before.
- **Old meta-analysis outputs** in `meta_analysis/<selection>/` stay
  readable by Figure Export's saved-run list (ADR-014).
- **Runs without a recorded significance rule** were made with the
  configuration value of their time. They are read with
  `lfc_threshold = 0.25`, the shipped default, and labelled *assumed*
  wherever the value is shown.
- A project mixing old and new layouts is valid; the run lists show each
  run once.

## Consequences

- Original-counts and DecontX analyses of one project can be run, kept and
  compared side by side, at the study and the meta-analysis level.
- The number of significant genes a run reports is the same on every page,
  in every file and in the Methods text, unless the user deliberately
  chooses a different gate for one figure or one enrichment run.
- About 27 modules read or write DE outputs and change with this ADR:
  `kosmic/paths.py`, `kosmic/de/batch.py`, `kosmic/de/de_analysis.py`,
  `kosmic/de/saved.py`, `kosmic/meta_analysis/io.py`,
  `kosmic/gui/meta_analysis/pseudobulk_discovery.py`, the DE workspace
  pages (setup, Gene DE, methods, enrichment), the Meta-Analysis workspace
  and its gene page, the Figure Export DE pages and their base, the home
  workspace's project overview and `interactive_top_de`. The Script
  Folder's `shared_universe.py` and the KOSMIC Web exporters read the old
  layout and are updated to read both.
- Paths get one level deeper; a folder of one cell type holds four files.
- Pathway scoring outputs (`results/de_analysis/pathway_scoring/`) keep
  their layout for now; moving them under `<counts>/<cell type>/` is a
  follow-up with the same rules.

## Stages

Each stage ends with `pytest tests/` and `ruff check kosmic/ main.py
tests/` passing.

1. **Paths and readers.** The path helper; `saved.py`, Meta discovery and
   Figure Export find runs in both layouts. Tests: a fixture with the old
   flat layout and one with the new layout list the same runs; the counts
   source is read from the folder or the record.
2. **Writers and summary.** The batch and Gene DE writers use the new
   layout; the new summary columns. Tests: files land in
   `<counts>/<cell type>/`; the summary's donors per arm and counts, with
   no stored paths; a copied project still lists its runs.
3. **Significance per run.** The Min |log2FC| controls, the record
   fields, and every consumer reading them; `lfc_threshold` removed from
   `config.toml`. Tests: Off gives FDR-only counts in the summary,
   `significant_genes.csv` and the Methods text; a gate of 0.5 is recorded
   and used; Figure Export's control starts from the run's value; an old
   run is read as 0.25, labelled assumed.
4. **Meta-analysis by counts.** Pooling one counts source; refusing a
   mixed selection; outputs under `meta_analysis/<counts>/<selection>/`.
   Tests: original and DecontX selections write to separate folders; a
   mixed selection is refused naming its studies.
5. **Documentation.** CODEBASE.md, the DE and Meta help pages and
   screenshots, `paths.py`'s layout docstring, and the Script Folder's
   `shared_universe.py`.
