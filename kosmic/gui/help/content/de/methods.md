# Methods

The Methods step shows the values and parameters that describe the study's differential expression analysis, with references for the tools used, as text you can copy into a methods section. **Copy to clipboard** copies the whole sheet; **Refresh** rebuilds it after a setting changes.

## Which run it describes

The sheet describes the DE run chosen under **Results for** at the top of the page: the whole-study run or one cell type's run. Choosing another run there rebuilds the sheet for that run.

## Sections

| Section | What it gives |
|---|---|
| DATASET | The loaded file and accession, its cells and genes, and the cells and donors per condition |
| DE RUN | The run described: whole study or the cell type, the cell-type column it was taken from, the counts and the donors per condition |
| PREPROCESSING | The study's processing steps, from its provenance record |
| GENE DIFFERENTIAL EXPRESSION | The DE method and its settings |
| FILTERS | The gene and donor filters the DE used |
| RESULTS | Genes tested and significant under the significance rule, genome-wide and within the selected gene sets |
| PATHWAY DIFFERENTIAL EXPRESSION | The pathway scoring of this run: method, counts, gene detection and coverage filters, and the pathways scored and dropped |
| ENRICHMENT | Each enrichment run made on this run: the method, the library, the genes put in and the terms significant |
| OTHER RUNS | The study's other runs (cell types or counts) that have enrichment or pathway results of their own; choose one under **Results for** to describe it |
| REFERENCES, SOFTWARE | Citations for the methods used and the versions of the packages |

The pathway and enrichment sections read the study's provenance record, which keeps one entry per run: a rerun of a cell type's enrichment or pathway DE replaces only that cell type's entry. Runs recorded before KOSMIC recorded the cell type are not listed under OTHER RUNS.
