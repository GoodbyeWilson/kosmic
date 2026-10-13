# Load Data

Load an h5ad file and configure sample and condition columns for DE analysis.

## Input format

KOSMIC expects an AnnData .h5ad file with raw counts in `layers['counts']` (what the scRNA QC step writes), or in `.raw` or `X` for a file from elsewhere. The file should already be QC-filtered and annotated with metadata columns in adata.obs.

## Counts

The **Counts** button chooses which counts the analysis reads:

- **Original counts** — the counts the study was imported with (the default).
- **DecontX-corrected counts** — `layers['decontX_counts']`, written by the scRNA Decontaminate step. Offered when the study's file has that layer.

KOSMIC loads only the chosen counts, so a study costs the same memory either way, and every step reads them: Gene DE, the per-cell-type run, pathway scoring and the pseudobulk exported for the meta-analysis. Changing the choice reloads the study. The choice stays for the session, so each study loaded after it uses it; a study without DecontX-corrected counts is loaded with its original counts, and the **Counts** row of the DATASET panel says so. The counts used are recorded in provenance with each result and stated on the Methods page.

## Sample column

Select the column in adata.obs that identifies biological replicates (e.g. 'patient', 'donor', 'sample_id'). Each unique value becomes one pseudobulk sample. You need at least 2 samples per condition.

## Condition column

Select the column that defines your experimental groups (e.g. 'disease', 'treatment', 'genotype'). Then assign which value is the control (reference) and which is the disease/test group. **Combine Conditions** adds a column that merges several values into one group; it is saved into the study's file, and nothing else in the file is changed.

## Minimum cells

Samples with fewer cells than this threshold are excluded from pseudobulk aggregation. The default is 10 (`[de] min_cells` in `config.toml`). Lower values include more samples but with noisier estimates.
