# Load Data

Load an h5ad file and configure sample and condition columns for DE analysis.

## Input format

KOSMIC expects an AnnData .h5ad file with raw counts in `layers['counts']` (what the scRNA QC step writes), or in `.raw` or `X` for a file from elsewhere. The file should already be QC-filtered and annotated with metadata columns in adata.obs.

## Sample column

Select the column in adata.obs that identifies biological replicates (e.g. 'patient', 'donor', 'sample_id'). Each unique value becomes one pseudobulk sample. You need at least 2 samples per condition.

## Condition column

Select the column that defines your experimental groups (e.g. 'disease', 'treatment', 'genotype'). Then assign which value is the control (reference) and which is the disease/test group.

## Minimum cells

Samples with fewer cells than this threshold are excluded from pseudobulk aggregation. Default is 100. Lower values include more samples but with noisier estimates.
