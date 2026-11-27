# Pathway Meta-Analysis

Pool each study's Pathway DE results across studies.

## What it does

In each study, the Pathway DE page (DE workspace) gives every donor one score per pathway: DESeq2-normalised log2 expression averaged over the pathway's genes. Welch's t-test compares that study's disease and control donors, giving each pathway an effect (the difference in mean score, in log2) and its standard error. This page pools those per-study effects across studies with a random-effects model, as the gene meta-analysis pools genes, and corrects the pooled p-values with Benjamini–Hochberg across pathways.

Run Pathway DE in each study first, with the same gene sets.

## Which studies are pooled

Only tables made with the current Pathway DE are pooled: their score and test are recorded beside each table. A table saved by an earlier version of KOSMIC is not, because it may have been made with one of the scoring methods since removed (CPM, VST and per-cell scores among them), and pooling tables made different ways would combine different quantities. Such studies are named under **Loaded Studies** and in the log; re-run Pathway DE in them.

## Pooling method

REML (restricted maximum likelihood) random effects, the same code as the gene meta-analysis. **Min studies** sets how many studies must have a pathway for it to be pooled.

Across studies, the pathway question can also be asked of the pooled gene results: the meta-analysis **Enrichment** step runs preranked GSEA on the genes ranked by their pooled z-statistic.

## Output columns

| Column | Meaning |
|---|---|
| Pathway | Pathway name |
| log2FC, SE | Pooled effect (difference in mean pathway score, log2) and its standard error |
| CI low, CI high | 95% confidence interval of the pooled effect |
| Pooled P, FDR | Pooled p-value and its Benjamini–Hochberg FDR across pathways |
| k | Number of studies pooled for the pathway |
| I^2 | Proportion of the variation between studies attributed to heterogeneity |
| Up, Down | Studies in which the pathway is significant (study FDR < 0.05) with a positive or a negative effect |
| Conflict | Marked when the pathway is significantly up in one study and significantly down in another |

Up, Down and Conflict use each study's adjusted p-values from its Pathway DE results. See [Statistical approach](statistics.md#direction-of-effect-across-studies).
