# Pathway Meta-Analysis

Pool per-study pathway scores using DESeq2 + VIF (CAMERA-style). Single pooling method; default REML.

## What it does

*TODO.*

## Pooling method

*TODO.*

## Output columns

| Column | Meaning |
|---|---|
| Pathway | Pathway name |
| log2FC, SE | Pooled log2 fold change and its standard error |
| CI low, CI high | 95% confidence interval of the pooled log2 fold change |
| Pooled P, FDR | Pooled p-value and its Benjamini–Hochberg FDR across pathways |
| k | Number of studies pooled for the pathway |
| I^2 | Proportion of the variation between studies attributed to heterogeneity |
| Up, Down | Studies in which the pathway is significant (study FDR < 0.05) with a positive or a negative log2 fold change |
| Conflict | Marked when the pathway is significantly up in one study and significantly down in another |

Up, Down and Conflict use each study's adjusted p-values. With CPM +
Welch scoring these come from the pathway DE results. With DESeq2 + VIF
scoring each study's pathway p-value is a Wald test of the pooled log2
fold change against its VIF-adjusted standard error, adjusted by
Benjamini–Hochberg across that study's pathways. See
[Statistical approach](statistics.md#direction-of-effect-across-studies).
