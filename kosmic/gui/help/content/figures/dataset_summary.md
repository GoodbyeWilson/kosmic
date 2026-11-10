# Dataset Summary

A table of the gene sets with at least one significant DE gene, grouped
by direction. A gene is significant at FDR < 0.05; **Min |log2FC|**
also requires its log2 fold change to exceed a value. It starts from the value the DE run was made with, so the figure counts the same genes as the run's files; a run saved before KOSMIC recorded the value is shown with 0.25, the setting of the time, and the tooltip says the value is assumed.

## Columns

| Column | Meaning |
|---|---|
| **Mean log2FC (sig. genes)** | Mean log2 fold change of the gene set's significant genes. Its sign places the gene set under *Upregulated* or *Downregulated*. |
| **Sig. up / Sig. down / Total** | Numbers of significant genes up, down, and in all. |

These summarise gene-level DE results; they are not a pathway-level test.
For that, use Pathway DE.
