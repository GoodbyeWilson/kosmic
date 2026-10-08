# Dataset Summary

A table of the gene sets with at least one significant DE gene, grouped
by direction. A gene is significant at FDR < 0.05; **Min |log2FC|**
(default *Off*) also requires its log2 fold change to exceed a value.

## Columns

| Column | Meaning |
|---|---|
| **Mean log2FC (sig. genes)** | Mean log2 fold change of the gene set's significant genes. Its sign places the gene set under *Upregulated* or *Downregulated*. |
| **Sig. up / Sig. down / Total** | Numbers of significant genes up, down, and in all. |

These summarise gene-level DE results; they are not a pathway-level test.
For that, use Pathway DE.
