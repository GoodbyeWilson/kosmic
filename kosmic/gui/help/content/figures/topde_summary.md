# DE Genes by Pathway

For each gene set, the number of DE genes up in the disease group (right
of zero) and down (left of zero), ordered by their total.

**Significant genes only** (default on) counts genes with FDR < 0.05;
**Min |log2FC|** also requires the log2 fold change to exceed a value.
It starts from the value the DE run was made with, so the figure counts the same genes as the run's files; a run saved before KOSMIC recorded the value is shown with 0.25, the setting of the time, and the tooltip says the value is assumed. Off, the first counts every tested gene by the sign of its
fold change.

The numbers after each gene set's name are its coverage: genes measured in
the dataset / genes in the set.
