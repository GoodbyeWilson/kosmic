# Meta-analysis volcano

Volcano plot of the latest meta-analysis in the Meta-Analysis workspace:
each point is a gene (**Gene Meta-Analysis Volcano**) or a pathway
(**Pathway Meta-Analysis Volcano**), with its pooled log2 fold change on
the x axis and -log10 FDR on the y axis. The second line of the title
names the selection (for example the cell type) and the number of
studies pooled.

Red and blue mark items with FDR < 0.05 (up and down), the same
definition the Meta-Analysis workspace uses when it counts significant
results. **Min |log2FC|** (default *Off*) also requires the pooled log2
fold change to exceed a value, e.g. 0.25. The dotted line marks the FDR
threshold, and the most significant items are labelled without overlap,
as on the [Volcano](volcano.md) page.

The figure shows whatever the Meta-Analysis workspace currently holds. To
draw another selection, run or reload its meta-analysis first.
