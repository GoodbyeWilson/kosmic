# Leading-edge heatmap

The genes that drive a preranked GSEA result, and how each of them
changes in every dataset.

Preranked GSEA reports, for each gene set, its *leading edge*: the
members ranked before the running enrichment score reaches its peak,
which carry the enrichment. This figure takes the strongest leading-edge
genes of the gene sets you tick and shows each gene's log2 fold change as
a heatmap, one block of rows per gene set.

## Source

The figure draws on the last library GSEA run. It is listed in two
sections of Figure Export, and the section decides the run:

- **Meta section** (Pooled Enrichment), from the Meta-Analysis
  Enrichment step — the leading edge comes from the
  GSEA of the pooled gene ranking. There is one column per study of the
  meta-analysis, showing that study's own log2 fold change, and a
  **Pooled** column, separated by a gap, showing the pooled estimate.
  Each study column is labelled with the study's name; the cell type is
  in the figure title, and the exported table keeps the full name. This
  shows whether the studies agree on the genes behind the pooled
  result.
- **DE section**, from the DE workspace's Enrichment page — the leading
  edge comes from that study's GSEA, and there is one column, the
  study's log2 fold change.

## Reading the figure

- **Rows** are the leading-edge genes of each ticked gene set, strongest
  first in the direction of the set's enrichment: the most increased
  first for a set with a positive NES, the most decreased first for a
  negative one. Genes are ordered by the statistic the GSEA ranked on:
  the pooled z for a meta-analysis, the log2 fold change for a study.
- **Colour** is the log2 fold change, red for higher in the disease
  group and blue for lower. Changes beyond the colour range take the end
  colour.
- **\*** marks a gene that is significant in that column: adjusted p
  below 0.05 for a study, FDR below 0.05 for the pooled estimate.
- **Grey** marks a gene that column did not test, for example a gene
  filtered out of a study with few cells of this type.

GSEA matches gene symbols in upper case, so for mouse data the leading
edge is mapped back to the data's own symbols (`XBP1` to `Xbp1`) before
it is drawn.

## Controls

| Control | Effect |
|---|---|
| **Gene sets** | The gene sets of the run, most significant first; a set's tooltip gives its NES and FDR. GSEA estimates the FDR from 1,000 permutations, so an FDR of 0 (no permutation reached the observed score) is shown as FDR < 0.001, in the tooltip and in the block's title. Each ticked set is one block of rows. The two most significant are ticked when a run is first shown. |
| **Genes per set** | How many leading-edge genes to show for each set (default 15). |
| **Colour range** | The colour scale runs from minus to plus this log2 fold change. **Auto** (the default) uses the largest change shown, rounded up to the next 0.5. |

Exporting also writes a CSV next to the figure with one row per gene set,
gene and column: the log2 fold change and adjusted p-value behind each
cell.
