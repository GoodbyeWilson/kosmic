# Enriched term links

Shows which DE genes are behind the top enriched terms, as a chord diagram
or a grid. It reads the last enrichment run in the DE workspace (GO, or an
Enrichr library) and the genes' log2 fold changes from Gene DE.

## Which terms and genes

The top terms are those passing the enrichment's significance threshold
(FDR, or elim P for GO elim), most significant first; **Terms** sets how
many (default 8). If no term passes, the most significant terms are shown
and the title says so. For each term, **Genes per term** sets how many of
its genes from the query are shown: those with the largest |log2 fold
change| (default 8).

**Exclude mitochondrial genes (MT-)** leaves MT- genes out of the figure.
In single-nucleus data their fold changes are largely technical (ambient
or carry-over RNA), and they can otherwise take most of the places among
the largest fold changes. Only the genes drawn change; the enrichment
result is not recomputed. To exclude them from the analysis itself, use
the same option on the Gene DE page.

## Chord

Genes sit around the left of the circle and terms around the right. Each
gene's arc is coloured by its log2 fold change (blue down, red up); each
term's arc is sized by its number of genes. A ribbon, in the term's
colour, joins every gene to each term it belongs to, so genes shared
between terms show as fans of ribbons.

## Grid

One row per term and one column per gene, with a dot coloured by log2 fold
change where the gene belongs to the term.

## Order genes

*By fold change* orders genes from most down to most up, so the chord's
gene ring is a colour gradient and each grid row reads as a direction
profile. *By term* groups genes by the first term they belong to, which
keeps the chord's ribbons from crossing.

The colour scale is symmetric and capped at the 90th percentile of
|log2 fold change| among the genes shown; the arrows on the colour bar
mark the cap.

GO elim removes genes shared between parent and child terms, so its terms
share few genes. Fisher + BH and Enrichr libraries keep them.
