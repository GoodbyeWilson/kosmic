# Enriched terms

The top terms of the last enrichment run in the DE workspace (GO or an
Enrichr library), most significant first. **Terms** sets how many are
shown (default 20); terms passing the 0.05 threshold are shown first.

## Dot view (default)

- **Position:** the term's fold enrichment: the share of your query genes
  that are in the term, divided by the share expected from the background
  genes. Values above 1 mean the term holds more of your genes than
  chance would give.
- **Dot size:** the number of query genes in the term.
- **Colour:** -log10 of the significance value, darker for more
  significant.

## Bar view

One bar per term, its length -log10 of the significance value; the dotted
line marks the 0.05 threshold.

## Significance value

The legend and axis name the statistic. GO results from the elim method
are not adjusted for multiple testing, following topGO, and are labelled
"elim P (uncorrected)". Use Fisher + BH, or an Enrichr library, for
FDR-adjusted results. If no term passes the threshold, the most
significant terms are shown and the title says that none passed.
