# Pathway Bubble

One bubble per gene set from Pathway DE, ordered by effect size (disease
vs control).

- **Position:** the effect size; left of the vertical line is lower in
  the disease group.
- **Size:** -log10 FDR; the legend gives the scale.
- **Colour:** red or blue when the gene set passes FDR < 0.05 (up or down
  in the disease group), grey otherwise.

The effect size is on the scale of the Pathway DE scoring method; for the
DESeq2 method it is the mean per-gene log2 fold change.
