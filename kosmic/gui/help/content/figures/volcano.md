# Volcano

Volcano plot of the current Gene DE result: each point is a gene, with its
log2 fold change on the x axis and -log10 FDR (or nominal P, if chosen)
on the y axis.

## Colours and threshold

Genes with FDR < 0.05 are red (up in the disease group) or blue (down);
all other genes are grey. **Min |log2FC|** also requires the log2 fold
change to exceed a value. It starts from the value the DE run was made with, so the figure counts the same genes as the run's files; a run saved before KOSMIC recorded the value is shown with 0.25, the setting of the time, and the tooltip says the value is assumed. The dotted line marks the FDR
threshold.

## Labels

The most significant genes are labelled, split between up and down
(Max labels sets how many; 0 turns labels off). **Label named genes only**
(on by default) passes over genes with no descriptive symbol when choosing
the labels: Ensembl IDs, clone names such as AC004936.1, numbered lncRNAs
(LINC00670) and antisense or divergent transcripts (CDC42EP3-AS1, CRIM1-DT).
The next named gene is labelled instead; the skipped genes are still
plotted and counted. Each label is placed at
the nearest position that does not overlap another label or a labelled
point, with a thin line to its point when it has been moved away. A label
with no free position is left out. Labels can be moved further in a
vector graphics editor after exporting as PDF or SVG.
