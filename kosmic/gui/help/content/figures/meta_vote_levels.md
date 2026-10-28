# Study agreement (vote levels)

How many studies agree on each gene of a gene-level meta-analysis. Each
bar is a level: the number of studies that call a gene significant
(study FDR < 0.05) in the same direction, from 0 of k to k of k, split
into up-regulated (red) and down-regulated (blue). Genes called up in
one study and down in another form a separate grey **conflict** bar, as
do genes no study calls (0 of k).

**Overlay** hatches, at each level, the genes that one method of the
run calls significant (FDR < 0.05), for example REML. Read together, the
bars show what the pooled analysis confirms among the genes the studies
already call, and what it adds among genes few or no studies call on
their own. Choose **None** to leave the overlay out.

The figure is the same as the **Votes** tab of the Discovery
Meta-Analysis step. Exporting also writes a CSV next to the figure with
the counts per level (up, down, total and, with an overlay, the
significant genes).

See [Vote counting](../meta/statistics.md#vote-counting) for how the
levels are counted.
