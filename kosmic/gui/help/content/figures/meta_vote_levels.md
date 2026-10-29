# Study agreement (vote levels)

How many studies agree on each gene of a gene-level meta-analysis. Each
bar is a level: the number of studies that call a gene significant
(study FDR < 0.05) in the same direction, from 0 of k to k of k, split
into up-regulated (red) and down-regulated (blue). Genes called up in one study and down in another form a
separate grey **conflict** bar, as do genes no study calls (0 of k).

The levels come from the per-study calls that every gene-level run
records, so the figure is available after any gene-level run, with or
without vote counting.

**Significant by** hatches, in each bar, the genes that one pooling
method of the run calls significant (FDR < 0.05), for example REML. It
lists the methods of the run currently open in the Meta-Analysis
workspace; run REML, or REML with vote counting, to choose it. It is
**None** by default.

Exporting also writes a CSV next to the figure with the counts per
level (up, down, total and, with a method chosen, the significant
genes).

See [Vote counting](../meta/statistics.md#vote-counting) for how the
levels are counted.
