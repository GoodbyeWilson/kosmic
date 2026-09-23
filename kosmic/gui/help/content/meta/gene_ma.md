# Gene Meta-Analysis (Consensus)

Pool per-study gene-level DE results with one or more methods. Single mode = one method; consensus mode = several methods, intersection of significance.

## Family checkboxes (Single / Consensus)

*TODO.*

## Top 50% and HKSJ modifiers

*TODO.*

## Calibration (analytical / CC permutation)

*TODO.*

## Tabs: Plots / Pathway Forest / Results / Venn / Enrichment / Cross-Dataset Reproducibility / LOO Validation

*TODO.*

### Direction conflicts in the Results tab

The Results table has three columns that count, for each gene, the
studies in which it is significant (study FDR < 0.05) with a positive
log2 fold change (**Up**) and with a negative one (**Down**).
**Conflict** is marked when both counts are above zero: the gene is
significantly up in at least one study and significantly down in
another. The line above the table reports how many of the significant
genes have such a conflict, and the same count is written to the output
panel and the Methods page.

Tick **Show only genes with opposite-direction effects** to list only
these genes. They occur mainly with Fisher's method, SumRank and gwOP,
which do not use the sign of the effect; see
[Statistical approach](statistics.md#direction-of-effect-across-studies).

Studies whose DE results have no adjusted p-values are not counted; the
output panel names them when the run starts.

## Olink panel overlay

*TODO.*
