# Gene Heatmap

Expression of one gene set's genes in every sample (donor), as a heatmap
with genes as rows and samples as columns.

## Values

Each cell is the sample's mean log-normalised expression of the gene
(counts scaled to 10,000 per cell, then log1p), z-scored per gene across
samples. Only genes tested in Gene DE are shown, and samples with fewer
cells than the DE minimum are left out.

## Sample order and colours

Samples are ordered control first, then disease; the grey strip under the
columns marks the group. The default colour map is blue-white-red, centred
on zero, so white is the gene's average across samples. **Min / Max
percentile** clip the colour scale; other colour maps can be chosen.
