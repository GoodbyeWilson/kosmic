# Marker score UMAPs

One UMAP per cell type of a marker panel, each cell coloured by that
panel's module score. A good panel lights up one population of the
embedding and leaves the rest grey.

## Controls

| Control | Effect |
|---|---|
| **Markers** | The marker panel, as on the [marker dot plot](marker_dotplot.md). |
| **Embedding** | UMAP or t-SNE, where the dataset has both. |
| **Cells shown** | Cells drawn, sampled at random (default 150,000). Scoring and drawing every cell of a large atlas is slow and looks the same. |
| **Point size** | Size of each cell's point. |
| **Columns** | UMAPs per row. |

## The score

The module score is the mean expression of the panel's genes in a cell
minus the mean of control genes drawn from the same expression range
(scanpy's `score_genes`, after Tirosh et al. 2016). Subtracting the
controls stops a panel of highly expressed genes from scoring high in
every cell. Panels with fewer than two of their genes in the dataset are
left out; the genes used are written under each UMAP.

For display, each panel's score is shifted so the median cell is at zero
(scores below it are drawn as zero) and divided by its 99.9th
percentile, so every panel runs from low to high on the same colour
scale. Anchoring the scale on the top of the signal lets a type that is
only 0.5% of the cells reach the top colour.

> **Warning:** This scaling cannot tell a rare population from an
> absent one. If a dataset has no cells of a panel's type, stray
> (ambient) counts of its genes are what reach the top of the scale, and
> the UMAP shows scattered speckle instead of one bright cluster. Check
> such a panel on the marker dot plot before reading it as a population.
