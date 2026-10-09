# UMAP

Cell-level UMAP (or t-SNE) of the loaded dataset, coloured by an obs
column such as cell type, cluster, condition or sample.

## Controls

| Control | Effect |
|---|---|
| **Color by** | The column to colour by. A cell-type column is chosen by default. |
| **Embedding** | UMAP or t-SNE, where the dataset has both. |
| **Point size** | Size of each cell's point; large datasets need values below 1. |
| **Colours** | The colour scheme for the categories. |
| **Labels** | *On plot* writes each category's name in bold at the centre of its cells, with no legend; names are placed so they do not overlap. *Legend* lists the categories, with their cell counts, beside the plot. Use *Legend* for columns with many categories. |
| **Highlight** | Colours one category and draws every other cell grey. |

The title gives the number of cells and of categories shown.
