# Variable Genes

Each point is a gene: its mean expression (log scale) against its
normalised variance, or its normalised dispersion when the genes were
picked with the `seurat` method. The genes selected as highly variable
are shown in red; the title gives how many were selected out of how many
genes.

## How HVGs were picked

The scRNA workflow's HVG step selects them with scanpy's
`highly_variable_genes`: the `seurat_v3` method on raw counts, or the
`seurat` method on normalised data.
