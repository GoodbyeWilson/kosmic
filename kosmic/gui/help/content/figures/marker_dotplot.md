# Marker dot plot

Shows how well the dataset's cell-type labels agree with known marker
genes. Each column is a marker gene, grouped under the cell type it marks
(brackets, top); each row is a cell type of the loaded dataset. A correct
annotation reads as a diagonal.

## Controls

| Control | Effect |
|---|---|
| **Cell types** | The obs column holding the cell-type labels. |
| **Markers** | The marker panel. *Cardiac (human, curated)* is 71 genes for 14 cardiac cell types. *General (built-in)* is a short list for common cell types across tissues. *From file...* reads a GMT file (set name, description, then genes, tab-separated) with one set per cell type. |

## Reading the plot

Dot area is proportional to the percentage of the type's cells with
non-zero expression of the gene, so a dot half the area of another marks
half the share of cells. Colour is the mean expression over all cells of
the type. When the dataset's expression is log1p of counts per 10,000
(the scRNA workflow's normalisation), the colour bar says so. The colour
scale is capped at the 99.5th percentile of the values shown, so one very
high gene does not wash out the rest; an arrow on the colour bar marks
the cap.

Rows are placed under the panel whose genes they express most, in the
panel's order, rather than matched by name. A dataset's labels and a
panel rarely name types the same way ("LYVE1+ Macrophage" against
"Myeloid"), and placing by expression still gives a diagonal. A row that
sits under an unexpected panel is worth checking.

Every cell of each type is used, so the numbers describe the whole
dataset. Cells labelled Unknown or Unassigned are left out.

Exporting also writes a CSV beside the figure with the percentage and
mean expression for every gene and cell type.
