# Cell-type composition

The percentage of cells of each type, as 100% stacked bars.

## Controls

| Control | Effect |
|---|---|
| **Cell types** | The obs column holding the cell-type labels. |
| **Bars** | *By condition* draws one bar for control and one for disease. *By sample* draws one bar per sample, controls first, with a bracket over each condition's samples. |
| **Colours** | The colour scheme. With the same scheme and column as the UMAP page, each cell type has the same colour in both figures. |

Control and disease are the roles set on the Inspect tab, so studies
that name their conditions differently (NF, Donor, normal) are grouped
together. Cells with neither role, and cells labelled Unknown or
Unassigned, are left out.

The most abundant type sits at the bottom of each bar, and the legend
lists the types in the order they are stacked. *By condition* pools the
cells of every sample in a condition, so a sample with many cells
counts for more; *By sample* shows the spread between samples, which
is what a test of a composition difference would use.

Exporting also writes a CSV beside the figure with the percentages, one
column per bar.
