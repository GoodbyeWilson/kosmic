# Select Studies

Pick which DE results in the project to pool. Each row is one result:
a whole study, or one cell type of a study when DE was run per cell
type.

## What gets discovered

KOSMIC looks in every study folder of the project for DE results in
`results/de_analysis/statistics/`. A study's own result is named
`{study}_DE_{method}.csv`; a per-cell-type run writes one file per cell
type, `{study}_{cell type}_DE_{method}.csv`, and each of those appears
here as its own row (for example `Koenig_Endothelial_Cell`).

The list is refreshed when you open the Meta-Analysis workspace, if DE
results were added or rewritten since the last scan, so results from a
DE run in the same session appear without restarting. **Rescan** looks
again on demand.

## Choosing a cell type

Pooling is only meaningful within one cell type: a disease effect in
endothelium and one in cardiomyocytes are different contrasts. The
**Cell type** box lists every cell type found, with the number of
studies that have it. Choosing one selects that cell type in every
study and shows only those rows. **All entries** shows every row and
leaves the ticks as they are.

A cell type found in fewer studies than the minimum for pooling (3 by
default) is listed but cannot be chosen. **Whole study** collects the
results of DE run on a whole study, and **External imports** the files
added with **Import external result...**.

For the cell types to line up across studies, run DE by cell type on
the same label column in every study, such as `cell_type_atlas` after
propagating labels from a shared atlas.

## DE method filter

When a result exists for more than one DE method, the method chosen in
**View → Meta-Analysis Settings** is used (DESeq2 by default). A result
without that method uses the one it has, and the output panel says so.

## Skipping a problematic study

Untick its row after choosing the cell type. The other studies of that
cell type stay selected, and the Discovery step pools only the ticked
rows.
