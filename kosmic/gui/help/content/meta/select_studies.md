# Select Studies

Pick which DE results in the project to pool. Each row is one result:
a whole study, or one cell type of a study when DE was run per cell
type.

## What gets discovered

KOSMIC looks in every study folder of the project for DE results in
`results/de_analysis/statistics/`, in the folder of the counts each run
used (`original_counts/`, `decontX_counts/`) and, for a per-cell-type
run, the cell type's folder inside it. A study's own result is named
`{study}_DE_{method}.csv`; a per-cell-type run writes one file per cell
type, `{study}_{cell type}_DE_{method}.csv`, and each of those appears
here as its own row (for example `Koenig_Endothelial_Cell`). Results
saved directly in `statistics/` by earlier versions are found too.

DE results of the combined atlas (`_master`) are not listed. The atlas
contains the cells of every study, so pooling its result with the
per-study results would count each study twice.

The list is refreshed when you open the Meta-Analysis workspace, if DE
results were added or rewritten since the last scan, so results from a
DE run in the same session appear without restarting. **Rescan** looks
again on demand.

## Choosing the counts

A study can hold DE results on its original counts and on corrected
counts such as the DecontX counts. These are separate analyses, so a
meta-analysis pools results made on one of them. The **Counts** box
lists the counts the project's results were made on, with the number of
studies that have results on each, and only results on the chosen counts
are listed. The counts that most results were made on are chosen first, and the original counts when two are equal. A
study with no result on the chosen counts is named in the output panel.

Results saved before KOSMIC recorded the counts of a run are listed
under **Counts not recorded**.

## Choosing a cell type

Pooling is only meaningful within one cell type: a disease effect in
endothelium and one in cardiomyocytes are different contrasts. The
**Cell type** box lists every cell type found, with the number of
studies that have it. Choosing one selects that cell type in every
study and shows only those rows. **All entries** shows every row and
leaves the ticks as they are.

When the project has more than four results, as it does after DE is run
per cell type, the box starts at **Choose a cell type...** and no rows
are listed or selected until you choose one.

A cell type found in fewer studies than the minimum for pooling (3 by
default) is listed but cannot be chosen. **Whole study** collects the
results of DE run on a whole study, and **External imports** the files
added with **Import external result...**.

For the cell types to line up across studies, run DE by cell type on
the same label column in every study, such as `cell_type_atlas` after
propagating labels from a shared atlas.

## Where the results are saved

Each selection saves its results in its own folder under the project's
`meta_analysis/` folder: a folder for the counts the results were made
on, and inside it a folder named after what was pooled. Running the
same method on another cell type, or on the other counts, therefore
does not replace earlier results:

- one cell type across studies:
  `meta_analysis/original_counts/Endothelial_Cell/` or
  `meta_analysis/decontX_counts/Endothelial_Cell/`;
- whole-study results: `meta_analysis/original_counts/all_cells/`;
- more than one cell type: `meta_analysis/original_counts/mixed/`, with
  a warning in the output panel.

A selection that includes an externally imported result, or results
whose counts were not recorded, is saved under
`meta_analysis/unspecified_counts/`, because the counts behind those
results are not known. Externally imported results do not change the
cell-type folder name.

The folder holds the result tables, their settings files, and the
record the Methods step renders for that selection; the settings file
and the record name the counts. The line under the page title shows the
folder for the current selection. Results saved by earlier versions in
`meta_analysis/<cell type>/` stay readable by the Figures workspace.

## DE method filter

When a result exists for more than one DE method, the method chosen in
**View → Meta-Analysis Settings** is used (DESeq2 by default). A result
without that method uses the one it has, and the output panel says so.

## Skipping a problematic study

Untick its row after choosing the cell type. The other studies of that
cell type stay selected, and the Discovery step pools only the ticked
rows.
