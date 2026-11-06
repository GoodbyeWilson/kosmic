# Cell Filtering

A flow chart of the cells left after each step that removed cells,
from the loaded dataset to the cells in the study now. The boxes down
the left hold the cells remaining; between two of them, a box to the
right names the step and the number of cells it removed.

The counts come from the study's provenance record
(`processed_data/provenance.json`), which each scRNA step writes when
it runs. The steps are:

| Step | Removes |
|---|---|
| Loaded | Nothing: the cells in the imported dataset. |
| Pre-analysis filter | Cells outside the samples or conditions kept for analysis. |
| QC limits | Cells outside the QC limits (genes and counts per cell, mitochondrial %). |
| Doublet removal | Cells Scrublet predicted to be doublets, when they were filtered. |
| QC limits on decontaminated counts | Cells that fall below the QC minimums once DecontX has removed ambient RNA from their counts. |
| Cell-type subset | Cells outside the cell types kept, for a study made by subsetting. |

A step that ran again without removing cells is left out. Studies
processed before doublet removal was recorded show it as *Doublet
removal* when the drop equals the number of doublets detected just
before it. A drop that no recorded step accounts for is shown as *Not
recorded*.

## Controls

| Control | Effect |
|---|---|
| **Start at** | *Loaded cells* starts from the imported dataset. *After the pre-analysis filter* starts after the last pre-analysis filter, which often removes whole samples or conditions; the QC steps are then drawn at a readable scale. |

Exporting also writes a CSV beside the figure with each step, the cells
remaining and the cells removed.
