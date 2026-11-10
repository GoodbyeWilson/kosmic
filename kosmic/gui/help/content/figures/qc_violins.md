# QC Violins

Distributions of cell-level QC metrics for the study's dataset: genes
detected per cell, counts per cell, and the percentage of counts from
mitochondrial (MT-) genes. The line in each violin marks the median.

The figure shows the cells that are in the dataset now, so the cells a
QC step removed are not in it. The [Cell Filtering](cell_filtering.md)
figure shows how many cells each step removed.

## Controls

| Control | Effect |
|---|---|
| **Violins** | *Pooled* draws one violin per metric for all cells. *By condition* draws control and disease side by side. *By sample* draws one violin per sample, controls first, with the metrics stacked in rows. |
| **Log scale for genes and counts** | On (the default), genes and counts per cell are drawn on a log axis, where their spread is visible. |
| **Show the QC limits** | On (the default), dotted lines mark the limits of the study's last recorded QC step, read from its provenance record. A limit that was not set has no line. |
| **Add DecontX contamination and doublet score** | Adds a panel for each of these when the dataset has them: the share of each cell's counts that DecontX estimated as ambient RNA, and Scrublet's doublet score. |
| **Y-axis range** (one row per metric) | The lowest and highest value on each metric's axis, in the metric's own units (cells for genes and counts, % for mitochondrial counts and DecontX). A box left empty shows *Auto* and sets that end from the data. Type a number and press Enter, or move to another box, to redraw. The DecontX and doublet rows apply when their panels are drawn. |

To compare studies, type the same ranges in each study's figure before
exporting: the exported figures then share their scales and can be
placed side by side. With *Auto*, each study's axes fit its own data,
so violins of the same height can stand for different values.

Control and disease are the roles set on the Inspect tab. In *By
condition* and *By sample*, controls are blue and disease red; a sample
is drawn under the role most of its cells have.

*By sample* is the view to check before analysis: a sample whose
violins sit well below the others has fewer genes or counts per cell,
and a sample with high mitochondrial or ambient RNA counts stands out
in its row.

Exporting also writes a CSV beside the figure with the number of cells
and the median of each metric for each violin.
