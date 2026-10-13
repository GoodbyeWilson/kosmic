# ADR-011: A shared atlas may be built from the studies' DecontX counts

- Status: Accepted by Kirk (2026-10-13).

## Problem

Combine builds the shared atlas from each study's original counts
(`layers['counts']`). Ambient RNA in those counts misled the atlas's
annotation on the DCM project:

- About 1,000 fibroblast nuclei carrying cardiomyocyte ambient RNA fell
  in a Leiden cluster of cardiomyocytes and were labelled cardiomyocyte.
  In one Youness sample they were 38% of the cells called cardiomyocyte,
  and that study's cardiomyocyte DE reported fibroblast genes (LUM,
  COL1A2, DCN) as strongly up in DCM.
- Subclustering the cluster and annotating it with CellTypist did not
  separate them. The contaminated cardiomyocytes and the contaminated
  fibroblasts each carry both programmes, they share subclusters at any
  resolution, and CellTypist's calls for them are split even cell by
  cell. Majority voting then labelled about 1,100 cardiomyocytes as
  fibroblast, adipocyte or smooth muscle.

Each study has DecontX-corrected counts already (`layers['decontX_counts']`),
and the DE workspace can run on them. DecontX on the atlas itself does
not run on the machine the project is analysed on, so an atlas, and a
merged pseudobulk DE on it, could only use original counts. Comparing
corrected with uncorrected counts across both the per-study and the
merged analysis needs an atlas built from the corrected counts.

DecontX estimates contamination from the cell-type labels it is given,
so its output follows those labels: a cell given the wrong label has its
own programme removed as contamination. Corrected counts are only as
good as the labels the studies were decontaminated with.

## Decision

**Combine can build the atlas from a named layer of every study. The
Combine dialog offers the DecontX layer when every selected study has
it.**

1. `concat_studies(counts_layer='decontX_counts')` reads that layer of
   each study where it would read the original counts. A study without
   the layer stops the build, naming the study.
2. The atlas holds those counts as its own. After QC its
   `layers['counts']` are the corrected counts: the data contract that
   `layers['counts']` holds a dataset's original counts does not hold for
   such an atlas. The original counts remain in the study files.
3. The `combine` provenance stage records `counts_layer`.
   `provenance.combined_counts_layer()` reads it back. The DE Methods
   page uses it to describe a DE on such an atlas as run on DecontX
   counts and to cite DecontX, and the Decontaminate step warns before
   correcting such an atlas a second time.
4. The default is unchanged: without the option, the atlas is built from
   the original counts.

## Consequences

- An atlas's counts are original or corrected, and only the provenance
  says which. Code that needs to know reads `combined_counts_layer()`;
  nothing in the h5ad itself marks it.
- An atlas built from corrected counts cannot be used to compare
  corrected with uncorrected counts on its own. That comparison needs
  two atlases, built separately; Combine writes to the project's
  `_master` folder, so the other atlas has to be moved out of the
  project first.
- Label propagation writes the corrected atlas's labels to the studies
  as before; the studies' own counts are not touched.
