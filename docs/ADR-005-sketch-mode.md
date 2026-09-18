# ADR-005: Sketch mode for large atlases

- Status: Proposed (2026-09-18). Awaiting Calum and Kirk. Depends on
  ADR-003 for the memory figures; independent of ADR-004.

## Problem

Embedding and clustering cost time and memory in proportion to the
number of cells, and two of the steps are single-threaded: Harmony
(25 minutes on the 1.09 million-cell DCM atlas) and Leiden. UMAP's
spectral initialisation is another. A ten-million-cell atlas is out of
reach on any desktop this way, and even the DCM atlas is an hour on a
64 GB workstation.

But the labels do not need every cell to be embedded. A few hundred
thousand representative cells find the same broad cell types, and the
remaining cells take their label from their nearest neighbours in the
embedding space. This is how large atlases are labelled in practice
(Seurat's sketch-based integration, the CELLxGENE Census workflows).
What must not be sketched is the differential expression: the
pseudobulk per donor per cell type has to sum every cell, or the
mega-analysis is a comparison against a subsample.

KOSMIC already has both halves in a form that was built for a
different reason: the per-study cell cap in Combine, and the kNN
projection in propagation for cells the cap left out. The cap was
withdrawn from the DCM run because the mega-analysis DE would have run
on the capped atlas. Sketch mode is the cap done properly.

## Decision

**An atlas may be embedded, clustered and annotated on a sketch of its
cells; every cell is labelled by projection; and every downstream
computation over the counts -- pseudobulk, marker means, annotation
means -- runs over all cells, from the file, in chunks.**

1. The atlas file holds every included cell's counts (ADR-003: one
   matrix). Nothing about the file changes; the sketch is an obs flag,
   `obs['in_sketch']`.
2. The sketch is drawn per study, stratified by sample so each donor is
   represented, to a configurable total (default 250,000). Leverage-score
   sketching (as in Seurat v5) is preferred to uniform sampling because
   it keeps rare types; uniform is the fallback.
3. HVG selection, PCA, Harmony, the neighbour graph, Leiden, marker
   ranking and UMAP run on the sketch. The reference annotation labels
   the sketch's clusters. This is the existing cluster-step code on a
   masked working subset (`workset.py`); the mask is `in_sketch`.
4. Every cell gets `leiden` and `cell_type` by projection: its
   coordinates in the sketch's PCA space (the stored loadings and gene
   means) and the majority label of its nearest sketch cells. This is
   the propagation path that already exists for capped cells, applied
   within the atlas before propagation to studies. The projection is
   chunked over the file; the full matrix is never loaded.
5. Pseudobulk for the mega-analysis DE sums counts over all cells of a
   donor and type in one chunked pass over the file. Annotation means
   for Marker Check are computed the same way. The UMAP shows the
   sketch, and says so in its caption.
6. Provenance records the sketch (size, method, seed) and the
   projection (k, agreement of projected labels with sketch labels on a
   held-out portion of the sketch, which is the quality check).

Sketch mode is off by default and offered when a dataset exceeds a
configured size (`sketch_above_cells`, default 500,000). A study small
enough is processed whole, as now.

## Consequences

- The DCM atlas: embedding and clustering on 250,000 cells takes about
  five minutes and a few GB; projection of the remaining 840,000 is a
  chunked pass of a few minutes. The whole atlas step fits a 16 GB
  laptop. Ten million cells is the same procedure, longer.
- The mega-analysis DE is unchanged in what it computes: every cell,
  every donor. Only the way labels were obtained differs, and that is
  recorded.
- The comparison the paper makes -- labels from a sketch against labels
  from embedding every cell -- is available for the DCM atlas, where
  both have been done.
- The cell cap in Combine is retired; the sketch replaces it.
- Effort: about a week, most of it the chunked passes over the file and
  the projection quality check. It builds on `workset.py`,
  `propagate_labels` and `h5ad_meta.py`.
