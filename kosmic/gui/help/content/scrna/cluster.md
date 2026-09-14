# Cluster

Turns a normalised matrix into cell types, in three stages the sidebar
shows as cards: **PCA**, **Clustering**, **Annotate**. The UMAP fills
the screen; the *Elbow Plot* tab sits beside it, and a strip of
resolution thumbnails runs underneath. Each card says what it did once
it has run; *View full report* gives the numbers.

![Cluster: PCA, Clustering and Annotate cards; the UMAP coloured by cell type; Elbow Plot and Annotate tabs above it](img/cluster.png)

## PCA

**HVG count** (2,000) — how many highly variable genes to keep. The
method follows the matrix: on raw counts KOSMIC uses Seurat v3's
mean–variance model; on log-normalised data it uses the dispersion-based
method, because v3 gives nonsense on log values.

**Compute PCs** (50) — how many components to compute. You will use
fewer; computing more costs little and lets the elbow plot show the
tail.

Two options that matter more than they look:

- **Ignore cells marked 'exclude'** — cells whose condition you set to
  Exclude in Inspect are left out of HVG selection, PCA, the neighbour
  graph, clustering and the UMAP. They stay in the file with their
  metadata; they get no coordinates and no cluster. Untick to embed
  every cell (Exclude then applies to the DE contrast only).
- **Scale genes before PCA** — z-scores each gene so PCA weighs a
  lowly-expressed transcription factor as heavily as a ribosomal gene.
  Off by default because it makes the matrix dense, which on a large
  dataset costs more memory than everything else together. HVG
  selection has already removed most of what scaling would help with.

**Run PCA**, then look at the **Elbow Plot**. KOSMIC suggests a number
of PCs to use — the largest of the classic elbow (max distance to the
line), the point where 80% of variance is explained, and a floor of 10.
Under-estimating is the common mistake, so it errs high. The suggestion
lands in *PCs to use*; you can change it.

## Batch correction

**Harmony batch correction** adjusts the PCA space so cells from
different batches overlap where they are the same type. **Batch key**
is the grouping to correct for:

- On a **shared atlas** it is `study`. Without it the clusters will
  substantially *be* the studies. This is the one place Harmony is not
  optional.
- Within **one study** it is the sample or donor column — useful when
  samples were processed on different days and the UMAP shows it.

**Merge across** — an optional second variable, so its groups land in
the same clusters. Set it to the condition column when you want one
cardiomyocyte cluster rather than a diseased one and a healthy one.
Leave it blank to keep disease-specific cell states separate.

This affects clustering and the UMAP only. Differential expression runs
on the raw counts and never sees the corrected embedding — so Harmony
cannot create or remove a DE result.

## Clustering

**PCs to use** (from the elbow), **Neighbors (k)** (15) and
**Resolution** (0.5) → Leiden clustering on the neighbour graph
(igraph flavour, two iterations), then a UMAP. **Keep existing
embeddings** skips recomputing the UMAP when you only changed the
resolution.

Resolution is the knob you will actually turn. Higher splits clusters,
lower merges them; there is no right value, only the one that matches
the cell-type level you want to annotate at. Two tools for choosing:

- **Resolution Sweep** — runs Leiden from *Min* to *Max* in *Step*s
  (0.1 to 0.8 by 0.1) and shows a thumbnail UMAP for each. Click one to
  adopt it.
- **Stability Analysis** — a clustree: how clusters split as resolution
  rises. A cluster that keeps its members across several resolutions is
  real; one that appears at one resolution and dissolves at the next is
  noise.

Colour the UMAP by `study` or `sample` before you go further. If the
clusters are the batches, fix that first.

## Annotate

Gives each cluster a cell-type name, written to `obs['cell_type']`.
Four methods, picked from *Method*:

| Method | What it is | When |
|---|---|---|
| **CellTypist** | Trained classifiers per tissue; per-cell labels, *Majority voting* gives one label per cluster, *Broad labels* collapses subtypes | the default for human tissues CellTypist has a model for |
| **Reference Atlas** | Label transfer from an annotated h5ad you supply (*Build from h5ad…*) | you have a trusted atlas of the same tissue |
| **CellMarker 2.0** | Marker-gene scoring against the CellMarker database, by species and tissue | a tissue without a good classifier |
| **PanglaoDB** | Marker-gene scoring against PanglaoDB; *Canonical markers only* tightens it | the same; the two databases disagree usefully |

*Run ORA* on the marker methods adds an over-representation test per
cluster, which is more robust than a raw score when a cluster is small.

The cluster table underneath shows each cluster's assigned type and
lets you change it by hand — pick from the dropdown, **Save Changes**.
**Use existing cell_type** keeps labels that came with the dataset and
skips annotation.

Whatever the method, go to **Marker Check** next and look. Automated
annotation is right often enough to be trusted and wrong often enough
that you must not.

## What is saved

`obs['leiden']`, `obs['cell_type']`, the PCA and UMAP coordinates, the
neighbour graph, and the settings for all of it in provenance. The
counts are untouched.
