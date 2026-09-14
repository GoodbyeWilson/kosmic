# The scRNA workspace

The scRNA workspace is the single-dataset pipeline: it takes **one
study** from an imported file to a QC'd, clustered, annotated dataset
ready for differential expression. For a single-dataset analysis it
is most of the work; for a meta-analysis you run it once per study.
It opens on the study selected in the Project table, and everything
it does is written back to that study's `processed_data/` h5ad — so
when you return to Project, the ticks have moved.

Eight steps, in the order you meet them:

| Step | What it does | Leaves behind |
|---|---|---|
| [Load Data](load_data.md) | Shows the study's dataset and its basic properties; converts a raw import if one is waiting | the working h5ad |
| [Gene Names](gene_names.md) | Renames genes to current HGNC symbols so studies can be joined | `var_names` harmonised |
| [Inspect](inspect.md) | Tells KOSMIC which column is the sample, which is the condition, and what each condition value means | `sample`, `condition`, roles |
| [Quality Control](qc.md) | Removes low-quality cells and doublets, optionally ambient RNA, then normalises | filtered cells, `.raw` snapshot |
| [Cluster](cluster.md) | HVGs → PCA → (Harmony) → neighbours → Leiden → UMAP, then cell-type annotation | `leiden`, `cell_type`, embeddings |
| [Marker Check](marker_check.md) | Sense-checks the annotation against known markers | nothing — it is a check |
| [Decontaminate](decontx.md) *(optional)* | Removes ambient RNA per sample using the cell-type labels | `layers['decontX_counts']` |
| [Subset](subset.md) | Cuts one cell type out as a new study | a new study folder |

Steps are gated: each needs the one before it, and the STATUS panel in
the explorer says what is missing.

![The scRNA workspace: eight steps in the explorer, the step's analysis summary cards on the left, and the main view — here the UMAP after clustering and annotation](img/cluster.png)

## What order the data goes through

```
raw import  →  harmonised names  →  roles set  →  QC'd + normalised
            →  clustered + annotated  →  (decontaminated)  →  subset
```

Two things are worth holding in mind:

- **The original counts are preserved.** QC keeps them in `.raw` and
  `layers['counts']`; normalisation, clustering and annotation all work
  on copies or add columns. Differential expression later runs on the
  counts, not the normalised values.
- **Everything is recorded.** Each step writes its settings to the
  study's `provenance.json`. The DE workspace's Methods step turns that
  into a paragraph you can paste into a paper — so you do not need to
  note parameters as you go.

## Doing this for several studies

A meta-analysis wants the studies processed *the same way*. Do them one
after another, keeping the QC thresholds, the annotation reference and
the cell-type level the same. Where a study genuinely needs a different
choice — a different mitochondrial cut for nuclei, say — the provenance
record will show the difference, and so will the Methods text.

If you want cell types defined once across all studies, do steps 1–4
for each study, then build the **shared atlas** from the Project
workspace and come back here to cluster and annotate it. See
[Shared atlas](../project/atlas.md).
