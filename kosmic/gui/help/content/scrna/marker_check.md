# Marker Check

Sense-check an annotation before you build on it. Pick a few of your cell
types; the plot shows known markers for those types against those types. A
correct annotation reads as a diagonal — each type bright on its own
markers, dark on everyone else's.

## Reading the plot

The plot is laid out like **Figures → Marker Dot Plot**, for the types
you ticked; use that page to export a figure. Columns are marker genes, kept together for each cell
type they mark; the **Marks** column of the table below names it. Rows are the cell types you ticked, placed under
the marker set they express most. Dot area is proportional to the
percentage of cells expressing the gene; colour is mean log-normalised
expression (a dataset holding raw counts is normalised for these genes
only, without changing it). The
key to both is beside the plot, and the colour scale is capped at the
99.5th percentile (a `+` on its top label marks the cap). Hover over a dot
for its values.

- **Bright block on the diagonal** — that type expresses what it should.
- **A column bright everywhere** — the gene is not specific in this tissue.
  Common for `DES`, `ENG`, and anything mitochondrial.
- **A row dark on its own markers** — the interesting one. That type is
  not expressing the genes it is named after.
- **Two rows bright on the same block** — those types are not separated
  by these markers. Often real: pericytes and smooth muscle sit on a
  continuum, and endocardium is specialised endothelium.

The table below lists the same genes in the plot's left-to-right order. **Marks** is the
type the gene is a marker for; **Peak group** is the type that actually
expresses it most. They should agree. Sorting is deliberately off, so the
rows stay aligned with the plot.

## Picking cell types

The list shows *your* cell types, not the marker database's. Those
vocabularies rarely match — CellTypist says "Ventricular Cardiomyocyte"
where a curated set says "Cardiomyocyte" — so each of your types is matched
to a marker set by whole words. Plurals count as the singular, and a few
abbreviations are expanded ("SMC" matches "Smooth Muscle Cell", "NK"
matches "Natural Killer"). Hover over a type to see which marker set it
was matched to.

A type with no match is greyed out and cannot be ticked. That is not a
failure of the annotation; it means this marker source has nothing for it.
Try another source, or supply your own.

Four types are ticked by default. That is not a limitation, it is the
readable size: everything at once is dozens of rows against a dozen columns
and tells you nothing.

## Marker sources

- **Cardiac (human, curated)** — three to six markers for each of the
  broad cell types of the human heart, named so that the CellTypist heart
  model's broad labels match.
- **General (built-in)** — around a dozen broad types from several
  tissues at five markers each. Small, canonical, no lookup.
- **PanglaoDB** — by organ. Broader and more citable; panels run long, so
  the per-type cap matters.
- **CellMarker2** — by tissue and species.
- **From file…** — any GMT, which is just `name<TAB>description<TAB>genes`.
  Use this for a panel from a paper, or for the bundled gene-family sets in
  `kosmic/reference/gene_sets/`.

When a dataset is first opened, the tab selects whichever of the two
bundled sources above has markers for more of its cell types (the general
list when they tie). Once you pick a source yourself, the tab keeps it.

## Is this a real check?

It depends on how the labels were made.

**CellTypist** is a model trained on reference atlases, weighting thousands
of genes. Canonical markers had no special role in producing its calls, so
plotting them is a genuinely independent check.

**Marker scoring** assigns each cluster the type whose markers it scores
highest on. Plotting those same markers gives a perfect diagonal by
construction — it confirms the arithmetic ran, not that the biology is
right. Use a different source than the one you annotated with.

## Independent vs atlas labels

Once labels have been propagated from an atlas, the sidebar reports how far
the two annotations agree. Both are kept: `cell_type` is what the study said
on its own, `cell_type_atlas` is what the combined atlas said.

Switch **Group by** between them to check each. Agreement in the high 90s
with disagreement concentrated at genuinely ambiguous boundaries — pericyte
against smooth muscle, endocardial against endothelial — is the expected
result. Cells moving between unrelated lineages, like macrophage to smooth
muscle, are more likely doublets than a disagreement worth resolving.
