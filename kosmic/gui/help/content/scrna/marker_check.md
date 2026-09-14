# Marker Check

Sense-check an annotation before you build on it. Pick a few of your cell
types; the plot shows known markers for those types against those types. A
correct annotation reads as a diagonal — each type bright on its own
markers, dark on everyone else's.

## Reading the plot

Rows are marker genes, grouped by the cell type they mark. Columns are the
cell types you ticked. Dot **size** is the fraction of cells expressing the
gene; **colour** is mean expression.

- **Bright block on the diagonal** — that type expresses what it should.
- **A row bright everywhere** — the gene is not specific in this tissue.
  Common for `DES`, `ENG`, and anything mitochondrial.
- **A column dark on its own markers** — the interesting one. That type is
  not expressing the genes it is named after.
- **Two columns bright on the same block** — those types are not separated
  by these markers. Often real: pericytes and smooth muscle sit on a
  continuum, and endocardium is specialised endothelium.

The table below lists the same genes in the same order. **Marks** is the
type the gene is a marker for; **Peak group** is the type that actually
expresses it most. They should agree. Sorting is deliberately off, so the
rows stay aligned with the plot.

## Picking cell types

The list shows *your* cell types, not the marker database's. Those
vocabularies rarely match — CellTypist says "Ventricular Cardiomyocyte"
where a curated set says "Cardiomyocyte" — so each of your types is matched
to a marker set by whole words, and the match is shown in brackets.

A type with no match is greyed out and cannot be ticked. That is not a
failure of the annotation; it means this marker source has nothing for it.
Try another source, or supply your own.

Four types are ticked by default. That is not a limitation, it is the
readable size: everything at once is dozens of rows against a dozen columns
and tells you nothing.

## Marker sources

- **Built-in curated** — around a dozen broad types at five markers each.
  Small, canonical, no lookup.
- **PanglaoDB** — by organ. Broader and more citable; panels run long, so
  the per-type cap matters.
- **CellMarker2** — by tissue and species.
- **From file…** — any GMT, which is just `name<TAB>description<TAB>genes`.
  Use this for a panel from a paper, or for the bundled gene-family sets in
  `kosmic/reference/gene_sets/`.

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
