# Figures

Figures drawn from the analyses in the other workspaces, with a preview, size controls and export to PNG, PDF or SVG.

## Available figures

The figure list is grouped by the analysis each figure draws on. A figure whose results are not available shows a message naming the step to run.

| Group | Figures | Draws on |
|---|---|---|
| Embeddings | [UMAP / t-SNE](umap.md) | The dataset loaded in the scRNA workspace |
| Cell Types | [Marker Dot Plot](marker_dotplot.md), [Marker Score UMAPs](marker_scores.md), [Cell-Type Composition](composition.md) | The loaded dataset and its cell-type labels |
| Single-cell QC | [QC Violins](qc_violins.md), [Variable Genes](variable_genes.md) | The loaded dataset |
| Clustering | [PCA Elbow](elbow.md) | The loaded dataset, after PCA |
| Differential Expression | [Volcano Plot](volcano.md), [Gene Heatmap](gene_heatmap.md), [DE Genes by Pathway](topde_summary.md), [Dataset Summary](dataset_summary.md) | Gene DE in the DE workspace |
| Pathway Activity | [Pathway Bubble](pathway_bubble.md), [Patient Dotplot](patient_dotplot.md) | Pathway DE (scoring) in the DE workspace |
| Enrichment | [Enriched Terms](enriched_terms.md), [Enriched Term Links](term_links.md) | An over-representation run (GO or an Enrichr library) from the DE workspace or the Meta-Analysis Enrichment step |
| Enrichment | Pathway Enrichment (GSEA) | Preranked GSEA on the Pathway DE page of the DE workspace |
| Enrichment | Library Enrichment (GSEA), Library GSEA (running enrichment), [Leading-Edge Heatmap](leading_edge.md) | A preranked GSEA run against a library, from the DE workspace's Enrichment page or the Meta-Analysis Enrichment step |
| Meta-Analysis | [Gene and Pathway Meta-Analysis Volcano](meta_volcano.md), [Gene and Pathway Forest Plot](meta_forest.md), [Leave-One-Out Sensitivity](meta_loo.md), [Study Agreement (Vote Levels)](meta_vote_levels.md) | A gene or pathway meta-analysis in the Meta-Analysis workspace |

## Exporting

**Export figure…** saves the current figure as PNG, PDF or SVG. Choose the format in the **Format** list above the button; the save dialog opens with that format selected, and a different file type chosen in the dialog takes precedence. The format applies to every figure page and is remembered between sessions. Use the size controls to adjust dimensions before export.

Exported figures have a white background in both the light and dark themes. Text in PDF and SVG files stays editable in a vector graphics editor. Scatter plots with 20,000 points or more, such as UMAPs, are embedded as an image inside PDF and SVG files so that the files remain small; axes, labels and other text stay vector. The resolution of PNG files and of these embedded images is the export DPI in Plot Settings (default 300).
