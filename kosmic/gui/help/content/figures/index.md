# Figures

Figures drawn from the analyses in the other workspaces, with a preview, size controls and export to PNG, PDF or SVG.

## Which study and which results

The scRNA and DE figures show one study, named above the switch (*Study: Chaffin*). It is the active study when Figure Export is opened from the left-hand rail, or the study you chose with **Open in Figures** on the Project page. Open in Figures does not change the active study and does not load the study into any workspace.

Each figure draws on the workspace that holds the study when one does, so a change not yet saved, such as an open subclustering, is shown. Otherwise it reads the study's folder: the processed dataset file and the results the DE steps saved. You do not need to open the scRNA or DE workspace, or keep a study loaded, to export its figures. A large study's file takes a few seconds to read; until then its scRNA figures are greyed with *Reading* followed by the file name.

When a study has more than one saved DE run, for example one per cell type, a **Run** list appears above the DE figures. It lists the run in the DE workspace, when it holds one for this study, and then each saved run by cell type and method. The Meta figures have their own **Run** list: the run in the Meta-Analysis workspace, then each meta-analysis run saved in the project, by selection and methods.

A DE figure drawn from saved results that are older than the study's dataset file carries the note *These results were made before the study's file was last saved*, with the date. A later step, such as relabelling cell types, may have changed the data since; run the DE step again to update the figure.

Some results are saved only by DE runs made with this version or later: the gene set selection of a gene DE run (needed by DE Genes by Pathway, Dataset Summary and the Gene Heatmap), over-representation results and the running-enrichment curves of a library GSEA run. For an earlier run these figures stay greyed until the step is run again.

## Available figures

The switch above the figure list shows the figures of one workspace at a time: **scRNA**, **DE** or **Meta**. Figure Export opens on the section of the workspace you came from. Within a section the figures are grouped by the analysis they draw on, as in the table below; the Enrichment figures are in the DE section and also draw an enrichment run from the Meta-Analysis Enrichment step.

A figure whose results are not available is greyed in the list, and its tooltip names the step to run, or says why the saved results cannot be used; selecting it shows the same message. Tick **Show available only** to hide these figures.

| Group | Figures | Draws on |
|---|---|---|
| Embeddings | [UMAP / t-SNE](umap.md) | The study's dataset, with an embedding |
| Cell Types | [Marker Dot Plot](marker_dotplot.md), [Marker Score UMAPs](marker_scores.md), [Cell-Type Composition](composition.md) | The study's dataset and its cell-type labels |
| Single-cell QC | [QC Violins](qc_violins.md), [Variable Genes](variable_genes.md) | The study's dataset |
| Clustering | [PCA Elbow](elbow.md) | The study's dataset, after PCA |
| Differential Expression | [Volcano Plot](volcano.md), [Gene Heatmap](gene_heatmap.md), [DE Genes by Pathway](topde_summary.md), [Dataset Summary](dataset_summary.md) | A gene DE run of the study |
| Pathway Activity | [Pathway Bubble](pathway_bubble.md), [Patient Dotplot](patient_dotplot.md) | The study's pathway DE (scoring) |
| Enrichment | [Enriched Terms](enriched_terms.md), [Enriched Term Links](term_links.md) | An over-representation run (GO or an Enrichr library) from the DE workspace or the Meta-Analysis Enrichment step |
| Enrichment | Pathway Enrichment (GSEA) | Preranked GSEA on the Pathway DE page of the DE workspace |
| Enrichment | Library Enrichment (GSEA), Library GSEA (running enrichment), [Leading-Edge Heatmap](leading_edge.md) | A preranked GSEA run against a library, from the DE workspace's Enrichment page or the Meta-Analysis Enrichment step |
| Meta-Analysis | [Gene and Pathway Meta-Analysis Volcano](meta_volcano.md), [Gene and Pathway Forest Plot](meta_forest.md), [Leave-One-Out Sensitivity](meta_loo.md), [Study Agreement (Vote Levels)](meta_vote_levels.md) | A gene or pathway meta-analysis, in the Meta-Analysis workspace or saved in the project |

The forest and leave-one-out figures of a saved meta-analysis run read each pooled study's DE table from the project. When one of those tables is no longer in the project, these two figures are greyed and the tooltip names the study.

## Making figures on another computer

The figures can be made on a computer other than the one that ran the analysis, for example a laptop at home after the analysis ran on a workstation in the lab. Copy the project folder to the other computer and open it in KOSMIC; nothing is computed again. The `raw_data/` folders are not needed and can be left out of the copy. The figures read:

| Section | Reads |
|---|---|
| scRNA | The study's `processed_data/` (its dataset file) and `study.json` |
| DE | The study's `results/de_analysis/` and its dataset file (the Gene Heatmap reads the counts of its gene sets' genes) |
| Meta | The project's `meta_analysis/` and each pooled study's `results/de_analysis/statistics/` |

Only the parts of a dataset file that a figure needs are read, so a large atlas can be exported on a computer with 16 GB of memory.

## Exporting

**Export figure…** saves the current figure as PNG, PDF or SVG. The save dialog opens in the study's `figures/` folder for the scRNA figures, its `results/figures/` folder for the DE figures, and the selection's folder under `meta_analysis/` for the Meta figures. Choose the format in the **Format** list above the button; the save dialog opens with that format selected, and a different file type chosen in the dialog takes precedence. The format applies to every figure page and is remembered between sessions. Use the size controls to adjust dimensions before export.

Exported figures have a white background in both the light and dark themes. Text in PDF and SVG files stays editable in a vector graphics editor. Scatter plots with 20,000 points or more, such as UMAPs, are embedded as an image inside PDF and SVG files so that the files remain small; axes, labels and other text stay vector. The resolution of PNG files and of these embedded images is the export DPI in Plot Settings (default 300).
