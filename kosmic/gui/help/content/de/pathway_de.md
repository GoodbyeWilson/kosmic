# Pathway Differential Expression

Test the pathways you chose for a change between disease and control donors.

## What it does

Each donor gets one score per pathway. The donor's counts are summed over its cells, divided by DESeq2's median-of-ratios size factors (computed by pydeseq2, from every gene) and log2(x + 1) transformed, which is DESeq2's normTransform, the same normalisation as the gene DE; the score is the average of the pathway's genes. Disease and control donors are compared with Welch's t-test. The difference in mean score is in log2 units: it is the average of the pathway's genes' log2 differences.

This is a self-contained test: it asks whether the pathway's genes are higher or lower in disease. The **GSEA of these pathways** tab asks a different, competitive question: whether the pathway's genes stand out from the rest of the transcriptome. Both are standard (Goeman & Bühlmann, *Bioinformatics* 2007); when they agree, the result is solid.

The test does not adjust for covariates such as sex, which the gene DE may do; the methods say so.

## Settings

| Setting | Effect |
|---|---|
| **Gene detection ≥** (5%) | A gene counts towards a pathway's score only if it is detected (non-zero count) in at least this percentage of cells in disease or in control. |
| **Min pathway coverage ≥** (30%) | A pathway is scored only if at least this percentage of its genes pass the detection filter, and at least 3 genes do. **Genes measured** lists each pathway's numbers and why any was not scored. |

A pathway that is not scored has no score: no other score is put in its place.

## Above the results

**Results for** chooses the gene DE run, and **Pathway** the pathway the tabs show. Each pathway in the list gives its direction and FDR, most significant first; when the gene sets have parent categories, these follow.

## Tabs

| Tab | Shows |
|---|---|
| **Pathway scores** | The picked pathway's score for each donor, a box per group, with its difference, p-value and FDR. |
| **GSEA of these pathways** | Preranked GSEA of the chosen pathways: every tested gene ranked by its DESeq2 Wald statistic, and a pathway tested when 15-500 of its genes are in the ranking (GSEA's defaults, as on the Enrichment page). The bar chart gives each pathway's NES; the curve is the picked pathway's. Pathways with too few genes in the ranking are listed as not tested. |
| **Ambient RNA check** | Whether the picked pathway's change could come from ambient RNA of other cell types: each donor's score against its contamination (a QC contamination panel, or DecontX's estimate), and the disease effect before and after adjusting for contamination. |
| **Genes measured** | For every chosen pathway, how many of its genes are annotated, present and detected, and whether it was scored. |

**Export CSV** saves every pathway's result as a table.

## Which cells are scored

Pathway DE uses the DE run chosen under **Results for**. For the whole-study run it scores every cell of the loaded dataset. For a cell type's run it scores only that cell type's cells, taken from the cell-type column the run was made with, and keeps the donors that run kept (its minimum cells and transcripts per donor). The plots and the ambient RNA check use the same cells.

The results are saved in that run's folder, `results/de_analysis/pathway_scoring/<counts>/` for the whole-study run and `.../<counts>/<cell type>/` for a cell type, so each cell type keeps its own pathway results. Choosing another run clears the results shown and shows that run's saved results, if it has any.

## Earlier versions

Earlier versions offered seven other scores (DESeq2 VST, z-scored VST, CPM, and four per-cell scores tested after normality and variance pre-tests). They are removed; the methods still name them for runs made with them.
