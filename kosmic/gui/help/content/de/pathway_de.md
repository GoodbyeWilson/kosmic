# Pathway Differential Expression

Test entire pathways for coordinated expression changes between conditions.

## What it does

Pathway DE scores each pathway as a whole, testing whether its member genes show coordinated up- or down-regulation between conditions. This complements single-gene DE by detecting subtle shifts across many genes that might not reach significance individually.

## Pathway scoring

KOSMIC computes a per-cell pathway activity score (mean expression of pathway genes), then tests for differences between conditions using pseudobulk aggregation.

## Which cells are scored

Pathway DE uses the DE run chosen under **Results for** at the top of the page. For the whole-study run it scores every cell of the loaded dataset. For a cell type's run it scores only that cell type's cells, taken from the cell-type column the run was made with, and keeps the donors that run kept (its minimum cells and transcripts per donor). The pathway gene plots and the contamination check use the same cells.

The results are saved in that run's folder, `results/de_analysis/pathway_scoring/<counts>/` for the whole-study run and `.../<counts>/<cell type>/` for a cell type, so each cell type keeps its own pathway results. Choosing another run clears the results shown and shows that run's saved results, if it has any.

## Cascade plot

The cascade plot ranks pathways by their test statistic, showing which pathways are most strongly up- or down-regulated. Error bars indicate uncertainty.
