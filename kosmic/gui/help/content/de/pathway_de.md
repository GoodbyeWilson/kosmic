# Pathway Differential Expression

Test entire pathways for coordinated expression changes between conditions.

## What it does

Pathway DE scores each pathway as a whole, testing whether its member genes show coordinated up- or down-regulation between conditions. This complements single-gene DE by detecting subtle shifts across many genes that might not reach significance individually.

## Pathway scoring

KOSMIC computes a per-cell pathway activity score (mean expression of pathway genes), then tests for differences between conditions using pseudobulk aggregation.

## Cascade plot

The cascade plot ranks pathways by their test statistic, showing which pathways are most strongly up- or down-regulated. Error bars indicate uncertainty.
