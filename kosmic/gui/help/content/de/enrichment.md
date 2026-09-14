# Pathway Enrichment

Identify over-represented pathways among your significant DE genes.

## What it does

Enrichment analysis tests whether your significant DE genes are enriched for membership in known biological pathways or gene sets. This is the interpretation endpoint of discovery mode.

## GO enrichment

Gene Ontology enrichment uses local GO annotations. The 'elim' algorithm is DAG-aware and removes redundant parent terms. 'Fisher + BH' is the standard approach with FDR correction.

## Enrichr libraries

Over-representation analysis against pathway databases downloaded from Enrichr (KEGG, Reactome, MSigDB Hallmark, WikiPathways, etc.). Uses Fisher's exact test with BH correction.

## Thresholds

FDR and |log2FC| cutoffs define which genes from Gene DE are considered 'significant' for enrichment. The P threshold controls which enrichment results are highlighted as significant.
