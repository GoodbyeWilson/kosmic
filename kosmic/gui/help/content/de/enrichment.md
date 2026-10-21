# Pathway Enrichment

Identify over-represented pathways among your significant DE genes.

## What it does

Enrichment analysis tests whether your significant DE genes are enriched for membership in known biological pathways or gene sets. This is the interpretation endpoint of discovery mode.

## GO enrichment

Gene Ontology enrichment uses local GO annotations. The 'elim' algorithm is DAG-aware and removes redundant parent terms. 'Fisher + BH' is the standard approach with FDR correction.

## Enrichr libraries

Over-representation analysis against pathway databases downloaded from Enrichr (KEGG, Reactome, MSigDB Hallmark, WikiPathways, etc.). Uses Fisher's exact test with BH correction.

Gene symbols are matched without regard to case, because the Enrichr libraries use human upper-case symbols. Mouse genes therefore match their human counterparts when the two symbols are spelled the same (for example, Col4a1 and COL4A1). Mouse genes whose human ortholog has a different symbol are not matched. The Genes column lists the genes as they are named in your data.

## Thresholds

FDR and |log2FC| cutoffs define which genes from Gene DE are considered 'significant' for enrichment. The P threshold controls which enrichment results are highlighted as significant.
