# Pathway Enrichment

Identify over-represented pathways among your significant DE genes.

## What it does

Enrichment analysis tests whether your significant DE genes are enriched for membership in known biological pathways or gene sets. This is the interpretation endpoint of discovery mode.

## GO enrichment

Gene Ontology enrichment uses local GO annotations. The 'elim' algorithm is DAG-aware and removes redundant parent terms. 'Fisher + BH' is the standard approach with FDR correction.

## Enrichr libraries

Over-representation analysis against pathway databases downloaded from Enrichr (KEGG, Reactome, MSigDB Hallmark, WikiPathways, etc.). Uses Fisher's exact test with BH correction.

Enrichr libraries are downloaded each time the enrichment runs. Enrichr keeps no version history, so a library served under the same name could change between two runs.

## Bundled libraries

**MSigDB: Hallmark 2026.2 (human)** is the MSigDB Hallmark collection, release 2026.2.Hs (`h.all.v2026.2.Hs.symbols.gmt`, 50 gene sets), bundled with KOSMIC and used exactly as published; no symbol is updated. Tables and figures show its set names without MSigDB's `HALLMARK_` prefix and with spaces for underscores (*OXIDATIVE PHOSPHORYLATION*); the exported results and the provenance keep the names as published. Because the file is part of KOSMIC, every run uses the same gene sets, with or without an internet connection. Cite Liberzon et al., *Cell Systems* 2015.

**MSigDB: Hallmark 2026.2 (mouse)** is MSigDB's own mouse edition of the same release, 2026.2.Mm (`mh.all.v2026.2.Mm.symbols.gmt`), with each set's genes given as mouse symbols mapped by MSigDB. Use it for mouse data. Matching a mouse dataset against the human file works for genes whose mouse symbol is the human one in different case, but misses those whose mouse ortholog is named differently, such as the MHC genes (*HLA-DRA* against *H2-Ea*), which sit in the interferon and inflammatory sets.

## What is recorded

Each run is recorded in the provenance, so the methods state which gene sets were used: the method, the library, where it came from (the bundled file and its release, or the Enrichr address), a SHA-1 fingerprint of its content, the number of gene sets and genes, the time it was loaded, the thresholds, and how many terms were tested and significant. A SHA-1 that differs between two runs of the same Enrichr library shows that Enrichr changed it in between. A DE-workspace run is recorded in the study's provenance and a meta-analysis run in the selection's; each keeps one entry per method and library, so running a second library adds an entry rather than replacing the first.

An enrichment run describes the DE result or pooled gene list it was run on. When Gene DE, or the gene meta-analysis, is run again, the enrichment runs made before it no longer describe the current result: the methods leave them out and say how many were left out. They remain in the provenance file. Run the enrichment again to include it.

## Gene symbol matching

Gene symbols are matched without regard to case, because the libraries use human upper-case symbols. Mouse genes therefore match their human counterparts when the two symbols are spelled the same (for example, Col4a1 and COL4A1). Mouse genes whose human ortholog has a different symbol are not matched. The Genes column lists the genes as they are named in your data.

## Thresholds

FDR and |log2FC| cutoffs define which genes from Gene DE are considered 'significant' for enrichment. The P threshold controls which enrichment results are highlighted as significant.
