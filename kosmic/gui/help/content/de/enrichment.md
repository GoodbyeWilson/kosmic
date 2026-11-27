# Pathway Enrichment

Which pathways change in disease, tested with preranked GSEA.

## What it does

Every tested gene is ranked by its test statistic: in the DE workspace DESeq2's Wald statistic, the log2 fold change divided by its standard error; in the meta-analysis Enrichment step the pooled z-statistic. For each gene set in the chosen library, preranked GSEA (Subramanian et al., *PNAS* 2005, as implemented by GSEApy) asks whether its genes gather towards the top of the ranking (up in disease) or the bottom (down), and gives a normalised enrichment score (NES), a p-value and an FDR. Every gene is used, so no significance cut-off is needed.

Ranking by the test statistic, rather than the fold change alone, is the usual choice for DESeq2 results: a gene measured with little precision, such as one with few counts, has a large standard error, so a noisy fold change does not carry it to either end of the ranking. A DE result saved by a version of KOSMIC that did not keep the test statistic cannot be ranked; run the DE again.

The page has three settings: the library, the gene-set size, and the FDR below which a set is called significant.

## Libraries

**MSigDB: Hallmark 2026.2 (human)**, the default, is the MSigDB Hallmark collection, release 2026.2.Hs (`h.all.v2026.2.Hs.symbols.gmt`, 50 gene sets), bundled with KOSMIC and used exactly as published; no symbol is updated. Tables and figures show its set names without MSigDB's `HALLMARK_` prefix and with spaces for underscores (*OXIDATIVE PHOSPHORYLATION*); the exported results and the provenance keep the names as published. Because the file is part of KOSMIC, every run uses the same gene sets, with or without an internet connection. Cite Liberzon et al., *Cell Systems* 2015.

**MSigDB: Hallmark 2026.2 (mouse)** is MSigDB's own mouse edition of the same release, 2026.2.Mm (`mh.all.v2026.2.Mm.symbols.gmt`), with each set's genes given as mouse symbols mapped by MSigDB. Use it for mouse data. Matching a mouse dataset against the human file works for genes whose mouse symbol is the human one in different case, but misses those whose mouse ortholog is named differently, such as the MHC genes (*HLA-DRA* against *H2-Ea*), which sit in the interferon and inflammatory sets.

**Enrichr libraries** (KEGG, Reactome, WikiPathways, BioPlanet, the Gene Ontology libraries and an older Hallmark release) are downloaded from Enrichr each time the analysis runs. Enrichr keeps no version history, so a library served under the same name could change between two runs; the SHA-1 recorded with each run shows whether it did.

## Which gene sets are tested

A gene set is tested only when enough of its genes are in the ranking. **Min genes per set** and **Max genes per set** set the bounds; the defaults, 15 and 500, are those of the Broad Institute's GSEA and of GSEApy. A smaller set gives an unstable enrichment score; a larger one is too broad to interpret. When few genes are tested, as in a small cell population, many sets fall below 15 and are not tested: the status line gives how many of the library's sets were tested. The FDR is computed over the sets tested.

The a-priori panel on the Pathway DE page tests the pathways you chose whatever their size.

## P-values and FDR

Each gene set's p-value and FDR are estimated from 1,000 permutations of the gene sets. When none of the permutations reaches the gene set's observed enrichment score, the estimate is 0, which means the value is below 0.001 rather than zero. The results table, the bar plot and the figures show such a value as < 0.001; the saved results table keeps the 0. **FDR threshold** decides which sets are shown as significant.

## What is recorded

Each run is recorded in the provenance, so the methods state which gene sets were used: the library, where it came from (the bundled file and its release, or the Enrichr address), a SHA-1 fingerprint of its content, the number of gene sets and genes, the time it was loaded, what the genes were ranked by, the gene-set size bounds, the FDR threshold, and how many sets were tested and significant. A DE-workspace run is recorded in the study's provenance, with the cell type and counts of the DE run it used, and a meta-analysis run in the selection's. Each keeps one entry per library, and in the DE workspace per cell type and counts, so running a second library or another cell type adds an entry rather than replacing the first.

A run describes the DE result or pooled gene list it was run on. When Gene DE, or the gene meta-analysis, is run again, the runs made before it no longer describe the current result: the methods leave them out and say how many were left out. They remain in the provenance file. Run the enrichment again to include it.

## What is saved

Each run's results table is saved with the running-enrichment curves beside it, so the Figures workspace can draw it without the run in memory. A DE-workspace run is saved in the folder of the DE run chosen under **Results for**: `results/de_analysis/pathway_scoring/<counts>/` for a whole-study run and `results/de_analysis/pathway_scoring/<counts>/<cell type>/` for a cell type's run, so each cell type keeps its own results; each run replaces the previous one in its folder. A meta-analysis run is saved in the selection's folder under `meta_analysis/`, per pooling run and library, as `<run>_gsea_<library>.csv` beside the pooling run it ranked (`consensus_reml_analytical_gsea_MSigDB_Hallmark_2026.2_Hs.csv`), so a GSEA of REML's pooled ranking and one of Stouffer's are kept apart.

Choosing another run under **Results for** clears the results shown and shows that run's saved GSEA, if it has one.

Earlier versions of KOSMIC also offered over-representation analysis (ORA) of the significant genes. A project's saved ORA results are still drawn in the Figures workspace.

## Gene symbol matching

Gene symbols are matched without regard to case, because the libraries use human upper-case symbols. Mouse genes therefore match their human counterparts when the two symbols are spelled the same (for example, Col4a1 and COL4A1). Mouse genes whose human ortholog has a different symbol are not matched. The leading-edge genes are listed as they are named in your data.
