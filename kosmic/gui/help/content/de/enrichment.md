# Pathway Enrichment

Identify over-represented pathways among your significant DE genes.

## What it does

Enrichment analysis tests whether your significant DE genes are enriched for membership in known biological pathways or gene sets. This is the interpretation endpoint of discovery mode.

## Enrichr libraries

Over-representation analysis against pathway databases downloaded from Enrichr (KEGG, Reactome, MSigDB Hallmark, WikiPathways, the Gene Ontology libraries, etc.). Uses Fisher's exact test with BH correction.

Enrichr libraries are downloaded each time the enrichment runs. Enrichr keeps no version history, so a library served under the same name could change between two runs.

## Bundled libraries

**MSigDB: Hallmark 2026.2 (human)** is the MSigDB Hallmark collection, release 2026.2.Hs (`h.all.v2026.2.Hs.symbols.gmt`, 50 gene sets), bundled with KOSMIC and used exactly as published; no symbol is updated. Tables and figures show its set names without MSigDB's `HALLMARK_` prefix and with spaces for underscores (*OXIDATIVE PHOSPHORYLATION*); the exported results and the provenance keep the names as published. Because the file is part of KOSMIC, every run uses the same gene sets, with or without an internet connection. Cite Liberzon et al., *Cell Systems* 2015.

**MSigDB: Hallmark 2026.2 (mouse)** is MSigDB's own mouse edition of the same release, 2026.2.Mm (`mh.all.v2026.2.Mm.symbols.gmt`), with each set's genes given as mouse symbols mapped by MSigDB. Use it for mouse data. Matching a mouse dataset against the human file works for genes whose mouse symbol is the human one in different case, but misses those whose mouse ortholog is named differently, such as the MHC genes (*HLA-DRA* against *H2-Ea*), which sit in the interferon and inflammatory sets.

## What is recorded

Each run is recorded in the provenance, so the methods state which gene sets were used: the method, the library, where it came from (the bundled file and its release, or the Enrichr address), a SHA-1 fingerprint of its content, the number of gene sets and genes, the time it was loaded, the thresholds, and how many terms were tested and significant. A SHA-1 that differs between two runs of the same Enrichr library shows that Enrichr changed it in between. A DE-workspace run is recorded in the study's provenance, with the cell type and counts of the DE run it used, and a meta-analysis run in the selection's. Each keeps one entry per method and library, and in the DE workspace per cell type and counts, so running a second library or another cell type adds an entry rather than replacing the first.

An enrichment run describes the DE result or pooled gene list it was run on. When Gene DE, or the gene meta-analysis, is run again, the enrichment runs made before it no longer describe the current result: the methods leave them out and say how many were left out. They remain in the provenance file. Run the enrichment again to include it.

## Ranking genes for GSEA

Preranked GSEA orders every gene by one value. In the DE workspace, **Rank genes by** chooses it:

| Choice | What orders the genes |
|---|---|
| **Log2 fold change** (default) | The effect size alone. |
| **Test statistic** | The log2 fold change divided by its standard error: DESeq2's Wald statistic, or Welch's t for the t-test methods. A gene measured with little precision, such as one with few counts, has a large standard error, so a noisy fold change does not carry it to either end of the ranking. Preranked GSEA of DESeq2 results usually ranks by it. |

A DE result saved by a version of KOSMIC that did not keep the test statistic cannot be ranked by it; run the DE again. The meta-analysis Enrichment step has no choice: it ranks by the pooled z-statistic, the pooled effect divided by its pooled standard error. The ranking used is recorded with each run.

## Which gene sets GSEA tests

A gene set is tested only when enough of its genes are in the ranking. **Min genes per set** and **Max genes per set**, under Term filters, set the bounds; the defaults, 15 and 500, are those of the Broad Institute's GSEA and of GSEApy. A smaller set gives an unstable enrichment score; a larger one is too broad to interpret. When few genes are tested, as in a small cell population, many sets fall below 15 and are not tested: the status line gives how many of the library's sets were tested. The bounds used are recorded with each run. (**Min genes per term** applies to ORA only.)

The a-priori panel on the Pathway DE page tests the pathways you chose whatever their size.

## GSEA p-values and FDR

A library GSEA run estimates each gene set's p-value and FDR from 1,000 permutations of the gene sets. When none of the permutations reaches the gene set's observed enrichment score, the estimate is 0, which means the value is below 0.001 rather than zero. The results table, the bar plot and the figures show such a value as < 0.001; the saved results table keeps the 0.

## What is saved

The results of each run are saved for the Figures workspace, which can then draw them without the run in memory: an over-representation run as `enrichment_results.csv`, and a library GSEA run as its results table with the running-enrichment curves beside it. A DE-workspace run is saved in the folder of the DE run chosen under **Results for**: `results/de_analysis/pathway_scoring/<counts>/` for a whole-study run and `results/de_analysis/pathway_scoring/<counts>/<cell type>/` for a cell type's run, so each cell type keeps its own results. A meta-analysis run is saved in the selection's folder under `meta_analysis/`. Each run replaces the previous one of its kind in its folder, except that a meta-analysis library GSEA run is saved per pooling run and library, as `<run>_gsea_<library>.csv` beside the pooling run it ranked (`consensus_reml_analytical_gsea_MSigDB_Hallmark_2026.2_Hs.csv`), so a GSEA of REML's pooled ranking and one of Stouffer's are kept apart, and the Figures workspace draws each with its own run.

Choosing another run under **Results for** clears the results shown and shows that run's saved results of the chosen method, if it has any. Results that earlier versions saved directly in `pathway_scoring/` are shown for the whole-study run.

## Gene symbol matching

Gene symbols are matched without regard to case, because the libraries use human upper-case symbols. Mouse genes therefore match their human counterparts when the two symbols are spelled the same (for example, Col4a1 and COL4A1). Mouse genes whose human ortholog has a different symbol are not matched. The Genes column lists the genes as they are named in your data.

## Thresholds

FDR and |log2FC| cutoffs define which genes from Gene DE are considered 'significant' for enrichment. They start from the values the DE run shown was made with, each time a different run is shown; values typed here stay until then. The P threshold controls which enrichment results are highlighted as significant.
