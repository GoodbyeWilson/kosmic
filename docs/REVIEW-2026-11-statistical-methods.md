# KOSMIC statistical methods review, November 2026

A review of KOSMIC's statistical methods against their published definitions
and standard practice. The rule: a method in KOSMIC is the textbook method, and
any departure from it is either removed or shown to the user and documented.

## Already fixed

| Commit | What |
|---|---|
| 21c57a1 | Stouffer, Fisher and gwOP floored each study's p at 1/n genes and dropped p = 1 studies. Removed: each p is used as it is (1e-300 only to keep an underflow finite). SumRank's z now takes its sign from SumRank's own direction. |
| 168564c | DE-page GSEA can rank by the Wald statistic (log2FC / SE), the usual choice for DESeq2 results. Meta GSEA files are saved per pooling run, so methods no longer overwrite each other's. |
| 46e1f80 | REML I² is tau² / (tau² + s²) with the REML tau², as metafor computes it (was the DerSimonian-Laird (Q - df) / Q). Methods Comparison takes Stouffer's direction from its z. Preranked GSEA tests sets with 15-500 genes in the ranking (GSEA's defaults; was a hidden 3-5000), now a setting. `bh_fdr` counts a p of 0 instead of giving it FDR 1. |
| ee00a55 | SumRank top-50%: the Monte Carlo p is (count + 1) / (n + 1) (Phipson & Smyth 2010), never 0. Before, a gene could get p = 0, which BH, now counting p = 0, would have called significant. HKSJ and REML comments corrected (the truncation is metafor's `adhoc`). Help: GO elim no longer described as topGO's; the claim that DL and REML coincide for two studies removed. |
| 149f12b | Local GO enrichment removed: its "elim" was not topGO's (it removed a significant term's genes from every later term, not only its ancestors, and never shrank the background). GO stays available as the Enrichr GO libraries, with the same Fisher + BH ORA and preranked GSEA as every library. `kosmic/de/go_enrichment.py`, the 37 MB of bundled GO files and the goatools dependency are gone. |
| e14ef49 (#116) | Gene DE runs DESeq2 only: the log-CPM Welch t-test (no covariates, no weighting for low counts) and its home-made empirical Bayes moderation are removed; limma-voom is listed as coming. The Welch test kept for per-donor pathway scores reports the real Welch SE (NaN where undefined), with no SE cap, floor or substituted value. |
| (with this) | The Enrichment page runs preranked GSEA only (library, gene-set size, FDR threshold). Over-representation analysis is removed from it, which settles its background problem (it counted genes DESeq2 never tested); so are the log2FC ranking choice and the detection filter that recomputed the FDR after DE. |

## Still to fix (Kirk to decide how)

None of these is used in the DCM paper.

| # | Area | Where | Problem | Textbook version |
|---|---|---|---|---|
| 2 | Pathway meta-analysis, "DESeq2 + VIF (recommended)" | `kosmic/gui/meta_analysis/pages/pathway_ma_page.py:164-184`, `PATHWAY_VIF_RHO` | Home-made: a fixed inter-gene correlation of 0.1, not CAMERA's estimated one. If it cannot run it falls back to Welch while still recording "DESeq2 + VIF". Mouse gene names never match the human sets. Help page is a TODO. | CAMERA (Wu & Smyth 2012) with the correlation estimated from the data, or a per-donor pathway score. |
| 3 | HKSJ option | `kosmic/meta_analysis/pooling/pooling_core.py:256` | Truncates q at 1, which is metafor's `test="adhoc"`, while the code and labels say `knha`. | Either plain HKSJ (`knha`), or keep the truncation and call it the ad hoc variant (Knapp & Hartung 2003). |
| 4 | gwOP | `kosmic/meta_analysis/pooling/gwop.py:61-165` | The per-gene weights cancel out. The null CDF lookup grid (10,001 points, clipped at 1e-15) makes the strongest p-values about 10^17 times too small. The permutation engine's gwOP differs from the analytical one. | Li & Tseng 2011 adaptively weighted / ordered p as published, with an exact or permutation null. |
| 6 | Fast and Welch permutation engines | `kosmic/meta_analysis/cc_permutation.py` | Donors marked "exclude" are shuffled into the comparison (the exact engine leaves them out). | Permute disease / control labels among included donors only. |
| 7 | Cell-level pathway scoring | `kosmic/de/pathway_scoring.py` | Chooses the test from Shapiro-Wilk and Levene pre-tests (discouraged: the two-stage procedure does not keep its error rate). Writes Cohen's d into the `logfoldchanges` column. | One test chosen in advance (Welch); effect in its own column. |
| 9 | Imported tables | `methods_page.py:981-985`, `enrichment_step.py:53`, `reproducibility.py:49` | Missing raw p-values silently fall back to adjusted ones. edgeR (`PValue`, `FDR`, `logFC`) and limma (`P.Value`, `adj.P.Val`) column names are not recognised. | Map the column names; refuse a table with no raw p. |
| 10 | Consensus FDR scope | `kosmic/meta_analysis/pooling_dispatch.py:213-222` | Each method's consensus FDR is computed over the first method's genes; genes the method did not test count as p = 1. (The Methods Comparison has its own FDR over the genes each method tested, and is what the paper uses.) | BH over the genes each method tested. |
| 11 | GSEA a-priori panel | `kosmic/gui/de_analysis/workers.py` (`FgseaWorker`) | Tests the user's chosen pathways at 3-5000 genes, below GSEA's 15. Kept deliberately for chosen sets; decide whether to warn below 15. | Warn, or use 15 with a visible override. |
| 12 | Permutation calibration, DESeq2 engine | `kosmic/meta_analysis/cc_permutation.py` (around 1517) | When a permutation's PyDESeq2 fit fails it silently falls back to a Welch t-test on log-CPM, so the null can mix two different tests. Check before the paper's supplement calibration. | Skip and count the failed permutation, or stop; never substitute another test. |
| 13 | pandas 3 deprecations | across `kosmic/` (`inplace=`, `mode.chained_assignment`) | Work today but are removed in pandas 4; the test suite shows about 3,000 warnings under pandas 3.1. | Replace before pandas 4. |

## Meanwhile in the app

- Pathway meta-analysis: "DESeq2 + VIF (KOSMIC's own, under review)", no longer
  "recommended"; the tooltip gives the fixed correlation of 0.1.
- HKSJ: the tooltip says it is the ad hoc (truncated) variant.
- Empirical Bayes moderation (t-test): the tooltip says it is not limma's.
- gwOP: not relabelled (its label is also its lookup key in the method
  registry); covered by this list.
