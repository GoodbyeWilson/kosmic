# Meta-Analysis Pooling Methods -- Design Notes

## Overview

KOSMIC implements multiple meta-analysis pooling methods for combining
gene-level DE results across studies. This document records the rationale
behind each method, the permutation strategies used, and the benchmarking
results that informed the default choice.

---

## Two families of methods

### Effect-size methods (DerSimonian-Laird family)

These pool log2 fold changes across studies using inverse-variance weighting.
They produce a pooled effect size, SE, and confidence interval per gene.

- **DL Random Effects** -- standard DerSimonian-Laird. Z-test p-value.
- **DL REML** -- REML tau-squared estimator instead of DL moment estimator.
  Less biased with few studies. Matches R's metafor default.
- **Hartung-Knapp** -- DL estimate with t-distribution CI. More conservative.
  Unreliable with k < 10 studies (benchmarks confirm this).
- **Hunter-Schmidt** -- sample-size weighting instead of inverse-variance.
  Poor performance with few studies. SE formula was corrected (divide by k,
  not by total N).
- **Permutation-Calibrated RE** -- DL effect estimate with gene-label
  permutation-calibrated p-value. Best overall performer in benchmarks.
- **DL Top-50% Proportion** -- for each gene, selects the top 50% of studies
  by p-value, runs DL on those. Permutation-calibrated. Conservative.
- **DL GeneGate Top-50%** -- for each study, keeps only genes in the top 50%
  by |log2FC|, then runs DL across studies where the gene passed. A novel
  within-study quality gate approach. Permutation-calibrated.

### Rank-based methods (SumRank family)

These rank genes within each study by signed -log10(p-value), normalise
to [0,1], and sum ranks across studies. They do not produce effect size
estimates -- only consistency scores and p-values.

- **SumRank (Irwin-Hall)** -- analytical p-value from the Irwin-Hall CDF
  (sum of k uniform[0,1] variables). Includes k/2 cap per Nakatsuka et al.
- **SumRank (Permutation)** -- same statistic, permutation-calibrated p-value.
  Results are nearly identical to Irwin-Hall (validates the k/2 cap).
- **SumRank Top-50% (Permutation)** -- sums only the best ceil(k*0.5)
  dataset ranks per gene. Permutation-calibrated.

---

## Permutation strategies

### Gene-label permutation (used by DL methods)

**What it does:** Within each study's DE results, shuffle which gene name
is attached to which effect size + SE. Run the pooling method on the
shuffled data.

**Null hypothesis:** "Given the distribution of effect sizes in each study,
how often would this gene's pooled Z be this large if the gene-to-effect
mapping were random?"

**Why it works for DL:** The null preserves the marginal distribution of
effect sizes within each study. Some null genes randomly inherit large
effects from multiple studies, creating a realistic tail in the null |Z|
distribution. This gives proper calibration across the full p-value range.

**Why it doesn't work for SumRank:** Gene-label permutation shuffles which
gene gets which effect, but the *ranks* within each study stay the same
(the same set of values are just reassigned). The SumRank statistic only
depends on relative ranks, so the null distribution is identical to the
observed distribution -- no discrimination.

### Case-control permutation (used by SumRank, per Nakatsuka et al. 2025)

**What it does:** Shuffle case/control labels at the sample level within
each study, re-run DE from scratch (e.g. DESeq2), then run the pooling
method on the permuted DE results.

**Null hypothesis:** "If there were no true disease effect in any study,
what distribution of pooling statistics would I see?"

**Why it works for SumRank:** Under the null, no genes are systematically
top-ranked across studies. The null rank sums have a genuine null
distribution with structure (genes still get ranks 1 to N). This allows
proper calibration of the observed rank sums.

**Why it doesn't work well for DL:** Under CC permutation, ALL effects
are destroyed. Null log2FCs are ~1e-6, null |Z| values are near zero.
The observed |Z| for any real signal (3-10+) exceeds every null value,
giving p = 1/(n_perms+1) for all DE genes. The result is binary
(significant/not) with no useful gradation -- equivalent to FWER control
rather than per-gene calibration.

### Summary

| Permutation type | Good for | Bad for | Reason |
|------------------|----------|---------|--------|
| Gene-label | DL (effect-size methods) | SumRank | Preserves effect distribution, breaks gene mapping |
| Case-control | SumRank (rank methods) | DL | Destroys all signal, null |Z| degenerate for DL |

---

## Benchmark results (simulated data)

### Simulation setup

- 15,000 genes, 500 truly DE (3.3%)
- True effect sizes: 30% strong (|log2FC| 0.8-2.0), 40% moderate (0.3-0.8),
  30% weak (0.1-0.3), random direction
- 6 studies, 3-8 pseudobulk samples per condition
- 80-400 cells per sample, negative binomial counts
- Between-study heterogeneity: tau = 0.15
- DE via KOSMIC's actual DESeq2 pipeline
- Evaluation: BH FDR < 0.05, AUC-ROC, AUC-PR, Spearman with true log2FC

### Results

```
Method                   Sens   eFDR   Prec  AUC-ROC  AUC-PR  Spearman
-----------------------------------------------------------------------
Perm_Calib_RE            0.470  0.004  0.996  0.896   0.704   0.970
SumRank_Top50_Perm       0.404  0.029  0.971  0.887   0.622   0.936
DL_random                0.350  0.022  0.978  0.896   0.705   0.970
DL_REML                  0.336  0.018  0.982  0.896   0.694   0.970
SumRank_IH               0.310  0.055  0.945  0.866   0.567   0.936
SumRank_Perm             0.310  0.055  0.945  0.866   0.567   0.936
DL_Top50_Perm            0.218  0.000  1.000  0.872   0.659   0.956
DL_GeneGate_Top50        0.182  0.000  1.000  0.877   0.687   0.950
```

### Key findings

1. **Perm_Calibrated_RE is the best overall method.** Highest sensitivity
   (0.470) with near-zero FDR (0.4%), same effect size quality as DL Random
   (identical Spearman/AUC-ROC). It combines DL's inverse-variance-weighted
   effect estimates with gene-label permutation p-values.

2. **SumRank methods have worse effect size ranking.** Spearman 0.936 vs
   0.970 for DL methods. SumRank does not estimate effect sizes -- it only
   measures consistency. For downstream interpretation (forest plots, effect
   magnitudes), DL methods are superior.

3. **SumRank Top-50% has high sensitivity but higher FDR.** 0.404 sensitivity
   but 2.9% FDR vs 0.4% for Perm_Calibrated_RE.

4. **SumRank IH and SumRank Perm are identical.** The k/2 cap (matching
   Nakatsuka et al.) makes the analytical and permutation p-values converge.
   This validates the cap implementation.

5. **Conservative methods (Top50_Perm, GeneGate) achieve zero FDR** but at
   the cost of sensitivity (0.18-0.22). Useful when false positives must be
   eliminated entirely.

6. **Hartung-Knapp is unreliable with few studies.** 13.3% eFDR with k=6.
   Should only be used with k >= 10.

7. **Hunter-Schmidt SE formula was incorrect.** Dividing by total sample
   size N instead of number of studies k inflated precision. After fix,
   method is very conservative (sensitivity 0.150) and not recommended.

---

## Default recommendation

**Perm_Calibrated_RE** is the default pooling method in KOSMIC for gene-level
meta-analysis. It provides the best balance of sensitivity, FDR control,
and effect size estimation quality.

For users who want maximum confidence with zero false positives,
**DL GeneGate Top-50%** is recommended as a secondary analysis.

---

## Novel contributions

### DL GeneGate (within-study gene filtering)

A novel approach not present in the SumRank paper or existing meta-analysis
software. Instead of selecting the top 50% of *studies* per gene (which
cherry-picks which studies to use), it filters the top 50% of *genes* within
each study by |log2FC|. This is a quality gate: a gene only contributes
from studies where it has a strong signal.

Key properties:
- No study cherry-picking: each study contributes its strongest signals
- A gene consistently in the top half across all studies gets full k-study power
- Permutation-calibrated to account for the selection
- Zero empirical FDR in benchmarks

### Permutation-Calibrated RE (vectorized)

The permutation-calibrated random-effects approach itself is not novel
(see Follmann & Proschan 1999), but our implementation is:
- Vectorized numpy DL computation on (n_genes x n_datasets) matrices
- Null |Z| values grouped by k for proper per-gene calibration
- Parallel permutation batches via ProcessPoolExecutor
- ~1000x faster than naive per-gene-loop + per-permutation implementation

---

## Limitations of current benchmarks

- Single simulation replicate (no parameter sweep yet)
- Fixed k=6 studies (real-world ranges from 2-20+)
- Fixed tau=0.15 (moderate heterogeneity only)
- Fixed DE fraction 3.3% (real varies 1-10%+)
- Synthetic count generation (not semi-synthetic from real data)
- All cells treated as one type (real analysis is cell-type stratified)
- Balanced study designs (real data often unbalanced)

These would need to be addressed for a methods paper but are sufficient
for justifying default choices in the KOSMIC software paper.

---

## References

- DerSimonian R, Laird N (1986). Meta-analysis in clinical trials.
  Controlled Clinical Trials, 7(3), 177-188.
- Nakatsuka N et al. (2025). A reproducibility focused meta-analysis method
  for single-cell transcriptomic case-control studies uncovers robust
  differentially expressed genes. Nature Communications.
- Follmann DA, Proschan MA (1999). Valid inference in random effects
  meta-analysis. Biometrics, 55(3), 732-737.
- Viechtbauer W (2010). Conducting meta-analyses in R with the metafor
  package. Journal of Statistical Software, 36(3), 1-48.
