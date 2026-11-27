# ADR-019: A vectorised per-study fit as the fast case-control permutation engine

- Status: Proposed (Calum approved integrating the engine, 2026-11-23; this
  record is for review before the stages below are merged)
- Extends: ADR-018 (exact case-control permutation), section 5

## Problem

ADR-018 added the exact case-control (CC) permutation, which re-runs each
study's own DE (`fit_deseq2_wald`, PyDESeq2) on every shuffle. It is exact
and slow: about 300 s per permutation for five studies on one core, so 1,000
permutations take about 6 hours on 14 cores for one cell type. The default
fast engine, `deseq2_fast`, is quick (0.7 s per permutation) but is an
approximation: it fits a condition-only model and reuses dispersions fitted
on the observed labels. Against the exact permutation on the DCM endothelial
cells (five studies, 1,000 shuffles, the same shuffles for each engine) it
calls fewer genes at BH FDR 5%: DL 4,251 against 4,350, REML 4,158 against
4,287, top-50% SumRank 178 against 291, gwOP 5,492 against 5,951, Fisher
6,790 against 7,447, Stouffer 6,306 against 6,612.

ADR-018 section 5 allows the fast engine to be replaced by an engine that
re-estimates the dispersions on every permutation, if its null and its FDR
5% counts lie within the 95% bootstrap interval of the exact permutation.

## Decision

**A new engine, `deseq2_vec`, fits each study's model on every shuffle with
array code over genes, following PyDESeq2 0.5.4 step by step. It takes the
place of `deseq2_fast` as the *fast* choice in the Gene Meta-Analysis page. The
exact engine is kept as the reference. `deseq2_fast.py` stays, named
*approximate*, for the pathway CC permutation, leave-one-study-out and as the
fallback where the vectorised engine cannot run; the selector `deseq2_fast`
keeps its meaning.**

### 1. The engine

- `kosmic/meta_analysis/vec_deseq2.py`, `fit_vec`: size factors, rough and
  moment dispersions, gene-wise Cox-Reid dispersions, parametric trend with
  the mean-trend fallback, prior variance, MAP dispersions, fold changes by
  IRLS, Cook's distances with outlier replacement and refit, the Wald test
  and the Cook's p-value filter, for every gene at once. It returns the
  unshrunk log2 fold change, its standard error and the Wald p-value, as
  `fit_deseq2_wald` does.
- `kosmic/meta_analysis/vec_lbfgsb.py` runs scipy's own compiled L-BFGS-B in
  lockstep over genes, with the objective and gradient answered for all genes
  at once. PyDESeq2 runs the same optimiser gene by gene, so the optimiser's
  steps and stopping decisions are unchanged. It calls
  `scipy.optimize._lbfgsb.setulb`, a private function whose signature has
  changed between scipy releases.
- The engine uses the same shuffle, covariates, usable-covariate rule, gene
  filter, resume folder and manifest as the exact engine. It shares
  `exact_study_state` and differs only in the function that fits one study on
  one shuffled labelling.

### 2. Availability and drift

- `vec_lbfgsb.available()` runs a small problem through the lockstep
  optimiser once and compares the answer with `scipy.optimize.minimize`. If
  scipy's private interface is missing or has changed, the engine reports
  itself unavailable; the Gene Meta-Analysis page then offers the exact
  engine only and says why.
- `tests/test_vec_deseq2.py` compares `fit_vec` with `fit_deseq2_wald` on
  simulated counts with and without a covariate, so a PyDESeq2 release that
  changes its fit fails a test instead of changing results unnoticed. The
  provenance records the pydeseq2 and scipy versions.

### 3. Where it is offered

- Calibration choice *CC permutation (fast)* on the Gene Meta-Analysis page
  runs `deseq2_vec` for DESeq2 input, saving its permutations in
  `cc_permutation_vec/` and its results under the file tag `ccvec`. *CC
  permutation (exact, slow)* is unchanged. With Welch input the permutation
  keeps its Welch engine.
- A study whose pseudobulk file has no covariate columns that its DE run used
  is refused, as in the exact engine. `deseq2_fast` accepted such a study; the
  analytical p-values still do.
- Permutations default to 1,000 for `deseq2_vec` and 200 for the exact engine.

## Evidence

On the DCM endothelial cells (old data, five studies, 1,000 shuffles, the same
shuffles for every engine, drawn over all donors without stratification;
`dev/validation/`), genes at BH FDR 5% for the
exact and the vectorised engine were identical for all seven methods (DL
4,350, REML 4,287, SumRank 1,440, top-50% SumRank 291, gwOP 5,951, Fisher
7,447, Stouffer 6,612), with a difference of 0 in each 95% paired-bootstrap
interval, and the null tail fractions at p < 0.05, 0.01 and 0.001 equal to
within 1e-5. Time per permutation was 30 s for the vectorised engine (DE and
pooling, one core) and 307 s for the exact engine.

Stage 5 repeats this on the final DCM cardiomyocytes (143 donors) with the
integrated code before the engine is named in any result.

## Compatibility

- Saved calibrated runs are read as before; the provenance of a new run names
  `deseq2_vec`. Runs made with `deseq2_fast` keep that name and the note that
  it is an approximation.
- No file format changes.

## Consequences

- A calibrated p-value from the same model as the observed one takes minutes
  per permutation set instead of hours.
- A dependency on a private scipy function, guarded as in section 2.
- About 6 modules change: `kosmic/meta_analysis/cc_permutation.py`, two new
  modules, the Gene Meta-Analysis page and its worker, the help and Methods
  pages, CODEBASE.md.

## Stages

Each stage ends with `pytest tests/` and `ruff check kosmic/ main.py
tests/` passing.

1. **Engine modules.** `vec_lbfgsb.py` and `vec_deseq2.py`. Tests: lockstep
   L-BFGS-B equals scipy's `minimize`; `fit_vec` equals `fit_deseq2_wald` on
   simulated counts with and without a covariate; `available()` is false when
   the private function is replaced by a stub.
2. **Wiring.** `deseq2_vec` in `consensus_cc_permutation`, the engine selector
   and the Gene Meta-Analysis page. Tests: the engine on a small fixture
   equals the exact engine on the same shuffles; a resumed run equals an
   uninterrupted one; a missing covariate column stops with the message.
3. **Documentation.** CODEBASE.md, the Methods and Gene Meta-Analysis help
   pages, the glossary, `main.py` TODO item 12, this ADR's status.
4. **Reference validation script.** `dev/validation/cc_exact_vs_vec.py`
   reproduces the comparison above from a project's pseudobulk files.
5. **Final-data validation.** The exact and vectorised engines on the final
   DCM cardiomyocytes, recorded under `dev/validation/`.
