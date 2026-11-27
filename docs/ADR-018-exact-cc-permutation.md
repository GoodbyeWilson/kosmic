# ADR-018: Exact case-control permutation, and covariates in the pseudobulk file

- Status: Accepted (Calum, 2026-11-19)
- Extends: ADR-007 (meta output per cell type), ADR-015 (DE outputs per counts and cell type)

## Problem

The meta-analysis can calibrate its pooled p-values by case-control (CC)
permutation: the disease and control labels are shuffled within each study,
the per-study differential expression and the pooling are repeated, and each
observed pooled statistic is compared with the null distribution of that
statistic over all genes and permutations. The calibration was compared with
an exact permutation on the DCM endothelial cells (five studies, 200
shuffles, the same shuffles for every engine). The exact permutation re-runs
the real per-study DE on every shuffle. Four problems were found.

1. **The permuted model is not the observed model.** The per-study DE is
   `~sex_inferred + condition`. The fast permutation engine
   (`kosmic/meta_analysis/deseq2_fast.py`) fits `~condition` with a
   two-column design and has no way to take a covariate. Any covariate chosen
   in the DE workspace is ignored in the permutation. The pseudobulk file
   cannot supply covariates: it holds `condition`, `n_cells`, `total_counts`,
   `role` and the genes, and `load_pseudobulk` reads every column it does not
   recognise as a gene.
2. **The fast engine reuses dispersions fitted on the observed labels.**
   Its null is heavier-tailed than the exact null, so it calls fewer genes
   (genes at BH FDR 5%, exact against fast: REML 4,427 against 4,287;
   top-50% SumRank 305 against 193; Stouffer 6,789 against 6,462). Adding
   sex to the fast engine changed nothing measurable. Fitting the
   dispersions without the condition term moved the null the other way
   (REML 4,784; Stouffer 7,160). The exact permutation sits between the two
   because it re-estimates the dispersions with each shuffled label in the
   model.
3. **Study-level p-values differ between methods.** SumRank reads the raw
   per-study `pvals`. The matrix used by Fisher, Stouffer and gwOP
   (`pooling_core.py`) and the WOP and gwOP modules read `pvals_adj` when it
   exists. Adjusted p-values are not uniform under the null, so the
   analytical p-values of these methods are not valid, and the permutation
   null (which uses raw p-values) does not describe the observed statistic.
4. **The fast engine's name does not match its selector.**
   `consensus_cc_permutation` defaults to `de_method='deseq2_fast'`, but
   dispatches on `'deseq2'`, so that default falls through to the Welch
   t-test.

## Decision

**The exact permutation is added as an engine, run as an opt-in background
job that can be stopped and resumed. The pseudobulk file carries each run's
covariates. The shuffle respects the covariates. Methods that combine
per-study p-values use the raw p-values. The fast engine is kept and
labelled as an approximation.**

### 1. Covariates in the pseudobulk file

- `pseudobulk_frame` writes one column `cov:<name>` for every sample-level
  covariate of the DE run, after `role`. The prefix cannot occur in a gene
  symbol, so a loader that does not know about covariates still has to
  exclude them by name; `load_pseudobulk` does so and keeps its return
  signature.
- `load_pseudobulk_covariates(path)` returns the covariate table
  (index: donor id), or `None` for a file written before this change.
  `load_pseudobulk_set` adds it to each study's dict as `'covariates'`.
- The DE run's `.json` record already names the covariates and the design
  (`covariates`, `design`). The permutation reads the names from the
  pseudobulk file's `cov:` columns, so the two files agree.
- **Files without `cov:` columns.** A study whose DE record lists covariates
  but whose pseudobulk file has none cannot be permuted exactly. The exact
  engine stops with a message naming the study and asking for the DE run to
  be repeated. Nothing is guessed. The fast engine and the analytical
  p-values do not need the columns.

### 2. The exact engine

- `de_method='deseq2_exact'` in `consensus_cc_permutation`. On every
  permutation and every study the per-study model is fitted by the function
  the DE workspace uses. The model fit and Wald test in `_run_deseq2`
  (`kosmic/de/de_analysis.py`) are moved into one function,
  `fit_deseq2_wald`, which `_run_deseq2` and the permutation both call, so
  the two cannot differ. It applies the same expression filter, the same
  usable-covariate rule, the same design and PyDESeq2 with the same
  settings. It returns the unshrunk log2 fold change, its standard error and
  the Wald p-value. Shrinkage and the descriptive columns are not computed.
- Permutations run in parallel over processes, one core per fit. Each
  finished permutation is written to the selection's `meta_analysis/` folder
  (`cc_permutation/perm_<n>.npz`) so a stopped run resumes where it stopped.
  A progress line gives permutations done, time per permutation and the
  estimated time left.
- **Cost.** About 30–37 s per study on one core, 165–170 s per permutation
  for five studies, and about 40–75 minutes for 200 permutations on 14 cores
  (measured). The run reports the permutations done and the time left as it
  goes, and the Calibration box's tooltip states the cost.
- **Number of permutations.** The exact engine defaults to 200 (a null of
  about 2.7 million values, pooled over genes), the fast engine to 1,000.
  Both are recorded with the run, with the smallest calibrated p-value they
  can give, 1 / (null values + 1).
- **Stopping.** During an exact run the run button reads *Stop*. A stopped
  run keeps its finished permutations and is not reported as a failure;
  running again with the same settings continues from them.
- **Where it is offered.** On the Gene Meta-Analysis page, as the third
  Calibration choice, *CC permutation (exact, slow)*, for DESeq2 input.
  With Welch input the permutation keeps its Welch engine. Leave-one-study-out
  repeats the whole pipeline once per study, so it uses the fast engine
  whichever is chosen.
- **The calibrated p-value** is (1 + null values at least as extreme as the
  observed statistic) / (1 + null values), for every engine. Earlier versions
  counted only strictly more extreme values and floored the result at
  1 / (n + 1).

### 3. The shuffle

- Only donors with role `disease` or `control` are shuffled. Excluded donors
  are not in the analysis and are not shuffled into it.
- When the design has categorical covariates, roles are permuted within each
  combination of covariate levels, so the covariate-by-condition table of
  each study is the same on every permutation. A level combination with one
  donor stays fixed. If no donor can move, the engine stops with a message.
- The shuffle is independent per study. Its seed is `(seed, permutation,
  study)`, so a resumed or extended run uses the same shuffles.

### 4. Raw per-study p-values

- Fisher, Stouffer, gwOP, WOP and the permutation null use `pvals`.
  `pvals_adj` is used only when `pvals` is absent. Vote counting keeps
  `pvals_adj`, because it applies an FDR threshold to each study.
- Pooled results for these methods change relative to runs made before this
  ADR. The run's provenance records which p-value column each method used.

### 5. The fast engine

- It stays, named *Fast (approximate)* in the interface, the Methods page and
  the provenance, with the measured difference from the exact permutation
  stated on the Methods help page. The selector mismatch is corrected so that
  each name selects the engine it describes.
- A covariate-aware fast engine that re-estimates the dispersions on every
  permutation replaces it only if its null and its FDR 5% counts fall within
  the 95% bootstrap interval of the exact permutation on the DCM endothelial
  reference (stage 5).

### 6. Provenance

The meta-analysis `provenance.json` and `methods.md` record the engine, the
number of permutations, the seed, whether the shuffle was stratified, the
covariates, the smallest attainable calibrated p-value (1 / (n + 1) over the
pooled null) and the per-study DE package version.

## Compatibility

- Pseudobulk files written before this change load unchanged. They carry no
  covariates, so exact permutation of a study with covariates needs its DE
  run repeated.
- Saved meta-analysis runs are read as before. Runs made before section 4
  are marked as pooled with the earlier p-value rule in their provenance
  when it can be told from the saved settings, and are otherwise unmarked.
- The current DCM project (`New_QC_DecontX`) has no pseudobulk file with
  covariates; the studies are rerun after the pending sample removals.

## Consequences

- A calibrated p-value can be produced from the same model as the observed
  one, at the cost of hours rather than minutes.
- Pooled Fisher, Stouffer and gwOP p-values change. The default consensus
  includes Stouffer, so default consensus results change.
- The pseudobulk file format gains optional columns.
- About 8 modules change: `kosmic/de/batch.py`, `kosmic/de/de_analysis.py`,
  `kosmic/meta_analysis/cc_permutation.py`, `io.py`, `consensus_loo.py`,
  `pooling/pooling_core.py`, `pooling/wop.py`, `pooling/gwop.py`, the Gene
  Meta-Analysis page and its worker, and the help and Methods pages.

## Stages

Each stage ends with `pytest tests/` and `ruff check kosmic/ main.py
tests/` passing.

1. **Covariates in the pseudobulk file.** Writer, loader,
   `load_pseudobulk_covariates`, `load_pseudobulk_set`; CODEBASE.md data
   contract. Tests: a file written with covariates loads with the same genes
   as one written without; covariates round-trip; a file without them gives
   `None`.
2. **One per-study model function.** `fit_deseq2_wald` extracted from
   `_run_deseq2`; the DE workspace's results are unchanged. Tests: the
   existing DE tests pass; `_run_deseq2` and `fit_deseq2_wald` give the same
   fold changes and p-values on a fixture.
3. **Raw p-values.** `pooling_core.py`, `wop.py`, `gwop.py`; provenance
   records the column. Tests: Fisher and Stouffer on a fixture with both
   columns equal the values computed from raw p; the permutation null and
   the observed statistic use the same column.
4. **The exact engine and the shuffle.** `deseq2_exact`, the stratified
   shuffle, resume, progress, the engine selector and its estimate on the
   Gene Meta-Analysis page, the corrected selector for the fast engine.
   Tests: the shuffle keeps the covariate-by-condition table; a resumed run
   equals an uninterrupted one; the exact engine on a small fixture equals
   a direct `fit_deseq2_wald` call on the same shuffled labels; a missing
   covariate column stops with the message.
5. **Validation, fast engine, documentation.** The 200-permutation exact
   reference on the DCM endothelial cells recorded under `dev/`; the
   covariate-aware fast engine tried against it (section 5); the Methods and
   help pages, CODEBASE.md, and the paper Methods text.
