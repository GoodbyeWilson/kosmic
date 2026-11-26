# ADR-017: A Methods Comparison workspace for merged pseudobulk and meta-analysis

- Status: Accepted (Kirk, 2026-11-26)
- Extends: ADR-010 (Figures reads a workspace through `figure_inputs()`),
  ADR-015 (DE outputs per counts and cell type)

## Problem

A project's studies can be analysed in two ways that answer the same
question. In a **meta-analysis**, each study is analysed separately and
the study-level results are combined. In a **merged pseudobulk**
analysis, one DESeq2 model is fitted to the donors of every study in the
shared atlas (`_master`), with study as a covariate. Comparing the genes
that each approach calls is how a user decides which result to report.
KOSMIC cannot make that comparison.

1. **Merged pseudobulk can be run, but nothing compares it.** Running the
   DE workspace on `_master` with `study` among the covariates writes
   `_master_<cell type>_DE_deseq2.csv` (ADR-015). The Meta-Analysis
   workspace skips `_master` on purpose (`discover_de_results`), so that
   a study is never pooled with the atlas built from it.
2. **The existing comparison covers meta-analysis methods only.** The
   Meta-Analysis workspace's Methods Comparison mode compares DL, REML,
   Stouffer, Fisher, SumRank, gwOP and vote counting, using leave-one-out
   replication and held-out AUC. A merged pseudobulk analysis is not a
   meta-analysis, so it does not belong among those methods.
3. **The designs are recorded under different names.** The two
   approaches are comparable only when the merged model adjusts for the
   same covariates as the per-study models, plus study. The studies and
   the atlas name the same variable differently. In the DCM project,
   each study's run records `~sex_inferred + condition`, and the atlas's
   run records `~sex + study + condition`. This is because Combine writes
   each study's designated sex column (`uns['sex_column']`) to the
   atlas's `sex` column. One study identifies its donors by `donor_id`,
   while the atlas uses `sample`. A string comparison of the recorded
   formulas would treat matching designs as different. It would also
   miss a study whose designated sex column is not the one its model
   used.

## Decision

**A new workspace, Methods Comparison, is placed between Meta-Analysis
and Figures. For one counts source and cell type, it collects the
results of merged pseudobulk and of four meta-analysis methods. It
converts them into one table, compares the genes each method calls, and
saves the comparison in the project, where Figures can read it.**

### 1. Packages and navigation

- `kosmic/comparison/` is the analysis package. It does not import Qt.
- `kosmic/gui/comparison/` is the workspace. It has one `SidebarPage`.
- The navigation rail order becomes Project, scRNA, Differential
  Expression, Meta, Methods Comparison, Figure Export. No other entry
  moves.
- Like Meta-Analysis, the workspace operates on the project and needs an
  open project.
- `tests/test_architecture.py` adds `comparison` to its lists of
  workspaces, page directories and GUI packages.

### 2. Methods compared

| ID | Method | Source |
|---|---|---|
| `pb_merged` | Merged pseudobulk | DESeq2 on `_master`: the per-study covariates, then `study`, then condition |
| `ma_reml` | REML random effects | `get_analytical_pool_fn('reml')` |
| `ma_sumrank` | SumRank | `get_analytical_pool_fn('sumrank')` |
| `ma_stouffer` | Stouffer | `get_analytical_pool_fn('stouffer')` |
| `ma_vote` | Vote counting | `get_analytical_pool_fn('vote_count')` |

The workspace calls these existing functions. It does not change how any
of them computes its result.

### 3. Where each result comes from

- **Meta-analysis methods.**
  - The workspace reads a saved single-method run for the selection
    (`consensus_<key>_analytical.csv`) when its recorded `study_names`
    are the selection's studies.
  - Otherwise it pools the method again from the saved per-study tables
    (`saved.study_tables`). It uses the expression filter and minimum
    number of studies recorded for the selection's latest run, or the
    Meta-Analysis defaults when there is none.
  - Results are adjusted with the same Benjamini–Hochberg procedure as
    in Meta-Analysis.
- **Merged pseudobulk.**
  - The workspace reads the saved `_master` run for the selection's
    counts and cell type when its recorded design matches (§4).
  - When there is no saved run, it fits the model with
    `run_de_by_cell_type` on `master.h5ad`, using the minimum cells and
    counts per donor, the gene filters and the scope of the FDR
    correction that the studies' runs recorded. Gene filters are applied
    within each study.
  - DE run records gain `filters` and `fdr_scope` for this purpose. A
    saved merged run matches only if these equal the studies' values.
    Runs saved before these fields existed do not have them. For those
    runs, the filters are not compared, a refit uses KOSMIC's default
    gene filters, and the page says so. When a saved run's design does not match, the user chooses
    what happens in a warning dialog (§4).
  - The run is written to `_master`'s usual statistics folder, so the DE
    workspace lists it as well.
- The page states, for each method, whether its result was read from a
  saved run or computed, and why.

### 4. Matching the merged design to the per-study designs

The per-study design is read from each pooled study's DE run record:
covariates, sample column, filter settings and counts. Each covariate is
mapped to an atlas column as follows:

1. the column of the same name, if the atlas has one;
2. otherwise `sex`, if the covariate is a sex column (the study's
   designated `uns['sex_column']`, or a recognised recorded sex column);
3. otherwise the covariate is reported as missing from the atlas.

The merged design is the mapped covariates, then `study`, then
condition. Study is a fixed effect, because pydeseq2 fits only fixed
effects.

A saved `_master` run matches when all of the following hold:
- its covariates, after the same mapping, equal the merged design's;
- its counts equal the studies';
- its filter settings equal the studies'.

**When the saved merged run does not match.** Pressing *Run comparison*
opens a warning dialog before anything runs. The dialog shows the
per-study design, the saved merged design and the differences between
them, and offers three choices:

- **Resolve.** Fit merged pseudobulk again on `_master` with the matched
  design (§3), then run the comparison. The new fit replaces the saved
  `_master` run for that counts source and cell type. The dialog says
  so. Resolve cannot be chosen when a per-study covariate has no column
  in the atlas, and the dialog names that covariate.
- **Run anyway.** Use the saved merged run as it is. The mismatch is
  recorded in the comparison's parameters and in its provenance stage,
  and the page shows it beside the results.
- **Cancel.** Run nothing.

When there is no saved merged run, the model is fitted with the matched
design and no dialog opens.

The page shows the resolved formula before the run. It shows a warning
in each of the cases below, and the warnings are also recorded in the
comparison's parameters. No covariate is removed without a warning.

- The studies used different covariates.
- A study's designated sex column is not the covariate its model used,
  so the atlas's `sex` column holds a different variable.
- The studies' sample columns differ from the atlas's.
- A covariate is missing from the atlas, varies within one study only,
  or is confounded with study or condition.
- The atlas holds donors of studies the meta-analysis does not pool, or a
  pooled study has no donors of the cell type in the atlas. Merged
  pseudobulk is fitted to the atlas's studies.

### 5. The comparison

- **One table.** Each method's output is converted to a long table with
  the columns `gene, method, effect, direction, stat, pvalue, padj,
  rank, significant`.
  - Values that a method does not define are NA.
  - SumRank's effect is the mean log2 fold change, and its statistic is
    the rank sum.
  - Vote counting has no p-value.
- **Genes compared.** The comparison uses the genes tested by every
  selected method.
- **A call.** Under the default rule, a method calls a gene when its
  adjusted p-value is below 0.05. Vote counting calls a gene when its
  existing `vote_pass` is true.
  - Two other rules can be chosen: the top N genes by rank, or the
    significance rule plus the same direction in every method compared.
  - Direction concordance among shared calls, and the Spearman
    correlation of ranks, are reported under every rule.
- **Pairwise measures.** For each pair of methods, the comparison
  reports the number of genes both call, the Jaccard index, direction
  concordance and rank correlation.

### 6. Saved output, provenance and hand-off to Figures

**Saved output.** The comparison is saved in
`project/methods_comparison/<counts>/<selection>/`, a path built by
`comparison_output_dir` in `paths.py`. The folder holds:
- `results.csv`, the long table;
- `pairwise.csv`;
- `params.json`, which records the methods, the rule and thresholds, the
  merged design formula, where each result came from, and the warnings.

**Provenance.** Each run records a `methods_comparison` stage in that
folder's `provenance.json`. The stage replaces the previous entry, as
meta-analysis stages do, because a new run overwrites the saved
comparison.

**Hand-off to Figures.**
- Figures reads the comparison through `ComparisonWorkspace.figure_inputs()`
  when the workspace holds one. Otherwise it reads the newest saved
  comparison, following the pattern of ADR-010.
- A new Figures section, Methods Comparison, holds a Venn diagram of
  the calls of two or three chosen methods. Exporting it also writes a
  CSV of each region's genes beside the figure.
- An UpSet plot of all five methods is planned. The saved table already
  holds the data it needs.

### 7. The Meta-Analysis Methods Comparison mode remains

Amendment (2026-11-26): that mode is renamed Consensus Overview, so the
two places no longer share a name. Its code and help file names are
unchanged.

That mode compares meta-analysis methods with each other, using
leave-one-out replication, held-out AUC and combinations of methods. The
new workspace answers a different question. The two overlap for the
four shared methods. Whether to retire the older mode is left to a later
decision.

## Consequences

- Users can compare merged pseudobulk with meta-analysis, and export the
  figure, within KOSMIC rather than in separate scripts.
- A project gains a `methods_comparison/` folder beside `meta_analysis/`.
- Fitting the merged model reads `master.h5ad`, which can take minutes
  for a large atlas. This happens only when no matching saved run exists.
- The covariate mapping in §4 handles only sex columns. A future
  covariate that Combine renames would need its own rule; until then, it
  is reported as missing.
- Plasmode simulation is not part of this workspace.
