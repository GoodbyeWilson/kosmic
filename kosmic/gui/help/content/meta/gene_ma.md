# Gene Meta-Analysis (Consensus)

Pool each study's DESeq2 gene results across studies with one or more
methods. With one method ticked, the run reports that method's significant
genes. With several, a gene is a consensus gene when every ticked method
calls it.

## Methods

The **Methods** box has one tick box per method:

- **REML**: random-effects pooling of each gene's log2 fold change, with
  the between-study variance estimated by restricted maximum likelihood.
  It gives the pooled effect, its confidence interval and I².
- **SumRank**: pools each gene's rank within every study, so it does not
  depend on the scale of the effects. **best N% of studies** sets how many
  of each gene's study ranks are summed (see below).
- **Stouffer**: combines each study's signed z-score.
- **Vote counting**: counts the studies that call the gene significant
  (see below).

REML, SumRank and Stouffer are ticked by default; vote counting is not.
Stouffer and vote counting give no pooled effect of their own, so their
volcano and results table show REML's pooled log2 fold change. See
[Statistical approach](statistics.md#random-effects-meta-analysis-across-studies).

### SumRank: best N% of studies

At 100% (the default), SumRank sums all of a gene's study ranks. Below
100%, each gene's score is the sum of its best r study ranks, r = ceil(k ×
N/100) of its k studies (at 50%, the best 3 of 5), tested against the null
distribution of that partial sum. It asks whether a gene is highly ranked
in at least that share of the studies, which suits an effect present in a
subset of cohorts. When the cohorts agree it is stricter than the full sum,
because the other studies' evidence is not used: at 50% on five DCM
cardiomyocyte studies it called 378 genes against 1,348 for the full sum.
It does not remove genes before pooling. The Methods page states the
setting in words.

### Vote counting

**Vote counting** passes a gene when at least *N* studies
call it significant (study FDR < 0.05) in the same direction and none
calls it significant in the opposite direction. **Vote threshold** (*N*
of *k*, all studies by default) appears when it is ticked.
Run with REML, the consensus genes are those with a pooled effect that
at least *N* studies also call on their own. It ignores effect sizes and
gives no pooled estimate, so use it beside REML, not in place of it. Run on its
own, the page counts the genes that pass the rule rather than
significant genes. It has no p-value to draw a volcano from, so the
**Plots** tab shows the genes at each level of study agreement, split
up and down: the figure that Figure Export's
[Study Agreement](../figures/meta_vote_levels.md) page exports.

**Vote threshold** decides which genes the run calls significant: the
significant-gene count, the Venn diagram and the gene set the later
steps use. Every gene's level is kept whatever the threshold: the
results file has a `vote_level` column (for example `3/5` or
`conflict`), with `n_up`, `n_down` and `votes`.

In a run that includes vote counting, the **Results** tab shows, in
place of an FDR for vote counting, a column named after the threshold
(for example **≥ 5 of 5 agree**) that is ticked for the genes that
pass; the gene's
**Level** beside **k** (the studies in which the gene was tested) and a
**Studies agreeing** row with one box per level, 0/*k* to *k*/*k* and
*conflict*, each with its number of genes; untick levels to hide their
genes from the table. Figure Export's
[Study Agreement](../figures/meta_vote_levels.md) page draws the genes
at each level, split up and down, and can hatch the genes REML or
another method of the run calls significant. See
[Statistical approach](statistics.md#vote-counting).

## Calibration (analytical / CC permutation)

The Calibration box sets how the p-value of each pooling method is
obtained.

- **Analytical** uses each method's own p-value. Nothing is shuffled, and
  the result is immediate.
- **CC permutation (fast)** shuffles the disease and control labels within
  each study and fits each study's own DE model, with its covariates, on
  every shuffle by array code. It gives the exact permutation's results in
  a fifth to a tenth of the time and runs 1,000 permutations. Five studies of a
  large cell type take about an hour and a half on 14 cores. It relies on a private
  scipy function; where that function is not usable, KOSMIC says so and does
  not run, and the exact engine is the one to choose.
- **CC permutation (exact, slow)** re-runs each study's own DE (PyDESeq2) on
  every shuffle. It runs 200 permutations. A permutation takes about
  30–40 s per study on one core for a small cell type, and several
  minutes for a large one.

During a fast or exact run the run button reads *Stop*. A stopped run keeps
its finished permutations in the `cc_permutation_vec/` or `cc_permutation/`
folder beside the results; running again with the same settings continues
from them, and settings or data that differ start the permutations again.
Both need each study's pseudobulk file to carry the
covariates of its DE run; a study that does not says so and asks for its
DE run to be repeated. Leave-one-study-out validation always uses the
approximate engine.

Two settings under the box apply to a permutation only. Both are off by
default.

- **Null from the pooled gene set only.** Each shuffle still fits every
  study on all its genes, so size factors and the dispersion trend are as
  in the real fit. Only the genes this run pools, after the gene-set
  settings below, enter the null distribution and are calibrated. Off,
  every gene pooled from enough studies enters the null.
- **Don't re-apply the CPM filter on shuffles.** Fast and exact engines
  only. A study's saved pseudobulk holds only the genes its DE tested, so
  its library sizes are smaller than in the real fit, and re-applying the
  CPM rule on every shuffle drops some genes the real fit tested. Ticked,
  the CPM rule is off on shuffles; the summed-count rule still applies,
  since a shuffle cannot change it.

The run records the engine, the number of permutations, the seed, whether
the shuffle was within covariate groups, and the smallest calibrated
p-value the null can give. How the permutation works and how the
engines differ is on [Statistical approach](statistics.md#permutation-calibration).

## Gene set

These settings choose the genes the run pools. Each is off by default, so
a run that sets none pools every gene the studies report. The log gives
the number of genes each setting removed, and the run's settings file
(`.json`) records them.

- **Independent filter** drops genes whose detection rate across studies
  is below the threshold.
- **Current gene symbols only** keeps genes named by a current approved
  symbol (HGNC for human, MGI for mouse). Withdrawn or renamed symbols,
  Ensembl IDs and clone names are left out.
- **Leave out genes merged inconsistently across studies.** Gene-name
  harmonisation sums two features into one symbol where HGNC has merged
  the loci, but only in studies whose annotation carried both. Such a gene
  is a sum of two features in one study and a single feature in another.
  This is the same check the Combine step offers when building an atlas,
  run on the studies of this meta-analysis. Studies with no harmonisation
  record, and tables imported from outside the project, cannot be checked;
  the log names them.

The genes every pooling method tested are saved beside the results as
`consensus_<methods>_<calibration>_genes.csv`. A per-cell-type DE run on
the shared atlas can correct its FDR over that list (see
[Gene DE](../de/gene_de.md#one-gene-universe-for-both-arms)), so the
merged analysis tests the same genes as the meta-analysis.

## Tabs: Plots / Pathway Forest / Results / Venn / Enrichment / Cross-Dataset Reproducibility / LOO Validation

*TODO.*

### Direction conflicts in the Results tab

The Results table has three columns that count, for each gene, the
studies in which it is significant (study FDR < 0.05) with a positive
log2 fold change (**Up**) and with a negative one (**Down**).
**Conflict** is marked when both counts are above zero: the gene is
significantly up in at least one study and significantly down in
another. The line above the table reports how many of the significant
genes have such a conflict, and the same count is written to the output
panel and the Methods page.

Tick **Show only genes with opposite-direction effects** to list only
these genes. They occur mainly with SumRank over the best N% of
studies, which does not use the studies it leaves out; see
[Statistical approach](statistics.md#direction-of-effect-across-studies).

Studies whose DE results have no adjusted p-values are not counted; the
output panel names them when the run starts.

## Olink panel overlay

*TODO.*
