# Statistical approach

This page states the positions KOSMIC takes on how single-cell
differential expression and cross-study meta-analysis should be done,
and why. The settings a particular run used are recorded in that
project's provenance and rendered by the Methods step; this page is the
reasoning behind the defaults.

## Pseudobulk differential expression, per study

Within each study, KOSMIC aggregates counts to one value per donor per
cell type and runs differential expression on those pseudobulk profiles
with PyDESeq2, a Python implementation of DESeq2. The biological donor,
not the individual cell, is the unit of replication. This is the
position established by Squair et al. (*Nature Communications* 2021),
who showed that methods treating each cell as an independent replicate
(Wilcoxon, MAST, edgeR on cells) inflate the effective sample size and
produce large numbers of false positives. Welch's t-test on
CPM-normalised or raw pseudobulk profiles is available as an
alternative.

The DESeq2 design is `~ covariates + condition`, where the covariates
are sample-level variables you choose to adjust for (age, sex, batch);
covariates that cannot be estimated — constant, unique per sample, or
collinear with condition — are dropped and the fitted design is
recorded with the result.

## Uncorrected counts, not corrected matrices

Differential expression always runs on the original counts. Batch
correction in KOSMIC (Harmony) acts on the PCA embedding used for
clustering and never produces a corrected expression matrix, so DE
cannot be run on corrected values. This follows Nguyen et al. (*Nature
Communications* 2023), who found that batch-effect-correction methods
which return corrected expression matrices distort log fold-changes
when DE is run on their output, and that modelling batch as a
covariate on uncorrected counts performs better — more clearly so as
batch effects grow, which is the cross-study regime KOSMIC operates in.

When several studies have been combined into a shared atlas and DE is
run on that combined dataset — the "mega-analysis" arm, as opposed to
the meta-analysis — the study is included as a covariate in the same
way.

## Random-effects meta-analysis across studies

Independently published studies differ in protocol, platform,
dissociation, and cohort, so between-study heterogeneity (τ²) is
structurally large. A random-effects model is the appropriate one for
that regime and is KOSMIC's primary pooling method: DerSimonian–Laird
(DerSimonian & Laird, *Controlled Clinical Trials* 1986), with REML as
an alternative estimator and the Hartung–Knapp–Sidik–Jonkman
adjustment (IntHout et al., *BMC Medical Research Methodology* 2014)
available for small numbers of studies.

I² is the share of the total variation that comes from between-study
heterogeneity (Higgins & Thompson, *Statistics in Medicine* 2002):
τ² / (τ² + s²), where s² is the typical within-study variance,
(k − 1) Σw / ((Σw)² − Σw²) with w = 1 / SE². This is how metafor
computes it for every τ² estimator. For DerSimonian–Laird it equals
(Q − df) / Q; for REML it uses the REML τ².

Two other families are implemented alongside it, because they make
different assumptions and fail in different ways:

- **Rank aggregation** — SumRank, and gwOP — which pool the *ranking*
  of genes within each study and are insensitive to effect-size scale.
- **p-value combination** — Fisher's method (unsigned) and Stouffer's
  (signed) — which pool evidence without an effect size. Both are the
  textbook methods: each study's p-value is used as it is, with no floor,
  and a study with p = 1 counts (it adds z = 0 to Stouffer's sum). The
  only limit is numerical: a p-value that underflowed to 0 is taken as
  10⁻³⁰⁰. With p-values as small as single-cell studies give, the pooled
  p-value of the strongest genes can therefore be very small; their order
  is still informative.
- **Vote counting** — the number of studies that call a gene
  significant in the same direction. See
  [Vote counting](#vote-counting).

The Discovery Meta-Analysis step can run any combination and show the
overlap. Genes found by every family are a **high-confidence consensus
set**: a gene that survives effect-size pooling, rank aggregation and
p-value combination is not depending on the assumptions of any one of
them. It is a consensus of statistical methods on the same data, not an
independent validation; the Validation step's proteomics (Olink) and
GWAS comparisons are the external evidence.

With two studies, τ² rests on a single difference between them and is
very uncertain whichever estimator is used, Hartung–Knapp has one degree of freedom
in its t-distribution and finds nothing, and the rank methods have
almost no resolution. Two studies is the minimum that runs; five or
more is where the method choice starts to matter.

### Vote counting

Vote counting asks in how many studies a gene is significant (study
FDR < 0.05) with the same sign of log2 fold change. A study that did
not measure the gene, or whose results have no adjusted p-values, is
left out. A gene **passes** when at least *N* studies call it
significant in the same direction and no study calls it significant in
the opposite direction; *N* is the **Vote threshold** (*N* of *k*) in
the Discovery Meta-Analysis settings, all studies by default. It is a
rule, not a test: there is no p-value, and CC permutation calibrates
only the other methods of a run. In a consensus run, a gene that fails
the rule fails the consensus.

Every gene also has a **level**: 0/*k* to *k*/*k* for the studies
agreeing, or *conflict* when it is significant up in one study and down
in another. The Results tab can show any set of levels, and the Votes
tab and the Study Agreement figure count the genes at each level, with
the genes another method of the run calls significant overlaid.

Vote counting ignores the size of the effects and their standard
errors, and gives no pooled estimate or confidence interval: the log2
fold change shown for it in the volcano and results table is the
inverse-variance mean, as for Fisher's method. It is a robustness check
beside DerSimonian–Laird or REML, not a replacement: a gene with a
consistent but modest effect that no single study calls significant
gets no votes, however many studies there are. In the DCM ventricular
cardiomyocyte meta-analysis (five studies, 9,263 genes), 1,131 of the
3,324 genes REML called significant were called by no more than one
study, and 571 by four or five.

### Direction of effect across studies

Fisher's method, SumRank and gwOP combine per-study p-values or ranks
without using the sign of the fold change. A gene that is significantly
up-regulated in one study and significantly down-regulated in another
can therefore be called significant by these methods, although the
studies disagree about what happens to it. Stouffer's method uses the
sign, and DerSimonian–Laird and REML pool the signed effect, so such a
gene is less likely to reach significance with them, but it is not
excluded.

Every pooled table (gene and pathway) therefore carries three columns:

- `n_up` — the number of studies in which the gene has a study-level
  FDR below 0.05 and a positive log2 fold change;
- `n_down` — the same for a negative log2 fold change;
- `direction_conflict` — true when both counts are above zero.

The counts use each study's adjusted p-value (`pvals_adj`) as written
by the study's DE step. The DESeq2 + VIF pathway summaries carry only a
pooled log2 fold change and a VIF-adjusted standard error per study, so
for them each study's p-value is a two-sided Wald test of log2FC / SE,
adjusted by Benjamini–Hochberg across that study's pathways. A study
whose results table has no adjusted p-value is not counted, and the
output panel names it when the run starts. When no study has adjusted
p-values the columns are omitted rather than reported as zero. A significant gene with `direction_conflict` set should be read
as a disagreement between studies, not as a consistent change.

## Permutation calibration

The analytical p-value of a pooling method rests on assumptions about
the per-study tests that pseudobulk DE does not always meet. KOSMIC can
calibrate the p-values by case–control permutation. On each permutation
the disease and control labels are shuffled among the donors of each
study, the study-level DE is repeated on the shuffled labels, and the
studies are pooled by each method. The pooled statistics of all genes and
all permutations form one null distribution. Each observed pooled
statistic then receives the permutation p-value

(1 + *n*<sub>extreme</sub>) / (1 + *n*<sub>null</sub>),

where *n*<sub>extreme</sub> is the number of null values at least as
extreme as the observed one and *n*<sub>null</sub> is the number of null
values. The observed statistic counts as one draw of its own null, so a
calibrated p-value is never zero.

Three engines can produce the study-level DE on the shuffled labels.

**Exact.** Every permutation re-runs the study's own DE: the same
function the DE workspace uses, the same design including the covariates
the run adjusted for (for example `sex_inferred`), the same expression
filter, and dispersions estimated again from the shuffled labels. Donors
with role *exclude* are not in the analysis and do not take part. When
the design has covariates, the roles are shuffled within each combination
of covariate levels, so the covariate-by-condition table of every study
is the same on every permutation. The covariates come from the
`cov:<name>` columns of the pseudobulk file; a study whose DE run
adjusted for covariates but whose pseudobulk file has none must have its
DE run repeated. A permutation takes about 30–40 s per study on one core
for a small cell type and several minutes for five studies of a large
one, so the default is 200 permutations. Each finished permutation is
saved in the `cc_permutation/` folder beside the results, and a stopped
run continues from them.

**Fast.** The same model, the same shuffles and the same saved
permutations as the exact engine, with the fit of each study written as
array code over genes (`kosmic/meta_analysis/vec_deseq2.py`). It follows
PyDESeq2 step by step, including the dispersion fits, the Cook's outlier
refit and the Wald test, and returns the same log2 fold changes, standard
errors and p-values. It takes a fifth to a tenth of the time (4.9 s against 31.9 s per permutation of
five studies on 14 cores for DCM cardiomyocytes; 30 s against 307 s per
permutation on one core for endothelial cells), so the default is
1,000 permutations, saved in `cc_permutation_vec/`. It uses scipy's own
L-BFGS-B through a private scipy function; KOSMIC checks at run time that
the function behaves as expected and, if it does not, uses the
approximation below and says so. On the DCM endothelial cells (five
studies, 1,000 permutations, the same shuffles for every engine) the exact
and fast engines gave identical numbers of genes at a Benjamini–Hochberg
FDR below 0.05 for every method, and null tail fractions equal to within
0.00001.

**Approximate.** Used only where the fast engine cannot run. The
dispersions are estimated once on the observed labels, and each
permutation refits a vectorised two-group model (condition only, no
covariates) with them. A permutation of five studies takes under a second.
Its null has heavier tails than the exact null, so it calls fewer genes. On
the same cells:

| Method | Exact | Fast | Approximate |
|---|---|---|---|
| DL | 4,350 | 4,350 | 4,251 |
| REML | 4,287 | 4,287 | 4,158 |
| SumRank, top 50% | 291 | 291 | 178 |
| gwOP | 5,951 | 5,951 | 5,492 |
| Fisher | 7,447 | 7,447 | 6,790 |
| Stouffer | 6,612 | 6,612 | 6,306 |

This is one cell type, so the size of the difference in another data set
is not known. Leave-one-study-out validation repeats the pipeline once
per study and always uses the approximate engine.

## Pathway-level testing

Hypothesis mode tests pathways you choose *a priori*. For each pathway,
KOSMIC computes a per-donor pathway score — the mean log2 normalised
expression of the pathway's genes in that donor's pseudobulk — and
tests disease against control across donors with Welch's t-test, giving
a per-study pathway effect that is pooled across studies exactly as a
gene would be. Collapsing to one value per donor absorbs the
correlation between genes in the set, so no inter-gene correlation
correction is needed, and the test is powered by donors rather than by
genes.

Two things follow, and they are the reason hypothesis mode can find
what genome-wide discovery misses:

- **Multiple testing is over the pathways you chose**, not over the
  transcriptome. Testing ten pathways carries a ten-fold correction;
  testing 20,000 genes carries a 20,000-fold one.
- **A pathway score accumulates small, consistent shifts.** Genes that
  individually fall short of significance can move a pathway's score
  together; this is what a coordinated biological programme looks like
  and what per-gene testing is least able to see.

Over-representation and gene-set enrichment on discovery results use
Fisher's exact test against the tested-gene background and preranked
GSEA respectively; Gene Ontology enrichment uses Fisher's exact test
with BH correction by default (KOSMIC's *elim* option is its own
variant, not topGO's, and is under review).

## Gene identity across studies

Everything above pools *by gene name*. Studies annotated against
different releases name the same gene differently, so KOSMIC harmonises
symbols to the current nomenclature (HGNC for human, MGI for mouse)
before pooling (see
[Gene Names](../scrna/gene_names.md)), and counts how many genes are
shared across the studies before an atlas is built. A gene absent from
one study is absent from the pooled result; harmonisation is what keeps
that number small.
