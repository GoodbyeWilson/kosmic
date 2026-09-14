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

Two other families are implemented alongside it, because they make
different assumptions and fail in different ways:

- **Rank aggregation** — SumRank, and gwOP — which pool the *ranking*
  of genes within each study and are insensitive to effect-size scale.
- **p-value combination** — Fisher's method (unsigned) and Stouffer's
  (signed) — which pool evidence without an effect size.

The Discovery Meta-Analysis step can run any combination and show the
overlap. Genes found by every family are a **high-confidence consensus
set**: a gene that survives effect-size pooling, rank aggregation and
p-value combination is not depending on the assumptions of any one of
them. It is a consensus of statistical methods on the same data, not an
independent validation; the Validation step's proteomics (Olink) and
GWAS comparisons are the external evidence.

With two studies, the DerSimonian–Laird and REML estimates coincide
(τ² has one degree of freedom), Hartung–Knapp has one degree of freedom
in its t-distribution and finds nothing, and the rank methods have
almost no resolution. Two studies is the minimum that runs; five or
more is where the method choice starts to matter.

## Permutation calibration

Analytical p-values from pooled single-cell DE are anticonservative
when the per-study tests share structure. KOSMIC calibrates them by
case–control permutation: shuffling the disease/control labels within
each study, re-running DE and the pooling, and building an empirical
null for each gene. Re-fitting DESeq2 per permutation is not tractable,
so the permutation loop reuses each gene's dispersion estimated once on
the observed data and refits only the coefficients — an accelerator
validated head-to-head against the full model
(`dev/simulations/validate_deseq2_fast.py`). The assumption this makes,
that dispersions are fixed across label shuffles, is recorded here so
that it can be stated in a methods section.

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
GSEA respectively; Gene Ontology enrichment uses the topGO-style *elim*
algorithm, which removes the genes of a significant child term before
testing its parents.

## Gene identity across studies

Everything above pools *by gene name*. Studies annotated against
different releases name the same gene differently, so KOSMIC harmonises
symbols to the current HGNC nomenclature before pooling (see
[Gene Names](../scrna/gene_names.md)), and counts how many genes are
shared across the studies before an atlas is built. A gene absent from
one study is absent from the pooled result; harmonisation is what keeps
that number small.
