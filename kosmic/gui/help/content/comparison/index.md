# Methods Comparison

The Methods Comparison workspace compares the genes called by merged
pseudobulk with the genes called by four meta-analysis methods, for one
counts source and one cell type of the project.

- **Merged pseudobulk** fits one DESeq2 model to the donors of every study
  in the shared atlas (`_master`). It adjusts for the covariates the
  per-study models used, plus study.
- **REML**, **SumRank**, **Stouffer** and **vote counting** analyse each
  study separately and combine the study-level results. They are the
  methods of the Meta-Analysis workspace, computed in the same way; see
  [Statistical approach](../meta/statistics.md).

The workspace needs an open project with a shared atlas, made with
*Combine* on the Project page, and a DESeq2 run per cell type in two or more
studies, made on the same counts.

The Meta-Analysis workspace's [Consensus Overview](../meta/methods_comparison.md)
mode is different: it compares the meta-analysis methods with each other,
using leave-one-out replication and held-out AUC, and does not include
merged pseudobulk.

## Choosing what to compare

| Control | What it does |
|---|---|
| **Selection** | A cell type and counts source with a DESeq2 run in two or more studies, for example *Ventricular Cardiomyocyte, decontX counts (5 studies)*. |
| Studies | The studies of the selection, all ticked. Untick a study to leave it out of the meta-analysis methods. Merged pseudobulk is fitted to every study in the atlas, so the page warns when the atlas holds a study that is not ticked. |
| **Methods** | The methods compared. At least two must be ticked. |
| **Calls** | The rule that decides which genes each method calls (below), with its threshold, and **FDR over the genes compared** (ticked by default; below). |
| **Run comparison** | Collects each method's result, compares them and saves the comparison. |

## The merged pseudobulk design

When a selection is chosen, the page reads each study's DE run record and
the atlas, and shows the design merged pseudobulk uses, for example
`~sex + study + condition`.

The studies and the atlas can name the same variable differently. A study's
model may adjust for `sex_inferred`, while the atlas holds `sex`, which
Combine writes from each study's designated sex column. The page maps each
per-study covariate to an atlas column and shows the mapping
(`sex_inferred → sex`). It then checks the values, not the names:

- a sex covariate is accepted when it agrees, for every donor, with the
  column the atlas's `sex` was written from;
- a study that identifies its donors by another column (`donor_id`) is
  accepted when that column names the same donors as the atlas's `sample`.

A warning is shown under the design in each of these cases. No covariate is
removed from the design because of a warning.

- The studies used different covariates.
- A study's sex covariate differs from the atlas's `sex` for some donors.
- A study's donors are not the atlas's donors.
- A covariate has no column in the atlas, varies within one study only, does
  not vary within any study (it is confounded with study), or is confounded
  with condition.
- The studies' runs differ in their counts, gene filters or minimum cells or
  counts per donor.
- The atlas holds donors of studies that the meta-analysis methods do not
  pool (an unticked study, or a study with no run for this cell type), or a
  pooled study has no donors of the cell type in the atlas. Merged
  pseudobulk is fitted to the atlas's studies, so the two sides then do not
  compare the same studies.

The page also says whether the atlas has a saved merged pseudobulk run for
the selection, and whether its design matches. A saved run matches when it
adjusts for the mapped covariates and study, on the same counts, with the
same minimum cells and counts per donor and, where the runs recorded them,
the same gene filters and multiple-testing scope. DE runs made before
KOSMIC recorded their gene filters cannot be checked on them; the page says
so.

## Running the comparison

**Run comparison** reads or computes each method's result:

| Method | Read from | Computed when |
|---|---|---|
| Merged pseudobulk | The atlas's saved DESeq2 run for the selection, when its design matches | There is no saved run, or you choose **Resolve** (below). It is fitted with the DE workspace's per-cell-type run and saved in the atlas's DE results, so the DE workspace lists it. |
| REML, SumRank, Stouffer, vote counting | A run of that method alone, saved by the Meta-Analysis workspace for the same studies with the same settings | There is no such run. The method is pooled again from the studies' saved DE tables with the settings of the selection's latest meta-analysis run, or the Meta-Analysis defaults: the minimum number of studies, the expression filter, and the gene-set settings **Current gene symbols only** and **Leave out genes merged inconsistently across studies** (see [Gene Meta-Analysis](../meta/gene_ma.md#gene-set)). |

When the atlas's saved merged run does not match the design, a dialog
opens before anything runs. It shows the two designs and their differences,
and offers three choices:

- **Resolve** fits merged pseudobulk again with the matched design and
  replaces the saved merged run for this counts source and cell type.
  Resolve is not offered when a per-study covariate has no column in the
  atlas.
- **Run anyway** uses the saved merged run as it is. The differences are
  recorded with the comparison and listed under the results.
- **Cancel** runs nothing.

## Calls and the summary

The comparison uses only the genes that every selected method tested. A
gene that DESeq2 filtering removed from one method's results is therefore
not counted as a gene that method did not call.

**FDR over the genes compared** (ticked by default) corrects every method
again, with Benjamini-Hochberg, from its raw p-values, over the genes
compared. Every method's FDR then has the same denominator, so a
difference between two methods' calls is not a difference in how many
genes each corrected over. A method has tested a gene when it gave a raw
p-value for it: DESeq2's independent filter withholds only the adjusted
p-value, and the correction is made again here, as the DE workspace does
when it corrects over a committed gene list. Merged pseudobulk fitted over
every gene is then corrected over the same genes as the meta-analysis
without fitting it again; the fit, and so its raw p-values, do not depend
on the genes corrected over. Unticked, each method's adjusted p-value is
used as the method gave it, and the genes compared are those with one in
every method. Vote counting has no p-value and is not affected.

Unticked, the comparison gives what you get by comparing the tables
yourself, when merged pseudobulk was run on the DE page with **Cell-type
batch: each cell type's meta-analysis genes only** (see
[Gene DE](../de/gene_de.md#one-gene-universe-for-both-arms)): merged
pseudobulk corrected over the meta-analysis genes it tested, each
meta-analysis method over its own genes. The two settings differ only when
one method tests genes another does not, for example a gene the merged
model's expression filter removes.

| Rule | A method calls a gene when |
|---|---|
| **Adjusted p-value** (default) | its adjusted p-value is below the threshold (0.05). Vote counting calls the genes its vote rule passes: significant in the same direction in the minimum number of studies, all of them by default. |
| **Top N by rank** | the gene is among the method's N strongest genes, ranked within the genes compared. Every method then calls the same number of genes. |
| **Adjusted p-value, same direction** | it is called under the first rule, and the gene changes in the same direction in every method compared. |

Changing the rule, its threshold or the methods after a run updates the
tables and the saved comparison from the results already collected; nothing
is fitted again.

**Genes called** gives, for each method, where its result came from, the
genes compared, and the genes called, split into up and down in disease.

**Agreement between methods** gives, for each pair of methods:

| Column | Meaning |
|---|---|
| Both call | Genes both methods call. |
| Jaccard | Genes both call divided by genes either calls. |
| Same direction | Of the genes both call, the percentage that change in the same direction in both. |
| Rank correlation | Spearman correlation of the two methods' ranks over the genes compared. |

**Genes by the methods that call them** groups the genes called by at
least one method by exactly which methods call them (all of them, merged
pseudobulk and REML only, and so on). For each group:

| Column | Meaning |
|---|---|
| Genes | Genes in the group. |
| Same direction | The percentage whose log2 fold change has the same sign in every study's own DE result. |
| Median \|log2FC\| (REML) | The median absolute REML pooled log2 fold change. |
| Median I² | The median REML I², the share of the variation between studies not explained by chance. |
| Prediction interval excludes 0 | The percentage whose REML 95% prediction interval, pooled ± t(k−2)·√(SE² + τ²) with k the studies pooled, lies wholly above or below 0. |
| No single study | The percentage that no study's own analysis calls (study FDR below the threshold). |

The REML columns use the REML fit with the comparison's meta-analysis
settings; REML is pooled for them when it is not among the methods
compared.

The notes under the tables state where each method's result came from.

## Saved output

Each comparison is saved in the project in
`methods_comparison/<counts>/<cell type>/`:

| File | Contents |
|---|---|
| `results.csv` | One row per gene and method: effect, direction, statistic, p-value, adjusted p-value, rank, whether the method tested and called the gene, and whether the gene was compared. |
| `pairwise.csv` | The agreement between each pair of methods. |
| `genes.csv` | Per gene: the number of studies, same direction in every study, the number of studies calling it, REML's effect, I² and prediction interval. |
| `method_groups.csv` | The table *Genes by the methods that call them*. |
| `params.json` | The methods, the call rule, the merged design and its warnings, where each result came from, and the studies. |
| `provenance.json`, `methods.md` | The comparison's methods record. A new run replaces it. |

Choosing a selection with a saved comparison shows it. The
[Venn diagram](../figures/methods_venn.md) in Figure Export draws the
comparison.
