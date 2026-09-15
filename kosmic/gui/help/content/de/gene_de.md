# Gene DE

Pseudobulk differential expression on the active dataset. Aggregates
cells to one value per (sample, gene) before testing, so the unit of
analysis is the donor / sample rather than the cell.

> **Tip:** Set disease / control roles on the Inspect tab before
> running. KOSMIC reads ``adata.obs['_role']`` here, not the raw
> condition label.

## Settings (sidebar)

| Setting | Effect |
|---|---|
| **DE method** | Welch's t-test or DESeq2 (default). |
| **EB moderation** | Limma-style variance shrinkage (t-test only). |
| **Min cells per donor** | A donor contributing fewer cells is dropped. |
| **Min transcripts per donor** | A donor whose cells sum to fewer transcripts is dropped. Off (0) by default. |
| **Filter genes within each study** | Atlas only -- see below. |

## The two filters

Pseudobulk gives a **donors x genes** table, and there is one filter per
axis. Both are at the values a standard DESeq2 pseudobulk analysis uses,
and only the first is worth changing.

**Min cells per donor (drops donors).** A donor who contributed three
cells produces a profile summed from three cells, and it then enters the
model weighted like a donor with two thousand. DESeq2 sees counts, not
how many cells made them, so a thin donor inflates the gene-wise
dispersion -- which is shared across samples -- and costs power on
*every* gene. 10 is the convention.

The count is over whatever is loaded. On a whole dataset that is all of
a donor's cells; in a per-cell-type run it is only their cells of that
type. That is why this setting decides which cell types are testable: a
donor with 5,000 cells but 30 adipocytes fails it for the adipocyte run.

**Min transcripts per donor (drops donors).** The cell-count filter
does not guard depth. Ten shallow nuclei clear it and can still sum to a
few thousand transcripts, which is too thin a profile for a count model
for the same reason a three-cell profile is. This setting drops a donor
whose cells sum, over every gene, to fewer transcripts than the value.
It is off by default; Gao et al. use 50,000 for heart snRNA-seq. The two
filters count the same cells, so in a per-cell-type run both apply to
the donor's cells of that type.

Before setting a floor, look at where it would fall. Each donor's summed
depth is written to the pseudobulk file (`{accession}_pseudobulk.csv`,
column `total_counts`) beside the cell count, and donors dropped by
either filter are named in the output panel when the analysis runs. In
shallow datasets the floor mostly removes donors from rare cell types
(mast cells, neurons, lymphatics); it rarely touches the abundant ones.
In a small study losing one donor can leave an arm with a single
replicate, which stops that study being testable at all.

**The gene filter (drops genes).** A gene is kept when it reaches roughly
10 counts in a donor of median depth, in at least the smaller arm's
number of donors. The threshold is converted to CPM once, using the
median donor library, so a deeply sequenced donor needs more raw counts
than a shallow one to clear the same bar.

The donor count is set to the *smaller* arm deliberately: a gene
expressed only in disease should still be testable, and a higher bar
would filter out exactly the on/off genes worth finding. The count is
then taken across all donors, so any six of fourteen will do -- it does
not have to be six within an arm.

## Filtering within each study

This appears only when the dataset holds more than one study, and it is a
methodological choice rather than a tuning knob.

**Off** is the conventional setting: the gene filter is applied once to
the pooled table. That is what a standalone mega-analysis normally does,
and it keeps genes that only one cohort measures well.

**On** (the default here) runs the filter inside each study and requires
a gene to pass in all of them.

Turn it on when the pooled result will be compared against per-study DE.
A pooled cohort hides a problem: with 12 disease and 12 control donors
the filter asks for 11.4 samples, so a gene present in exactly one
study's 12 donors clears it and gets tested as though both cohorts could
see it. On the DCM endothelial data that was 592 genes -- genes the
pooled arm tested and the per-study arm was never offered, which would
have shown up as findings unique to pooling.

Turn it off for a standalone atlas analysis, where losing
cohort-specific biology costs more than the comparison is worth.

Either setting is recorded in the run's provenance as
`gene_filter_per_study`, so the methods text says which was used.

The **Run DE Analysis** button starts the worker. Progress shows in the
output panel; results render across the tabs on the right when it
finishes.

## Reading the results

### Volcano

Each dot is a gene. X = pooled log2 fold change (disease vs control),
Y = -log10 FDR. Dots above the horizontal line and beyond the vertical
lines pass the current significance thresholds. Click a dot to select
the gene in the table and other tabs.

### DE Results table

Sortable, filterable. The sidebar's display filters apply to this tab
*and* the volcano simultaneously.

> **Note:** Pathway / pct-cells filters trigger a BH FDR recompute
> over the visible scope (Bourgon Independent Filtering). The
> displayed FDR therefore matches the visible gene universe, not the
> genome-wide one. The find-gene field is a display lookup and does
> **not** redefine the scope.

### Gene Heatmap

Per-sample expression heatmap for the genes in a chosen pathway.

### Gene Expression

Per-condition mean expression with per-sample dots, for the gene(s)
typed in the input field. Comma-separate to plot a panel.

### Pathway Genes

Same per-gene bar render as **Gene Expression**, but driven by a
pathway dropdown rather than a free-text gene list.

### Top DE Genes

Stacked summary of the strongest hits, optionally grouped by pathway.

## Methods

### Welch's t-test (default)

Welch's t-test on mean-per-cell pseudobulk values. Fold changes
directly reflect per-cell expression differences, so they match what
you see when plotting means. Fast, and makes no assumption about the
proportion of genes that are unchanged.

### DESeq2

Negative-binomial GLM on sum-count pseudobulk with median-of-ratios
size-factor normalisation. Powerful for whole-transcriptome analyses
where most genes are unchanged. KOSMIC always runs DESeq2 on the full
transcriptome to keep size factors unbiased, then filters to your
gene set if applicable.

> **Warning:** In single-cell pseudobulk, library-size differences can
> be biological (disease cells expressing less per cell), not just
> technical. DESeq2's median-of-ratios normalisation may interpret
> that as a technical effect and normalise it away -- so a gene with
> near-identical per-cell expression (0.13 vs 0.12 counts/cell) can be
> reported as up-regulated because it "resisted" the global trend.
> The t-test avoids this because it normalises at the aggregation
> step. Default to t-test when testing curated gene sets.

### Empirical Bayes moderation (t-test only)

Optional. Shrinks per-gene variance toward a common prior (Smyth 2004 /
limma's `eBayes`). Stabilises estimates with few replicates and reduces
false positives from genes with unusually low variance.

## Multiple testing

BH FDR at 5%. The displayed scope depends on which display filters are
active -- see the **DE Results table** note above.

## Exporting

Use the export buttons in the sidebar to save the full table, the
significant-only subset, and the analysis settings JSON. The CSVs are
the inputs the meta-analysis workspace looks for.
