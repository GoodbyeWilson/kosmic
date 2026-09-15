# Inspect

Tell KOSMIC how to read this dataset. Nothing downstream — QC, clustering,
differential expression, meta-analysis — can run until it knows which column
identifies a **sample**, which identifies a **condition**, and what each
condition value *means*.

The sidebar shows dataset-level context that stays in view across all three
tabs: **Data Summary** (cells, genes, metadata columns, whether a raw-counts
snapshot exists, sparsity) and **Current Mapping** (what you have designated
so far, and which condition values are excluded).

## The three tabs

They are a sequence, not alternatives.

| Tab | What it is for |
|---|---|
| **Setup** | Configure. Designate the columns, then assign a role to each condition value. **Save & continue** writes both to the h5ad |
| **Samples** | Verify. One row per sample with QC metrics, so you can see the mapping produced the groups you expected |
| **Tools** | Occasional. Filter the dataset down to a subset of cells, or merge in external metadata from a CSV |

## Setup tab

### 1. Dataset mapping

Three dropdowns — Sample, Condition, Cell type — each picking an `obs` column.

**Sample** must identify the *donor or biological replicate*, not the
sequencing run. Pseudobulk aggregates cells within a sample and every
downstream test treats samples as the unit of replication, so getting this
wrong inflates your n and your significance with it.

**Browse all N columns…** opens the full column list: every `obs` column with
its type, how many distinct values it holds, and example values. This is one
row per *column*, unlike the Samples tab's one row per *sample*. Use it when
you are not sure which column holds what — the `Designate as` dropdown in its
last column does the same job as the three combos.

### 2. Condition roles

One row per distinct value of the condition column, each assigned a role:

- **Control** — the reference group
- **Disease** — the comparison group
- **Exclude** — kept in the dataset, left out of the contrast

You need at least one Control and one Disease before Save will do anything
useful. KOSMIC pre-fills a guess from the value names; check it, because the
guess is a pattern match and cannot know your study's conventions.

**Renaming a value**: double-click it. The new name replaces the old one
throughout the condition column when you press Save & continue. Renaming onto
a name another row already holds is refused — merging two conditions is a
different operation and is not what this does.

**Set condition per sample…** opens a dialog for labelling each sample
individually. You need it only when the dataset arrived with no condition
column; the labels you enter there *become* the values this card assigns roles
to. When a condition column already exists the dropdown offers all of its
values, and changing one rewrites that column for every cell of that sample.

## What Exclude actually does

The cells stay in the dataset. They are simply never part of a comparison.

```
Inspect        obs['_role'] = disease | control | exclude
   ↓
Pseudobulk     one row per sample, carrying its role
   ↓
DE             disease_mask = role == 'disease'
               control_mask = role == 'control'
               excluded samples are in neither
   ↓
Meta-analysis  reads per-study DE results, so exclusions are already applied
```

Use it for an arm that is not part of your question — a HCM cohort in a study
you are mining for DCM, or a treatment group you are not comparing.

Two things worth knowing:

- Exclude works on a **condition value**, not a sample. To drop a single bad
  donor you would have to give it its own label first.
- Excluded samples still influence the **gene filter** — the prevalence
  test that decides which genes are worth testing runs over the whole
  pseudobulk matrix. They never affect a fold change or a p-value, but they
  can affect which genes get tested at all.

## Samples tab

One row per sample: cell count, condition, median genes detected per cell,
mitochondrial percentage, doublet rate once Scrublet has run, and the
sample's sex.

**Sex** is read from the data, not the metadata: XIST is expressed in every
female cell and no male cell, and the Y-chromosome genes (DDX3Y, UTY,
KDM5D, RPS4Y1, EIF1AY, USP9Y) the reverse, so summed over a sample's cells
the two signals separate completely. Saving the setup writes the call to
`obs['sex_inferred']`, one value per cell, so it is available as a
covariate in differential expression and for the methods text. If the
file has a sex column, designate it as **Sex** in the column mapping (the
usual names are picked up automatically); the study then keeps that
column, and the **Sex check** column reports any sample whose call
disagrees with it — a mislabelled or swapped sample. With no recorded
sex, the inferred column becomes the study's sex column.

Two other things the check can show. **Mixed signal**: a sample carrying
both XIST and Y-gene counts well above zero, which means cells of both
sexes are in it — ambient contamination, or a multiplexed library that
was not demultiplexed. The call is still made (from the larger signal),
but the sample deserves a look. **No signal**: neither gene set is
detected, usually because the matrix does not contain them. Hover the
cell for the counts per million behind the flag.

This is the verification step. Look for a sample with far fewer cells than the
rest, an outlying mitochondrial fraction, or — most importantly — a sample in
the wrong group.

Once roles are set the condition column shows `control` / `disease` /
`exclude` rather than the original values, because that is what the analysis
will actually use. It is the fastest way to confirm an exclusion landed where
you meant it.

## Tools tab

**Filter dataset** keeps a subset of cells — for example only the author's
endothelial cells — while preserving existing embeddings and annotations. This
is not the same as the Subset workflow step, which derives a new analysis
dataset from clusters you have computed.

**Merge external metadata** joins a CSV of per-sample information (age, sex,
aetiology) onto `obs` by sample ID. Useful when the deposited object omits
clinical variables the paper describes.

## Why roles rather than labels

Studies name the same thing differently — `NF`, `Donor`, `Healthy`, `control`
all mean the reference group. Roles normalise that once, here, so DE and
meta-analysis compare like with like across studies without pattern-matching
condition strings at every step.

It also makes the third state explicit. A study with three arms has a group
that belongs to neither side of your contrast, and `exclude` says so, rather
than leaving it to be silently dropped or silently included.

## When you are done

Save & continue writes `obs['_role']`, `uns['role_map']` and the standardised
`sample` / `condition` / `cell_type` columns to the h5ad, records the step in
the study's provenance, and lands you on the Samples tab to check the result.

Nothing is written until you press it, so a mapping or rename you think better
of is abandoned by leaving the tab.
