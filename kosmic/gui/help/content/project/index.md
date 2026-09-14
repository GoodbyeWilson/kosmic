# The Project workspace

The Project workspace is where datasets are organised and prepared for
analysis. A KOSMIC **project** is one cross-study analysis — for
example, a meta-analysis of human heart single-nucleus datasets
comparing dilated cardiomyopathy with non-failing donors.

Within a project, each published dataset is a separate **study**.
Studies are processed independently through the scRNA and Differential
Expression workspaces, and only then are their results pooled in
Meta-Analysis. Keeping them separate is deliberate: a meta-analysis
needs per-study results, and combining raw data first would throw that
away.

The Project workspace is the control centre for that. From here you
add studies, get data into them, see how far each has got, and check
they are ready to pool.

## The five steps

The steps are listed in the left panel; the one you are on decides what
the sidebar offers. Ticks are computed from what is on disk, so they
stay right if you change files outside KOSMIC.

| Step | What you do | Ticked when |
|---|---|---|
| [1. Project](project.md) | Open or create the folder for this analysis | a project is open |
| [2. Studies](studies.md) | Add each dataset as a study, named by its accession | at least two studies |
| [3. Data](data.md) | Import or download data for each study | every study has a processed dataset |
| [4. Shared atlas](atlas.md) *(optional)* | Combine studies so cell types are annotated consistently | the atlas exists |
| [5. Review](review.md) | Check every study is ready for analysis | at least two studies have DE results and pseudobulk counts |

When a project opens, KOSMIC selects the first step that is not done.
The steps are a map, not a gate — you can click any of them at any time.

![The Project workspace on the Review step: workflow steps and files on the left, the step's sidebar, the study table and the selected study's details](img/step5_review.png)

## What is in a project folder

```
my_project/
  GSE183852/            a study
    raw_data/           files as imported — h5ad, Seurat .rds, matrices
    processed_data/     the working h5ad KOSMIC produces and analyses
    results/            DE tables, pathway scoring, figures
  GSE298023/            another study
  _master/              the shared atlas, if you built one
  GSE183852_ECs/        a subset — endothelial cells cut from GSE183852
  meta_analysis/        cross-study output
```

Three kinds of folder appear in the study table, and they look alike:

- A **study** is one published dataset. It is what you add, import into,
  and analyse.
- The **shared atlas** (`_master`) is every study combined, built by
  step 4. It exists to be clustered and annotated once; it is never
  itself pooled.
- A **subset** is one cell type carved out of a study (or of the atlas)
  by the scRNA *Subset* step. It becomes a study in its own right — with
  its own DE results — and it is usually the subsets, not the full
  studies, that go into the meta-analysis.

Selecting a row tells you which kind it is: the details panel says
*Shared atlas* or *Subset of GSE183852* next to the name. KOSMIC reads
this from the study's provenance record, not from the folder name.

## The selected study

The row highlighted in the table is the selected study. The Project
workspace's actions apply to it: *Add Data* imports into it, and
*Open in scRNA* and *Open in DE* open it, as does clicking scRNA
Analysis or Differential Expression in the left-hand rail while this
page is showing. No separate step is needed to make a study active.

The **Active** chip marks the study the scRNA and DE workspaces are
currently holding. It follows your selection when you open one of them.

## Where to go next

The details panel's **Next step** card says what the selected study
needs. Broadly: a study with data but no conditions goes to **scRNA →
Inspect**; a study that has been through QC goes to **Differential
Expression**; once two or more studies have DE results, **Review** offers
*Go to Meta-Analysis*.
