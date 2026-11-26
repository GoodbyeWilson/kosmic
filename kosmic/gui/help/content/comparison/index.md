# Methods Comparison

The Methods Comparison workspace compares the genes called by merged
pseudobulk with the genes called by four meta-analysis methods (REML,
SumRank, Stouffer and vote counting), for one counts source and one cell
type of the project.

- **Merged pseudobulk** is one DESeq2 model fitted to the donors of every
  study in the shared atlas (`_master`). It adjusts for the covariates the
  per-study models used, plus study.
- **Meta-analysis** analyses each study separately and combines the
  study-level results. These are the same methods as in the
  Meta-Analysis workspace.

The workspace needs an open project that has a shared atlas, made with
*Combine* on the Project page, and saved per-study DE runs for the cell
type.

> **Note:** This workspace is being built. Its controls, the design check
> and the summary will be described here as they are added.

The Meta-Analysis workspace also has a mode named Methods Comparison. That
mode compares meta-analysis methods with each other, using leave-one-out
replication and held-out AUC, and does not include merged pseudobulk.
