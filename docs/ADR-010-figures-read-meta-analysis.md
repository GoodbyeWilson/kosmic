# ADR-010: Figures reads the Meta-Analysis workspace

- Status: Proposed (2026-10-08). Agreed by Kirk; awaiting Calum.

## Problem

The Figures workspace reads only the scRNA and Differential Expression
workspaces, so meta-analysis results cannot be exported as figures. A
legacy copy of KOSMIC had meta-analysis figure pages (volcano, forest,
leave-one-out) that read saved result files from the top of
`meta_analysis/`. Since ADR-007 those files sit in one folder per cell
type, and the saved tables hold only pooled values: a forest or
leave-one-out plot also needs each study's effect and standard error.

## Decision

1. The Figures workspace holds a reference to the Meta-Analysis
   workspace, set by the application window in the same way as the scRNA
   and DE references, and reads the results the Meta-Analysis workspace
   currently holds. It does not read meta-analysis files from disk. To
   export figures for another selection, the user runs or reloads that
   meta-analysis first.
2. The Meta-Analysis workspace exposes these results through one method,
   `figure_inputs()`, returning the gene- and pathway-level pooled tables,
   the studies' own DE tables (as loaded from each study's DE file for
   the current selection) and the selection's name. Figure pages do not
   reach into the workspace's pages.
3. Per-study effects and standard errors come from those study DE tables.
   The meta-analysis output format is unchanged.
4. Leave-one-out estimates are re-pooled for the one gene or pathway
   shown, with the pooling function the meta-analysis used
   (`kosmic.meta_analysis.pooling_dispatch`).
5. Forest and leave-one-out exports write a text file of the numbers
   behind the figure next to it (pooled estimate, heterogeneity and the
   per-study or per-fold rows). `FigurePage` gains a hook,
   `_export_companions`, that pages override to write such files.

## Consequences

- Meta-analysis figures are available only after a meta-analysis has run
  (or been reloaded) in this session.
- The figures and the Meta-Analysis workspace always describe the same
  run and the same study selection.
- CODEBASE.md's description of the Figures workspace names the
  Meta-Analysis workspace as a source.
- The batch-correction page from the legacy copy is not ported: KOSMIC
  does not keep an uncorrected embedding.
