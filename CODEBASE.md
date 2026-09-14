# KOSMIC codebase

KOSMIC is a PyQt6 desktop application for the analysis and comparison
of single-cell and single-nucleus RNA-seq datasets. The codebase is
divided into two layers: a Qt-independent analysis package and a
graphical user interface.

The core `kosmic/` package implements the analysis functionality —
data loading, quality control, clustering, cell-type annotation,
differential expression, meta-analysis and dataset integration. These
modules operate on standard Python objects such as AnnData and pandas
DataFrames and do not depend on the GUI. They can be used
programmatically.

The graphical interface is implemented in `kosmic/gui/`. It organises
KOSMIC's functionality into a set of workspaces for managing projects
and carrying out the stages of an analysis. GUI components call
functions from the core package rather than implementing analysis logic
themselves.

A KOSMIC project can contain multiple studies. Each study has its own
data, metadata, provenance and analysis results, while project-level
operations can act across several studies. The distinction between a
project, a study and the current workspace is central to the structure
of the application and is defined below.

Documentation for individual modules is in their docstrings and is
rendered in the Developer Reference section of the documentation site.

## Core concepts

**Project.** The top-level working directory for an analysis. A project
contains one or more study folders and a `meta_analysis/` folder for
project-level outputs. `AppWindow` holds the open project as
`current_project_dir`.

**Study.** One dataset within a project, stored as a folder named by
its accession. A study owns its raw import (`raw_data/`), its processed
AnnData file and provenance record (`processed_data/`), and its
analysis results (`results/`). Two further kinds of folder appear in a
project and are treated as studies by the application: the shared atlas
(`_master/`), produced by combining studies, and subsets, produced by
extracting cell types from a study or from the atlas. A subset's
provenance record points at its parent dataset; the kind of a folder is
determined from that record, not from its name.

**Workspace.** One major area of the application: Project, scRNA
Analysis, Differential Expression, Meta-Analysis or Figures, plus the
Home page. A workspace operates either on a single study (scRNA,
Differential Expression) or on the project as a whole (Project,
Meta-Analysis). Figures reads from the scRNA and Differential
Expression workspaces.

**Step.** One stage within a workspace. Steps are listed in the
Explorer pane, have completion state, and are what F1 resolves to a
help page.

**Selected and active study.** The study highlighted in the Project
workspace's table is the selected study; it is what the Project
workspace's actions apply to. When a per-study workspace is opened from
the Project page, the selected study becomes the active study — the one
that scRNA Analysis and Differential Expression are currently holding
(`current_accession`). There is no separate user action for making a
study active.

## Repository structure

```
main.py                 AppWindow: navigation rail, Explorer pane, page stack, output panel
config.toml             analysis defaults, bound to constants in kosmic/__init__.py
kosmic/                 the installable package
  paths.py              project and study folder layout
  manifest.py           per-study study.json (dataset identity and semantics)
  provenance.py         per-study provenance.json and methods.md
  project_workflow.py   the Project workspace's steps and their completion rules
  scrna/                loading, gene names, roles, QC, clustering, annotation
  de/                   pseudobulk DE, pathway scoring, enrichment, per-cell-type batch runs
  meta_analysis/        pooling methods, consensus, leave-one-out, permutation calibration
  combine/              shared atlas: concatenation, gene merging, label propagation
  visualisation/        matplotlib figures
  reference/            bundled reference data (HGNC, gene sets, markers, GO, GWAS, Olink)
  gui/                  the PyQt6 layer
    shared/             theme, widgets, page archetypes, Explorer pane, plot widgets
    project/            Project workspace
    scrna/              scRNA Analysis workspace
    de_analysis/        Differential Expression workspace
    meta_analysis/      Meta-Analysis workspace
    figure_export/      Figures workspace
    home/               Home page
    intake/             the Add Data dialog and import workers (opened from Project)
    combine/            the Combine and Propagate dialogs (opened from Project)
    help/content/       the user documentation, served by F1 and built into the website
tests/                  pytest suite
dev/                    benchmarks, simulations, screenshot capture, smoke drive (not shipped)
docs/                   architecture decision records
```

## Architectural boundaries

These rules are enforced by `tests/test_architecture.py`.

1. **The analysis package does not import Qt.** Nothing under `kosmic/`
   outside `kosmic/gui/` may import PyQt6.
2. **Shared GUI code does not import workspaces.** `kosmic/gui/shared/`
   provides widgets, page archetypes and the theme, and must not depend
   on any workspace package.
3. **Every page subclasses a page archetype.** Concrete pages inherit
   from `SimplePage`, `TabbedPage`, `SidebarPage`, `SidebarTabbedPage`
   or `ModeChooserPage` in `kosmic/gui/shared/widgets/` (Figures pages
   inherit from their own `FigurePage`).

Two further conventions are not enforced by tests but are relied on
throughout:

- All project and study paths are constructed through
  `kosmic/paths.py`. Other modules should not build these paths
  themselves.
- Every analysis stage records the parameters it ran with through
  `kosmic/provenance.py`. The provenance record is the only description
  of what was run; the Methods steps are rendered from it.

## Data flow between workspaces

The user selects a study in the Project workspace. Opening a per-study
workspace makes that study active, and the workspace loads the study's
processed AnnData file from `processed_data/`.

**scRNA Analysis** processes one study. Each step writes the updated
AnnData back to the study's processed file and records a provenance
stage.

**Differential Expression** operates on one study. When the workspace
is opened it first checks whether the scRNA workspace holds the same
study and whether that dataset is ready for differential expression
(sample and condition columns and roles set, normalised). If so, it uses
the in-memory AnnData object; otherwise it loads the study's processed
file from disk. The hand-off is one-way: Differential Expression never
modifies the scRNA workspace's state. Results are written to the study's
`results/de_analysis/statistics/` folder as
`{accession}_DE_{method}.csv`, with the pseudobulk count matrix beside
it as `{accession}_pseudobulk.csv`.

**Meta-Analysis** operates on the project. It discovers each study's
differential expression results and pseudobulk files on disk, pools
them, and writes its outputs to the project's `meta_analysis/` folder.
It does not read anything from the per-study workspaces in memory.

**Figures** renders figures from the datasets and results currently
held by the scRNA and Differential Expression workspaces.

**Combine** is a dialog opened from the Project workspace. It reads the
selected studies' processed files, writes the shared atlas to
`_master/processed_data/`, and the atlas then appears in the project as
a study. Label propagation is a second dialog that writes atlas labels
into each study's processed file.

## Data contracts

**Raw counts are preserved.** Quality control stores the original
counts in `adata.raw` and `adata.layers['counts']` before normalising.
Differential expression is computed from those counts, never from the
normalised matrix. DecontX writes its output to
`adata.layers['decontX_counts']` and leaves the originals unchanged.

**The processed file holds the study's state.** Every workspace reads and writes
the study's `processed_data/*.h5ad`. Metadata that other workspaces
depend on — the sample and condition columns, the role of each
condition value, cell-type labels — is stored in that file and in the
study's `study.json` manifest, not in application settings.

**Differential expression output.** Meta-Analysis requires, for each
study it pools, a `{accession}_DE_{method}.csv` and a
`{accession}_pseudobulk.csv` containing a `role` column. The role
column is how Meta-Analysis identifies disease and control samples.

**Subsets are studies.** Extracting a cell type in the scRNA Subset
step creates a new study folder whose provenance record's `source`
field points at the parent's processed file. Embeddings and clusterings
are not carried into the subset, because they were computed over the
full dataset.

**Provenance.** Study stages append to the study's record, so the
record is a lineage. Meta-analysis stages replace the previous entry
for the same stage, because a re-run overwrites its output file.

## GUI state and navigation

- **Workspace state persists during navigation.** Each workspace owns
  its own dataset, results and step-completion state. Switching to
  another workspace must not reset the state of the workspace being
  left.
- **Views refresh on dataset version.** A workspace increments a
  version counter whenever its AnnData object changes. Each page
  compares that counter against the version it last rendered when it is
  activated, and refreshes if they differ. A page that does not follow
  this pattern will display stale data.
- **Step completion.** In the Project workspace, completion is derived
  from what exists on disk each time the project is scanned. In the
  scRNA and Differential Expression workspaces, completion is session
  state and resets when the active study changes.
- **The Explorer pane reflects the visible workspace.** Its steps,
  completion indicators and status line are re-installed whenever a
  workspace is activated. Inactive workspaces must not modify it.
- **Selection versus active study.** See *Core concepts*. Navigating
  from the Project page to a per-study workspace adopts the selected
  study; navigating between two per-study workspaces keeps whatever
  study they already hold.

## Where to change what

| To change… | Look in |
|---|---|
| where a file lives within a project or study | `kosmic/paths.py` |
| an analysis default | `config.toml`, then the constant bound in `kosmic/__init__.py` |
| what is recorded for a stage | the `record_stage` call in that stage's worker or page |
| what a Project step considers complete | `kosmic/project_workflow.py` |
| what a scRNA or DE step considers complete | the workspace's `mark_step_complete` calls |
| the page F1 opens for a screen | the page's `help_id` and the matching file under `kosmic/gui/help/content/` |
| a screenshot in the documentation | re-run `dev/scripts/capture_help_screenshots.py` |
| a user-facing explanation | `kosmic/gui/help/content/` (serves both F1 and the website) |
| colours, fonts, widget styling | `kosmic/gui/shared/theme.py` |
| the layout of a dataset overview screen | `kosmic/gui/shared/widgets/dataset_overview.py` |

## Further documentation

- `CONTRIBUTING.md` — setting up, development conventions, Qt pitfalls, the PR process.
- The documentation site's Developer Reference — one page per module,
  generated from docstrings on every build.
- `docs/` — architecture decision records.
- `kosmic/gui/help/content/` — the user documentation.
