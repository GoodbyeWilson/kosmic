# ADR-002: Project workspace as a five-step workflow

- Status: Accepted (Calum, 2026-09-11)

## Problem

Every analysis workspace presents itself as a numbered workflow in the
explorer pane (scRNA: eight steps; DE and Meta: per mode). The Project
workspace — the one a new user meets first — passes an empty step list,
so the explorer shows only a file tree, and the page's own sidebar is an
accretion of unrelated things: project name, Switch Project, Add Study,
Add Data, Combine, an asset filter list, and a health card.

Nothing on the screen says what to do next. A user with an empty project
sees a table with no rows and a column of buttons with no order. The
help page (`project/overview`) describes the buttons; it cannot supply
the sequence the interface lacks.

The Home workspace already orients ("which route applies to me"); the
Project workspace is where the work starts, and it has no route.

## Decision

Give the Project workspace the same shape as the other workspaces:
workflow steps in the explorer pane, a step-specific sidebar on the
page, project health in the explorer's STATUS panel.

Five steps, each with a completion rule computed from what is on disk:

| # | step | sidebar holds | complete when |
|---|---|---|---|
| 1 | **Project** — choose the folder for this analysis | open / create / recent projects; current path | a project folder is open (reopened from the last session automatically) |
| 2 | **Studies** — add the datasets to compare | + Add Study; remove selected | at least two studies (real datasets, not the atlas or subsets) |
| 3 | **Data** — import or download data for each study | Add Data… routes for the selected study; which studies still lack data | every study has a processed dataset |
| 4 | **Shared Atlas** *(optional)* — combine studies for consistent cell-type annotation | Create shared atlas; shared-gene count; Propagate labels | the `_master` atlas exists |
| 5 | **Review** — check studies are ready for analysis | readiness filters; Go to Meta-Analysis | at least two studies have DE results and pseudobulk counts |

The three columns keep their existing jobs:

- **Explorer pane** (owned by `AppWindow`): WORKFLOW lists the five
  steps; FILES is the project tree as today; STATUS carries project
  health, which is what `StepStatusPanel` is for in every other
  workspace.
- **Page sidebar** (`SidebarPage`): shows only the controls for the
  selected step. Clicking a step swaps the sidebar. The old sidebar's
  name/path block is dropped (the explorer's directory indicator has
  it), Switch Project moves into step 1, the health card moves to
  STATUS, the asset filter list becomes step 5's readiness view.
- **Content**: the studies table and the details panel, unchanged
  across steps. With no project open, the content shows the
  recent-projects view instead of an empty table. The details panel's
  head names what kind of folder the row is — study, atlas, or subset
  (derived from `provenance.json`'s `source` pointer, not the name).

The empty-state page and the `QStackedWidget` that swapped it for the
loaded page go away; "no project" is simply step 1 being the only step
that can be active.

## Consequences

- `ProjectWorkspace` gains `WORKFLOW_STEPS`, `on_sidebar_step()`,
  `iter_completed_steps()` and `compute_sidebar_status()`, and
  `AppWindow._activate_project` goes through
  `_activate_workspace_common` like every other workspace.
- Step completion lives in a Qt-free module (`kosmic/project_workflow.py`)
  so it is unit-testable and the docs can state the rules exactly.
- Subsets and the atlas are counted honestly: step 2 counts real
  datasets; step 5 counts anything with DE + pseudobulk, since the
  meta-analysis pools whatever was analysed (typically cell-type
  subsets). `_master` is never "ready for meta-analysis".
- The help pages for the workspace follow the five steps, one page per
  step, and F1 opens the page for the selected step.
- "Skipped" is not observable for the optional step 4; it reads as
  incomplete until an atlas exists. The label says optional, and step 5
  does not depend on it.
