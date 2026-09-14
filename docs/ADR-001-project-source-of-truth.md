# ADR-001: Project as source of truth (study manifest + per-module selection)

- Status: Accepted (Calum, 2026-08-14)
- Phases: A, C, D done (2026-08-14); E (intake extraction) done
  (2026-08-16): GEO + Assemble tools moved to `kosmic/gui/intake/`,
  hosted by the Project "Import External..." menu via the Add Data
  dialog; scRNA Load Data is a single Dataset page. B remains
  (per-module selection; the riskiest phase).

## Problem

KOSMIC has a Project concept, but the project acts as a folder browser.
Each workspace then recreates the analytical context on its own:

- Project discovers studies and their Raw / Processed / DE / Pseudobulk
  status, then forgets why it matters.
- scRNA asks which dataset to load.
- DE asks which processed h5ad to load, again, and asks for the
  sample / condition columns, again.
- Meta asks the user to Add Files for DE results the project already
  knows about.

The user tells KOSMIC the same facts several times, through file
pickers -- which are exactly where non-computational users get lost
(see the presumed-knowledge audit: every dead end is a file-choosing
moment).

## Decision

**Import and configure once at Project level; modules choose project
studies, not files.**

Three rules:

1. Import / register data once, at Project level (scRNA Load Data is
   the import surface; what it produces is registered to the study).
2. Dataset semantics (sample column, condition column, role levels,
   cell-type column) are stored with the study and pre-fill every
   downstream surface. Editable, never re-asked from scratch.
3. Modules select *studies already in the project*. Direct file
   loading remains as an explicit escape hatch ("Import external
   dataset...") -- the exception, not the workflow.

**Project membership is distinct from analytical selection.** There is
no single global "active study": scRNA and DE each hold their own
current study (inheriting sensibly on first visit); Meta holds a study
*set* (default: all studies with DE results). Every module shows a
persistent study header ("GSE183852_Endothelium · Change study") so
inherited selection is always visible, never silent.

## The manifest (`<study>/study.json`)

A small JSON cache at the study root:

```json
{
  "version": 1,
  "updated": "2026-08-14T10:30:00",
  "dataset": {
    "file": "GSE183852_Endothelium.h5ad",
    "n_cells": 38029,
    "n_genes": 29420
  },
  "semantics": {
    "sample_column": "sample",
    "condition_column": "condition",
    "cell_type_column": "cell_type",
    "role_map": {"Donor": "control", "DCM": "disease"}
  }
}
```

**Hard rules (the drift guard):**

- The manifest is a *cache*, never a second source of truth. The h5ad
  (obs columns, `uns['role_map']`, `uns['role_condition_col']`) and the
  pseudobulk `role` column remain authoritative.
- **One writer:** the Inspect tab's save flow (the moment semantics are
  committed to the h5ad) regenerates the manifest. Nothing else writes
  it.
- **Absence is always tolerated.** No manifest means "not configured
  yet", never an error. Old projects work unchanged.
- Asset flags (raw / processed / DE / pseudobulk) are NEVER stored --
  always derived from disk via `paths.study_status()`, so they cannot
  go stale.
- On conflict (manifest names a column the h5ad lacks), the h5ad wins
  and the consumer treats the field as unset.

Module: `kosmic/manifest.py` (Qt-free, tested).

## Phases

- **A. Manifest layer + Project as Study Manager.** `kosmic/manifest.py`
  + tests; Inspect save writes the manifest; Project workspace gains a
  per-study details view (assets from disk, semantics from manifest).
  "Set as active" stays until Phase B replaces it (removing it first
  would leave no selection mechanism).
- **B. Per-module selection.** Workspace-local current study with
  inheritance; persistent study header card per module; Project row
  action becomes "Open study"; `AppWindow.study_changed` contract
  rewired. Riskiest phase: touches the state contract every workspace
  listens to; manual regression pass required (auto-load on switch,
  scRNA→DE `set_adata` handoff, step-completion resets).
- **C. DE setup consumes the manifest.** Study header card + "Change
  study"; column configuration pre-filled from semantics; "Import
  external dataset..." escape hatch replaces the file-loading UI.
- **D. Meta select-studies simplification.** Eligible studies
  (DE results present) listed with checkboxes, all selected by
  default; "Import external result..." escape hatch replaces Add
  Files.

## Consequences

- Project becomes a Study Manager (the table: Study | Raw | Processed |
  Metadata | DE | Pseudobulk, click to inspect/edit).
- The module rail reads as: Project = what data do I have; scRNA =
  analyse a study; DE = compare groups within a study; Meta = combine
  studies; Figure Export = export.
- Help content and the tutorial are written AFTER this lands, so they
  document the new flow once.
- The DE→Meta filesystem contract and the scRNA→DE in-memory handoff
  (CODEBASE.md "When modifying the software") are unchanged by Phase A
  and C/D; only Phase B touches state flow.
