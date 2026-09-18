# ADR-004: Headless execution

- Status: Proposed (2026-09-18). Awaiting Calum and Kirk.

## Problem

Every analysis step in KOSMIC runs inside the desktop application, on
the machine the window is open on. Three things follow from that:

- Work is bounded by that machine. A million-cell atlas needs 35 GB for
  its cluster step today (20 GB after ADR-003); the University's cluster
  (ARCHIE-WeSt) has nodes with 3 TB, and no way to use them for KOSMIC.
- A run is not reproducible from a script. The provenance record says
  what was done, but re-doing it means clicking through the same tabs.
  A software paper reviewer will ask for the command that reproduces
  the worked example.
- Long steps hold the window. The atlas's Harmony step took 25 minutes
  with the application unresponsive to anything else.

The analysis package already has no Qt in it -- `tests/test_architecture.py`
enforces that -- and the DCM atlas's Reichart-LV study was rebuilt on
17 September by calling the same functions the tabs call, from a
script, with no window. The capability exists; it has no entry point.

## Decision

**A command-line entry point runs any pipeline step on a project
directory, with the same code, the same files and the same provenance
as the desktop application.**

```
kosmic run --project "C:\dcm raw" --study reichart_LV --step qc
kosmic run --project "C:\dcm raw" --study reichart_LV --step cluster --harmony sample --resolution 0.3
kosmic run --project "C:\dcm raw" --study _master --step de --cell-type-column cell_type --types Cardiomyocyte
kosmic run --project "C:\dcm raw" --step meta --cell-type Cardiomyocyte
kosmic status --project "C:\dcm raw"
```

Rules:

1. One implementation. The GUI workers and the command line call the
   same functions in `kosmic/`; the workers' `_run` bodies move into
   plain functions the workers then call. A step run headless writes
   the same files and the same provenance stage as the tab does.
2. Parameters have the same names and defaults as the tabs, and are
   recorded in provenance the same way, so a methods text generated
   from a headless run reads identically.
3. The command reads and writes the project directory and nothing else.
   It needs no display and no Qt (`kosmic.gui` is never imported).
   Progress goes to stdout as the lines the output panel would show,
   including the memory line, so a Slurm log is a complete record.
4. `kosmic status` prints what the Project page shows: each study's
   stages, and what is stale.
5. A step that needs a decision the tabs ask for interactively
   (roles, relabels, the Unknown clusters) takes it from a file
   (`study.json` semantics, a JSON of manual labels) or refuses with a
   message naming what is missing. Nothing is guessed.

The desktop application stays the place where results are looked at
and decisions are made; the command line is where the compute can go.

## Consequences

- KOSMIC runs on ARCHIE-WeSt: a conda environment in the project space,
  `pip install -e .`, one Slurm job per step. The 3 TB nodes make any
  atlas size a matter of time, not memory.
- The worked example in the software paper becomes a script that
  reproduces every result from the deposited files.
- Kirk can run the large steps of the DCM project on the cluster and
  open the same project directory on his laptop to look at them.
- A later step (not this ADR): the GUI submitting a step to a remote
  machine and watching its progress. This ADR is its prerequisite.
- Effort: two days. Most of it is moving the worker bodies into
  functions and adding the argument parsing; the functions exist.
