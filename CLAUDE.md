# CLAUDE.md — instructions for AI coding assistants

## Start here

1. Read `CODEBASE.md` — the current architecture, core concepts and
   data flow. It is the one authoritative description of how the
   application is structured.
2. Read `CONTRIBUTING.md` — development conventions and Qt pitfalls.
3. Check the TODO block at the top of `main.py` for in-flight work.

Responsibilities of the documentation files:

- `CLAUDE.md` (this file) — how AI assistants are expected to work in
  this repository.
- `CODEBASE.md` — the current architecture and codebase map.
- `CONTRIBUTING.md` — human development conventions.
- Module docstrings and headers — implementation and API detail,
  rendered into the documentation site's Developer Reference.
- `docs/` — architecture decision records.

Do not duplicate content between these files.

## Project

KOSMIC (Kirk's Open-source Single-cell Meta-analysis Integration and
Comparison) is a PyQt6 desktop application for single-cell and
single-nucleus RNA-seq analysis. It supports study-level processing and
differential expression, cross-study meta-analysis, dataset integration
and publication-ready figure export.

The application is intended for distribution as a standalone desktop
application via GitHub Releases and for publication as a JOSS software
paper. See `CODEBASE.md` for the workspace structure and data flow.

## Rules

- **Preserve the architectural boundaries** in `CODEBASE.md`: the
  analysis package does not import Qt; shared GUI code does not import
  workspaces; every page subclasses a page archetype.
  `tests/test_architecture.py` enforces these and will fail otherwise.
- **Do not make architectural or structural changes without
  approval.** If a task appears to require changing an established
  boundary, data contract or project structure, explain the proposed
  change and wait for confirmation before implementing it. Approved
  changes of this kind require an ADR in `docs/`.
- **Verify APIs before using them.** Open the definition before writing
  code against a function, class or attribute. Signatures, attribute
  names and defaults in this codebase have all been wrong when guessed.
- **Test empirical assumptions.** When a question can be answered from
  the code, the data or the running application, verify it rather than
  estimating.
- **Work incrementally.** Complete and test each coherent change before
  moving to the next part of a larger task.
- **Report inconsistencies.** Broken code, a stale file or a wrong
  statement in these documents should be surfaced, not worked around.
- **Tests and lint must pass** before a task is reported as done:
  `pytest tests/` and `ruff check kosmic/ main.py tests/`.
- **Update documentation with the change.** A changed boundary or
  hand-off updates `CODEBASE.md`; a changed screen updates its page
  under `kosmic/gui/help/content/` and, if the screen looks different,
  the screenshots (`dev/scripts/capture_help_screenshots.py`); a new
  module documents itself in its header and docstrings.

## Development

```sh
python -m venv .venv && .venv\Scripts\activate   # or source .venv/bin/activate
pip install -e ".[dev,docs]"
python main.py
pytest tests/
ruff check kosmic/ main.py tests/
mkdocs serve --livereload                        # documentation site
```

A clean virtual environment supports the full test suite without
additional setup. R is required only for importing Seurat objects.

To inspect a GUI change without a display, set
`QT_QPA_PLATFORM=offscreen`, build a `QApplication`, apply the theme,
instantiate the widget and call `widget.grab().save(path)`. Layout
structure and colours are accurate; text renders as boxes and font
metrics are roughly twice their real width, so never tune pixel sizes
against an offscreen render. `dev/scripts/capture_help_screenshots.py`
drives a real window for real-font checks.

When relaunching the app, kill the existing process first (match on
the command line, not the process name) or two windows appear.

## KOSMIC Web

A separate static educational explorer reuses this package: its
exporters import from `kosmic.de`, so changes to KOSMIC's analysis
code change what that site shows. Keep the boundary clean — analysis
logic belongs in `kosmic/`; the website consumes its output. A local
checkout may exist at `C:/00_Code/kosmic-web` on the primary
development machine; do not assume that path exists elsewhere.
