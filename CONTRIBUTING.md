# Contributing to KOSMIC

This file covers how to get set up, the development conventions, and
what to do before opening a PR.

For what KOSMIC *is* and how to use it, see [`README.md`](README.md).
For how the code is currently laid out, see [`CODEBASE.md`](CODEBASE.md).

---

## Getting set up

```sh
git clone https://github.com/GoodbyeWilson/kosmic.git
cd kosmic
python -m venv .venv
.venv\Scripts\activate          # Windows
source .venv/bin/activate       # macOS / Linux
pip install -e ".[dev]"
python main.py
```

Python 3.11+. No compiler needed — every dependency has a prebuilt
wheel on all three platforms. `pyproject.toml` is the single source of
truth for dependencies. Conda users: `conda create -n kosmic
python=3.11` and then the same `pip install`; there is no separate
conda path to maintain.

Two things are deliberately **not** required, because neither is on
PyPI: `soupx` and `decontx`. Both are imported lazily and the app
reports their absence rather than failing; the DecontX test skips when
the package is missing.

### Running the app

```sh
python main.py
```

With the venv activated, that is all there is to it.

### Tests and linting

```sh
pytest tests/
ruff check kosmic/ main.py tests/
```

Both must pass before a PR. There is no separate formatter — `ruff`'s
line length is 100 and that is the whole style rule.

### Documentation

The user docs are Markdown under `kosmic/gui/help/content/`, one page
per screen, named by the screen's `help_id`. `mkdocs.yml` builds them
into the website, and a push to `main` publishes it. So a docs change
is an ordinary PR.

In the app, **F1 opens the website at the page for the screen in
front** (`HelpManager.open`). The same Markdown ships inside the
package as the offline fallback: Help → *Offline Help Browser*, and
where F1 lands when no site is configured or it cannot be reached.
The site is `https://docs.scmetaanalysis.com/` (`site_url` in `mkdocs.yml`
and `[help] site_url` in `config.toml`; keep them the same). When the
site is not reachable — which includes before it exists — F1 in a
source checkout starts `mkdocs serve` on `http://127.0.0.1:8001/`
itself and opens the page there, stopping the server when KOSMIC
exits; so every F1 press is the live preview with no setup.

```sh
pip install -e ".[docs]"
mkdocs serve --livereload   # live preview at http://127.0.0.1:8000
                            # (--livereload is needed with click >= 8.2, or edits never rebuild)
mkdocs build --strict   # what CI runs; fails on a broken link
```

The **Developer reference** section of the site is generated from the
package at build time (`dev/docs/gen_ref_pages.py`): one page per
module, showing the module's leading comment block or docstring and
the API from its docstrings (numpy style). To document a module, write
its header and docstrings; nothing else is needed.

Rules that keep the two renderings in step:

- Plain Markdown only: headings, lists, tables, fenced code, and
  `> **Note:**` blockquotes. The in-app renderer is markdown-it with
  no extensions, so Material admonitions and tabs will not show there.
- A new screen gets a page: `tests/test_help_content.py` fails on a
  `help_id` without one, and on a page missing from the `nav` in
  `mkdocs.yml`. Pages not yet written are listed in `KNOWN_MISSING`.
- **Never hand-capture screenshots.** Run
  `python dev/scripts/capture_help_screenshots.py <project>` against a
  realistic project; it drives the app on a real display and writes
  every image under `<section>/img/`. Add a capture entry when a page
  needs a new picture. Re-run after any UI change.

---

## Repository layout

- `kosmic/` — the application. This is what ships.
- `tests/` — the suite; must pass before a PR.
- `docs/` — architecture decision records.
- `kosmic/gui/help/content/` — the user documentation. It is both the
  in-app F1 help and the website (`mkdocs.yml`); see below.
- `dev/` — our benchmarking, validation and one-off tooling. Not
  shipped, not linted, may be rough. See `dev/README.md`.

---

## Development conventions

The architecture — the two-layer split, the page archetypes, the data
contracts between workspaces — is described in `CODEBASE.md`. These
are the conventions for working within it. If you think one is wrong,
open an issue and argue for it rather than working around it.

**1. Respect the enforced boundaries.** `tests/test_architecture.py`
fails a PR that imports PyQt6 outside `kosmic/gui/`, imports a workspace
from `kosmic/gui/shared/`, or adds a page that does not subclass one of
the page archetypes. Subclassing an existing page is fine; the check
follows inheritance.

**2. It stays PyQt6.** Not Streamlit, not a web app, not Electron. The
target is a one-click desktop installer for people who do not want to
run a server.

**3. Use the page archetypes.** New pages subclass `SimplePage`,
`TabbedPage`, `SidebarPage`, `SidebarTabbedPage` or `ModeChooserPage`
rather than assembling their own splitter and sidebar geometry; two
pages that should look alike share a widget (`DatasetOverview` is the
example) rather than copying dimensions.

**4. Structural changes get an ADR.** New workspace, new page archetype,
change to the provenance model, anything that moves a boundary: write it
up in `docs/` as `ADR-NNN-short-name.md` before implementing, and agree
it before the code. Ordinary changes — fixes, new analysis functions,
rework inside a page — do not need one.

An ADR records a *decision* someone might otherwise reopen: the problem,
what was decided, and what it costs (`ADR-002` is the shape). It is not
a change log; finished work is described in `CODEBASE.md` and git. Number
them consecutively from the last one in `docs/`. When a decision is
reversed, add `Superseded by ADR-NNN` to its status line rather than
deleting it, so the numbering and the reasoning both survive.

**5. New analysis modules need tests.** Anything under `kosmic/` that
isn't `kosmic/gui/` should have tests. GUI code is harder to test and
the bar is lower, but layout regressions and widget-lifetime bugs are
testable — see `tests/test_meta_workflow_steps.py`.

**6. Document in the module.** A module's leading docstring or comment
block says what it is for and why it is the way it is; functions have
numpy-style docstrings. The documentation site's Developer Reference
is generated from these on every build. There is no separate place to
describe a module.

---

## Building the desktop package

`packaging/kosmic.spec` builds a one-directory PyInstaller bundle
(`dist/KOSMIC/`, ~700 MB): `pip install -e ".[build]"` then
`pyinstaller packaging/kosmic.spec`. Data files keep their package
paths inside the bundle, so code locates them from `__file__` exactly
as in a checkout; never special-case `sys.frozen`. Set
`KOSMIC_BUILD_CONSOLE=1` to keep a console window and see a startup
traceback. On Windows, build from a conda environment or a python.org
Python — a venv created from a conda Python is missing the DLLs
PyInstaller needs to find (`pyexpat` fails to load).

---

## Qt pitfalls

**Qt object lifetime.** A widget with no parent and no Python reference
is garbage-collected, and Qt then deletes its children. Symptoms are
`RuntimeError: wrapped C/C++ object ... has been deleted`, a control
that silently stops working, or a hard interpreter crash with no
traceback. Every widget needs a parent or a layout. In tests, hold the
`QApplication` in a module-scoped fixture — `QApplication.instance() or
QApplication([])` as a bare statement collects it immediately and
crashes mid-test.

**A parentless `QWidget` is a window.** Calling `setVisible(True)` on
one before adding it to a layout flashes an empty frame on screen.

**`QLabel` does not wrap by default.** A non-wrapping label reports its
entire string as its minimum width, which forces its whole container
that wide. Use `HintLabel` (wraps by default) or set `setWordWrap(True)`.

**`isVisible()` is False whenever the page is not the active tab.**
Don't use it to decide what a widget's state is; check the underlying
setting instead.

**Headless font metrics are wrong.** Under `QT_QPA_PLATFORM=offscreen`
Qt substitutes a fallback font roughly twice the width of a real UI
font. Offscreen rendering is fine for checking layout structure and
colours, but do not tune pixel widths against it. Real-font checks come
from `dev/scripts/capture_help_screenshots.py`, which drives the app
on a real display.

**A stacked widget is as tall as its tallest page.** Put
`setAlignment(Qt.AlignmentFlag.AlignTop)` (or a trailing stretch) on
each page's layout, or shorter pages spread their widgets down the
stack's full height.

**An ID selector beats a pseudo-state.** `QPushButton#primary_button`
outranks `QPushButton:disabled`; a themed widget needs its own
`:disabled` rule or it never looks disabled.

**One screen shape, one widget.** When two pages should look alike
(the scRNA and DE dataset screens), give them a shared widget that
owns the geometry (`DatasetOverview`) rather than copying numbers
between them; copies drift within weeks.

---

## Provenance

Every analysis stage writes what it did to a `provenance.json` sidecar,
and the Methods pages render those records. If you add an analysis step,
record it:

- Per-study stages **append** (`provenance.record_stage`) — each one
  transforms the study's h5ad, so the chain is a real lineage.
- Meta-analysis stages **replace** (`record_meta_stage`, which passes
  `replace=True`) — a run overwrites its output file, so appending one
  entry per click builds a history describing results that no longer
  exist.

Add a title for any new stage key to `_STAGE_TITLES` in
`kosmic/provenance.py`, or it renders as a raw key.

---

## Opening a PR

1. Branch off `main`.
2. Make sure `pytest tests/` and `ruff check kosmic/ main.py tests/` pass.
   For a change that touches more than one workspace, also run
   `python dev/scripts/smoke_drive.py <project>` against a real project.
3. Update `CODEBASE.md` if you changed a boundary or a hand-off between
   workspaces. Per-module documentation is the module's own header
   and docstrings — the docs site's Developer reference is generated
   from them.
4. Open an issue first for anything non-trivial, so we can agree the
   approach before you spend time on it.
5. No AI attribution in commit messages or PR descriptions: no
   assistant co-author trailer, no "generated with" line. GitHub turns
   such trailers into a listed contributor. Turn the byline off in your
   own assistant's settings rather than stripping it by hand.

Small, focused PRs get reviewed faster than large ones.
