# ADR-009: One export path for figures

- Status: Proposed (2026-10-08). Agreed by Kirk; awaiting Calum.

## Problem

Figures exported from the Figures workspace were not reliably usable in
a vector graphics editor or in a manuscript:

1. Text in a PDF stays editable only with `pdf.fonttype = 42`. Eight
   plotting modules set it as a side effect inside their own functions;
   the others, including the meta-analysis forest plots, did not. Whether
   an export had editable text depended on which figure had been drawn
   earlier in the session.
2. `svg.fonttype` was not set anywhere, so SVG text was always written
   as outlines.
3. The save dialog offered PNG, PDF and SVG but always wrote a PDF next
   to a PNG and ignored the choice otherwise.
4. The figure's own background colour was exported, so an export made
   under the dark theme had a dark background.

## Decision

`kosmic/visualisation/style.py` holds the export settings and the one
function that writes a figure to disk:

- `apply_export_rcparams()` sets `pdf.fonttype = 42`,
  `ps.fonttype = 42` and `svg.fonttype = 'none'` for the process. It
  runs when `kosmic.visualisation` is imported and again in
  `save_figure`; plotting modules no longer set these values.
- `save_figure(fig, path, dpi)` writes exactly the format given by the
  file suffix (PNG, PDF or SVG), on a white background. Scatter
  collections with 20,000 points or more are rasterised first, so a
  UMAP of many cells stays a small file while its axes and text stay
  vector. Smaller collections, such as a volcano plot's genes, stay
  vector so each point can be selected.
- The Figures workspace's export button calls `save_figure` with the
  export DPI from Plot Settings. The default export DPI is 300.

The module sits in `kosmic/visualisation/`, not `kosmic/gui/`, so that
analysis code can call it without importing Qt
(`tests/test_architecture.py`). The DPI is therefore passed in by the
caller rather than read from the GUI settings.

## Consequences

- An exported PDF or SVG has editable text whatever was drawn before it.
- A PNG export writes one file. A user who wants a PDF as well exports
  twice.
- No figure changes appearance in the preview. Exported files change:
  white background, and SVG text is text.
- A user who has already saved Plot Settings keeps their stored DPI;
  the new default applies only where no value is stored.
- Out of scope, for later decisions: physical sizing in mm with journal
  column presets, a shared house style for forest plots, volcano plots
  and embeddings, and writing a statistics file next to a figure.
