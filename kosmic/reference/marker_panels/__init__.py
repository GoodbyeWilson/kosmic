"""Curated marker panels for cell-type figures.

A marker panel lists a few genes for each cell type of a tissue, chosen
because each gene is expressed in that type and little elsewhere. The
Figures workspace draws a panel as a marker dot plot and as one UMAP per
cell type coloured by the panel's module score.

Each panel is a JSON file in this directory with a ``name``, a
``species`` and ``panels`` (``{cell_type: [gene, ...]}``, in the order the
figures show them). ``builtin_panels()`` globs them, so adding a panel
means adding a JSON file here.
"""

from __future__ import annotations

import json
from importlib.resources import files


def builtin_panels() -> dict:
    """Every bundled panel as ``{name: {cell_type: [gene, ...]}}``."""
    panels = {}
    root = files("kosmic.reference.marker_panels")
    for entry in sorted(root.iterdir(), key=lambda p: p.name):
        if entry.name.endswith(".json"):
            spec = json.loads(entry.read_text(encoding="utf-8"))
            panels[spec["name"]] = spec["panels"]
    return panels
