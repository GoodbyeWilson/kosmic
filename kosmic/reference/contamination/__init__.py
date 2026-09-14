"""Ambient-contamination signature panels.

A contamination panel is a set of genes *specific to a dominant, RNA-rich
cell type* (e.g. cardiomyocytes in heart) that leak into other cell types'
droplets/nuclei as ambient RNA. Scoring the panel inside a *non*-source
cell type gives a ``percent.mito``-style per-cell contamination fraction
(see ``kosmic.scrna.qc.signatures``), which can then be used to check
whether a differential-expression signal is driven by ambient spill-over
from the dominant cell type rather than by biology of the cell type under
study.

These are distinct from the pathway gene sets in ``gene_sets/`` (which are
biological programmes for scoring) and must NOT appear in the pathway
picker -- hence a separate package. Each panel is a plain JSON file in
this directory; ``builtin_panels()`` globs them. Adding a panel means
dropping a JSON file here; nothing in this module changes.
"""

from __future__ import annotations

import json
from importlib.resources import files


def builtin_panels() -> dict:
    """Load every bundled contamination panel.

    Returns
    -------
    dict
        ``{key: {"display_name", "dominant_cell_type", "tissue",
        "species", "provenance", "genes": [...]}}``. The ``key`` is the
        JSON's ``key`` field (falling back to the filename stem).
    """
    panels = {}
    root = files("kosmic.reference.contamination")
    for entry in sorted(root.iterdir(), key=lambda p: p.name):
        if not entry.name.endswith(".json"):
            continue
        spec = json.loads(entry.read_text(encoding="utf-8"))
        key = spec.get("key") or entry.name[: -len(".json")]
        spec["key"] = key
        panels[key] = spec
    return panels


def get_panel(key: str) -> dict | None:
    """Return a single bundled panel spec by key, or ``None``."""
    return builtin_panels().get(key)


__all__ = ["builtin_panels", "get_panel"]
