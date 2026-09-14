"""Pathway / gene-set infrastructure.

`PathwayRegistry` (in-memory store), GMT/JSON/CSV I/O (`formats`),
parent-grouping rules (`hierarchy`), and the Enrichr web client
(`enrichr`).

Curated gene-set data ships as plain GMT files under
`kosmic/reference/gene_sets/`. `builtin_registry()` returns a lazily-loaded
`PathwayRegistry` populated by globbing every `.gmt` / `.json` / `.csv`
in that directory. Adding a new collection means dropping a file there;
nothing in this package needs to change.
"""

from importlib.resources import files

from kosmic.reference.pathways.hierarchy import (
    derive_hierarchy,
    build_parent_gene_sets,
)
from kosmic.reference.pathways.registry import PathwayRegistry
from kosmic.reference.pathways import formats


_BUILTIN_REGISTRY: PathwayRegistry | None = None


def builtin_registry() -> PathwayRegistry:
    """Return the singleton registry of built-in gene-set collections.

    Loads every `.gmt`/`.json`/`.csv` under `kosmic/gene_sets/` on first
    call. Each filename stem becomes the collection name.
    """
    global _BUILTIN_REGISTRY
    if _BUILTIN_REGISTRY is None:
        reg = PathwayRegistry()
        reg.load_directory(files("kosmic.reference.gene_sets"))
        _BUILTIN_REGISTRY = reg
    return _BUILTIN_REGISTRY


__all__ = [
    "derive_hierarchy",
    "build_parent_gene_sets",
    "PathwayRegistry",
    "builtin_registry",
    "formats",
]
