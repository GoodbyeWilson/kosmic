"""Curated reference data and the Python that loads it.

Two flavours of reference, both pure-Python (no Qt):

- `pathways/` -- pathway / gene-set infrastructure (`PathwayRegistry`,
  GMT/JSON/CSV I/O, parent-grouping rules, Enrichr web client). Backed
  by GMT files in `gene_sets/`.
- `markers/` -- cell-type marker DBs (PanglaoDB, CellMarker 2.0).
  Backed by TSV/XLSX files in `marker_dbs/`.

The two have superficially similar shapes (named gene lists from
biomedical sources) but different schemas and access patterns. Pathways
get a flat `{name: [genes]}` registry; markers expose a query interface
that filters by organ, cell type, species, specificity, and canonical
flag.
"""
