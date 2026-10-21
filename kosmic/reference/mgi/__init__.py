"""Bundled MGI mouse symbol table.

`mgi_symbols.tsv.gz` -- one row per current mouse gene or pseudogene:
its MGI symbol, its locus group (MGI's feature type reduced to the four
HGNC locus groups) and its synonyms, which include its earlier symbols.
Built from MGI's `MRK_List2.rpt` (informatics.jax.org) by
`dev/scripts/build_mgi_symbols.py`, which also writes
`mgi_symbols.json`: when the report was downloaded and the table built.
Consumed by `kosmic.scrna.inspect.gene_names` to harmonise mouse gene
names, as `kosmic.reference.hgnc` is for human. Python package only so
`importlib.resources.files("kosmic.reference.mgi")` can locate the data
under pip-install / PyInstaller bundles.
"""
