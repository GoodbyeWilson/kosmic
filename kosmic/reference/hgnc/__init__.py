"""Bundled HGNC approved-symbol table.

`hgnc_symbols.tsv` -- previous-symbol / alias -> approved-symbol map,
used for gene-name harmonisation across datasets. Consumed by
`kosmic.scrna.inspect.gene_names`. Python package only so
`importlib.resources.files("kosmic.reference.hgnc")` can locate the
data under pip-install / PyInstaller bundles.
"""
