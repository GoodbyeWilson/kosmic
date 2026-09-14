"""Bundled GWAS reference data.

`gene_positions_GRCh37.tsv` -- HGNC gene loci on GRCh37, used for
nearest-gene assignment in `dev.bench.gwas_overlap`. Python
package only so `importlib.resources.files("kosmic.reference.gwas")`
can locate the data under pip-install / PyInstaller bundles.
"""
