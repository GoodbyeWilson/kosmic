"""Bundled MSigDB gene-set collections.

`h.all.v2026.2.Hs.symbols.gmt` -- the MSigDB Hallmark collection, release
2026.2.Hs, human gene symbols: 50 gene sets, 4,384 genes. Downloaded
unaltered from the Broad Institute's public release folder
(https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.2.Hs/).
`mh.all.v2026.2.Mm.symbols.gmt` -- MSigDB's mouse edition of the same
release, 2026.2.Mm, in mouse gene symbols mapped by MSigDB: 50 gene sets,
4,291 genes. Same source folder for Mm.

CC BY 4.0; cite Liberzon et al., Cell Systems 2015.

The file is used exactly as published -- symbols are not updated -- so a
run can be repeated with the same gene sets. Consumed by
`kosmic.reference.pathways.libraries`. Python package only so
`importlib.resources.files("kosmic.reference.msigdb")` can locate the data
under pip-install / PyInstaller bundles.
"""
