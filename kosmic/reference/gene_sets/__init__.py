"""Bundled gene-set data files.

Only `.gmt` files live here. The Python registry / loaders that consume
them live in `kosmic.reference.pathways`. This module is a Python package
solely so `importlib.resources.files("kosmic.reference.gene_sets")` can
locate the data when the project is pip-installed or PyInstaller-bundled.
"""
