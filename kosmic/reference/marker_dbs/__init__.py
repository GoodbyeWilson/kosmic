"""Bundled cell-type marker DB files (PanglaoDB TSV, CellMarker 2.0 XLSX).

Only data files live here. The Python that consumes them is in
`kosmic.reference.markers`. This module is a Python package solely so
`importlib.resources.files("kosmic.reference.marker_dbs")` can locate
the bundled data under pip-install / PyInstaller bundles.
"""
