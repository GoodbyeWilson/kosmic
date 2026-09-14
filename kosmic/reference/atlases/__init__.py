"""Bundled cell-type reference atlases.

Compressed JSON snapshots used by `kosmic.scrna.annotate.reference`
to score query datasets against curated atlases. Currently ships
`heart_atlas.json.gz` (Human Heart Cell Atlas, 704K cells, 12 cell
types). Python package only so
`importlib.resources.files("kosmic.reference.atlases")` can locate
the data under pip-install / PyInstaller bundles.
"""
