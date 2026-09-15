"""Bundled cell-type reference atlases.

Compressed JSON snapshots used by `kosmic.scrna.annotate.reference`
to score query datasets against curated atlases. Ships two human
left-ventricle references built with `dev/scripts/build_lv_reference.py`:
`heartmap_lv_broad.json.gz` (HeartMap, Datar et al. 2026) and
`gao_lv_broad.json.gz` (Gao/Wu, Gao et al. 2026), 13 broad cell types
each; every filter and count is in the file's `description`. Any other
`.json.gz` placed here (built with `dev/scripts/build_reference.py` from
an annotated atlas) is discovered automatically. Python package only so
`importlib.resources.files("kosmic.reference.atlases")` can locate
the data under pip-install / PyInstaller bundles.
"""
