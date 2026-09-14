"""Bundled MGI human-mouse 1:1 ortholog table.

`mgi_human_mouse_1to1.tsv` -- human HGNC symbol -> mouse MGI symbol,
restricted to homology classes with exactly one gene per species (no
ambiguous many-to-many families). Built from MGI's
`HOM_MouseHumanSequence.rpt` (informatics.jax.org). Consumed by
`kosmic.scrna.inspect.detection.format_gene_for_species` so that
matching a human-authored gene set (pathways, markers) against mouse
data uses real orthology instead of a naive upper/lower-case guess.
Python package only so `importlib.resources.files(...)` can locate the
data under pip-install / PyInstaller bundles.
"""
