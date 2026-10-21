"""Bundled MGI human-mouse 1:1 ortholog table.

`mgi_human_mouse_1to1.tsv` -- human HGNC symbol -> mouse MGI symbol,
one row per human symbol and one row per mouse symbol. Built from MGI's
`HOM_MouseHumanSequence.rpt` (informatics.jax.org). Each mouse symbol in
that report belongs to one homology class, but some classes hold one
human gene and several mouse genes (for example ALDH1A1 -> Aldh1a1,
Aldh1a7; AGTR1 -> Agtr1a, Agtr1b). For such a human symbol the table
keeps only the row whose mouse symbol, upper-cased, equals the human
symbol (ALDH1A1 -> Aldh1a1); if no row matches in that way, the human
symbol is left out of the table (AGTR1, SCD, ABCB1). Consumed by
`kosmic.scrna.inspect.detection.format_gene_for_species` so that
matching a human-authored gene set (pathways, markers) against mouse
data uses MGI orthology instead of a naive upper/lower-case guess.
Python package only so `importlib.resources.files(...)` can locate the
data under pip-install / PyInstaller bundles.
"""
