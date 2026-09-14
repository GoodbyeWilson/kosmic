"""Bundled Gene Ontology data + writable cache.

Ships with `go-basic.obo`, `gene2go_human.gz`, `hgnc_entrez.tsv`. The
"Update GO Data" feature in the GUI re-downloads these from NCBI / GO
into the same directory. Consumed by `kosmic.de.go_enrichment`.

Note: under PyInstaller bundles `importlib.resources` is read-only, so
the refresh path needs a writable cache dir (e.g. user app-data) once
the project ships as a frozen bundle.
"""
