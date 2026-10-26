"""
Enrichment Libraries
=====================
Load a gene-set library for enrichment, with a record of where it came from.

A library is either bundled with KOSMIC (a published file, used exactly as
released) or downloaded from Enrichr when the enrichment runs. Enrichr
serves a library under a fixed name but keeps no version history, so a
download is fingerprinted as well as named: the record says what was
used, from where and when, and a SHA-1 of the content shows whether a
later run used the same gene sets.

Pure Python -- no GUI imports.
"""

from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from importlib.resources import files

#: Libraries bundled with KOSMIC: library id -> description.
BUNDLED = {
    'MSigDB_Hallmark_2026.2_Hs': {
        'package': 'kosmic.reference.msigdb',
        'file': 'h.all.v2026.2.Hs.symbols.gmt',
        'provider': 'MSigDB (file bundled with KOSMIC)',
        'version': '2026.2.Hs',
        'url': ('https://data.broadinstitute.org/gsea-msigdb/msigdb/release/'
                '2026.2.Hs/h.all.v2026.2.Hs.symbols.gmt'),
        'citation': 'Liberzon et al., Cell Systems 2015',
    },
}


def is_bundled(library_id: str) -> bool:
    return library_id in BUNDLED


def _parse_gmt_text(text: str) -> dict:
    pathways = {}
    for line in text.splitlines():
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 3:
            continue
        genes = [g.strip() for g in parts[2:] if g.strip()]
        if parts[0].strip() and genes:
            pathways[parts[0].strip()] = genes
    return pathways


def _source(library_id, content: bytes, gene_sets: dict, **fields) -> dict:
    """The record of one library load."""
    return {
        'library': library_id,
        **fields,
        'sha1': hashlib.sha1(content).hexdigest()[:12],
        'n_sets': len(gene_sets),
        'n_genes': len({g for genes in gene_sets.values() for g in genes}),
        'accessed': datetime.now(timezone.utc).isoformat(timespec='seconds'),
    }


def load_library(library_id: str, timeout: int = 60):
    """Load a library for enrichment.

    Parameters
    ----------
    library_id : str
        A key of :data:`BUNDLED`, or an Enrichr library name
        ('KEGG_2026', 'MSigDB_Hallmark_2020', ...).

    Returns
    -------
    (gene_sets, source)
        ``gene_sets`` is ``{set_name: [genes]}``, unaltered. ``source``
        records the load: 'library', 'provider', 'version' (bundled files),
        'url', 'file' (bundled files), 'sha1' (of the file, or of the text
        Enrichr served), 'n_sets', 'n_genes' and 'accessed' (UTC time of
        the load).
    """
    if library_id in BUNDLED:
        spec = BUNDLED[library_id]
        content = (files(spec['package']) / spec['file']).read_bytes()
        gene_sets = _parse_gmt_text(content.decode('utf-8'))
        return gene_sets, _source(
            library_id, content, gene_sets, provider=spec['provider'],
            version=spec['version'], file=spec['file'], url=spec['url'],
            citation=spec['citation'])

    from kosmic.reference.pathways.enrichr import (
        fetch_library_text, library_url, parse_library_text,
    )
    text = fetch_library_text(library_id, timeout)
    gene_sets = parse_library_text(text)
    return gene_sets, _source(
        library_id, text.encode('utf-8'), gene_sets,
        provider='Enrichr (downloaded when run)', url=library_url(library_id))
