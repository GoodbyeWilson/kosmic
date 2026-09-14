"""
Enrichr Provider
=================
Fetch gene set libraries from Enrichr (Ma'ayan Lab).
No authentication required. Returns {pathway_name: [gene_list]}.

API: https://maayanlab.cloud/Enrichr/
Pure Python — no GUI imports.
"""

import urllib.request

BASE_URL = "https://maayanlab.cloud/Enrichr"

# Curated list of the most useful pathway libraries for researchers.
# Users can also fetch the full catalogue via list_libraries().
RECOMMENDED_LIBRARIES = [
    "KEGG_2026",
    "Reactome_Pathways_2024",
    "WikiPathways_2024_Human",
    "GO_Biological_Process_2025",
    "GO_Cellular_Component_2025",
    "GO_Molecular_Function_2025",
    "MSigDB_Hallmark_2020",
    "BioPlanet_2019",
    "Reactome_2022",
    "KEGG_2019_Human",
    "Panther_2016",
    "WikiPathway_2023_Human",
    "TRRUST_Transcription_Factors_2019",
    "ChEA_2022",
]


def fetch_library(library_name, timeout=60):
    """Download a complete gene set library from Enrichr.

    Parameters
    ----------
    library_name : str
        Exact library name (e.g. 'KEGG_2026', 'MSigDB_Hallmark_2020').
    timeout : int
        Request timeout in seconds.

    Returns
    -------
    dict
        {pathway_name: [gene_list]}.

    Raises
    ------
    urllib.error.URLError
        If the request fails.
    ValueError
        If the response is empty or unparseable.
    """
    url = f"{BASE_URL}/geneSetLibrary?mode=text&libraryName={library_name}"
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        text = resp.read().decode("utf-8")

    if not text.strip():
        raise ValueError(f"Empty response for library '{library_name}'")

    pathways = {}
    for line in text.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        name = parts[0].strip()
        # parts[1] is often empty or a description
        genes = [g.strip() for g in parts[2:] if g.strip()]
        if name and genes:
            pathways[name] = genes

    return pathways
