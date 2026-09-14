# How many genes the studies in a project actually share.
#
# Combining studies is an inner join on gene names, so the master dataset
# carries only the genes present in *every* study. That number is not
# predictable from the per-study gene counts: two studies with 33,000 genes
# each can share barely 20,000 if one of them was annotated against an older
# GENCODE build and still uses clone-based names (RP11-34P13.3), or if a
# depositor left bare Ensembl IDs on the index where no symbol existed.
#
# The loss is silent -- concat just returns a narrower object -- so this
# module measures it up front, and measures what harmonising gene names to
# current HGNC symbols would recover, before anyone commits to a combine.
from __future__ import annotations

from functools import reduce
from pathlib import Path
from typing import Iterable, Optional


def read_var_names(h5ad_path) -> list[str]:
    """Read an h5ad's gene names without loading any expression data.

    Reads the ``var`` index alone -- a few hundred kilobytes of strings --
    so this stays cheap on files whose matrices run to tens of gigabytes.
    """
    import h5py
    from anndata.io import read_elem

    with h5py.File(str(h5ad_path), "r") as f:
        if "var" not in f:
            return []
        var = f["var"]
        index_key = var.attrs.get("_index", "_index")
        if isinstance(index_key, bytes):
            index_key = index_key.decode()
        if index_key not in var:
            return []
        return [str(name) for name in read_elem(var[index_key])]


def gene_overlap(datasets: dict, lookup: Optional[dict] = None,
                 protein_coding: Optional[Iterable[str]] = None) -> dict:
    """Summarise the gene sets of several studies and what they share.

    Parameters
    ----------
    datasets : dict
        ``{study_name: h5ad_path}``. A study whose file is missing or
        unreadable is reported with ``n_genes`` of ``None`` and left out
        of the shared set, so one broken file does not make the whole
        summary read as zero.
    lookup : dict, optional
        Old-symbol -> approved-symbol map from
        :func:`kosmic.scrna.inspect.gene_names.load_hgnc_lookup`. When
        given, a second set of figures is reported under
        ``'harmonised'`` showing what the join would yield if every
        study's names were harmonised first.
    protein_coding : iterable of str, optional
        Approved symbols of protein-coding genes, from
        :func:`kosmic.scrna.inspect.gene_names.load_locus_groups`. When
        given, the protein-coding subset of each figure is reported too.

    Returns
    -------
    dict
        ``{'studies': {name: {'n_genes', 'n_protein_coding', 'error'}},
        'shared': int, 'shared_protein_coding': int|None, 'union': int,
        'harmonised': {'shared', 'shared_protein_coding'}|None}``.
        ``shared`` is the size of the intersection across every readable
        study -- the gene count a combined dataset would have.
    """
    pc = set(protein_coding) if protein_coding is not None else None

    studies: dict[str, dict] = {}
    sets: dict[str, set] = {}
    for name, path in datasets.items():
        entry = {"n_genes": None, "n_protein_coding": None, "error": None}
        try:
            names = read_var_names(path)
        except Exception as exc:                      # unreadable / not h5ad
            entry["error"] = str(exc)
            studies[name] = entry
            continue
        if not names:
            entry["error"] = "no gene names found"
            studies[name] = entry
            continue
        genes = set(names)
        sets[name] = genes
        entry["n_genes"] = len(genes)
        if pc is not None:
            entry["n_protein_coding"] = len(genes & pc)
        studies[name] = entry

    out = {"studies": studies, "shared": 0, "shared_protein_coding": None,
           "union": 0, "harmonised": None}
    if not sets:
        return out

    shared = reduce(set.intersection, sets.values())
    out["shared"] = len(shared)
    out["union"] = len(reduce(set.union, sets.values()))
    if pc is not None:
        out["shared_protein_coding"] = len(shared & pc)

    if lookup:
        harmonised = [{lookup.get(n, n) for n in genes}
                      for genes in sets.values()]
        h_shared = reduce(set.intersection, harmonised)
        out["harmonised"] = {
            "shared": len(h_shared),
            "shared_protein_coding": (len(h_shared & pc)
                                      if pc is not None else None),
        }
    return out


def project_gene_overlap(project_dir, accessions: Optional[Iterable[str]] = None,
                         **kwargs) -> dict:
    """Run :func:`gene_overlap` over a project's studies.

    Each study contributes its newest ``processed_data/*.h5ad`` -- the
    same file the rest of the app treats as the study's dataset. Studies
    with nothing processed yet are skipped entirely rather than reported
    as errors: they simply have no genes to contribute.
    """
    from kosmic.paths import list_studies, processed_data_dir, study_dir

    root = Path(project_dir)
    names = list(accessions) if accessions is not None else list_studies(root)

    datasets = {}
    for accession in names:
        processed = processed_data_dir(study_dir(root, accession))
        if not processed.is_dir():
            continue
        h5ads = sorted(processed.glob("*.h5ad"),
                       key=lambda p: p.stat().st_mtime, reverse=True)
        if h5ads:
            datasets[accession] = h5ads[0]
    return gene_overlap(datasets, **kwargs)
