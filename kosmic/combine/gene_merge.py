# Genes that mean different things in different studies.
#
# Harmonising gene names to current HGNC symbols sometimes sums two columns
# into one, because HGNC has since decided the two loci are a single gene --
# HNRNPU-AS1 is a previous symbol of HNRNPU, CTAGE5 of MIA2, C10orf113 of
# NEBL. Summing is the right arithmetic: the counts are UMIs, each assigned
# to exactly one feature, so adding them recovers the locus total.
#
# The problem is that it only happens where a study's source annotation
# carried both features. A study built on a newer GENCODE release already
# has the merged locus as one feature and nothing is summed. So after
# harmonisation the column called 'HNRNPU' is a sum of two features in one
# study and a single feature in another, and the inner join at combine time
# stacks them as if they were the same measurement.
#
# Within a study this is harmless -- every cell gets the same treatment, so
# a DCM-vs-donor contrast is unaffected. Across studies it adds a constant
# offset that inflates between-study heterogeneity, and it makes the gene
# mean something slightly different in each cohort.
#
# It cannot be seen from the names alone: after harmonisation both studies
# just have a column called 'HNRNPU'. The only record of what was summed is
# the harmonisation stage each study wrote to its provenance sidecar, which
# is what this module reads.
from __future__ import annotations

from pathlib import Path
from typing import Optional

HARMONISE_STAGE = "gene_names"


def merge_map_from_provenance(sidecar_dir) -> Optional[dict]:
    """Return ``{approved_symbol: [source names summed into it]}``.

    Reads the most recent gene-name harmonisation stage from the provenance
    sidecar. ``sidecar_dir`` is the directory holding ``provenance.json``,
    which sits beside the working h5ad -- i.e. the study's
    ``processed_data/``, not the study folder itself.

    Returns ``None`` -- distinct from ``{}`` -- when there is no
    harmonisation record, so callers can tell "nothing was merged" apart
    from "we do not know what happened here".
    """
    from kosmic.provenance import load

    record = load(sidecar_dir)
    if not record:
        return None
    merged = None
    for stage in record.get("stages", []):
        if stage.get("stage") == HARMONISE_STAGE:
            merged = stage.get("params", {}).get("merged")
    if merged is None:
        return None
    return {target: list(sources) for target, sources in merged.items()}


def find_merge_asymmetry(studies: dict) -> dict:
    """Find genes built from different numbers of features in each study.

    Parameters
    ----------
    studies : dict
        ``{study_name: (gene_names, merge_map)}``. ``gene_names`` is the
        study's current (post-harmonisation) ``var_names``; ``merge_map``
        is from :func:`merge_map_from_provenance`, or ``None`` when the
        study has no harmonisation record.

    Returns
    -------
    dict
        ``{'genes': {gene: {study: n_features_summed}}, 'unknown':
        [study names with no harmonisation record]}``. ``genes`` holds only
        those present in every study -- a gene missing anywhere is dropped
        by the inner join regardless, so it is not an asymmetry to fix.

    Notes
    -----
    Asymmetry is judged on the *number* of source features summed into a
    gene, which is what changes the magnitude of its counts. A study that
    merely renamed one feature (``AARS`` -> ``AARS1``) contributes one
    feature like everyone else and is correctly not flagged.
    """
    unknown = sorted(name for name, (_, merges) in studies.items()
                     if merges is None)
    known = {name: (set(genes), merges or {})
             for name, (genes, merges) in studies.items()}
    if len(known) < 2:
        return {"genes": {}, "unknown": unknown}

    # Only genes every study carries can survive the inner join.
    shared = set.intersection(*(genes for genes, _ in known.values()))

    out: dict[str, dict] = {}
    for gene in shared:
        counts = {name: len(merges.get(gene, [gene]))
                  for name, (_, merges) in known.items()}
        if len(set(counts.values())) > 1:
            out[gene] = counts
    return {"genes": out, "unknown": unknown}


def project_merge_asymmetry(project_dir, accessions=None) -> dict:
    """Run :func:`find_merge_asymmetry` over a project's studies on disk.

    Each study contributes its newest ``processed_data/*.h5ad`` and the
    merge map recorded beside it. Studies with nothing processed are
    skipped.
    """
    from kosmic.paths import list_studies, processed_data_dir, study_dir
    from kosmic.scrna.load.gene_overlap import read_var_names

    root = Path(project_dir)
    names = list(accessions) if accessions is not None else list_studies(root)

    studies = {}
    for accession in names:
        sdir = study_dir(root, accession)
        processed = processed_data_dir(sdir)
        if not processed.is_dir():
            continue
        h5ads = sorted(processed.glob("*.h5ad"),
                       key=lambda p: p.stat().st_mtime, reverse=True)
        if not h5ads:
            continue
        try:
            genes = read_var_names(h5ads[0])
        except Exception:                       # unreadable: nothing to compare
            continue
        # The sidecar sits beside the h5ad, in processed_data/.
        studies[accession] = (genes,
                              merge_map_from_provenance(h5ads[0].parent))
    return find_merge_asymmetry(studies)
