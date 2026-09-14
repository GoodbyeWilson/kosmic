# Species Detection
# Heuristic species detection from gene naming conventions in scRNA-seq data.
# Human genes are uppercase (e.g. HK1, MT-ND1), mouse genes are mixed case (e.g. Hk1, mt-Nd1).
#
# Gene-name *conversion* between species (format_gene_for_species) has two
# tiers: a real ortholog table (kosmic.reference.orthologs, MGI 1:1 human-
# mouse homology) tried first, falling back to the naive upper/Title-case
# guess for anything not in the table. The naive guess alone is wrong for
# any gene family whose mouse symbol isn't just the human symbol Title-cased
# -- common for older Ensembl/CellRanger references that still carry
# pre-2016 nomenclature (e.g. the ATP synthase subunits) even though MGI's
# current symbol has moved on too. A short hand-verified patch for that
# specific family lives in _LEGACY_MOUSE_SYNONYMS below.

from functools import lru_cache
from importlib.resources import files
from typing import Dict, List


def detect_species(gene_names: List[str], sample_size: int = 100) -> str:
    """
    Detect species from gene naming convention.

    Samples the first 'sample_size' genes and checks the fraction that are
    fully uppercase (human convention). If >=50% are uppercase, returns 'human',
    otherwise 'mouse'.

    Parameters
    ----------
    gene_names : list of str
        Gene names from adata.var_names or similar.
    sample_size : int
        Number of genes to sample from the start of the list.

    Returns
    -------
    str
        'human' or 'mouse'
    """
    sample = list(gene_names[:sample_size])
    if not sample:
        return 'human'
    uppercase_count = sum(1 for g in sample if g.isupper())
    return 'human' if uppercase_count >= len(sample) / 2 else 'mouse'


# Pre-2016 HGNC/MGI symbols for the ATP synthase F1/F0 subunits. HGNC's 2016
# rename (ATP5A1->ATP5F1A, ATP5G1->ATP5MC1, etc.) is what current gene sets
# ship under, and MGI's mouse symbol table has followed suit -- but several
# widely-used Ensembl/CellRanger mouse references (this project's included)
# still annotate genes under the old symbol. Hand-verified against MGI;
# narrow and explicit on purpose rather than a broad legacy-synonym table.
_LEGACY_MOUSE_SYNONYMS = {
    'ATP5F1A': 'Atp5a1', 'ATP5F1B': 'Atp5b', 'ATP5F1C': 'Atp5c1',
    'ATP5F1D': 'Atp5d', 'ATP5F1E': 'Atp5e',
    'ATP5MC1': 'Atp5g1', 'ATP5MC2': 'Atp5g2', 'ATP5MC3': 'Atp5g3',
    'ATP5ME': 'Atp5k', 'ATP5MF': 'Atp5j2', 'ATP5MG': 'Atp5l',
    'ATP5PD': 'Atp5h', 'ATP5PF': 'Atp5j', 'ATP5PO': 'Atp5o',
}

# MGI's homology classing is occasionally wrong for short/generic human
# symbols (a coincidental 1:1 class with an unrelated mouse gene rather
# than a real ortholog). Verified case: human CS (citrate synthase) is
# classed against mouse 'Csl', which is not citrate synthase. Block known
# bad entries here rather than trust the table blindly; falls through to
# the naive guess ('Cs', correct).
_ORTHOLOG_TABLE_BLOCKLIST = {'CS'}


@lru_cache(maxsize=1)
def _load_ortholog_table() -> Dict[str, str]:
    """Human HGNC symbol (upper) -> mouse MGI symbol, 1:1 homology classes only.

    Loaded once from the bundled MGI table (see
    kosmic.reference.orthologs). Missing/unreadable file degrades to an
    empty table -- callers fall back to the naive heuristic, they don't
    fail.
    """
    try:
        path = files("kosmic.reference.orthologs") / "mgi_human_mouse_1to1.tsv"
        table: Dict[str, str] = {}
        with path.open("r", encoding="utf-8") as f:
            next(f)  # header: human\tmouse
            for line in f:
                human, _, mouse = line.rstrip("\n").partition("\t")
                if human and mouse:
                    table[human.upper()] = mouse
        return table
    except (OSError, FileNotFoundError):
        return {}


def format_gene_for_species(gene: str, species: str) -> str:
    """
    Format/convert a gene name for the given species namespace.

    Parameters
    ----------
    gene : str
        Gene name, typically in human HGNC (upper-case) format.
    species : str
        'human' or 'mouse'.

    Returns
    -------
    str
        Gene name in the target species' namespace: for 'mouse', the real
        MGI ortholog when 'gene' is a known human symbol with an
        unambiguous 1:1 mouse ortholog (checked via the hand-verified
        legacy-nomenclature patch first, then the bundled MGI table),
        otherwise the naive Title-case guess (first letter upper, rest
        lower; 'MT-' prefix genes go to lower-case 'mt-'). For 'human',
        the naive upper-case guess.
    """
    if species == 'mouse':
        gene_upper = gene.upper()
        if gene_upper in _LEGACY_MOUSE_SYNONYMS:
            return _LEGACY_MOUSE_SYNONYMS[gene_upper]
        if gene_upper not in _ORTHOLOG_TABLE_BLOCKLIST:
            ortholog = _load_ortholog_table().get(gene_upper)
            if ortholog:
                return ortholog
        # Mouse convention: first letter uppercase, rest lowercase
        # Special case for MT- prefix genes
        if gene_upper.startswith('MT-'):
            return 'mt-' + gene[3:].capitalize()
        return gene.capitalize()
    return gene.upper()
