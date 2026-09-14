# GO Term Enrichment
# topGO-style GO enrichment using goatools + custom elim algorithm.
#
# Ships with pre-built data files in kosmic/reference/go/:
#   - go-basic.obo        (GO DAG)
#   - gene2go_human.gz    (NCBI gene-to-GO, human only, ~5 MB)
#   - hgnc_entrez.tsv     (HGNC symbol -> Entrez Gene ID)
#
# An "Update GO Data" function can refresh these from source.
# All gene inputs/outputs use HGNC approved symbols.
# Pure Python -- no GUI imports.

import contextlib
import gzip
import urllib.request
from collections import defaultdict
from importlib.resources import files
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set, Tuple

import pandas as pd
from kosmic import DEFAULT_FDR

# URLs

OBO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"
GENE2GO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
HGNC_URL = (
    "https://www.genenames.org/cgi-bin/download/custom?"
    "col=gd_app_sym&col=gd_pub_eg_id&"
    "status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&"
    "format=text&submit=submit"
)

HUMAN_TAXON = 9606
DEFAULT_CACHE_DIR = Path(str(files("kosmic.reference.go")))

# Download helpers


def _download(url: str, dest: Path, progress_callback: Optional[Callable] = None,
              label: str = "") -> Path:
    """Download a URL to *dest*, showing progress via callback."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if progress_callback:
        progress_callback(f"Downloading {label or dest.name}...")

    req = urllib.request.Request(url, headers={"User-Agent": "KOSMIC/1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        total = int(resp.headers.get("Content-Length", 0))
        downloaded = 0
        chunk_size = 256 * 1024
        with open(dest, "wb") as fh:
            while True:
                chunk = resp.read(chunk_size)
                if not chunk:
                    break
                fh.write(chunk)
                downloaded += len(chunk)
                if progress_callback and total > 0:
                    pct = int(downloaded / total * 100)
                    progress_callback(f"Downloading {label or dest.name}... {pct}%")

    if progress_callback:
        progress_callback(f"Downloaded {label or dest.name}")
    return dest


def download_obo(cache_dir: Path = DEFAULT_CACHE_DIR,
                 progress_callback: Optional[Callable] = None) -> Path:
    """Download go-basic.obo if not cached. Return path."""
    dest = cache_dir / "go-basic.obo"
    if dest.exists():
        return dest
    return _download(OBO_URL, dest, progress_callback, "GO ontology (go-basic.obo)")


def download_gene2go(cache_dir: Path = DEFAULT_CACHE_DIR,
                     progress_callback: Optional[Callable] = None) -> Path:
    """Download gene2go and filter to human. Return path to gene2go_human.gz.

    Downloads the full gene2go.gz (~500 MB), filters to taxon 9606,
    saves as gene2go_human.gz (~5 MB), and deletes the full file.
    """
    dest = cache_dir / "gene2go_human.gz"
    if dest.exists():
        return dest

    # Download full file
    full_path = cache_dir / "gene2go_full.gz"
    _download(GENE2GO_URL, full_path, progress_callback, "NCBI gene2go annotations")

    # Filter to human only
    if progress_callback:
        progress_callback("Filtering to human annotations...")

    taxon_str = str(HUMAN_TAXON)
    with gzip.open(full_path, "rt", encoding="utf-8") as fin, \
         gzip.open(dest, "wt", encoding="utf-8") as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            if line.startswith(taxon_str + "\t"):
                fout.write(line)

    # Remove the large unfiltered file
    with contextlib.suppress(OSError):
        full_path.unlink()

    if progress_callback:
        progress_callback("Human GO annotations ready")
    return dest


def download_hgnc_entrez(cache_dir: Path = DEFAULT_CACHE_DIR,
                         progress_callback: Optional[Callable] = None) -> Path:
    """Download HGNC symbol->Entrez table if not cached. Return path."""
    dest = cache_dir / "hgnc_entrez.tsv"
    if dest.exists():
        return dest
    return _download(HGNC_URL, dest, progress_callback, "HGNC symbol-to-Entrez mapping")


def download_all(cache_dir: Path = DEFAULT_CACHE_DIR,
                 progress_callback: Optional[Callable] = None) -> Tuple[Path, Path, Path]:
    """Download all three data files. Returns (obo_path, gene2go_path, hgnc_path)."""
    obo = download_obo(cache_dir, progress_callback)
    g2g = download_gene2go(cache_dir, progress_callback)
    hgnc = download_hgnc_entrez(cache_dir, progress_callback)
    return obo, g2g, hgnc


def update_data(cache_dir: Path = DEFAULT_CACHE_DIR,
                progress_callback: Optional[Callable] = None) -> Tuple[Path, Path, Path]:
    """Force re-download of all GO data files (replaces existing).

    Use this to refresh the shipped data files with the latest versions
    from NCBI and HGNC.
    """
    # Remove existing files so download functions re-fetch
    for name in ("go-basic.obo", "gene2go_human.gz", "gene2go_full.gz", "hgnc_entrez.tsv"):
        f = cache_dir / name
        if f.exists():
            with contextlib.suppress(OSError):
                f.unlink()

    return download_all(cache_dir, progress_callback)


# Loading / mapping


def load_hgnc_entrez_map(cache_dir: Path = DEFAULT_CACHE_DIR) -> Dict[str, int]:
    """Parse hgnc_entrez.tsv. Return {HGNC_symbol: entrez_id}.

    Skips rows where the Entrez ID is missing or non-numeric.
    """
    path = cache_dir / "hgnc_entrez.tsv"
    if not path.exists():
        raise FileNotFoundError(
            f"HGNC-Entrez mapping not found at {path}. "
            "Call download_hgnc_entrez() first."
        )

    mapping = {}
    with open(path, "r", encoding="utf-8") as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            symbol = parts[0].strip()
            entrez_str = parts[1].strip()
            if symbol and entrez_str:
                try:
                    mapping[symbol] = int(entrez_str)
                except ValueError:
                    continue
    return mapping


def load_gene2go(cache_dir: Path = DEFAULT_CACHE_DIR,
                 taxon: int = HUMAN_TAXON) -> Dict[int, Set[str]]:
    """Parse gene2go_human.gz. Return {entrez_id: set(GO_IDs)}.

    Expects a pre-filtered file (human only). Falls back to gene2go.gz
    with taxon filtering if gene2go_human.gz is not found.

    gene2go format: tax_id, GeneID, GO_ID, Evidence, Qualifier, GO_term, PubMed, Category
    """
    # Prefer the pre-filtered human-only file
    path = cache_dir / "gene2go_human.gz"
    need_filter = False
    if not path.exists():
        path = cache_dir / "gene2go.gz"
        need_filter = True
    if not path.exists():
        raise FileNotFoundError(
            f"gene2go not found in {cache_dir}. "
            "Expected gene2go_human.gz or gene2go.gz."
        )

    taxon_str = str(taxon)
    entrez_to_go: Dict[int, Set[str]] = defaultdict(set)

    with gzip.open(path, "rt", encoding="utf-8") as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            if need_filter and parts[0].strip() != taxon_str:
                continue
            try:
                entrez_id = int(parts[1].strip())
            except ValueError:
                continue
            go_id = parts[2].strip()
            if go_id.startswith("GO:"):
                entrez_to_go[entrez_id].add(go_id)

    return dict(entrez_to_go)


def build_symbol_to_go(cache_dir: Path = DEFAULT_CACHE_DIR,
                       progress_callback: Optional[Callable] = None,
                       ) -> Dict[str, Set[str]]:
    """Build {HGNC_symbol: set(GO_IDs)} by chaining HGNC->Entrez->gene2go.

    Uses pre-shipped data files in kosmic/reference/go/. Raises FileNotFoundError
    if the files are missing (run update_data() to download them).
    """

    if progress_callback:
        progress_callback("Loading HGNC-Entrez mapping...")
    hgnc_map = load_hgnc_entrez_map(cache_dir)

    if progress_callback:
        progress_callback("Loading gene2go annotations (human)...")
    entrez_to_go = load_gene2go(cache_dir)

    if progress_callback:
        progress_callback("Building symbol-to-GO associations...")

    symbol_to_go: Dict[str, Set[str]] = {}
    n_mapped = 0
    for symbol, entrez_id in hgnc_map.items():
        go_terms = entrez_to_go.get(entrez_id)
        if go_terms:
            symbol_to_go[symbol] = go_terms
            n_mapped += 1

    if progress_callback:
        progress_callback(
            f"GO associations: {n_mapped:,} genes mapped, "
            f"{len(hgnc_map) - n_mapped:,} without GO annotations"
        )

    return symbol_to_go


# Elim algorithm


def _run_fisher(study_count: int, study_n: int,
                pop_count: int, pop_n: int) -> float:
    """Fisher's exact test (one-sided, enrichment)."""
    from scipy.stats import fisher_exact

    # 2x2 contingency table:
    #                  In term    Not in term
    # Study:           a          b
    # Background:      c          d
    a = study_count
    b = study_n - study_count
    c = pop_count - study_count
    d = (pop_n - study_n) - c

    # Clamp negatives (can happen with propagation edge cases)
    a, b, c, d = max(a, 0), max(b, 0), max(c, 0), max(d, 0)

    _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
    return pval


def run_elim(
    study_genes: Set[str],
    background_genes: Set[str],
    symbol_to_go: Dict[str, Set[str]],
    godag,
    ontology: str = "BP",
    alpha: float = DEFAULT_FDR,
    min_genes: int = 5,
    progress_callback: Optional[Callable] = None,
) -> pd.DataFrame:
    """Run the elim enrichment algorithm.

    The elim algorithm (Alexa et al. 2006, as implemented in R's topGO)
    processes GO terms bottom-up (most specific first). When a term is
    significant, its genes are removed from all ancestor terms before
    testing them, reducing redundancy in the results.

    Parameters
    ----------
    study_genes : set of str
        Significant DE genes (HGNC symbols).
    background_genes : set of str
        All tested genes (HGNC symbols).
    symbol_to_go : dict
        {HGNC_symbol: set(GO_IDs)} from build_symbol_to_go().
    godag : GODag
        Loaded GO DAG from goatools.
    ontology : str
        'BP', 'MF', or 'CC'.
    alpha : float
        Significance threshold for the elim step.
    min_genes : int
        Minimum genes annotated to a term to test it.
    progress_callback : callable, optional
        Called with progress messages.

    Returns
    -------
    pd.DataFrame
        Columns: GO_ID, Term, Namespace, P_value, Study_Count, Study_Total,
        Background_Count, Background_Total, Fold_Enrichment, Genes.
        Sorted by P_value ascending.
    """
    ns_map = {"BP": "biological_process", "MF": "molecular_function",
              "CC": "cellular_component"}
    target_ns = ns_map.get(ontology, ontology)

    if progress_callback:
        progress_callback(f"Building GO term associations ({ontology})...")

    # Build term -> genes mappings (propagate annotations up the DAG)
    # For each gene, all its direct GO terms + all ancestor terms get the gene
    term_to_bg_genes: Dict[str, Set[str]] = defaultdict(set)
    term_to_study_genes: Dict[str, Set[str]] = defaultdict(set)

    for gene in background_genes:
        go_terms = symbol_to_go.get(gene, set())
        for go_id in go_terms:
            if go_id not in godag:
                continue
            term = godag[go_id]
            if term.namespace != target_ns:
                continue
            # Direct annotation
            term_to_bg_genes[go_id].add(gene)
            # Propagate to all ancestors
            for ancestor_id in term.get_all_parents():
                if ancestor_id in godag and godag[ancestor_id].namespace == target_ns:
                    term_to_bg_genes[ancestor_id].add(gene)

    for gene in study_genes:
        go_terms = symbol_to_go.get(gene, set())
        for go_id in go_terms:
            if go_id not in godag:
                continue
            term = godag[go_id]
            if term.namespace != target_ns:
                continue
            term_to_study_genes[go_id].add(gene)
            for ancestor_id in term.get_all_parents():
                if ancestor_id in godag and godag[ancestor_id].namespace == target_ns:
                    term_to_study_genes[ancestor_id].add(gene)

    # Get all testable terms (enough genes in background)
    testable_terms = {
        go_id for go_id, genes in term_to_bg_genes.items()
        if len(genes) >= min_genes
    }

    if progress_callback:
        progress_callback(
            f"Testing {len(testable_terms):,} GO terms ({ontology}) "
            f"with elim algorithm..."
        )

    # Sort terms by depth (deepest first = most specific first)
    term_depths = {}
    for go_id in testable_terms:
        if go_id in godag:
            term_depths[go_id] = godag[go_id].depth
    sorted_terms = sorted(testable_terms, key=lambda t: term_depths.get(t, 0),
                          reverse=True)

    # Elim: process bottom-up, removing significant genes from ancestors
    study_n = len(study_genes & background_genes)
    pop_n = len(background_genes)
    eliminated_genes: Set[str] = set()
    results = []

    for i, go_id in enumerate(sorted_terms):
        if progress_callback and i % 500 == 0:
            progress_callback(
                f"Elim: testing term {i + 1:,}/{len(sorted_terms):,}..."
            )

        # Current study genes for this term, minus eliminated genes
        current_study = term_to_study_genes.get(go_id, set()) - eliminated_genes
        current_bg = term_to_bg_genes.get(go_id, set())

        study_count = len(current_study & study_genes)
        pop_count = len(current_bg)

        if study_count == 0:
            continue

        pval = _run_fisher(study_count, study_n, pop_count, pop_n)

        # Fold enrichment
        expected = (pop_count / pop_n) * study_n if pop_n > 0 else 0
        fold_enrichment = study_count / expected if expected > 0 else 0.0

        term_obj = godag[go_id]
        # Genes column: the DE genes that drove this term's significance
        # (after elim removal — for display and understanding)
        # All_Genes column: ALL genes annotated to this term in the background
        # (for pathway scoring — use the full biology, not the elim-reduced set)
        full_term_genes = term_to_bg_genes.get(go_id, set())
        results.append({
            "GO_ID": go_id,
            "Term": term_obj.name,
            "Namespace": ontology,
            "P_value": pval,
            "Study_Count": study_count,
            "Study_Total": study_n,
            "Background_Count": pop_count,
            "Background_Total": pop_n,
            "Fold_Enrichment": fold_enrichment,
            "Genes": sorted(current_study & study_genes),
            "All_Genes": sorted(full_term_genes),
        })

        # Elim step: if significant, remove these genes from ancestor terms
        if pval < alpha:
            eliminated_genes |= current_study

    if not results:
        return pd.DataFrame(columns=[
            "GO_ID", "Term", "Namespace", "P_value", "Study_Count",
            "Study_Total", "Background_Count", "Background_Total",
            "Fold_Enrichment", "Genes", "All_Genes",
        ])

    df = pd.DataFrame(results).sort_values("P_value").reset_index(drop=True)
    return df


# High-level API


def run_enrichment(
    study_genes: List[str],
    background_genes: List[str],
    ontology: str = "BP",
    method: str = "elim",
    alpha: float = DEFAULT_FDR,
    min_genes: int = 5,
    cache_dir: Optional[Path] = None,
    progress_callback: Optional[Callable] = None,
) -> pd.DataFrame:
    """Run GO enrichment analysis on HGNC gene symbols.

    Parameters
    ----------
    study_genes : list of str
        Significant DE genes (HGNC symbols).
    background_genes : list of str
        All tested genes (HGNC symbols). Should be the full set of genes
        that were tested for DE (not just the significant ones).
    ontology : str
        'BP' (Biological Process), 'MF' (Molecular Function), or
        'CC' (Cellular Component). Default: 'BP'.
    method : str
        'elim' (default) or 'fisher' (standard Fisher + BH correction).
    alpha : float
        Significance threshold. For elim, controls which terms trigger
        gene removal. For fisher, used for BH correction.
    min_genes : int
        Minimum genes annotated to a term to include it.
    cache_dir : Path, optional
        Directory for cached data files. Default: kosmic/reference/go/.
    progress_callback : callable, optional
        Called with (message: str) for progress updates.

    Returns
    -------
    pd.DataFrame
        Columns: GO_ID, Term, Namespace, P_value, Study_Count, Study_Total,
        Background_Count, Background_Total, Fold_Enrichment, Genes.
        Sorted by P_value ascending.
    """
    from goatools.obo_parser import GODag

    if cache_dir is None:
        cache_dir = DEFAULT_CACHE_DIR

    # Load pre-shipped data files
    obo_path = cache_dir / "go-basic.obo"
    if not obo_path.exists():
        raise FileNotFoundError(
            f"GO ontology not found at {obo_path}. "
            "Run update_data() to download GO reference files."
        )

    if progress_callback:
        progress_callback("Loading GO DAG...")
    godag = GODag(str(obo_path))

    # Build associations
    symbol_to_go = build_symbol_to_go(cache_dir, progress_callback)

    # Convert to sets
    study_set = set(study_genes) & set(symbol_to_go.keys())
    bg_set = set(background_genes) & set(symbol_to_go.keys())

    if not study_set:
        if progress_callback:
            progress_callback("No study genes mapped to GO annotations")
        return pd.DataFrame(columns=[
            "GO_ID", "Term", "Namespace", "P_value", "Study_Count",
            "Study_Total", "Background_Count", "Background_Total",
            "Fold_Enrichment", "Genes",
        ])

    n_study_mapped = len(study_set)
    n_study_total = len(study_genes)
    n_bg_mapped = len(bg_set)
    n_bg_total = len(background_genes)

    if progress_callback:
        progress_callback(
            f"Study genes: {n_study_mapped:,}/{n_study_total:,} mapped to GO. "
            f"Background: {n_bg_mapped:,}/{n_bg_total:,} mapped."
        )

    if method == "elim":
        df = run_elim(
            study_set, bg_set, symbol_to_go, godag,
            ontology=ontology, alpha=alpha, min_genes=min_genes,
            progress_callback=progress_callback,
        )
    elif method == "fisher":
        # Standard Fisher + BH correction (no DAG-aware filtering)
        df = run_elim(
            study_set, bg_set, symbol_to_go, godag,
            ontology=ontology,
            alpha=0.0,  # alpha=0 means no genes are eliminated
            min_genes=min_genes,
            progress_callback=progress_callback,
        )
        # Apply BH correction
        if not df.empty:
            from kosmic.numerical import bh_fdr
            df["P_adjusted"] = bh_fdr(df["P_value"].values)
    else:
        raise ValueError(f"Unknown method: {method!r}. Use 'elim' or 'fisher'.")

    if progress_callback:
        n_sig = (df["P_value"] < alpha).sum() if not df.empty else 0
        progress_callback(
            f"GO enrichment complete: {len(df):,} terms tested, "
            f"{n_sig:,} significant at alpha={alpha}"
        )

    return df

