"""
GWAS Overlap (Naive Nearest-Gene)
=================================
Annotates a consensus gene list with the closest GWAS-significant
variants from a user-supplied summary statistics file (IEU OpenGWAS
VCF, GWAS Catalog GWAS-SSF v1.0, or METAL tab-separated).

Pipeline:
    1. Stream-parse the GWAS file, retaining only variants with
       -log10(p) >= LP_THRESHOLD (default 7.30, ~5e-8).
    2. Map each significant variant to its nearest HGNC gene whose
       gene body extended by +/- WINDOW_KB lies within reach.
    3. For each consensus gene, report the nearest GWAS-significant
       variant (rsID, p-value, distance) within window.
    4. Set-level enrichment tests, all gene-length-aware:
         fisher_enrichment             -- naive, **BIASED**
         length_matched_fisher_enrichment -- Mantel-Haenszel stratified
         length_matched_permutation_test  -- non-parametric empirical

Set-level enrichment limitations:
    The naive Fisher test systematically overstates enrichment because
    consensus genes are typically longer than background genes (more
    SNPs per gene, more chances to have a nearby hit). Length-matched
    Mantel-Haenszel removes most of this bias but does NOT fully
    control for locus density or LD structure. Empirical floor: a
    height GWAS (zero biological connection to an EC consensus)
    produces a residual MH OR ~1.3 with p ~1e-5; CAD and HF GWAS on
    the same consensus produce MH ORs in the 1.16-1.25 range,
    indistinguishable from the height "null" floor. The enrichment
    tests in this module cannot reliably distinguish trait-specific
    biology from methodological bias on well-powered GWAS.

    For a publication-grade set-level enrichment claim, use MAGMA
    (gene-length + gene-density + LD corrections + proper V2G). The
    per-gene annotation outputs here remain useful (which consensus
    genes overlap known GWAS loci); the set-level p-values should
    not be quoted as validation.

Pure logic, no PyQt6 imports.
"""

from __future__ import annotations

import bisect
import contextlib
import gzip
from importlib.resources import files
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd


from kosmic import GWAS_LP_THRESHOLD, GWAS_WINDOW_KB as DEFAULT_WINDOW_KB

DEFAULT_GENE_POSITIONS = files("kosmic.reference.gwas") / "gene_positions_GRCh37.tsv"


# ----------------------------------------------------------------------
# VCF parsing
# ----------------------------------------------------------------------

def parse_gwas(path: str | Path,
               lp_threshold: float = GWAS_LP_THRESHOLD,
               progress_cb=None) -> pd.DataFrame:
    """Auto-detect GWAS sumstats format and parse.

    Supported formats:
      * IEU OpenGWAS VCF (file starts with ``##fileformat=VCF...``)
      * GWAS Catalog GWAS-SSF v1.0 tab-separated (header has
        ``p_value`` or ``hm_chromosome``)

    Returns a DataFrame in the same shape as `parse_gwas_vcf` /
    `parse_gwas_ssf`: columns chrom (str), pos (int), rsid (str),
    ref, alt, lp (float), beta (float), se (float).
    """
    path = Path(path)
    is_gz = str(path).endswith(".gz")
    opener = (lambda p: gzip.open(p, "rt", encoding="utf-8")) if is_gz \
        else (lambda p: open(p, "rt", encoding="utf-8"))

    with opener(path) as f:
        first = f.readline()
    first = first.lstrip("\ufeff")  # handle BOM if present

    if first.startswith("##fileformat=VCF") or first.startswith("##"):
        return parse_gwas_vcf(path, lp_threshold=lp_threshold,
                              progress_cb=progress_cb)

    # TSV header: route to GWAS-SSF parser. Sniff for any of the
    # common GWAS column-name variants -- canonical GWAS-SSF, METAL,
    # CARDIoGRAM, plain "CHR/BP/P" style.
    first_lower = first.lower()
    if "\t" in first and any(tok in first_lower for tok in (
            "p_value", "p-value", "pvalue", "\tp\t", "\tp\n",
            "hm_chromosome", "chromosome", "\tchr\t", "\tbp\t",
            "markername", "rs_number")):
        return parse_gwas_ssf(path, lp_threshold=lp_threshold,
                              progress_cb=progress_cb)

    raise ValueError(
        f"Could not auto-detect GWAS sumstats format for {path}.\n"
        f"Header looked like: {first[:200]!r}\n"
        "Supported: IEU OpenGWAS VCF (## headers) or GWAS Catalog "
        "GWAS-SSF v1.0 (TSV with p_value column).")


def parse_gwas_vcf(path: str | Path,
                   lp_threshold: float = GWAS_LP_THRESHOLD,
                   progress_cb=None) -> pd.DataFrame:
    """Stream-parse an IEU OpenGWAS-format VCF and return only variants
    passing the genome-wide significance threshold.

    The IEU format encodes per-sample (per-study) statistics in the
    FORMAT/sample columns. Per row we extract:

        chrom, pos, rsid, ref, alt, ES (effect size), SE, LP

    Parameters
    ----------
    path : str or Path
        VCF (optionally gzipped).
    lp_threshold : float
        Minimum -log10(p) to retain. Default 7.30103 = -log10(5e-8).
    progress_cb : callable(str), optional
        Called with status messages.

    Returns
    -------
    DataFrame with columns: chrom (str), pos (int), rsid (str), ref,
    alt, lp (float), beta (float, may be NaN), se (float, may be NaN).
    """
    path = Path(path)
    is_gz = str(path).endswith(".gz")
    opener = (lambda p: gzip.open(p, "rt", encoding="utf-8")) if is_gz \
        else (lambda p: open(p, "rt", encoding="utf-8"))

    chroms, poss, rsids, refs, alts = [], [], [], [], []
    lps, betas, ses = [], [], []

    line_n = 0
    n_total = 0
    n_kept = 0
    with opener(path) as f:
        for line in f:
            line_n += 1
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                # Header line -- could read sample names but not needed
                continue
            parts = line.rstrip("\n").split("\t")
            # VCF: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE...
            if len(parts) < 10:
                continue
            n_total += 1

            fmt = parts[8].split(":")
            sample_vals = parts[9].split(":")
            # Build a dict from FORMAT keys to sample values
            sd = dict(zip(fmt, sample_vals))

            lp_raw = sd.get("LP")
            if lp_raw is None or lp_raw in (".", ""):
                continue
            try:
                lp = float(lp_raw)
            except ValueError:
                continue
            if not np.isfinite(lp) or lp < lp_threshold:
                continue

            es_raw = sd.get("ES", ".")
            se_raw = sd.get("SE", ".")
            try:
                beta = float(es_raw) if es_raw not in (".", "") else float("nan")
            except ValueError:
                beta = float("nan")
            try:
                se = float(se_raw) if se_raw not in (".", "") else float("nan")
            except ValueError:
                se = float("nan")

            try:
                pos = int(parts[1])
            except ValueError:
                continue

            chroms.append(parts[0].replace("chr", ""))
            poss.append(pos)
            rsids.append(parts[2])
            refs.append(parts[3])
            alts.append(parts[4])
            lps.append(lp)
            betas.append(beta)
            ses.append(se)
            n_kept += 1

            if progress_cb is not None and (n_total % 500_000 == 0):
                progress_cb(
                    f"  parsed {n_total/1e6:.1f}M variants, "
                    f"{n_kept} above threshold...")

    if progress_cb is not None:
        progress_cb(
            f"VCF parsed: {n_total:,} variants total, "
            f"{n_kept:,} retained (LP >= {lp_threshold:.2f}).")

    return pd.DataFrame({
        "chrom": chroms,
        "pos":   poss,
        "rsid":  rsids,
        "ref":   refs,
        "alt":   alts,
        "lp":    lps,
        "beta":  betas,
        "se":    ses,
    })


def parse_gwas_ssf(path: str | Path,
                   lp_threshold: float = GWAS_LP_THRESHOLD,
                   progress_cb=None) -> pd.DataFrame:
    """Stream-parse a GWAS Catalog GWAS-SSF v1.0 tab-separated file
    and return only variants passing the genome-wide significance
    threshold.

    GWAS-SSF columns of interest (harmonised version uses ``hm_*``
    aliases for the canonical orientation):

        hm_chromosome    | chromosome
        hm_base_pair_location | base_pair_location
        hm_rsid          | rsid | variant_id
        hm_other_allele  | other_allele
        hm_effect_allele | effect_allele
        hm_beta          | beta
        standard_error
        p_value

    The parser prefers ``hm_*`` columns when present, falling back
    to the un-prefixed names. P-values in the file are RAW (not
    -log10), so we convert: ``lp = -log10(p)``.

    Returns the same DataFrame shape as `parse_gwas_vcf`.
    """
    path = Path(path)
    is_gz = str(path).endswith(".gz")
    opener = (lambda p: gzip.open(p, "rt", encoding="utf-8")) if is_gz \
        else (lambda p: open(p, "rt", encoding="utf-8"))

    p_threshold = 10 ** (-lp_threshold)

    chroms, poss, rsids, refs, alts = [], [], [], [], []
    lps, betas, ses = [], [], []

    n_total = 0
    n_kept = 0
    with opener(path) as f:
        header = f.readline().rstrip("\n").lstrip("\ufeff")
        cols = [c.strip() for c in header.split("\t")]
        # Resolve canonical column names (hm_* preferred)
        def _idx(*candidates):
            for c in candidates:
                if c in cols:
                    return cols.index(c)
            return None

        # Column lookups try canonical GWAS-SSF names first, then
        # common METAL / CARDIoGRAM / IEU variants. Case-sensitive.
        i_chrom = _idx("hm_chromosome", "chromosome", "chr",
                       "CHR", "Chromosome", "Chr")
        i_pos   = _idx("hm_base_pair_location", "base_pair_location",
                       "bp", "pos", "BP", "Position", "POS",
                       "base_pair", "position")
        i_rsid  = _idx("hm_rsid", "rsid", "variant_id", "snp", "snp_id",
                       "SNP", "rs_number", "RSID", "MarkerName",
                       "rs_id", "ID")
        i_ref   = _idx("hm_other_allele",  "other_allele",  "ref",
                       "Allele2", "non_effect_allele",
                       "reference_allele", "REF")
        i_alt   = _idx("hm_effect_allele", "effect_allele", "alt",
                       "Allele1", "ALT")
        i_beta  = _idx("hm_beta",  "beta", "Effect", "BETA",
                       "log_odds", "logOR")
        i_se    = _idx("standard_error", "se", "StdErr",
                       "standard_error_logOR", "SE")
        i_p     = _idx("p_value", "pvalue", "pval", "p",
                       "P-value", "Pvalue", "P.value", "P", "PVAL",
                       "P_value", "P_VAL")

        if i_chrom is None or i_pos is None or i_p is None:
            raise ValueError(
                f"GWAS-SSF parser: missing required column(s). "
                f"Need chromosome / base_pair_location / p_value. "
                f"Got header: {cols[:20]!r}...")

        for line in f:
            n_total += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= i_p:
                continue

            p_raw = parts[i_p]
            if not p_raw or p_raw in (".", "NA", "nan"):
                continue
            try:
                p = float(p_raw)
            except ValueError:
                continue
            if not np.isfinite(p) or p <= 0 or p > p_threshold:
                continue
            lp = -np.log10(p)

            try:
                pos = int(parts[i_pos])
            except (ValueError, IndexError):
                continue

            chrom = str(parts[i_chrom]).replace("chr", "")
            rsid  = parts[i_rsid] if i_rsid is not None and i_rsid < len(parts) else ""
            ref   = parts[i_ref]  if i_ref  is not None and i_ref  < len(parts) else ""
            alt   = parts[i_alt]  if i_alt  is not None and i_alt  < len(parts) else ""

            beta = float("nan")
            if i_beta is not None and i_beta < len(parts):
                with contextlib.suppress(ValueError):
                    beta = float(parts[i_beta])
            se = float("nan")
            if i_se is not None and i_se < len(parts):
                with contextlib.suppress(ValueError):
                    se = float(parts[i_se])

            chroms.append(chrom)
            poss.append(pos)
            rsids.append(rsid)
            refs.append(ref)
            alts.append(alt)
            lps.append(lp)
            betas.append(beta)
            ses.append(se)
            n_kept += 1

            if progress_cb is not None and (n_total % 1_000_000 == 0):
                progress_cb(
                    f"  parsed {n_total/1e6:.1f}M variants, "
                    f"{n_kept} above threshold...")

    if progress_cb is not None:
        progress_cb(
            f"GWAS-SSF parsed: {n_total:,} variants total, "
            f"{n_kept:,} retained (LP >= {lp_threshold:.2f}).")

    return pd.DataFrame({
        "chrom": chroms,
        "pos":   poss,
        "rsid":  rsids,
        "ref":   refs,
        "alt":   alts,
        "lp":    lps,
        "beta":  betas,
        "se":    ses,
    })


# ----------------------------------------------------------------------
# Gene positions
# ----------------------------------------------------------------------

def load_gene_positions(path: Optional[Path] = None) -> pd.DataFrame:
    """Load the bundled HGNC gene positions (GRCh37).

    Returns a DataFrame with columns: gene, chrom, start, end, strand.
    """
    p = Path(path) if path is not None else DEFAULT_GENE_POSITIONS
    if not p.exists():
        raise FileNotFoundError(
            f"Gene-position reference not found: {p}.\n"
            "Run scripts/generate_gene_positions.py to regenerate it.")
    df = pd.read_csv(p, sep="\t", dtype={"chrom": str})
    df["chrom"] = df["chrom"].astype(str).str.replace("chr", "", regex=False)
    return df


# ----------------------------------------------------------------------
# Variant -> gene mapping
# ----------------------------------------------------------------------

def _per_chrom_gene_index(genes_df: pd.DataFrame) -> dict:
    """Build a per-chromosome sorted index of (start, end, gene_idx)
    for fast nearest-gene lookups via bisect.

    Returns ``{chrom: (np.ndarray starts, np.ndarray ends, np.ndarray gene_indices)}``.
    """
    out = {}
    for chrom, grp in genes_df.groupby("chrom"):
        grp_sorted = grp.sort_values("start").reset_index(drop=False)
        # Original index in genes_df is preserved in 'index' col
        starts = grp_sorted["start"].to_numpy(dtype=np.int64)
        ends   = grp_sorted["end"].to_numpy(dtype=np.int64)
        idxs   = grp_sorted["index"].to_numpy(dtype=np.int64)
        out[str(chrom)] = (starts, ends, idxs)
    return out


def map_variants_to_nearest_gene(variants: pd.DataFrame,
                                 genes_df: pd.DataFrame,
                                 window_kb: int = DEFAULT_WINDOW_KB
                                 ) -> pd.DataFrame:
    """For each significant variant, find the nearest gene whose body
    sits within +/- window_kb of the variant.

    Parameters
    ----------
    variants : DataFrame
        Output of `parse_gwas_vcf`.
    genes_df : DataFrame
        Output of `load_gene_positions`.
    window_kb : int
        Maximum distance (kb) from variant to gene body.

    Returns
    -------
    DataFrame with columns:
        chrom, pos, rsid, lp, beta, se, gene, distance_bp
    Rows with no nearby gene are dropped. Distance is 0 when the
    variant lies inside the gene body.
    """
    window_bp = window_kb * 1_000
    index = _per_chrom_gene_index(genes_df)

    rows = []
    for _, v in variants.iterrows():
        chrom = str(v["chrom"]).replace("chr", "")
        if chrom not in index:
            continue
        starts, ends, idxs = index[chrom]
        pos = int(v["pos"])

        # Find candidate genes: those whose [start - window, end + window]
        # contains pos. Use bisect_right on starts to bracket the search.
        # Candidates = genes whose start <= pos + window AND end >= pos - window.
        right = bisect.bisect_right(starts, pos + window_bp)
        # We need genes ending at or after pos - window. With sorted
        # starts that's all i in [0, right) where ends[i] >= pos - window.
        if right == 0:
            continue
        cand_mask = ends[:right] >= (pos - window_bp)
        if not cand_mask.any():
            continue
        cand_idx = idxs[:right][cand_mask]
        cand_starts = starts[:right][cand_mask]
        cand_ends   = ends[:right][cand_mask]

        # Distance: 0 if pos is inside [start, end], else distance to
        # nearest edge.
        distances = np.where(
            (pos >= cand_starts) & (pos <= cand_ends), 0,
            np.minimum(np.abs(pos - cand_starts),
                       np.abs(pos - cand_ends)))
        best_i = int(np.argmin(distances))
        best_gene_row = genes_df.iloc[int(cand_idx[best_i])]
        rows.append({
            "chrom":       chrom,
            "pos":         pos,
            "rsid":        v["rsid"],
            "lp":          float(v["lp"]),
            "beta":        float(v["beta"]) if not pd.isna(v["beta"])
                           else float("nan"),
            "se":          float(v["se"]) if not pd.isna(v["se"])
                           else float("nan"),
            "gene":        str(best_gene_row["gene"]),
            "distance_bp": int(distances[best_i]),
        })

    if not rows:
        return pd.DataFrame(columns=[
            "chrom", "pos", "rsid", "lp", "beta", "se", "gene", "distance_bp"
        ])
    return pd.DataFrame(rows)


def gene_to_best_variant(variant_gene_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse a variant-to-gene table down to one row per gene
    (the most significant variant for each gene wins).

    Returns
    -------
    DataFrame indexed by gene, columns:
        rsid, chrom, pos, lp, beta, se, distance_bp, n_hits_in_window
    """
    if variant_gene_df is None or len(variant_gene_df) == 0:
        return pd.DataFrame(columns=[
            "rsid", "chrom", "pos", "lp", "beta", "se",
            "distance_bp", "n_hits_in_window"
        ])
    grp = variant_gene_df.sort_values("lp", ascending=False)
    # n_hits per gene
    counts = grp.groupby("gene").size().rename("n_hits_in_window")
    # Best variant per gene = first row (highest lp)
    best = grp.drop_duplicates(subset="gene", keep="first").set_index("gene")
    out = best[["rsid", "chrom", "pos", "lp", "beta", "se",
                "distance_bp"]].copy()
    out["n_hits_in_window"] = counts
    return out


# ----------------------------------------------------------------------
# Annotation + enrichment
# ----------------------------------------------------------------------

def annotate_consensus(consensus_df: pd.DataFrame,
                       gene_to_best: pd.DataFrame,
                       gene_col: str = "names") -> pd.DataFrame:
    """Add GWAS-overlap columns to a consensus DataFrame.

    Adds: ``gwas_hit_in_window`` (bool), ``gwas_nearest_rsid`` (str),
    ``gwas_nearest_lp`` (float), ``gwas_nearest_distance_kb`` (float),
    ``gwas_n_hits_in_window`` (int).
    """
    if gene_col not in consensus_df.columns:
        raise KeyError(
            f"gene_col '{gene_col}' not in consensus_df columns")

    out = consensus_df.copy()
    genes = out[gene_col].astype(str).fillna("").tolist()

    # Defensive: if gene_to_best is empty, set all columns to defaults
    if gene_to_best is None or len(gene_to_best) == 0:
        out["gwas_hit_in_window"]        = False
        out["gwas_nearest_rsid"]         = ""
        out["gwas_nearest_lp"]           = float("nan")
        out["gwas_nearest_distance_kb"]  = float("nan")
        out["gwas_n_hits_in_window"]     = 0
        return out

    in_window = []
    rsids = []
    lps = []
    dists_kb = []
    n_hits = []
    for g in genes:
        if g in gene_to_best.index:
            row = gene_to_best.loc[g]
            in_window.append(True)
            rsids.append(str(row["rsid"]))
            lps.append(float(row["lp"]))
            dists_kb.append(float(row["distance_bp"]) / 1000.0)
            n_hits.append(int(row["n_hits_in_window"]))
        else:
            in_window.append(False)
            rsids.append("")
            lps.append(float("nan"))
            dists_kb.append(float("nan"))
            n_hits.append(0)

    out["gwas_hit_in_window"]       = in_window
    out["gwas_nearest_rsid"]        = rsids
    out["gwas_nearest_lp"]          = lps
    out["gwas_nearest_distance_kb"] = dists_kb
    out["gwas_n_hits_in_window"]    = n_hits
    return out


def compute_length_diagnostic(consensus_genes: Iterable[str],
                              background_genes: Iterable[str],
                              genes_with_gwas_hit: Iterable[str],
                              gene_positions_df: pd.DataFrame) -> dict:
    """Quick gene-length bias diagnostic.

    Nearest-gene + Fisher's exact is known to be inflated when
    consensus genes are systematically longer than background genes
    (longer genes have more SNPs and are more likely to be near any
    GWAS hit). This function reports mean / median gene lengths for
    the consensus, background, and GWAS-hit gene sets so the bias
    magnitude can be eyeballed directly.

    Returns dict with:
        n_consensus_with_length, n_background_with_length,
        mean_length_consensus, mean_length_background,
        median_length_consensus, median_length_background,
        mean_length_hit_genes,
        length_ratio_consensus_over_background
            -- >1.5 implies material bias; >2 implies severe bias
    """
    gene_len = {
        str(row['gene']): int(row['end']) - int(row['start'])
        for _, row in gene_positions_df.iterrows()
    }

    consensus_set  = {str(g) for g in consensus_genes if g}
    all_set        = {str(g) for g in background_genes if g} | consensus_set
    bg_only        = all_set - consensus_set
    hit_set        = {str(g) for g in genes_with_gwas_hit if g}

    def _lens(gs):
        arr = np.array([gene_len[g] for g in gs if g in gene_len],
                       dtype=float)
        return arr if len(arr) else np.array([], dtype=float)

    cons_lens = _lens(consensus_set)
    bg_lens   = _lens(bg_only)
    hit_lens  = _lens(hit_set & all_set)

    mean_c = float(np.mean(cons_lens)) if len(cons_lens) else 0.0
    mean_b = float(np.mean(bg_lens))   if len(bg_lens)   else 0.0
    ratio = mean_c / mean_b if mean_b > 0 else float('nan')

    return {
        'n_consensus_with_length':  int(len(cons_lens)),
        'n_background_with_length': int(len(bg_lens)),
        'mean_length_consensus':    mean_c,
        'mean_length_background':   mean_b,
        'median_length_consensus':  (float(np.median(cons_lens))
                                     if len(cons_lens) else 0.0),
        'median_length_background': (float(np.median(bg_lens))
                                     if len(bg_lens) else 0.0),
        'mean_length_hit_genes':    (float(np.mean(hit_lens))
                                     if len(hit_lens) else 0.0),
        'length_ratio_consensus_over_background': float(ratio),
    }


def length_matched_fisher_enrichment(consensus_genes: Iterable[str],
                                     background_genes: Iterable[str],
                                     genes_with_gwas_hit: Iterable[str],
                                     gene_positions_df: pd.DataFrame,
                                     n_bins: int = 10) -> dict:
    """Gene-length-stratified enrichment test using Mantel-Haenszel.

    Bins all genes into length deciles (based on the combined
    consensus + background set). For each decile, computes a 2x2
    contingency table (consensus x has_hit). Pools via Mantel-Haenszel
    to return a length-corrected common odds ratio and its p-value.

    Unlike the naive `fisher_enrichment`, this test is robust to
    consensus genes being systematically longer than background.
    If the naive OR is 2.0 and the length-matched OR drops to 1.1
    with p > 0.1, the naive enrichment was gene-length bias.

    Parameters
    ----------
    consensus_genes, background_genes, genes_with_gwas_hit : iterable of str
    gene_positions_df : DataFrame
        Output of `load_gene_positions`; must have gene, start, end cols.
    n_bins : int
        Number of length deciles (default 10). Strata with zero
        row or column sums are dropped.

    Returns
    -------
    dict with:
        odds_ratio    : Mantel-Haenszel pooled OR
        pvalue        : MH test of null OR = 1 (two-sided by default;
                        we report the raw MH p-value)
        n_strata_used : number of deciles with non-degenerate tables
        n_consensus   : consensus genes with known length
        n_background  : background-only genes with known length
        n_consensus_with_hit : a-cell sum across strata
        strata        : list of 2x2 tables (for diagnostics)
        method        : human-readable string
    """
    from statsmodels.stats.contingency_tables import StratifiedTable

    gene_len = {
        str(row['gene']): int(row['end']) - int(row['start'])
        for _, row in gene_positions_df.iterrows()
    }

    consensus_set  = {str(g) for g in consensus_genes if g}
    all_set        = {str(g) for g in background_genes if g} | consensus_set
    bg_only        = all_set - consensus_set
    hit_set        = {str(g) for g in genes_with_gwas_hit if g}

    # Restrict to genes with known position (i.e. length known)
    all_with_len = [g for g in all_set if g in gene_len]
    if not all_with_len:
        return {
            'odds_ratio':           float('nan'),
            'pvalue':               float('nan'),
            'n_strata_used':        0,
            'n_consensus':          0,
            'n_background':         0,
            'n_consensus_with_hit': 0,
            'strata':               [],
            'method': 'mantel_haenszel (no genes with known length)',
        }

    lengths = np.array([gene_len[g] for g in all_with_len], dtype=float)

    # Decile edges from combined gene set. Handle ties by taking unique.
    edges = np.quantile(lengths, np.linspace(0, 1, n_bins + 1))
    edges = np.unique(edges)
    n_effective_bins = max(1, len(edges) - 1)

    def _bin(g):
        if g not in gene_len:
            return -1
        L = gene_len[g]
        b = int(np.searchsorted(edges, L, side='right') - 1)
        return max(0, min(n_effective_bins - 1, b))

    tables = []
    for b in range(n_effective_bins):
        a = b_cnt = c = d_cnt = 0
        for g in all_with_len:
            if _bin(g) != b:
                continue
            in_cons = g in consensus_set
            has_hit = g in hit_set
            if in_cons and has_hit:
                a += 1
            elif in_cons:
                b_cnt += 1
            elif has_hit:
                c += 1
            else:
                d_cnt += 1
        tables.append(np.array([[a, b_cnt], [c, d_cnt]]))

    # Keep strata with non-zero row + column sums so MH doesn't blow up
    valid = [t for t in tables
             if t[0].sum() > 0 and t[1].sum() > 0
             and t[:, 0].sum() > 0 and t[:, 1].sum() > 0]

    n_cons_len = sum(1 for g in consensus_set if g in gene_len)
    n_bg_len   = sum(1 for g in bg_only      if g in gene_len)
    a_total    = sum(int(t[0, 0]) for t in tables)

    if not valid:
        return {
            'odds_ratio':           float('nan'),
            'pvalue':               float('nan'),
            'n_strata_used':        0,
            'n_consensus':          n_cons_len,
            'n_background':         n_bg_len,
            'n_consensus_with_hit': a_total,
            'strata':               [t.tolist() for t in tables],
            'method':
                'mantel_haenszel (no usable strata; all degenerate)',
        }

    try:
        st = StratifiedTable(valid)
        or_mh = float(st.oddsratio_pooled)
        # two-sided MH test of null OR == 1
        res = st.test_null_odds(correction=False)
        p_mh = float(res.pvalue)
    except Exception as e:
        return {
            'odds_ratio':           float('nan'),
            'pvalue':               float('nan'),
            'n_strata_used':        len(valid),
            'n_consensus':          n_cons_len,
            'n_background':         n_bg_len,
            'n_consensus_with_hit': a_total,
            'strata':               [t.tolist() for t in tables],
            'method':               f'mantel_haenszel (error: {e})',
        }

    return {
        'odds_ratio':           or_mh,
        'pvalue':               p_mh,
        'n_strata_used':        len(valid),
        'n_consensus':          n_cons_len,
        'n_background':         n_bg_len,
        'n_consensus_with_hit': a_total,
        'strata':               [t.tolist() for t in tables],
        'method': f'mantel_haenszel (gene-length stratified, '
                  f'{len(valid)} deciles)',
    }


def length_matched_permutation_test(
        consensus_genes: Iterable[str],
        background_genes: Iterable[str],
        genes_with_gwas_hit: Iterable[str],
        gene_positions_df: pd.DataFrame,
        n_bins: int = 10,
        n_perms: int = 1000,
        seed: int = 42) -> dict:
    """Non-parametric gene-length-stratified permutation test.

    For each permutation, draws a null ``n_consensus_with_hit`` count
    by sampling from the hypergeometric distribution within each
    length decile (preserving per-decile consensus size and
    per-decile hit total). Sums across deciles and compares the null
    distribution of ``a`` (and derived OR) to the observed value.

    Conceptually equivalent to the Mantel-Haenszel test
    asymptotically, but:

    * provides an **empirical null distribution** that can be
      inspected directly rather than just a p-value
    * requires no parametric assumptions
    * gives a defensible p-value even in small-N strata
    * transparent for publication -- reviewers can see the null

    Fast: 1000 permutations with 10 strata run in well under a
    second (just hypergeometric draws, no gene-list reconstruction).

    Parameters
    ----------
    consensus_genes, background_genes, genes_with_gwas_hit : iterable of str
    gene_positions_df : DataFrame
    n_bins : int
        Length deciles (default 10).
    n_perms : int
        Permutation count (default 1000).
    seed : int
        RNG seed for reproducibility.

    Returns
    -------
    dict with:
        observed_n_consensus_with_hit : int
        observed_odds_ratio    : float (same as fisher_enrichment's OR)
        null_or_mean           : float (mean OR under null)
        null_or_median         : float
        null_or_ci_low / ci_high : float (2.5% and 97.5% quantiles)
        empirical_pvalue       : float (one-sided, observed >= null)
        n_perms_used           : int
        n_strata_used          : int
        null_or_samples        : list of float (sample of null ORs
                                 for inspection; capped at 500)
        method                 : str
    """
    gene_len = {
        str(row['gene']): int(row['end']) - int(row['start'])
        for _, row in gene_positions_df.iterrows()
    }

    consensus_set  = {str(g) for g in consensus_genes if g}
    all_set        = {str(g) for g in background_genes if g} | consensus_set
    hit_set        = {str(g) for g in genes_with_gwas_hit if g}

    all_with_len = [g for g in all_set if g in gene_len]
    if not all_with_len:
        return {
            'observed_n_consensus_with_hit': 0,
            'observed_odds_ratio': float('nan'),
            'null_or_mean':   float('nan'),
            'null_or_median': float('nan'),
            'null_or_ci_low': float('nan'),
            'null_or_ci_high': float('nan'),
            'empirical_pvalue': float('nan'),
            'n_perms_used':  0,
            'n_strata_used': 0,
            'null_or_samples': [],
            'method': 'permutation (no genes with known length)',
        }

    lengths = np.array([gene_len[g] for g in all_with_len], dtype=float)
    edges = np.unique(np.quantile(lengths, np.linspace(0, 1, n_bins + 1)))
    n_bins_eff = max(1, len(edges) - 1)

    def _bin(g):
        if g not in gene_len:
            return -1
        b = int(np.searchsorted(edges, gene_len[g], side='right') - 1)
        return max(0, min(n_bins_eff - 1, b))

    # Per-decile marginals
    strata = []  # list of (n_consensus_d, n_bg_only_d, n_hits_d, a_obs_d)
    for b in range(n_bins_eff):
        n_cons_d = n_bg_d = n_hit_d = a_obs_d = 0
        for g in all_with_len:
            if _bin(g) != b:
                continue
            in_cons = g in consensus_set
            has_hit = g in hit_set
            if in_cons:
                n_cons_d += 1
                if has_hit:
                    a_obs_d += 1
            else:
                n_bg_d += 1
            if has_hit:
                n_hit_d += 1
        strata.append((n_cons_d, n_bg_d, n_hit_d, a_obs_d))

    valid = [s for s in strata
             if s[0] + s[1] > 0 and s[2] > 0 and s[0] > 0 and s[1] > 0]

    if not valid:
        return {
            'observed_n_consensus_with_hit':
                sum(s[3] for s in strata),
            'observed_odds_ratio': float('nan'),
            'null_or_mean':   float('nan'),
            'null_or_median': float('nan'),
            'null_or_ci_low': float('nan'),
            'null_or_ci_high': float('nan'),
            'empirical_pvalue': float('nan'),
            'n_perms_used':  0,
            'n_strata_used': 0,
            'null_or_samples': [],
            'method': 'permutation (no usable strata)',
        }

    rng = np.random.default_rng(seed)

    # Observed totals
    a_obs = sum(s[3] for s in strata)
    n_cons_total = sum(s[0] for s in strata)
    n_bg_total   = sum(s[1] for s in strata)
    n_hits_total = sum(s[2] for s in strata)
    N_total = n_cons_total + n_bg_total

    def _or(a):
        """Fisher OR given per-total margins: a = consensus_hit,
        b = consensus_no_hit, c = background_hit, d = background_no_hit."""
        b_ = n_cons_total - a
        c_ = n_hits_total - a
        d_ = N_total - a - b_ - c_
        if b_ < 0 or c_ < 0 or d_ < 0:
            return float('nan')
        if b_ == 0 or c_ == 0:
            return float('inf')
        return float((a * d_) / (b_ * c_))

    or_observed = _or(a_obs)

    # Permutation: for each perm, sample a_null = sum over strata of
    # Hypergeometric(N_d = n_cons_d+n_bg_d, K_d = n_hit_d, n_d = n_cons_d)
    null_as = np.zeros(n_perms, dtype=np.int64)
    for p in range(n_perms):
        tot = 0
        for s in valid:
            n_cons_d, n_bg_d, n_hit_d, _ = s
            N_d = n_cons_d + n_bg_d
            # hypergeometric(ngood=n_hit_d, nbad=N_d - n_hit_d, nsample=n_cons_d)
            nbad_d = N_d - n_hit_d
            tot += int(rng.hypergeometric(n_hit_d, nbad_d, n_cons_d))
        null_as[p] = tot

    null_ors = np.array([_or(a) for a in null_as], dtype=float)
    finite = null_ors[np.isfinite(null_ors)]

    if len(finite) == 0:
        return {
            'observed_n_consensus_with_hit': a_obs,
            'observed_odds_ratio': or_observed,
            'null_or_mean':   float('nan'),
            'null_or_median': float('nan'),
            'null_or_ci_low': float('nan'),
            'null_or_ci_high': float('nan'),
            'empirical_pvalue': float('nan'),
            'n_perms_used':  n_perms,
            'n_strata_used': len(valid),
            'null_or_samples': [],
            'method': 'permutation (degenerate null)',
        }

    # One-sided empirical p: P(null OR >= observed OR)
    if np.isfinite(or_observed):
        emp_p = float((finite >= or_observed).sum() + 1) / (len(finite) + 1)
    else:
        emp_p = float('nan')

    ci_low, ci_high = np.quantile(finite, [0.025, 0.975])
    sample = null_ors.tolist()
    if len(sample) > 500:
        idx = rng.choice(len(sample), size=500, replace=False)
        sample = [sample[i] for i in sorted(idx)]

    return {
        'observed_n_consensus_with_hit': int(a_obs),
        'observed_odds_ratio':           float(or_observed),
        'null_or_mean':                  float(np.mean(finite)),
        'null_or_median':                float(np.median(finite)),
        'null_or_ci_low':                float(ci_low),
        'null_or_ci_high':               float(ci_high),
        'empirical_pvalue':              float(emp_p),
        'n_perms_used':                  int(n_perms),
        'n_strata_used':                 int(len(valid)),
        'null_or_samples':               sample,
        'method': f'permutation ({n_perms} draws, length-stratified, '
                  f'{len(valid)} deciles)',
    }


def fisher_enrichment(consensus_genes: Iterable[str],
                      background_genes: Iterable[str],
                      genes_with_gwas_hit: Iterable[str]) -> dict:
    """Fisher's exact enrichment of GWAS-hit-bearing genes in the
    consensus, against a background gene set.

    Parameters
    ----------
    consensus_genes : iterable of str
        Genes called significant by the consensus.
    background_genes : iterable of str
        All testable genes (typically the union of all studies'
        measured genes that survived gene filtering).
    genes_with_gwas_hit : iterable of str
        Genes that have at least one GWAS-significant variant within
        the chosen window (output of `gene_to_best_variant`).

    Returns
    -------
    dict with:
        ``n_consensus``, ``n_consensus_with_hit``,
        ``n_background``, ``n_background_with_hit``,
        ``odds_ratio``, ``pvalue``, ``contingency`` (2x2 list),
        ``method`` (string).

        Contingency table is::

            [[in_consensus AND has_hit, in_consensus AND no_hit],
             [in_background_only AND has_hit, in_background_only AND no_hit]]
    """
    from scipy.stats import fisher_exact

    consensus_set   = set(str(g) for g in consensus_genes if g)
    background_set  = set(str(g) for g in background_genes if g) | consensus_set
    hit_set         = set(str(g) for g in genes_with_gwas_hit if g)

    bg_only = background_set - consensus_set

    a = len(consensus_set & hit_set)
    b = len(consensus_set) - a
    c = len(bg_only & hit_set)
    d = len(bg_only) - c

    contingency = [[a, b], [c, d]]
    if a + b == 0 or c + d == 0:
        return {
            "n_consensus":            len(consensus_set),
            "n_consensus_with_hit":   a,
            "n_background":           len(bg_only),
            "n_background_with_hit":  c,
            "odds_ratio":             float("nan"),
            "pvalue":                 float("nan"),
            "contingency":            contingency,
            "method":                 "fisher_exact (degenerate)",
        }

    odds_ratio, pvalue = fisher_exact(contingency, alternative="greater")
    return {
        "n_consensus":            len(consensus_set),
        "n_consensus_with_hit":   a,
        "n_background":           len(bg_only),
        "n_background_with_hit":  c,
        "odds_ratio":             float(odds_ratio),
        "pvalue":                 float(pvalue),
        "contingency":            contingency,
        "method":                 "fisher_exact (one-sided, greater)",
    }
