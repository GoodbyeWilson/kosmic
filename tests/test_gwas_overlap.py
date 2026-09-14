"""Tests for src.exploration.gwas_overlap."""

from __future__ import annotations

import gzip
import textwrap

import numpy as np
import pandas as pd
import pytest

from dev.bench.gwas_overlap import (
    DEFAULT_GENE_POSITIONS,
    GWAS_LP_THRESHOLD,
    annotate_consensus,
    compute_length_diagnostic,
    fisher_enrichment,
    gene_to_best_variant,
    length_matched_fisher_enrichment,
    length_matched_permutation_test,
    load_gene_positions,
    map_variants_to_nearest_gene,
    parse_gwas,
    parse_gwas_ssf,
    parse_gwas_vcf,
)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

VCF_HEADER = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size">
    ##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error">
    ##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value">
    ##FORMAT=<ID=ID,Number=1,Type=String,Description="Variant ID">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
""")


def _write_vcf(tmp_path, body, gz=False):
    """Write a tiny VCF and return its path."""
    name = "fake.vcf.gz" if gz else "fake.vcf"
    path = tmp_path / name
    text = VCF_HEADER + body
    if gz:
        with gzip.open(path, "wt", encoding="utf-8") as f:
            f.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return path


# ----------------------------------------------------------------------
# parse_gwas_vcf
# ----------------------------------------------------------------------

def test_parse_vcf_filters_below_lp_threshold(tmp_path):
    body = (
        "1\t100\trs_a\tA\tG\t.\tPASS\t.\tES:SE:LP:ID\t0.1:0.05:8.0:rs_a\n"
        "1\t200\trs_b\tA\tT\t.\tPASS\t.\tES:SE:LP:ID\t0.05:0.04:5.0:rs_b\n"
        "2\t300\trs_c\tC\tG\t.\tPASS\t.\tES:SE:LP:ID\t0.2:0.07:10.0:rs_c\n"
    )
    path = _write_vcf(tmp_path, body)
    df = parse_gwas_vcf(path, lp_threshold=7.30)
    assert list(df["rsid"]) == ["rs_a", "rs_c"]
    assert df["lp"].tolist() == [8.0, 10.0]
    assert df["beta"].tolist() == [0.1, 0.2]


def test_parse_vcf_handles_gzip(tmp_path):
    body = "1\t100\trs_a\tA\tG\t.\tPASS\t.\tES:SE:LP:ID\t0.1:0.05:9.0:rs_a\n"
    path = _write_vcf(tmp_path, body, gz=True)
    df = parse_gwas_vcf(path, lp_threshold=7.30)
    assert len(df) == 1
    assert df.loc[0, "rsid"] == "rs_a"


def test_parse_vcf_strips_chr_prefix(tmp_path):
    body = "chr3\t999\trs_z\tA\tG\t.\tPASS\t.\tES:SE:LP:ID\t0.1:0.05:9.0:rs_z\n"
    path = _write_vcf(tmp_path, body)
    df = parse_gwas_vcf(path, lp_threshold=7.30)
    assert df.loc[0, "chrom"] == "3"


def test_parse_vcf_skips_rows_with_missing_lp(tmp_path):
    body = (
        "1\t100\trs_a\tA\tG\t.\tPASS\t.\tES:SE:LP:ID\t0.1:0.05:9.0:rs_a\n"
        "1\t200\trs_b\tA\tG\t.\tPASS\t.\tES:SE:LP:ID\t0.1:0.05:.:rs_b\n"
    )
    path = _write_vcf(tmp_path, body)
    df = parse_gwas_vcf(path, lp_threshold=7.30)
    assert list(df["rsid"]) == ["rs_a"]


def test_parse_vcf_handles_missing_es_se(tmp_path):
    """Missing ES/SE (".") should yield NaN, not crash."""
    body = "1\t100\trs_a\tA\tG\t.\tPASS\t.\tES:SE:LP:ID\t.:.:9.0:rs_a\n"
    path = _write_vcf(tmp_path, body)
    df = parse_gwas_vcf(path, lp_threshold=7.30)
    assert np.isnan(df.loc[0, "beta"])
    assert np.isnan(df.loc[0, "se"])


# ----------------------------------------------------------------------
# load_gene_positions (uses bundled data)
# ----------------------------------------------------------------------

def test_load_bundled_gene_positions():
    df = load_gene_positions()
    assert len(df) > 1000
    for col in ("gene", "chrom", "start", "end", "strand"):
        assert col in df.columns
    # Standard housekeeping genes should be present
    assert "GAPDH" in df["gene"].values
    assert "ACTB" in df["gene"].values
    # Chromosomes should be plain numbers / X / Y, no "chr" prefix
    assert all(not str(c).startswith("chr") for c in df["chrom"].head(50))


def test_load_gene_positions_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_gene_positions(tmp_path / "does_not_exist.tsv")


# ----------------------------------------------------------------------
# map_variants_to_nearest_gene
# ----------------------------------------------------------------------

@pytest.fixture
def small_genes():
    """Three genes on chrom 1, one on chrom 2."""
    return pd.DataFrame({
        "gene":   ["GENE_A", "GENE_B", "GENE_C", "GENE_X"],
        "chrom":  ["1", "1", "1", "2"],
        "start":  [1_000_000, 2_000_000, 3_000_000, 5_000_000],
        "end":    [1_100_000, 2_100_000, 3_100_000, 5_100_000],
        "strand": ["+", "-", "+", "+"],
    })


def test_map_variant_inside_gene_body_distance_zero(small_genes):
    variants = pd.DataFrame([
        {"chrom": "1", "pos": 1_050_000, "rsid": "rs1",
         "ref": "A", "alt": "G", "lp": 9.0, "beta": 0.1, "se": 0.05},
    ])
    out = map_variants_to_nearest_gene(variants, small_genes, window_kb=100)
    assert len(out) == 1
    assert out.loc[0, "gene"] == "GENE_A"
    assert out.loc[0, "distance_bp"] == 0


def test_map_variant_picks_nearest_within_window(small_genes):
    # Variant at 1.95M -> closer to GENE_B.start (2.0M, 50kb) than
    # GENE_A.end (1.1M, 850kb). With 100kb window only GENE_B reachable.
    variants = pd.DataFrame([
        {"chrom": "1", "pos": 1_950_000, "rsid": "rs2",
         "ref": "A", "alt": "G", "lp": 9.0, "beta": 0.1, "se": 0.05},
    ])
    out = map_variants_to_nearest_gene(variants, small_genes, window_kb=100)
    assert len(out) == 1
    assert out.loc[0, "gene"] == "GENE_B"
    assert out.loc[0, "distance_bp"] == 50_000


def test_map_variant_outside_window_dropped(small_genes):
    # Variant at 1.5M -> 400kb from GENE_A.end and 500kb from GENE_B.start.
    # With 100kb window, no gene is reachable.
    variants = pd.DataFrame([
        {"chrom": "1", "pos": 1_500_000, "rsid": "rs3",
         "ref": "A", "alt": "G", "lp": 9.0, "beta": 0.1, "se": 0.05},
    ])
    out = map_variants_to_nearest_gene(variants, small_genes, window_kb=100)
    assert len(out) == 0


def test_map_variant_on_unknown_chrom_dropped(small_genes):
    variants = pd.DataFrame([
        {"chrom": "99", "pos": 1_000, "rsid": "rs_alien",
         "ref": "A", "alt": "G", "lp": 9.0, "beta": 0.1, "se": 0.05},
    ])
    out = map_variants_to_nearest_gene(variants, small_genes, window_kb=100)
    assert len(out) == 0


def test_map_strips_chr_prefix_in_variants(small_genes):
    variants = pd.DataFrame([
        {"chrom": "chr2", "pos": 5_050_000, "rsid": "rs_x",
         "ref": "A", "alt": "G", "lp": 9.0, "beta": 0.1, "se": 0.05},
    ])
    out = map_variants_to_nearest_gene(variants, small_genes, window_kb=100)
    assert len(out) == 1
    assert out.loc[0, "gene"] == "GENE_X"


# ----------------------------------------------------------------------
# gene_to_best_variant
# ----------------------------------------------------------------------

def test_gene_to_best_keeps_most_significant_per_gene():
    vg = pd.DataFrame([
        {"chrom": "1", "pos": 1_010_000, "rsid": "rs_low",
         "lp": 8.0, "beta": 0.1, "se": 0.05,
         "gene": "GENE_A", "distance_bp": 0},
        {"chrom": "1", "pos": 1_050_000, "rsid": "rs_high",
         "lp": 12.0, "beta": 0.2, "se": 0.04,
         "gene": "GENE_A", "distance_bp": 0},
        {"chrom": "1", "pos": 2_010_000, "rsid": "rs_b",
         "lp": 9.5, "beta": 0.3, "se": 0.06,
         "gene": "GENE_B", "distance_bp": 0},
    ])
    out = gene_to_best_variant(vg)
    assert "GENE_A" in out.index and "GENE_B" in out.index
    assert out.loc["GENE_A", "rsid"] == "rs_high"
    assert out.loc["GENE_A", "lp"] == 12.0
    assert out.loc["GENE_A", "n_hits_in_window"] == 2
    assert out.loc["GENE_B", "n_hits_in_window"] == 1


def test_gene_to_best_empty_input_returns_empty_df():
    out = gene_to_best_variant(pd.DataFrame())
    assert len(out) == 0
    assert "rsid" in out.columns
    assert "n_hits_in_window" in out.columns


# ----------------------------------------------------------------------
# annotate_consensus
# ----------------------------------------------------------------------

def test_annotate_adds_columns_to_consensus():
    consensus = pd.DataFrame({
        "names": ["GENE_A", "GENE_B", "GENE_C", "BOGUS"],
        "logfoldchanges": [1.0, 2.0, 1.5, 0.5],
    })
    gtb = pd.DataFrame({
        "rsid": ["rs1", "rs2"],
        "chrom": ["1", "1"],
        "pos": [1_050_000, 2_050_000],
        "lp": [12.0, 9.5],
        "beta": [0.2, 0.3],
        "se": [0.04, 0.06],
        "distance_bp": [0, 50_000],
        "n_hits_in_window": [2, 1],
    }, index=pd.Index(["GENE_A", "GENE_B"], name="gene"))
    out = annotate_consensus(consensus, gtb, gene_col="names")
    assert list(out["gwas_hit_in_window"]) == [True, True, False, False]
    assert out.loc[0, "gwas_nearest_rsid"] == "rs1"
    assert out.loc[0, "gwas_nearest_lp"] == 12.0
    assert out.loc[1, "gwas_nearest_distance_kb"] == 50.0
    assert out.loc[3, "gwas_nearest_rsid"] == ""


def test_annotate_handles_empty_gene_to_best():
    consensus = pd.DataFrame({"names": ["A", "B"]})
    out = annotate_consensus(consensus, pd.DataFrame(), gene_col="names")
    assert not out["gwas_hit_in_window"].any()
    assert (out["gwas_n_hits_in_window"] == 0).all()


def test_annotate_missing_gene_col_raises():
    consensus = pd.DataFrame({"foo": ["A"]})
    with pytest.raises(KeyError):
        annotate_consensus(consensus, pd.DataFrame(), gene_col="names")


# ----------------------------------------------------------------------
# fisher_enrichment
# ----------------------------------------------------------------------

def test_fisher_enrichment_strong_positive():
    """All hit-bearing genes are in the consensus -> very small p."""
    background = ["G" + str(i) for i in range(100)]
    consensus  = ["G0", "G1", "G2", "G3", "G4"]
    hits       = ["G0", "G1", "G2"]   # 3/5 of consensus, 0 of others
    res = fisher_enrichment(consensus, background, hits)
    assert res["n_consensus_with_hit"] == 3
    assert res["n_background_with_hit"] == 0
    assert res["pvalue"] < 0.001
    assert res["odds_ratio"] == float("inf")


def test_fisher_enrichment_no_signal():
    """Hits distributed proportionally -> p ~ 1."""
    background = ["G" + str(i) for i in range(100)]
    consensus  = ["G0", "G1", "G2", "G3", "G4"]
    # 5% of background (5/100) carries a hit; 0 of consensus
    hits = ["G50", "G60", "G70", "G80", "G90"]
    res = fisher_enrichment(consensus, background, hits)
    assert res["pvalue"] > 0.5


def test_fisher_enrichment_consensus_subset_of_background():
    """Background should always include consensus genes (we union)."""
    background = ["G50"]               # tiny background; consensus extends it
    consensus  = ["G0", "G1", "G2"]
    hits       = ["G0"]
    res = fisher_enrichment(consensus, background, hits)
    # n_background = background_only = {G50} -> 1 (excludes consensus)
    assert res["n_background"] == 1
    assert res["n_consensus"] == 3


def test_fisher_enrichment_degenerate_returns_nan():
    """Empty consensus or empty background -> degenerate; no crash."""
    res = fisher_enrichment([], ["G1", "G2"], ["G1"])
    assert np.isnan(res["pvalue"])
    res = fisher_enrichment(["A"], [], [])
    assert np.isnan(res["pvalue"])


# ----------------------------------------------------------------------
# end-to-end: tiny VCF + tiny gene table -> annotated consensus
# ----------------------------------------------------------------------

def test_end_to_end_small(tmp_path, small_genes):
    body = (
        # Sig variant inside GENE_A
        "1\t1050000\trs_a\tA\tG\t.\tPASS\t.\tES:SE:LP:ID\t0.1:0.05:9.0:rs_a\n"
        # Sig variant 50kb outside GENE_B
        "1\t1950000\trs_b\tA\tG\t.\tPASS\t.\tES:SE:LP:ID\t0.1:0.05:8.5:rs_b\n"
        # Non-sig variant (should be dropped)
        "1\t3050000\trs_c\tA\tG\t.\tPASS\t.\tES:SE:LP:ID\t0.0:0.05:2.0:rs_c\n"
        # Sig variant on chrom 2 inside GENE_X
        "2\t5050000\trs_x\tA\tG\t.\tPASS\t.\tES:SE:LP:ID\t0.2:0.05:11.0:rs_x\n"
    )
    path = _write_vcf(tmp_path, body)

    sig = parse_gwas_vcf(path, lp_threshold=7.30)
    assert len(sig) == 3       # rs_a, rs_b, rs_x

    mapped = map_variants_to_nearest_gene(sig, small_genes, window_kb=100)
    assert len(mapped) == 3
    gtb = gene_to_best_variant(mapped)
    assert set(gtb.index) == {"GENE_A", "GENE_B", "GENE_X"}

    consensus = pd.DataFrame({"names": ["GENE_A", "GENE_C", "GENE_X"]})
    ann = annotate_consensus(consensus, gtb, gene_col="names")
    assert list(ann["gwas_hit_in_window"]) == [True, False, True]


# ----------------------------------------------------------------------
# parse_gwas_ssf (GWAS Catalog tab-separated format)
# ----------------------------------------------------------------------

def _write_ssf(tmp_path, header_cols, rows, gz=False):
    name = "fake.tsv.gz" if gz else "fake.tsv"
    path = tmp_path / name
    text = "\t".join(header_cols) + "\n"
    for r in rows:
        text += "\t".join(str(x) for x in r) + "\n"
    if gz:
        with gzip.open(path, "wt", encoding="utf-8") as f:
            f.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return path


def test_parse_ssf_basic_canonical_columns(tmp_path):
    cols = ["chromosome", "base_pair_location", "effect_allele",
            "other_allele", "beta", "standard_error",
            "p_value", "rsid"]
    rows = [
        ["1", 100, "A", "G", 0.1, 0.05, 1e-9, "rs_a"],   # sig
        ["1", 200, "A", "T", 0.05, 0.04, 0.5, "rs_b"],   # not sig
        ["2", 300, "C", "G", 0.2, 0.07, 1e-12, "rs_c"],  # sig
    ]
    path = _write_ssf(tmp_path, cols, rows)
    df = parse_gwas_ssf(path, lp_threshold=7.30)
    assert list(df["rsid"]) == ["rs_a", "rs_c"]
    assert abs(df.loc[0, "lp"] - 9.0) < 1e-6
    assert abs(df.loc[1, "lp"] - 12.0) < 1e-6
    assert df.loc[0, "beta"] == 0.1
    assert df.loc[1, "se"] == 0.07


def test_parse_ssf_prefers_hm_columns_when_present(tmp_path):
    """If hm_* columns are present they should be preferred over the
    un-prefixed canonical names."""
    cols = ["chromosome", "base_pair_location",
            "hm_chromosome", "hm_base_pair_location",
            "p_value", "rsid", "hm_rsid"]
    rows = [
        ["1", 100, "22", 999, 1e-9, "rs_orig", "rs_harmonised"],
    ]
    path = _write_ssf(tmp_path, cols, rows)
    df = parse_gwas_ssf(path, lp_threshold=7.30)
    assert df.loc[0, "chrom"] == "22"
    assert df.loc[0, "pos"] == 999
    assert df.loc[0, "rsid"] == "rs_harmonised"


def test_parse_ssf_handles_gzip(tmp_path):
    cols = ["chromosome", "base_pair_location", "p_value", "rsid"]
    rows = [["1", 100, 1e-9, "rs_a"]]
    path = _write_ssf(tmp_path, cols, rows, gz=True)
    df = parse_gwas_ssf(path, lp_threshold=7.30)
    assert len(df) == 1


def test_parse_ssf_skips_invalid_p_values(tmp_path):
    cols = ["chromosome", "base_pair_location", "p_value", "rsid"]
    rows = [
        ["1", 100, "NA", "rs_a"],
        ["1", 200, "0", "rs_b"],
        ["1", 300, "-0.1", "rs_c"],
        ["1", 400, "1e-12", "rs_d"],
    ]
    path = _write_ssf(tmp_path, cols, rows)
    df = parse_gwas_ssf(path, lp_threshold=7.30)
    assert list(df["rsid"]) == ["rs_d"]


def test_parse_ssf_strips_chr_prefix(tmp_path):
    cols = ["chromosome", "base_pair_location", "p_value", "rsid"]
    rows = [["chr3", 100, 1e-9, "rs_z"]]
    path = _write_ssf(tmp_path, cols, rows)
    df = parse_gwas_ssf(path, lp_threshold=7.30)
    assert df.loc[0, "chrom"] == "3"


def test_parse_ssf_missing_required_columns_raises(tmp_path):
    cols = ["foo", "bar", "rsid"]
    rows = [["x", "y", "rs"]]
    path = _write_ssf(tmp_path, cols, rows)
    with pytest.raises(ValueError, match="missing required column"):
        parse_gwas_ssf(path)


# ----------------------------------------------------------------------
# parse_gwas auto-detect dispatcher
# ----------------------------------------------------------------------

def test_parse_gwas_dispatches_to_vcf(tmp_path):
    body = "1\t100\trs_a\tA\tG\t.\tPASS\t.\tES:SE:LP:ID\t0.1:0.05:9.0:rs_a\n"
    path = _write_vcf(tmp_path, body)
    df = parse_gwas(path)
    assert len(df) == 1 and df.loc[0, "rsid"] == "rs_a"


def test_parse_gwas_dispatches_to_ssf(tmp_path):
    cols = ["chromosome", "base_pair_location", "p_value", "rsid"]
    rows = [["1", 100, 1e-9, "rs_a"]]
    path = _write_ssf(tmp_path, cols, rows)
    df = parse_gwas(path)
    assert len(df) == 1 and df.loc[0, "rsid"] == "rs_a"


def test_parse_gwas_dispatches_to_ssf_gzipped(tmp_path):
    cols = ["chromosome", "base_pair_location", "p_value", "rsid"]
    rows = [["1", 100, 1e-9, "rs_a"]]
    path = _write_ssf(tmp_path, cols, rows, gz=True)
    df = parse_gwas(path)
    assert len(df) == 1


def test_parse_gwas_unknown_format_raises(tmp_path):
    path = tmp_path / "weird.txt"
    path.write_text("just\nrandom\nlines\n")
    with pytest.raises(ValueError, match="auto-detect"):
        parse_gwas(path)


# ----------------------------------------------------------------------
# compute_length_diagnostic + length_matched_fisher_enrichment
# ----------------------------------------------------------------------

def _make_gene_positions(n_short=200, n_long=200, short_len=5_000,
                         long_len=500_000, seed=0):
    """Construct a synthetic gene-positions DataFrame with known
    length bimodality. First n_short genes are ~5 kb; next n_long are
    ~500 kb. All on chrom 1 with non-overlapping positions."""
    rng = np.random.default_rng(seed)
    rows = []
    cursor = 1_000_000
    # Short genes
    for i in range(n_short):
        L = short_len + int(rng.integers(-1000, 1000))
        rows.append({'gene': f'SHORT_{i:04d}', 'chrom': '1',
                     'start': cursor, 'end': cursor + L, 'strand': '+'})
        cursor += L + 200_000  # plenty of gap
    # Long genes
    for i in range(n_long):
        L = long_len + int(rng.integers(-10000, 10000))
        rows.append({'gene': f'LONG_{i:04d}', 'chrom': '1',
                     'start': cursor, 'end': cursor + L, 'strand': '+'})
        cursor += L + 200_000
    return pd.DataFrame(rows)


def test_length_diagnostic_detects_bias_ratio():
    """When consensus is drawn from long genes, the ratio should
    reflect the ~100x length difference between SHORT and LONG."""
    genes = _make_gene_positions()
    # Consensus = all 200 LONG genes; background = all 400 genes
    consensus = [f'LONG_{i:04d}' for i in range(200)]
    background = [f'LONG_{i:04d}' for i in range(200)] + \
                 [f'SHORT_{i:04d}' for i in range(200)]
    # Arbitrary hit set (irrelevant for the length ratio test)
    hits = [f'LONG_{i:04d}' for i in range(50)]

    d = compute_length_diagnostic(consensus, background, hits, genes)
    # Consensus is all long; bg-only is all short -> very large ratio
    assert d['length_ratio_consensus_over_background'] > 50
    assert d['n_consensus_with_length'] == 200
    # After bg -= consensus, the 200 LONG are removed; bg_only = 200 SHORT
    assert d['n_background_with_length'] == 200


def test_length_diagnostic_ratio_near_one_when_no_bias():
    """Random consensus drawn from the mixed background should have
    a length ratio near 1.0."""
    genes = _make_gene_positions()
    all_genes = [f'SHORT_{i:04d}' for i in range(200)] + \
                [f'LONG_{i:04d}'  for i in range(200)]
    rng = np.random.default_rng(7)
    consensus = list(rng.choice(all_genes, size=100, replace=False))
    background = all_genes
    hits = [f'LONG_{i:04d}' for i in range(10)]

    d = compute_length_diagnostic(consensus, background, hits, genes)
    # Random draw from mixed pool -> expect ratio ~1 (+/- 0.5 for noise
    # with this small sample)
    assert 0.3 < d['length_ratio_consensus_over_background'] < 3.0


def test_length_matched_strips_pure_length_bias():
    """Scenario: consensus is entirely LONG genes; hits are entirely
    LONG genes (long genes happen to be near GWAS hits because they
    are long). Naive Fisher reports huge enrichment. Length-matched
    MH should report OR ~1 / p > 0.05 -- because WITHIN the long
    decile, consensus is no more likely than background to have a hit."""
    genes = _make_gene_positions(n_short=200, n_long=200)
    # Consensus = 100 random LONG genes
    rng = np.random.default_rng(1)
    long_all = [f'LONG_{i:04d}' for i in range(200)]
    short_all = [f'SHORT_{i:04d}' for i in range(200)]
    consensus = list(rng.choice(long_all, size=100, replace=False))
    # Background = everyone
    background = long_all + short_all
    # Hits: half the LONG genes (regardless of consensus membership) +
    # none of the SHORT genes. This is the "only long genes have hits"
    # pure length-bias scenario.
    hits = list(rng.choice(long_all, size=100, replace=False))

    raw = fisher_enrichment(consensus, background, hits)
    mh  = length_matched_fisher_enrichment(
        consensus, background, hits, genes, n_bins=2)

    # Raw Fisher sees big enrichment (all hits are long, all consensus is long)
    assert raw['odds_ratio'] > 2 or not np.isfinite(raw['odds_ratio'])
    # Length-matched MH collapses it: OR ~1, not significant
    assert 0.5 < mh['odds_ratio'] < 2.0
    assert mh['pvalue'] > 0.05


def test_length_matched_preserves_real_enrichment():
    """Scenario: hits are enriched in consensus INDEPENDENTLY of length.
    Within each length decile, consensus genes have 3x the hit rate of
    background. Both raw Fisher and length-matched MH should report
    enrichment."""
    genes = _make_gene_positions(n_short=200, n_long=200)
    rng = np.random.default_rng(2)
    long_all = [f'LONG_{i:04d}' for i in range(200)]
    short_all = [f'SHORT_{i:04d}' for i in range(200)]
    # Consensus: 50 long + 50 short -> unbiased length-wise
    consensus = (list(rng.choice(long_all, size=50, replace=False))
                 + list(rng.choice(short_all, size=50, replace=False)))
    background = long_all + short_all
    # Hits: within each category, 60% of consensus genes are hits,
    # only 20% of non-consensus genes are hits -> real enrichment
    # independent of length
    bg_only_long = [g for g in long_all if g not in consensus]
    bg_only_short = [g for g in short_all if g not in consensus]
    hit_long = (list(rng.choice([g for g in consensus if g in long_all],
                                size=30, replace=False))
                + list(rng.choice(bg_only_long, size=30, replace=False)))
    hit_short = (list(rng.choice([g for g in consensus if g in short_all],
                                 size=30, replace=False))
                 + list(rng.choice(bg_only_short, size=30, replace=False)))
    hits = hit_long + hit_short

    raw = fisher_enrichment(consensus, background, hits)
    mh  = length_matched_fisher_enrichment(
        consensus, background, hits, genes, n_bins=2)

    assert raw['odds_ratio'] > 3
    # Length-matched MH should STILL see strong enrichment because
    # the signal is real and orthogonal to length
    assert mh['odds_ratio'] > 2
    assert mh['pvalue'] < 0.01


def test_length_matched_handles_empty_inputs():
    genes = _make_gene_positions()
    mh = length_matched_fisher_enrichment([], ['SHORT_0001'], [], genes)
    assert not np.isfinite(mh['odds_ratio']) or mh['n_strata_used'] == 0


# ----------------------------------------------------------------------
# length_matched_permutation_test
# ----------------------------------------------------------------------

def test_permutation_reports_nonsignificant_when_length_bias_only():
    """Length-bias-only scenario: the raw OR may be large (because
    all hits ARE concentrated in long genes), but under a
    length-matched null the observed OR sits in the MIDDLE of the
    null distribution -- empirical p should NOT be significant."""
    genes = _make_gene_positions(n_short=200, n_long=200)
    rng = np.random.default_rng(11)
    long_all = [f'LONG_{i:04d}' for i in range(200)]
    short_all = [f'SHORT_{i:04d}' for i in range(200)]
    consensus = long_all[:100]
    background = long_all + short_all
    hits = list(rng.choice(long_all, size=100, replace=False))

    perm = length_matched_permutation_test(
        consensus, background, hits, genes,
        n_bins=2, n_perms=500, seed=11)
    # Observed OR sits inside the null 95% CI: no real signal above
    # what length matching expects. (Empirical p is exactly the
    # CI test in distributional form, but stochastic from one seed
    # to the next; the CI containment check is the cleaner assertion.)
    assert (perm['null_or_ci_low']
            <= perm['observed_odds_ratio']
            <= perm['null_or_ci_high'])


def test_permutation_detects_real_signal():
    """When consensus genes are truly enriched for hits independently
    of length, both the observed OR and empirical p should reflect
    that."""
    genes = _make_gene_positions(n_short=200, n_long=200)
    rng = np.random.default_rng(22)
    long_all = [f'LONG_{i:04d}' for i in range(200)]
    short_all = [f'SHORT_{i:04d}' for i in range(200)]
    consensus = (list(rng.choice(long_all, size=50, replace=False))
                 + list(rng.choice(short_all, size=50, replace=False)))
    background = long_all + short_all
    bg_only_long = [g for g in long_all if g not in consensus]
    bg_only_short = [g for g in short_all if g not in consensus]
    # Hits: most of consensus + only a few of non-consensus -> real signal
    hits = (list(rng.choice([g for g in consensus if g in long_all],
                            size=40, replace=False))
            + list(rng.choice([g for g in consensus if g in short_all],
                              size=40, replace=False))
            + list(rng.choice(bg_only_long, size=10, replace=False))
            + list(rng.choice(bg_only_short, size=10, replace=False)))

    perm = length_matched_permutation_test(
        consensus, background, hits, genes, n_bins=2, n_perms=500)
    # Observed OR should clearly exceed the null 95% CI upper bound
    assert perm['observed_odds_ratio'] > perm['null_or_ci_high']
    assert perm['empirical_pvalue'] < 0.05


def test_permutation_is_reproducible_with_seed():
    """Given a scenario that produces finite ORs, two runs with the
    same seed must give identical null summaries."""
    genes = _make_gene_positions(n_short=100, n_long=100)
    rng = np.random.default_rng(5)
    # Mixed-length consensus + mixed-length hits so per-stratum
    # cells are non-degenerate
    consensus = ([f'LONG_{i:04d}' for i in range(40)]
                 + [f'SHORT_{i:04d}' for i in range(40)])
    background = ([f'LONG_{i:04d}' for i in range(100)]
                  + [f'SHORT_{i:04d}' for i in range(100)])
    hits = (list(rng.choice([f'LONG_{i:04d}' for i in range(100)],
                            size=40, replace=False))
            + list(rng.choice([f'SHORT_{i:04d}' for i in range(100)],
                              size=40, replace=False)))

    r1 = length_matched_permutation_test(
        consensus, background, hits, genes,
        n_bins=2, seed=99, n_perms=100)
    r2 = length_matched_permutation_test(
        consensus, background, hits, genes,
        n_bins=2, seed=99, n_perms=100)
    assert r1['null_or_mean'] == r2['null_or_mean']
    assert r1['empirical_pvalue'] == r2['empirical_pvalue']


def test_permutation_handles_empty_inputs():
    genes = _make_gene_positions()
    r = length_matched_permutation_test([], ['SHORT_0001'], [], genes,
                                        n_perms=50)
    assert 'observed_odds_ratio' in r


def test_length_matched_handles_genes_missing_from_positions():
    """A consensus gene not in the positions table shouldn't crash
    the test."""
    genes = _make_gene_positions(n_short=100, n_long=100)
    consensus = ['LONG_0001', 'LONG_0002', 'UNKNOWN_GENE']
    background = [f'LONG_{i:04d}' for i in range(100)] + \
                 [f'SHORT_{i:04d}' for i in range(100)] + \
                 ['UNKNOWN_GENE']
    hits = ['LONG_0050']
    mh = length_matched_fisher_enrichment(
        consensus, background, hits, genes, n_bins=2)
    # Doesn't crash; unknown gene is silently dropped
    assert 'odds_ratio' in mh


def test_lp_threshold_default_matches_5e_minus_8():
    # -log10(5e-8) = 7.30103
    assert abs(GWAS_LP_THRESHOLD - 7.30103) < 1e-4


def test_default_gene_positions_path_exists():
    assert DEFAULT_GENE_POSITIONS.exists()
