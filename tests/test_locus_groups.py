"""Locus-group classification of a dataset's gene names.

The point of the count is not the protein-coding number on its own but the
'unrecognised' bucket beside it: retired clone-based names, bare Ensembl IDs
and typos all fail to match the same gene in another study, and a study whose
bucket is large will lose genes in a cross-study join.
"""
import pytest

from kosmic.scrna.inspect.gene_names import (
    count_locus_groups,
    load_hgnc_lookup,
    load_locus_groups,
)


@pytest.fixture(scope="module")
def groups():
    return load_locus_groups()


def test_bundled_table_covers_the_whole_approved_symbol_set(groups):
    assert len(groups) > 40000
    assert groups["TTN"] == "protein-coding gene"
    assert groups["MALAT1"] == "non-coding RNA"


def test_counts_partition_the_gene_list(groups):
    names = ["TTN", "NPPA", "MALAT1", "ENSG00000239945", "RP11-34P13.3"]
    counts = count_locus_groups(names, locus_groups=groups)
    assert counts["total"] == 5
    assert sum(v for k, v in counts.items() if k != "total") == 5


def test_unmatched_names_land_in_unrecognised(groups):
    counts = count_locus_groups(
        ["ENSG00000239945", "RP11-34P13.3"], locus_groups=groups)
    assert counts["unrecognised"] == 2
    assert counts["protein-coding gene"] == 0


def test_protein_coding_and_noncoding_are_separated(groups):
    counts = count_locus_groups(["TTN", "MYH7", "MALAT1"], locus_groups=groups)
    assert counts["protein-coding gene"] == 2
    assert counts["non-coding RNA"] == 1


def test_lookup_recovers_retired_symbols(groups):
    """CARDIOMYOPATHY-era name -- harmonised it becomes a counted gene."""
    lookup, _ = load_hgnc_lookup()
    retired = next(old for old, new in lookup.items()
                   if groups.get(new) == "protein-coding gene")
    before = count_locus_groups([retired], locus_groups=groups)
    after = count_locus_groups([retired], locus_groups=groups, lookup=lookup)
    assert before["unrecognised"] == 1
    assert after["protein-coding gene"] == 1


def test_duplicate_names_are_counted_once(groups):
    """Harmonisation merges duplicates, so the count is of distinct names."""
    lookup = {"OLD1": "TTN", "OLD2": "TTN"}
    counts = count_locus_groups(["OLD1", "OLD2"], locus_groups=groups,
                                lookup=lookup)
    assert counts["total"] == 1
    assert counts["protein-coding gene"] == 1


def test_table_loads_from_the_bundle_without_an_explicit_path():
    counts = count_locus_groups(["TTN"])
    assert counts["protein-coding gene"] == 1
