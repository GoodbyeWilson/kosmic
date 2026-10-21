"""Human<->mouse gene-name conversion used to match human-authored gene
sets (pathways, markers) against mouse data.

format_gene_for_species has three tiers, tried in order: a hand-verified
legacy-nomenclature patch (Ensembl/CellRanger mouse references often still
carry pre-2016 ATP synthase symbols that neither current HGNC nor current
MGI use any more), the bundled MGI ortholog table, then a naive
Title-case guess. These tests pin each tier, the one confirmed-bad MGI
entry that must be blocked rather than trusted, and the reduction of
one-human-to-several-mouse homology classes to one row per human symbol.
"""
from __future__ import annotations

from importlib.resources import files

from kosmic.scrna.inspect.detection import detect_species, format_gene_for_species


def test_detect_species_human_vs_mouse():
    assert detect_species(["HK1", "PKM", "LDHA", "GAPDH"]) == "human"
    assert detect_species(["Hk1", "Pkm", "Ldha", "Gapdh"]) == "mouse"
    assert detect_species([]) == "human"


def test_format_gene_for_species_human_is_uppercase():
    assert format_gene_for_species("Hk1", "human") == "HK1"


def test_format_gene_for_species_naive_fallback():
    # No ortholog-table / legacy-patch entry -> plain Title-case guess.
    assert format_gene_for_species("GAPDH", "mouse") == "Gapdh"


def test_format_gene_for_species_mt_prefix():
    assert format_gene_for_species("MT-ND1", "mouse") == "mt-Nd1"


def test_format_gene_for_species_legacy_atp_synthase():
    # HGNC's 2016 ATP synthase rename; many Ensembl/CellRanger mouse
    # references still carry the pre-rename symbol under either naming.
    assert format_gene_for_species("ATP5F1A", "mouse") == "Atp5a1"
    assert format_gene_for_species("ATP5MC1", "mouse") == "Atp5g1"


def test_format_gene_for_species_uses_ortholog_table_when_naive_is_wrong():
    # Naive 'Gpi'/'Fh' are not real mouse gene symbols; the real
    # orthologs (from MGI) are Gpi1 / Fh1.
    assert format_gene_for_species("GPI", "mouse") == "Gpi1"
    assert format_gene_for_species("FH", "mouse") == "Fh1"


def test_format_gene_for_species_blocks_known_bad_ortholog_entry():
    # MGI classes human CS (citrate synthase) against mouse 'Csl', which
    # is not citrate synthase -- must fall through to the correct naive
    # guess instead of trusting the table.
    assert format_gene_for_species("CS", "mouse") == "Cs"


def test_ortholog_table_has_one_row_per_symbol():
    # The loader keeps the last row per human symbol, so a duplicate
    # would silently pick one paralog. Mouse symbols are unique as well.
    path = files("kosmic.reference.orthologs") / "mgi_human_mouse_1to1.tsv"
    rows = [line.split("\t")
            for line in path.read_text(encoding="utf-8").splitlines()[1:]]
    humans = [human.upper() for human, _ in rows]
    mice = [mouse for _, mouse in rows]
    assert len(set(humans)) == len(humans)
    assert len(set(mice)) == len(mice)


def test_format_gene_for_species_prefers_same_name_ortholog():
    # MGI classes each of these human genes with two mouse genes; the
    # same-name mouse gene is the one kept.
    assert format_gene_for_species("ALDH1A1", "mouse") == "Aldh1a1"  # not Aldh1a7
    assert format_gene_for_species("ABCG2", "mouse") == "Abcg2"  # not Abcg3
    assert format_gene_for_species("AHCY", "mouse") == "Ahcy"  # not Ahcyl


def test_format_gene_for_species_ambiguous_without_same_name_falls_back():
    # AGTR1 -> Agtr1a, Agtr1b has no same-name mouse gene, so it is not
    # in the table and gets the naive guess rather than one paralog.
    assert format_gene_for_species("AGTR1", "mouse") == "Agtr1"
