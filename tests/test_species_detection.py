"""Human<->mouse gene-name conversion used to match human-authored gene
sets (pathways, markers) against mouse data.

format_gene_for_species has three tiers, tried in order: a hand-verified
legacy-nomenclature patch (Ensembl/CellRanger mouse references often still
carry pre-2016 ATP synthase symbols that neither current HGNC nor current
MGI use any more), the bundled MGI 1:1 ortholog table, then a naive
Title-case guess. These tests pin each tier and the one confirmed-bad MGI
entry that must be blocked rather than trusted.
"""
from __future__ import annotations

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
    # Naive 'Gpi'/'Fh'/'Scd' are not real mouse gene symbols; the real
    # orthologs (from MGI) are Gpi1 / Fh1 / Scd1.
    assert format_gene_for_species("GPI", "mouse") == "Gpi1"
    assert format_gene_for_species("FH", "mouse") == "Fh1"


def test_format_gene_for_species_blocks_known_bad_ortholog_entry():
    # MGI classes human CS (citrate synthase) against mouse 'Csl', which
    # is not citrate synthase -- must fall through to the correct naive
    # guess instead of trusting the table.
    assert format_gene_for_species("CS", "mouse") == "Cs"
