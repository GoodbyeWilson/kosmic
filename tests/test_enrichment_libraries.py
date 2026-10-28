"""Enrichment libraries load with a record of where they came from.

A bundled library is a published file used exactly as released; an
Enrichr library is downloaded when the enrichment runs and has no version
history, so each load is fingerprinted. Enrichment runs are recorded one
entry per method and library.
"""
from __future__ import annotations

import pytest

from kosmic import provenance
from kosmic.reference.pathways import libraries


def test_bundled_hallmark_is_the_published_file_unaltered():
    sets, source = libraries.load_library('MSigDB_Hallmark_2026.2_Hs')
    assert len(sets) == 50 and 'HALLMARK_OXIDATIVE_PHOSPHORYLATION' in sets
    assert 'COXFA4' in sets['HALLMARK_OXIDATIVE_PHOSPHORYLATION']
    # SHA-1 of h.all.v2026.2.Hs.symbols.gmt as downloaded from the Broad.
    assert source['sha1'] == 'af8ba0e72b34'
    assert source['version'] == '2026.2.Hs' and source['n_genes'] == 4384
    assert source['accessed'].endswith('+00:00')


def test_bundled_mouse_hallmark_has_mouse_symbols():
    """MSigDB's mouse edition, including the MHC genes the human file's
    names cannot match ('HLA-DRA' against mouse 'H2-Ea')."""
    sets, source = libraries.load_library('MSigDB_Hallmark_2026.2_Mm')
    assert len(sets) == 50 and source['version'] == '2026.2.Mm'
    assert source['sha1'] == '6a3ff2442fa2' and source['n_genes'] == 4291
    assert 'H2-Ea' in {g for genes in sets.values() for g in genes}


def test_enrichment_of_an_earlier_pooling_is_left_out_of_the_methods(tmp_path):
    """Rerunning the gene meta-analysis makes the enrichment of the
    previous pooled list describe a result that no longer exists."""
    def record(stage, params):
        provenance.record_stage(tmp_path, 'S', stage, params,
                                replace=stage == 'meta_gene',
                                replace_matching=({'method': params['method'],
                                                   'library': params['library']}
                                                  if stage == 'meta_enrichment' else None))

    record('meta_gene', {'n_significant': 178})
    record('meta_enrichment', {'method': 'gsea', 'library': 'KEGG 2026'})
    record('meta_gene', {'n_significant': 451})          # pooling rerun
    record('meta_enrichment', {'method': 'gsea', 'library': 'Hallmark 2026.2'})
    text = provenance.render_methods(provenance.load(tmp_path))
    assert 'Hallmark 2026.2' in text and 'KEGG 2026' not in text
    assert '1 pathway enrichment of the pooled genes run(s) on an earlier run' in text
    assert len(provenance.load(tmp_path)['stages']) == 3     # the record keeps it


def test_set_names_are_shown_without_the_collection_prefix():
    assert libraries.term_label('HALLMARK_OXIDATIVE_PHOSPHORYLATION') == 'OXIDATIVE PHOSPHORYLATION'
    assert libraries.term_label('TNF-alpha Signaling via NF-kB') == 'TNF-alpha Signaling via NF-kB'
    assert libraries.term_label('CITRATE CYCLE (TCA CYCLE)') == 'CITRATE CYCLE (TCA CYCLE)'


def test_an_enrichr_download_is_fingerprinted(monkeypatch):
    from kosmic.reference.pathways import enrichr
    text = "SET A\t\tG1\tG2\nSET B\t\tG2\tG3\n"
    monkeypatch.setattr(enrichr, 'fetch_library_text', lambda name, timeout=60: text)
    sets, source = libraries.load_library('KEGG_2026')
    assert sets == {'SET A': ['G1', 'G2'], 'SET B': ['G2', 'G3']}
    assert source['provider'].startswith('Enrichr')
    assert source['url'].endswith('libraryName=KEGG_2026')
    assert source['n_sets'] == 2 and source['n_genes'] == 3
    assert len(source['sha1']) == 12


@pytest.mark.parametrize('meta', [False, True])
def test_one_entry_per_method_and_library(tmp_path, meta):
    """Rerunning a library replaces its entry; another library adds one,
    and the methods list both."""
    stage = 'meta_enrichment' if meta else 'enrichment'

    def run(library, n_sig):
        provenance.record_stage(
            tmp_path, 'S', stage, {'method': 'gsea', 'library': library,
                                   'n_terms_sig': n_sig},
            replace_matching={'method': 'gsea', 'library': library})

    run('KEGG 2026', 1)
    run('Hallmark 2026.2', 9)
    run('KEGG 2026', 3)
    stages = provenance.load(tmp_path)['stages']
    assert [(s['params']['library'], s['params']['n_terms_sig']) for s in stages] == [
        ('Hallmark 2026.2', 9), ('KEGG 2026', 3)]
    text = provenance.render_methods(provenance.load(tmp_path))
    assert 'Hallmark 2026.2' in text and 'KEGG 2026' in text


@pytest.fixture
def app():
    import os
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_a_study_enrichment_run_is_recorded_with_its_library(app, tmp_path):
    import pandas as pd

    from kosmic.gui.de_analysis.pages.enrichment_page import EnrichmentPage
    from kosmic.gui.meta_analysis.pages.enrichment_step import _MetaEnrichmentContext

    ws = _MetaEnrichmentContext()
    ws.set_meta_results(pd.DataFrame({'names': ['A', 'B'], 'logfoldchanges': [1.0, -1.0],
                                      'fdr': [0.01, 0.2], 'pvals_pooled': [0.001, 0.1]}))
    ws.h5ad_path = str(tmp_path / 'S.h5ad')          # a study, so the run is recorded
    page = EnrichmentPage(ws)
    page._results_df = pd.DataFrame({'Term': ['T1', 'T2'], 'FDR': [0.01, 0.5]})
    _, page._library_source = libraries.load_library('MSigDB_Hallmark_2026.2_Hs')

    page._record_study_enrichment('ora')
    (entry,) = provenance.load(tmp_path)['stages']
    assert entry['stage'] == 'enrichment'
    assert entry['params']['n_terms_sig'] == 1
    assert entry['params']['library_source']['sha1'] == 'af8ba0e72b34'
    assert 'Pathway enrichment' in provenance.render_methods(provenance.load(tmp_path))
