"""The meta-analysis must pool the studies that were ticked, not every
study it can find.

'discover_de_results' walks the whole project folder, so it also returns
studies the user unticked on the Select Studies page and any per-cell-type
results written beside a whole-dataset one. 'select_de_entries' is the gate
between discovery and pooling; these tests pin that gate open only for the
selected accessions.
"""
from __future__ import annotations

from kosmic.meta_analysis.io import select_de_entries


def _entry(accession, method, file_type='de_results'):
    return {
        'accession': accession,
        'file_type': file_type,
        'de_method': method,
        'variant': method,
        'path': f'/p/{accession}/{accession}_DE_{method}.csv',
    }


ENTRIES = [
    _entry('GSE1', 'deseq2'),
    _entry('GSE1', 'welch_cpm'),
    _entry('GSE2', 'deseq2'),
    _entry('GSE3', 'welch_cpm'),
]


def test_selection_gates_which_studies_are_pooled():
    chosen, _fallback, skipped = select_de_entries(
        ENTRIES, ['deseq2'], allowed=['GSE1'])
    assert [e['accession'] for e in chosen] == ['GSE1']
    assert skipped == ['GSE2', 'GSE3']


def test_no_selection_means_everything_discovered():
    """None is 'no gate set yet', not 'select nothing'."""
    chosen, _fallback, skipped = select_de_entries(ENTRIES, ['deseq2'])
    assert [e['accession'] for e in chosen] == ['GSE1', 'GSE2', 'GSE3']
    assert skipped == []


def test_empty_selection_is_treated_the_same_as_none():
    chosen, _f, skipped = select_de_entries(ENTRIES, ['deseq2'], allowed=[])
    assert len(chosen) == 3
    assert skipped == []


def test_preferred_method_wins_within_an_accession():
    chosen, _f, _s = select_de_entries(ENTRIES, ['deseq2'], allowed=['GSE1'])
    assert chosen[0]['de_method'] == 'deseq2'

    chosen, _f, _s = select_de_entries(ENTRIES, ['welch_cpm'], allowed=['GSE1'])
    assert chosen[0]['de_method'] == 'welch_cpm'


def test_alias_order_decides_the_fallback_within_an_accession():
    chosen, fallback, _s = select_de_entries(
        ENTRIES, ['welch_cpm_eb', 'welch_cpm'], allowed=['GSE1'])
    assert chosen[0]['de_method'] == 'welch_cpm'
    assert fallback == []


def test_missing_preferred_method_falls_back_and_is_reported():
    chosen, fallback, _s = select_de_entries(
        ENTRIES, ['deseq2'], allowed=['GSE3'])
    assert chosen[0]['de_method'] == 'welch_cpm'
    assert fallback == ['GSE3:welch_cpm']


def test_non_de_entries_are_ignored():
    entries = ENTRIES + [_entry('GSE9', 'Norm', file_type='gene')]
    chosen, _f, skipped = select_de_entries(entries, ['deseq2'])
    assert 'GSE9' not in [e['accession'] for e in chosen]
    assert 'GSE9' not in skipped


def test_per_cell_type_results_are_separate_accessions():
    """A per-cell-type run must not drag its parent study into the pool."""
    entries = [
        _entry('GSE1', 'deseq2'),
        _entry('GSE1_Endothelial', 'deseq2'),
        _entry('GSE2_Endothelial', 'deseq2'),
    ]
    chosen, _f, skipped = select_de_entries(
        entries, ['deseq2'], allowed=['GSE1_Endothelial', 'GSE2_Endothelial'])
    assert [e['accession'] for e in chosen] == [
        'GSE1_Endothelial', 'GSE2_Endothelial']
    assert skipped == ['GSE1']


def test_chosen_is_ordered_by_accession():
    entries = [_entry('B', 'deseq2'), _entry('A', 'deseq2'),
               _entry('C', 'deseq2')]
    chosen, _f, _s = select_de_entries(entries, ['deseq2'])
    assert [e['accession'] for e in chosen] == ['A', 'B', 'C']


def test_empty_input():
    assert select_de_entries([], ['deseq2']) == ([], [], [])
