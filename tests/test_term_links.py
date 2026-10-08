"""Enriched Term Links: term and gene selection, ordering and rendering."""
from __future__ import annotations

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import pytest

from kosmic.visualisation.de.term_links import (
    create_term_chord, create_term_grid, order_genes, select_term_genes,
)


def _inputs():
    enr = pd.DataFrame({
        'Term': ['GO:0006099 tricarboxylic acid cycle', 'GO:0006631 fatty acid metabolic process',
                 'GO:0000001 not significant'],
        'FDR': [0.001, 0.01, 0.5],
        'Genes': ['CS, IDH2, MDH2, SHARED', ['ACADL', 'CPT1A', 'SHARED'], 'XYZ'],
    })
    de = pd.DataFrame({'names': ['CS', 'IDH2', 'MDH2', 'SHARED', 'ACADL', 'CPT1A', 'XYZ'],
                       'logfoldchanges': [-0.3, -0.5, -2.0, -1.0, 1.5, -0.2, 0.1]})
    return enr, de


def test_selects_passing_terms_and_top_genes():
    enr, de = _inputs()
    members, lfc, n_passing = select_term_genes(enr, de, n_terms=8, genes_per_term=3)
    assert n_passing == 2
    assert list(members) == ['tricarboxylic acid cycle', 'fatty acid metabolic process']
    assert members['tricarboxylic acid cycle'] == ['IDH2', 'MDH2', 'SHARED']
    assert set(lfc.index) == {'IDH2', 'MDH2', 'SHARED', 'ACADL', 'CPT1A'}


def test_falls_back_to_top_terms_when_none_pass():
    enr, de = _inputs()
    enr['FDR'] = 0.9
    members, _, n_passing = select_term_genes(enr, de, n_terms=1)
    assert n_passing == 0 and len(members) == 1


def test_gene_order():
    enr, de = _inputs()
    members, lfc, _ = select_term_genes(enr, de)
    assert order_genes(members, lfc, 'lfc') == sorted(lfc.index, key=lfc.get)
    by_term = order_genes(members, lfc, 'term')
    assert set(by_term[:4]) == set(members['tricarboxylic acid cycle'])


@pytest.mark.parametrize('make', [create_term_chord, create_term_grid])
def test_renders(make):
    enr, de = _inputs()
    members, lfc, _ = select_term_genes(enr, de)
    fig = make(members, lfc, title='t')
    assert fig is not None
    plt.close(fig)
    assert make({}, lfc) is None


def test_excludes_mitochondrial_genes():
    enr = pd.DataFrame({'Term': ['electron transport chain'], 'FDR': [0.001],
                        'Genes': ['MT-ND1, mt-Co1, NDUFS1']})
    de = pd.DataFrame({'names': ['MT-ND1', 'mt-Co1', 'NDUFS1'],
                       'logfoldchanges': [-5.0, -4.0, -0.5]})
    members, lfc, _ = select_term_genes(enr, de, exclude_groups=('mitochondrial',))
    assert members == {'electron transport chain': ['NDUFS1']}
    assert list(lfc.index) == ['NDUFS1']
