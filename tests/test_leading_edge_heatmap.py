"""Leading-edge heatmap: the genes that drive a GSEA result, by dataset.

From the Meta-Analysis Enrichment step the columns are each study and the
pooled estimate; from the DE workspace the single column is that study.
The leading edge is the GSEA run's own (pooled for a meta-analysis), and
GSEA matches genes in upper case, so it must be mapped back to the data's
spelling.
"""
from __future__ import annotations

import os
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from kosmic.visualisation.de.leading_edge import (
    create_leading_edge_heatmap, leading_edge_genes, leading_edge_table,
)


@pytest.fixture
def qapp():
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


RANKING = pd.Series({'Xbp1': 5.0, 'Hspa5': 3.0, 'P4hb': 4.0, 'Axin2': -6.0,
                     'Lifr': -2.0, 'Tcf3': -4.0})


def test_leading_edge_is_mapped_to_the_data_spelling_and_ordered():
    """GSEA reports 'XBP1' for mouse 'Xbp1'; the strongest gene comes first
    in the direction of the set's NES."""
    assert leading_edge_genes('HSPA5;XBP1;P4HB', RANKING, nes=2.1) == ['Xbp1', 'P4hb', 'Hspa5']
    assert leading_edge_genes('LIFR;AXIN2;TCF3', RANKING, nes=-1.8) == ['Axin2', 'Tcf3', 'Lifr']
    assert leading_edge_genes('HSPA5;XBP1;P4HB;NOTINDATA', RANKING, nes=2.1, n=2) == ['Xbp1', 'P4hb']


def _de(lfc, p, names=('Xbp1', 'Hspa5', 'P4hb')):
    return pd.DataFrame({'names': list(names), 'logfoldchanges': lfc, 'pvals_adj': p})


def test_table_marks_untested_genes_as_missing():
    a = _de([1.0, 0.5, 0.2], [0.01, 0.2, 0.04])
    b = _de([0.8, 0.1], [0.03, 0.5], names=('Xbp1', 'Hspa5'))          # no P4hb
    pooled = pd.DataFrame({'names': ['Xbp1', 'Hspa5', 'P4hb'],
                           'logfoldchanges': [0.9, 0.3, 0.2], 'fdr': [0.001, 0.3, 0.2]})
    lfc, padj = leading_edge_table(['Xbp1', 'P4hb'],
                                   [('A', a, 'pvals_adj'), ('B', b, 'pvals_adj'),
                                    ('Pooled', pooled, 'fdr')])
    assert list(lfc.columns) == ['A', 'B', 'Pooled']
    assert np.isnan(lfc.loc['P4hb', 'B']) and lfc.loc['Xbp1', 'Pooled'] == 0.9
    assert padj.loc['P4hb', 'A'] == 0.04


def test_heatmap_draws_one_block_per_set():
    lfc, padj = leading_edge_table(['Xbp1', 'Hspa5'], [('A', _de([1.0, 0.5, 0.2], [0.01, 0.2, 0.04]), 'pvals_adj')])
    fig = create_leading_edge_heatmap(
        [{'title': 'SET ONE (NES +2.1)', 'lfc': lfc, 'padj': padj},
         {'title': 'SET TWO (NES -1.8)', 'lfc': lfc, 'padj': padj}])
    images = [ax for ax in fig.axes if ax.images]
    assert len(images) == 2
    plt.close(fig)
    assert create_leading_edge_heatmap([]) is None


# -- the Figure Export page --------------------------------------------------

def _gsea():
    return pd.DataFrame({'names': ['ER_PROCESSING', 'WNT'], 'nes': [2.1, -1.8],
                         'pvals_adj': [0.002, 0.08],
                         'leading_edge': ['XBP1;P4HB;HSPA5', 'AXIN2;LIFR;TCF3']})


def _workspace(study=True, meta=True):
    genes = list(RANKING.index)
    de = None
    if study:
        de = SimpleNamespace(
            de_results=pd.DataFrame({'names': genes, 'logfoldchanges': RANKING.values / 4,
                                     'pvals_adj': 0.01}),
            enrichment_results=None, gsea_library_results=_gsea(),
            gsea_library_details=None, gsea_library_name='KEGG',
            disease_label='MCAO', control_label='Sham', project_dir='C:/x/GSE174574',
            h5ad_path=None, current_adata=None)
    pooled = pd.DataFrame({'names': genes, 'logfoldchanges': RANKING.values / 5,
                           'z_stat': RANKING.values, 'fdr': 0.01})
    studies = [(f'S{i}', pd.DataFrame({'names': genes[:-1] if i == 2 else genes,
                                       'logfoldchanges': RANKING.values[:5 if i == 2 else 6] / 4,
                                       'pvals_adj': 0.2})) for i in (1, 2)]
    inputs = {'gene_results': pooled if meta else None, 'pathway_results': None,
              'gene_studies': studies if meta else [], 'pathway_studies': [],
              'selection': 'Endothelial_Cell', 'project_folder': None,
              'enrichment_results': None,
              'gsea_library_results': _gsea() if meta else None,
              'gsea_library_details': None, 'gsea_library_name': 'KEGG',
              'enrichment_ranking': pooled if meta else None}
    from tests.study_folder import figure_workspace
    return figure_workspace(de=de, meta=SimpleNamespace(figure_inputs=lambda: inputs),
                            name='GSE174574')


def _page(ws):
    from kosmic.gui.figure_export.pages.leading_edge import LeadingEdgeHeatmapPage
    page = LeadingEdgeHeatmapPage(ws)
    page._controls.blockSignals(True)
    page._sync_sources()
    page._sync_sets()
    return page


def test_meta_source_has_a_column_per_study_and_pooled(qapp):
    page = _page(_workspace(study=False))
    assert page._controls.ticked() == ['ER_PROCESSING', 'WNT']    # top two, by FDR
    blocks = page._blocks()
    assert list(blocks[0]['lfc'].columns) == ['S1', 'S2', 'Pooled']
    assert list(blocks[0]['lfc'].index) == ['Xbp1', 'P4hb', 'Hspa5']
    assert list(blocks[1]['lfc'].index) == ['Axin2', 'Tcf3', 'Lifr']
    assert np.isnan(blocks[1]['lfc'].loc['Tcf3', 'S2'])             # S2 did not test Tcf3
    plt.close(page._make_render_func()())


def test_study_source_has_one_column(qapp):
    page = _page(_workspace(meta=False))
    assert not page._is_meta()
    blocks = page._blocks()
    assert list(blocks[0]['lfc'].columns) == ['GSE174574']
    plt.close(page._make_render_func()())


def test_export_writes_the_numbers(qapp, tmp_path):
    page = _page(_workspace(study=False))
    out = page._export_companions(tmp_path / 'fig.png')
    table = pd.read_csv(out[0])
    assert set(table.columns) == {'gene_set', 'gene', 'column', 'log2fc', 'padj'}
    assert len(table) == 2 * 3 * 3            # sets x genes x columns


def test_no_gsea_run_means_unavailable(qapp):
    page = _page(_workspace(study=False, meta=False))
    assert not page.dependencies_met()
    assert page._make_render_func() is None
