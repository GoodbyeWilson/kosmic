"""Enrichment figure pages read the DE workspace or the Meta-Analysis
workspace's Enrichment step (ADR-010), defaulting to the meta-analysis."""
from __future__ import annotations

import os
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from kosmic.gui.figure_export.pages._enrichment import SOURCE_META, SOURCE_STUDY


@pytest.fixture
def qapp():
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _ora(term):
    return pd.DataFrame({'Term': [term, 'other term'], 'FDR': [0.001, 0.2],
                         'Fold_Enrichment': [3.0, 1.2], 'Gene_Count': [3, 2],
                         'Genes': ['A, B, C', 'D, E']})


def _gsea(term):
    res = pd.DataFrame({'names': [term, 'other set'], 'nes': [-2.1, 1.4],
                        'pvals_adj': [0.001, 0.04]})
    details = {term: {'nes': -2.1, 'RES': np.linspace(0, -0.6, 50).tolist(),
                      'hits': [3, 10, 20]},
               'other set': {'nes': 1.4, 'RES': np.linspace(0, 0.5, 50).tolist(),
                             'hits': [1, 5]}}
    return res, details


def _ranking(scale):
    return pd.DataFrame({'names': list('ABCDE'),
                         'logfoldchanges': scale * np.array([-1.0, 0.5, 2.0, -0.3, 1.0])})


def _workspace(study=True, meta=True):
    de = None
    if study:
        res, det = _gsea('study set')
        de = SimpleNamespace(
            enrichment_results=_ora('study term'), de_results=_ranking(1),
            gsea_library_results=res, gsea_library_details=det,
            gsea_library_name='Hallmark', disease_label='DCM', control_label='NF',
            project_dir=None, h5ad_path=None, current_adata=None)
    inputs = {'gene_results': None, 'pathway_results': None, 'gene_studies': [],
              'pathway_studies': [], 'selection': 'Ventricular_Cardiomyocyte',
              'project_folder': None, 'enrichment_results': None,
              'gsea_library_results': None, 'gsea_library_details': None,
              'gsea_library_name': None, 'enrichment_ranking': None}
    if meta:
        res, det = _gsea('meta set')
        inputs.update(enrichment_results=_ora('meta term'), enrichment_ranking=_ranking(2),
                      gsea_library_results=res, gsea_library_details=det,
                      gsea_library_name='Hallmark')
    return SimpleNamespace(de_ws=de, scrna_ws=None, data_version=lambda: 0,
                           meta_ws=SimpleNamespace(figure_inputs=lambda: inputs))


def _pages():
    from kosmic.gui.figure_export.pages.enriched_terms import EnrichedTermsPage
    from kosmic.gui.figure_export.pages.library_gsea import LibraryGseaPage
    from kosmic.gui.figure_export.pages.library_gsea_mountain import LibraryGseaMountainPage
    from kosmic.gui.figure_export.pages.term_links import TermLinksPage
    return [EnrichedTermsPage, TermLinksPage, LibraryGseaPage, LibraryGseaMountainPage]


def _ready(page):
    """Fill the Source (and gene-set) selectors as on_activated does, without
    starting the background preview render."""
    page._controls.blockSignals(True)
    page._sync_sources()
    if hasattr(page, '_sync_terms'):
        page._sync_terms()


@pytest.mark.parametrize('cls_index', range(4))
def test_both_sources_default_to_meta_and_switch(qapp, cls_index):
    page = _pages()[cls_index](_workspace())
    _ready(page)
    combo = page._controls.source
    assert [combo.itemText(i) for i in range(combo.count())] == [SOURCE_META, SOURCE_STUDY]
    assert not combo.isHidden()
    assert page._is_meta()
    plt.close(page._make_render_func()())

    combo.setCurrentText(SOURCE_STUDY)
    if hasattr(page, '_sync_terms'):
        page._sync_terms()
    assert not page._is_meta()
    plt.close(page._make_render_func()())


@pytest.mark.parametrize('study,meta', [(True, False), (False, True)])
def test_one_source_hides_the_selector(qapp, study, meta):
    for cls in _pages():
        page = cls(_workspace(study=study, meta=meta))
        _ready(page)
        assert page.dependencies_met()
        assert page._controls.source.isHidden()
        assert page._is_meta() == meta


def test_no_runs_means_unavailable(qapp):
    for cls in _pages():
        page = cls(_workspace(study=False, meta=False))
        _ready(page)
        assert not page.dependencies_met()
        assert page._make_render_func() is None


def test_meta_source_uses_pooled_fold_changes_and_titles(qapp):
    from kosmic.gui.figure_export.pages.term_links import TermLinksPage
    page = TermLinksPage(_workspace())
    _ready(page)
    src = page._source()
    assert src.de_results['logfoldchanges'].tolist() == _ranking(2)['logfoldchanges'].tolist()
    assert page._contrast(src) == 'Ventricular Cardiomyocyte (pooled)'
    assert page._default_export_basename().startswith('meta_')


def test_figure_inputs_carry_the_enrichment_step(qapp):
    from kosmic.gui.meta_analysis.workspace import MetaAnalysisWorkspace
    ws = MetaAnalysisWorkspace()
    ctx = ws._enrichment_page.context
    ctx.enrichment_results = _ora('t')
    ctx.de_results = _ranking(1)
    inputs = ws.figure_inputs()
    assert inputs['enrichment_results'] is ctx.enrichment_results
    assert inputs['enrichment_ranking'] is ctx.de_results
    assert inputs['gsea_library_results'] is None
