"""Figure Export's list: one section (scRNA, DE, Meta) at a time, with the
figures whose data is not available greyed and their step as tooltip."""
from __future__ import annotations

import os
from types import SimpleNamespace

import pandas as pd
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


@pytest.fixture(scope='module')
def qapp():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _rows(ws):
    """(title, greyed, tooltip) of the figure rows, headers left out."""
    from PyQt6.QtCore import Qt
    out = []
    for i in range(ws._list.count()):
        item = ws._list.item(i)
        if item.data(Qt.ItemDataRole.UserRole) is None:
            continue
        out.append((item.text().strip(), bool(item.toolTip()), item.toolTip()))
    return out


def test_every_figure_belongs_to_one_section():
    from kosmic.gui.figure_export.workspace import _PAGE_REGISTRY, _SECTION_OF_GROUP
    assert {group for group, _cls in _PAGE_REGISTRY} <= set(_SECTION_OF_GROUP)


def test_sections_list_their_own_figures_greyed_without_data(qapp):
    from kosmic.gui.figure_export.pages.meta_figures import VoteLevelsPage
    from kosmic.gui.figure_export.workspace import FigureExportWorkspace

    ws = FigureExportWorkspace()
    assert ws.current_section() == 'scrna'
    scrna = [t for t, _g, _tip in _rows(ws)]
    assert 'UMAP / t-SNE' in scrna and 'Study Agreement (Vote Levels)' not in scrna

    ws.show_section('meta')
    meta = _rows(ws)
    assert [t for t, _g, _tip in meta][-1] == VoteLevelsPage.TITLE
    assert all(greyed for _t, greyed, _tip in meta)     # no meta-analysis run
    assert 'Run a gene-level meta-analysis' in meta[0][2]

    ws._available_only.setChecked(True)
    assert _rows(ws) == []
    ws.show_section(None)                                # unknown: stays
    assert ws.current_section() == 'meta'


def test_a_figure_with_data_is_not_greyed_and_is_opened(qapp):
    from kosmic.gui.figure_export.pages import _base
    from kosmic.gui.figure_export.pages.meta_figures import VoteLevelsPage
    from kosmic.gui.figure_export.workspace import FigureExportWorkspace

    results = pd.DataFrame({'names': ['A', 'B'], 'n_up': [2, 0], 'n_down': [0, 0],
                            'logfoldchanges': [1.0, 0.1], 'fdr': [0.01, 0.5]})
    inputs = {'gene_results': results, 'gene_studies': [('S1', None), ('S2', None)],
              'pathway_results': None, 'pathway_studies': [],
              'selection': 'x', 'project_folder': None}
    ws = FigureExportWorkspace(meta_workspace=SimpleNamespace(figure_inputs=lambda: inputs))
    ws._available_only.setChecked(True)
    ws.show_section('meta')
    titles = [t for t, greyed, _tip in _rows(ws)]
    assert VoteLevelsPage.TITLE in titles
    assert not any(greyed for _t, greyed, _tip in _rows(ws))
    # The open figure is the first one that can draw.
    assert ws._stack.currentWidget().TITLE == titles[0]
    _base.wait_for_renders()
