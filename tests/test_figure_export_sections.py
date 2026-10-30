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
    assert VoteLevelsPage.TITLE in [t for t, _g, _tip in meta]
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


# -- Open in Figures (ADR-014) ---------------------------------------------------

def test_the_study_label_names_the_study(qapp, tmp_path):
    from kosmic.gui.figure_export.workspace import FigureExportWorkspace
    ws = FigureExportWorkspace()
    assert ws._study_label.text() == "No study selected"
    (tmp_path / 'Chaffin').mkdir()
    ws.set_study(tmp_path / 'Chaffin')
    assert ws._study_label.text() == "Study: Chaffin"
    assert ws._study_label.toolTip() == str(tmp_path / 'Chaffin')
    ws.source.wait()


def test_the_project_page_offers_open_in_figures(qapp):
    from kosmic.gui.project.workspace import ProjectWorkspace
    from kosmic.gui.shared.widgets import SecondaryButton
    ws = ProjectWorkspace()
    asked = []
    ws.open_in_requested.connect(lambda acc, mode: asked.append((acc, mode)))
    ws._details_accession = 'Chaffin'
    button = next(b for b in ws.findChildren(SecondaryButton) if b.text() == "Open in Figures")
    button.click()
    assert asked == [('Chaffin', 'figures')]
    ws.deleteLater()


_OPEN_IN_FIGURES = '''
import os, sys
from pathlib import Path
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
sys.path.insert(0, {root!r})
from PyQt6.QtCore import QSettings
from PyQt6.QtWidgets import QApplication
app = QApplication(sys.argv[:1])
import main
from tests.study_folder import make_study_folder, write_de_outputs

project = Path({project!r})
for name in ('StudyA', 'StudyB'):
    write_de_outputs(make_study_folder(project, name=name))
window = main.AppWindow()
# Shown: a hidden window defers every load, which would make the
# assertions below pass whatever Open in Figures did.
window.show()
# Not the user's settings: opening a project records it as recent.
window.settings = QSettings(str(project / 'settings.ini'), QSettings.Format.IniFormat)
window._on_project_open_requested(str(project))
window.set_active_study('StudyA')
app.processEvents()

window._on_open_study_in('StudyB', 'figures')
app.processEvents()
fw = window._figures_workspace
scrna, de = window._scrna_workspace, window._de_workspace
assert window.current_accession == 'StudyA', window.current_accession
assert window._stack.currentWidget() is fw
assert fw.source.study().name == 'StudyB' and fw.current_section() == 'scrna'
assert fw._study_label.text() == 'Study: StudyB'
fw.source.wait()
app.processEvents()
assert fw.source.scrna_adata() is not None          # read from the folder
assert scrna.current_adata is None and de.current_adata is None
assert scrna.download_tab._load_worker is None, 'scRNA started a load'
assert de.setup_page._load_worker is None, 'DE started a load'

# The navigation rail opens Figure Export on the active study again.
window._on_mode_selected('figures')
assert fw.source.study().name == 'StudyA'
fw.source.wait()
app.processEvents()
assert scrna.current_adata is None and de.current_adata is None
assert scrna.download_tab._load_worker is None and de.setup_page._load_worker is None
print('OPEN_OK', flush=True)
os._exit(0)
'''


def test_open_in_figures_loads_nothing_and_keeps_the_active_study(tmp_path):
    """Opening a project and then Figure Export, from the Project page's
    Open in Figures or the navigation rail, starts no scRNA or DE load."""
    import subprocess
    import sys
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script = _OPEN_IN_FIGURES.format(root=root, project=str(tmp_path / 'Project'))
    result = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True,
                            timeout=120, cwd=root)
    assert 'OPEN_OK' in result.stdout, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
