"""'Shared atlas genes only' on the Gene DE page restricts the BH
denominator and the results to the atlas's gene list."""
import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _project(tmp_path, atlas_genes):
    proj = tmp_path / 'proj'
    study = proj / 'S1' / 'processed_data'
    study.mkdir(parents=True)
    master = proj / '_master' / 'processed_data'
    master.mkdir(parents=True)
    a = ad.AnnData(X=sp.csr_matrix(np.ones((4, len(atlas_genes)), dtype=np.float32)),
                   obs=pd.DataFrame(index=[f'c{i}' for i in range(4)]),
                   var=pd.DataFrame(index=atlas_genes))
    a.write_h5ad(master / 'master.h5ad')
    return proj / 'S1'


@pytest.fixture
def workspace(app):
    """One DE workspace per test, destroyed before the next test starts.

    A DEWorkspace is hundreds of widgets and dozens of pyqtgraph plots.
    Left dangling, Python destroys it at an arbitrary later moment --
    on the Linux runner that was while a later test was building its own
    plots, and the process segfaulted.
    """
    from PyQt6.QtCore import QEvent
    from kosmic.gui.de_analysis.workspace import DEWorkspace
    ws = DEWorkspace()
    yield ws
    ws.close()
    ws.deleteLater()
    app.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    app.processEvents()


def _page(workspace, study_dir):
    workspace.project_dir = study_dir
    return workspace, workspace.gene_de_page


def test_option_reads_the_atlas_and_narrows_the_denominator(workspace, tmp_path):
    study_dir = _project(tmp_path, ['TTN', 'MYH7', 'DCN'])
    ws, page = _page(workspace, study_dir)
    page._refresh_atlas_genes_option()
    assert page.atlas_genes_check.isEnabled()
    assert '3' in page.atlas_genes_check.text()
    assert page._committed_fdr_genes() is None            # discovery, no universe
    page.atlas_genes_check.setChecked(True)
    assert ws.gene_universe == 'shared_atlas'
    assert page._committed_fdr_genes() == {'TTN', 'MYH7', 'DCN'}
    assert page._gene_universe_record() == 'shared atlas (3 genes)'
    page.atlas_genes_check.setChecked(False)
    assert ws.gene_universe is None and page._gene_universe_record() == 'all genes'


def test_option_is_disabled_without_an_atlas(workspace, tmp_path):
    study_dir = tmp_path / 'proj2' / 'S1'
    (study_dir / 'processed_data').mkdir(parents=True)
    ws, page = _page(workspace, study_dir)
    ws.gene_universe = 'shared_atlas'
    page._refresh_atlas_genes_option()
    assert not page.atlas_genes_check.isEnabled() and ws.gene_universe is None
