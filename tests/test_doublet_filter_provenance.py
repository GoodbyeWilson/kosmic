"""Filtering doublets is recorded in the provenance, not only detecting them.

Scrublet recorded a 'doublets' stage with the calls, but pressing Filter
Doublets recorded nothing, so a study's methods read as though doublets
were found and kept, and its data fingerprint still had the cell count
from before the filter.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PyQt6")


@pytest.fixture(scope="module")
def app():
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _study():
    obs = pd.DataFrame({
        'sample': ['s1'] * 4 + ['s2'] * 4,
        'predicted_doublet': [True, False, False, False, True, True, False, False],
    }, index=[f'c{i}' for i in range(8)])
    return ad.AnnData(X=np.ones((8, 3), dtype=np.float32), obs=obs)


def test_filtering_records_a_doublet_removal_stage(app, monkeypatch):
    from pathlib import Path

    from kosmic.gui.scrna.tabs import qc_tab
    monkeypatch.setattr(qc_tab.dialogs, 'info', lambda *a, **k: None)

    recorded = []

    class MW:
        current_adata = None
        current_h5ad_path = None
        _adata_version = 1

        def set_adata(self, adata, path=None, switch=False):
            self.current_adata = adata

        def record_provenance(self, stage, params, **kw):
            recorded.append((stage, params, kw))

        def mark_step_complete(self, i):
            pass

    tab = qc_tab.QCTab(MW())
    before = _study()
    tab.adata = before
    tab.h5ad_path = Path('study.h5ad')
    monkeypatch.setattr(tab, '_doublet_sample_col', lambda: 'sample')
    after = before[~before.obs['predicted_doublet']].copy()

    tab._on_filter_doublets_finished((after, 3))

    assert [s for s, _, _ in recorded] == ['doublet_filter']
    _, params, kw = recorded[0]
    assert params == {'n_cells_before': 8, 'n_cells_after': 5, 'n_removed': 3,
                      'removed_by_sample': {'s2': 2, 's1': 1}}
    # A data-changing stage: the fingerprint is refreshed.
    assert kw.get('data_changed', True) is True


def test_the_stage_has_a_title_in_the_methods():
    from kosmic.provenance import render_methods
    text = render_methods({'study': 'S', 'stages': [{
        'stage': 'doublet_filter', 'timestamp': '2026-10-26T09:00:00+00:00',
        'params': {'n_removed': 3}}]})
    assert 'Doublet removal' in text
