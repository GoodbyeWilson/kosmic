"""Re-applying the QC minimums to DecontX-corrected counts.

QC passes a near-empty droplet on its original counts; DecontX can then
leave it almost nothing. On the DCM atlas 1.3% of cells fell below the
studies' own limits after correction and formed UMAP islands and a
mixed "Mast" cluster. The re-check removes them as QC removes cells.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy.sparse import csr_matrix

from kosmic import provenance
from kosmic.scrna.qc.decontx import DECONTX_LAYER, qc_limits, recheck_qc

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


def _study():
    """Six cells; the corrected counts leave two of them near-empty."""
    raw = np.full((6, 10), 100.0, dtype=np.float32)          # all pass QC
    corrected = np.full((6, 10), 80.0, dtype=np.float32)
    corrected[4] = [3, 0, 0, 0, 0, 0, 0, 0, 0, 0]          # emptied
    corrected[5, :2] = 400.0                               # 2 genes, 800+ counts
    corrected[5, 2:] = 0.0
    a = AnnData(X=csr_matrix(raw), obs=pd.DataFrame(
        {'cell_type': ['CM'] * 3 + ['EC'] * 3},
        index=[f'c{i}' for i in range(6)]))
    a.layers['counts'] = csr_matrix(raw)
    a.layers[DECONTX_LAYER] = csr_matrix(corrected)
    a.obsm['X_umap'] = np.arange(12, dtype=float).reshape(6, 2)
    return a


def test_cells_below_the_minimums_are_removed_with_everything_aligned():
    a = _study()
    r = recheck_qc(a, min_genes=5, min_counts=500)
    assert list(a.obs_names) == ['c0', 'c1', 'c2', 'c3']
    assert a.layers[DECONTX_LAYER].shape == a.layers['counts'].shape == (4, 10)
    assert a.obsm['X_umap'].shape == (4, 2)
    assert r == {'layer': DECONTX_LAYER, 'min_genes': 5, 'min_counts': 500,
                 'n_cells_before': 6, 'n_cells_after': 4, 'n_removed': 2,
                 'n_below_min_genes': 2, 'n_below_min_counts': 1}


def test_a_limit_of_zero_is_off():
    a = _study()
    r = recheck_qc(a, min_genes=0, min_counts=500)
    assert r['n_removed'] == 1 and 'c5' in a.obs_names     # c5 has 800 counts


def test_nothing_to_remove_leaves_the_study_alone():
    a = _study()
    r = recheck_qc(a, min_genes=1, min_counts=1)
    assert r['n_removed'] == 0 and a.n_obs == 6


def test_without_corrected_counts_it_refuses():
    a = _study()
    del a.layers[DECONTX_LAYER]
    with pytest.raises(ValueError, match='run DecontX first'):
        recheck_qc(a, min_genes=5)


def test_the_limits_come_from_the_recorded_qc(tmp_path):
    assert qc_limits(None) is None
    provenance.record_stage(tmp_path, 'S', 'qc', {
        'min_genes': 300, 'max_genes': 5000, 'min_counts': 500, 'max_mt': 3})
    assert qc_limits(provenance.load(tmp_path)) == {'min_genes': 300, 'min_counts': 500}
    provenance.record_stage(tmp_path, 'S', 'qc', {'min_genes': 0, 'min_counts': 0})
    assert qc_limits(provenance.load(tmp_path)) is None     # every limit off


def test_the_methods_name_the_step(tmp_path):
    provenance.record_stage(tmp_path, 'S', 'qc_recheck', {
        'layer': DECONTX_LAYER, 'min_genes': 300, 'min_counts': 500,
        'n_cells_before': 1000, 'n_cells_after': 985, 'n_removed': 15})
    text = provenance.render_methods(provenance.load(tmp_path))
    assert 'QC limits re-applied to the decontaminated counts' in text
    assert 'n_removed' in text


# --- the Decontaminate screen ----------------------------------------------

@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_the_screen_offers_the_recheck_alone_on_a_decontaminated_study(app, tmp_path):
    from kosmic.gui.scrna.tabs.decontx_tab import DecontXTab

    h5ad = tmp_path / 'S.h5ad'
    provenance.record_stage(tmp_path, 'S', 'qc', {'min_genes': 5, 'min_counts': 500})

    class MW:
        current_adata = _study()
        current_h5ad_path = str(h5ad)
        _adata_version = 1

    tab = DecontXTab(MW())
    tab.on_tab_activated()
    assert 'at least 5 genes, 500 counts' in tab._recheck_check.text()
    assert tab._recheck_check.isEnabled() and tab._recheck_check.isChecked()
    assert not tab._recheck_btn.isHidden()

    del MW.current_adata.layers[DECONTX_LAYER]
    tab._refresh_gate()
    assert tab._recheck_btn.isHidden()


def test_without_recorded_qc_the_option_is_off(app, tmp_path):
    from kosmic.gui.scrna.tabs.decontx_tab import DecontXTab

    class MW:
        current_adata = _study()
        current_h5ad_path = str(tmp_path / 'S.h5ad')
        _adata_version = 1

    tab = DecontXTab(MW())
    tab.on_tab_activated()
    assert not tab._recheck_check.isEnabled()
    assert tab._recheck_btn.isHidden()


def test_the_worker_rechecks_without_running_decontx(tmp_path):
    from kosmic.gui.scrna.tabs.decontx_tab import DecontXWorker

    w = DecontXWorker(_study(), None, None, output_path=str(tmp_path / 'S.h5ad'),
                      limits={'min_genes': 5, 'min_counts': 500})
    adata, info, recheck = w._run()
    assert info is None and recheck['n_removed'] == 2 and adata.n_obs == 4
    assert (tmp_path / 'S.h5ad').exists()


def test_a_finished_run_records_decontx_then_the_recheck(app, monkeypatch):
    from kosmic.gui.scrna.tabs import decontx_tab
    monkeypatch.setattr(decontx_tab.dialogs, 'info', lambda *a, **k: None)

    recorded = []

    class MW:
        current_adata = _study()
        current_h5ad_path = None
        _adata_version = 1

        def set_adata(self, adata, path, switch=False):
            self.current_adata = adata

        def record_provenance(self, stage, params, **kw):
            recorded.append((stage, params))

        def mark_step_complete(self, i):
            pass

    tab = decontx_tab.DecontXTab(MW())
    recheck = {'layer': DECONTX_LAYER, 'min_genes': 5, 'min_counts': 500,
               'n_cells_before': 6, 'n_cells_after': 4, 'n_removed': 2,
               'n_below_min_genes': 2, 'n_below_min_counts': 1}
    tab._on_finished((MW.current_adata, {'mean_contamination': 0.1, 'n_samples': 2}, recheck))
    assert [s for s, _ in recorded] == ['decontx', 'qc_recheck']
    assert recorded[1][1] == recheck
    assert '2 of 6 cells were below 5 genes or 500 counts' in tab._result_label.text()

    recorded.clear()
    tab._on_finished((MW.current_adata, None, recheck))     # the re-check alone
    assert [s for s, _ in recorded] == ['qc_recheck']
