"""Sentinel: figure-export pages warn when the results they render predate the
current data. de_results_stale_note compares the provenance upstream token the
DE stage consumed against a live fingerprint of the loaded adata."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from kosmic import provenance
from kosmic.gui.figure_export.shared.staleness import de_results_stale_note


class _DE:
    pass


def _adata(seed, n_obs=8):
    rng = np.random.default_rng(seed)
    X = rng.integers(0, 30, size=(n_obs, 5)).astype(np.float32)
    return ad.AnnData(
        X=X, obs=pd.DataFrame(index=[f'c{i}' for i in range(n_obs)]),
        var=pd.DataFrame(index=[f'g{i}' for i in range(5)]))


def test_no_de_ws_or_no_results_is_not_stale():
    assert de_results_stale_note(None) is None
    de = _DE()
    de.h5ad_path = None
    de.current_adata = None
    assert de_results_stale_note(de) is None


def test_stale_when_data_changed_since_de_ran(tmp_path):
    study = tmp_path / 'study'
    study.mkdir()
    h5ad = study / 'study.h5ad'
    ran_on = _adata(0)
    provenance.record_stage(
        study, 'study', 'gene_de', {'method': 'deseq2'},
        consumed_token=provenance.compute_fingerprint(ran_on, h5ad_path=str(h5ad))['token'])

    de = _DE()
    de.h5ad_path = str(h5ad)
    # Current data is a different shape than gene_de consumed (e.g. re-QC
    # dropped cells) -> the structural fingerprint differs -> stale.
    de.current_adata = _adata(0, n_obs=6)
    de.de_results = pd.DataFrame({'names': ['g0']})
    note = de_results_stale_note(de)
    assert note and 'out of date' in note

    # Same data as the run -> not stale.
    de.current_adata = ran_on
    assert de_results_stale_note(de) is None
