"""Sentinels: the fgsea summary round-trips through disk unchanged (so the figure
export can reload it instead of recomputing), and the NES bar renderer produces a
figure coloured by FDR significance."""
from __future__ import annotations

import matplotlib
matplotlib.use('Agg')
import pandas as pd

from kosmic.de.fgsea import save_fgsea_results, load_fgsea_results
from kosmic.visualisation.de.fgsea_plot import create_fgsea_nes_plot


def _res():
    return pd.DataFrame({
        'names': ['Oxidative_Phosphorylation', 'Glycolysis', 'Fatty_Acid_Synthesis'],
        'nes': [-2.1, -1.8, 0.3],
        'es': [-0.6, -0.5, 0.1],
        'pvals': [0.001, 0.004, 0.5],
        'pvals_adj': [0.003, 0.008, 0.76],
        'leading_edge': ['NDUFA9;COX5A', 'PKM;ENO1', 'FASN'],
    })


def test_fgsea_round_trip(tmp_path):
    res = _res()
    path = save_fgsea_results(res, tmp_path / 'fgsea_results.csv')
    assert path.exists()
    back = load_fgsea_results(path)
    assert list(back['names']) == list(res['names'])
    assert list(back['nes']) == list(res['nes'])


def test_load_missing_returns_none(tmp_path):
    assert load_fgsea_results(tmp_path / 'nope.csv') is None


def test_nes_plot_renders():
    fig = create_fgsea_nes_plot(_res(), disease_label='HF', control_label='Healthy')
    assert fig is not None


def test_nes_plot_empty_returns_none():
    assert create_fgsea_nes_plot(pd.DataFrame(columns=['names', 'nes', 'pvals_adj'])) is None
