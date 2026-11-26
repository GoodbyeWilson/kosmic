"""The meta-analysis methods of the comparison and the standard table (ADR-017 §3, §5)."""
import json

import numpy as np
import pandas as pd
import pytest

from kosmic.comparison import ma
from kosmic.comparison.standard import COLUMNS, to_standard

SEL = 'decontX_counts/Ventricular_Cardiomyocyte'
NAMES = [f'S{i}_Ventricular_Cardiomyocyte' for i in range(1, 6)]


def _study_table(seed, n=40):
    rng = np.random.default_rng(seed)
    lfc = rng.normal(0, 0.05, n)
    lfc[:8] += 1.5                     # eight genes up in every study
    p = np.clip(np.exp(-np.abs(lfc) * 6), 1e-12, 1)
    return pd.DataFrame({'names': [f'G{i}' for i in range(n)], 'logfoldchanges': lfc,
                         'se': np.full(n, 0.2), 'pvals': p,
                         'pvals_adj': np.minimum(p * n / np.arange(1, n + 1), 1),
                         'pct_disease': 0.5, 'pct_control': 0.5})


@pytest.fixture
def tables():
    return [(n, _study_table(i)) for i, n in enumerate(NAMES)]


def _save_run(project, keys, params, df=None):
    folder = project / 'meta_analysis' / SEL
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"consensus_{'_'.join(keys)}_analytical.csv"
    (df if df is not None else pd.DataFrame({'names': ['G0'], 'pvals_pooled': [0.01]})
     ).to_csv(path, index=False)
    path.with_suffix('.json').write_text(json.dumps({'parameters': {
        'method_keys': keys, 'study_names': NAMES, **params}}), encoding='utf-8')
    return path


PARAMS = {'min_studies': 5, 'independent_filter': {'active': True, 'threshold': 0.05,
                                                   'stat': 'mean_max_pct'}}


@pytest.mark.parametrize('key', ['reml', 'sumrank', 'stouffer', 'vote_count'])
def test_each_method_pools_and_converts_to_the_standard_table(tables, key):
    settings = ma.MASettings(min_studies=5, filter_active=False, filter_threshold=None,
                             filter_stat=None)
    mid = {v: k for k, v in ma.MA_METHODS.items()}[key]
    std = to_standard(mid, ma.pool(tables, key, settings))
    assert list(std.columns) == COLUMNS and (std['method'] == mid).all()
    assert std['gene'].is_unique and std['tested'].all()
    top = std.sort_values('rank').head(8)
    assert set(top['gene']) == {f'G{i}' for i in range(8)} and (top['direction'] == 1).all()
    assert std.set_index('gene').loc[top['gene'], 'called'].all()
    assert std['rank'].min() == 1 and std['rank'].is_unique
    if key == 'vote_count':
        assert std['pvalue'].isna().all()
    else:
        top = std.sort_values('rank').iloc[0]
        assert top['pvalue'] == std['pvalue'].min()


def test_the_merged_adapter_leaves_filtered_genes_untested():
    df = pd.DataFrame({'names': ['A', 'B', 'C'], 'logfoldchanges': [1.0, -0.5, 2.0],
                       'pvals': [1e-6, 0.2, np.nan], 'pvals_adj': [1e-5, 0.3, np.nan],
                       'statistic': [5.0, -1.2, np.nan]})
    std = to_standard('pb_merged', df)
    row = std.set_index('gene')
    assert row.loc['A', 'called'] and not row.loc['B', 'called']
    assert not row.loc['C', 'tested'] and pd.isna(row.loc['C', 'rank'])
    assert row.loc['B', 'direction'] == -1 and row.loc['A', 'rank'] == 1


def test_sumrank_takes_its_named_direction():
    df = pd.DataFrame({'names': ['A', 'B'], 'logfoldchanges': [0.1, 0.1],
                       'direction': ['up', 'down'], 'rank_sum': [1.0, 2.0],
                       'pvals_pooled': [0.001, 0.002]})
    assert to_standard('ma_sumrank', df).set_index('gene')['direction'].to_dict() == \
        {'A': 1, 'B': -1}


def test_a_saved_single_method_run_with_the_same_settings_is_read(tmp_path, monkeypatch):
    saved = pd.DataFrame({'names': ['G0', 'G1'], 'logfoldchanges': [1.0, 0.1],
                          'z_stat': [4.0, 0.3], 'pvals_pooled': [1e-5, 0.7],
                          'fdr_reml': [2e-5, 0.7]})
    _save_run(tmp_path, ['reml'], PARAMS, saved)
    monkeypatch.setattr('kosmic.meta_analysis.saved.study_tables',
                        lambda *a, **k: pytest.fail('should not pool'))
    out = ma.load_ma(tmp_path, SEL, NAMES, ['ma_reml'])
    assert out['ma_reml'].source == 'saved'
    assert out['ma_reml'].table['fdr'].tolist() == [2e-5, 0.7]


def test_a_saved_run_of_other_studies_or_settings_is_pooled_again(tmp_path, tables,
                                                                   monkeypatch):
    monkeypatch.setattr('kosmic.meta_analysis.saved.study_tables',
                        lambda *a, **k: (tables, []))
    _save_run(tmp_path, ['reml'], {**PARAMS, 'min_studies': 3})        # older: other settings
    _save_run(tmp_path, ['reml', 'sumrank'], PARAMS)                    # newest: the settings
    import os
    import time
    later = time.time() + 5
    os.utime(tmp_path / 'meta_analysis' / SEL / 'consensus_reml_sumrank_analytical.csv',
             (later, later))
    out = ma.load_ma(tmp_path, SEL, NAMES, ['ma_reml', 'ma_sumrank'])
    assert out['ma_reml'].source == 'pooled' and 'other settings' in out['ma_reml'].note
    assert out['ma_sumrank'].source == 'pooled'
    assert out['ma_reml'].settings['min_studies'] == 5
    # A run of other studies is not one of this selection's.
    assert ma.saved_runs(tmp_path, SEL, NAMES[:4]) == []


def test_no_minimum_votes_equals_every_study(tmp_path, monkeypatch):
    _save_run(tmp_path, ['vote_count'], {**PARAMS, 'min_votes': 5},
              pd.DataFrame({'names': ['G0'], 'vote_pass': [True], 'z_stat': [5],
                            'n_tested': [5], 'logfoldchanges': [1.0]}))
    monkeypatch.setattr('kosmic.meta_analysis.saved.study_tables',
                        lambda *a, **k: pytest.fail('should not pool'))
    out = ma.load_ma(tmp_path, SEL, NAMES, ['ma_vote'])
    assert out['ma_vote'].source == 'saved'
    assert to_standard('ma_vote', out['ma_vote'].table)['called'].tolist() == [True]


def test_a_missing_study_table_is_an_error(tmp_path, monkeypatch):
    monkeypatch.setattr('kosmic.meta_analysis.saved.study_tables',
                        lambda *a, **k: ([], ['S1_Ventricular_Cardiomyocyte']))
    with pytest.raises(ValueError, match='S1_Ventricular'):
        ma.load_ma(tmp_path, SEL, NAMES, ['ma_reml'])
