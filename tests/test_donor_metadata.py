"""Per-donor metadata import (kosmic.scrna.inspect.donor_metadata).

A paper's supplementary table joined onto obs by sample id: the key
column is found by coverage, the match is reported rather than assumed,
only the chosen columns are written, unmatched samples get blanks, and
sex- and age-like columns are the ones suggested by default.
"""
import anndata as ad
import numpy as np
import pandas as pd

from kosmic.scrna.inspect.donor_metadata import (
    apply_metadata, column_summary, match_column, match_report, read_table,
    safe_obs_name, suggested_columns,
)


def _adata():
    obs = pd.DataFrame({'sample': ['D1'] * 3 + ['D2'] * 3 + ['D3'] * 2},
                       index=[f'c{i}' for i in range(8)])
    return ad.AnnData(X=np.zeros((8, 2), dtype=np.float32), obs=obs,
                      var=pd.DataFrame(index=['A', 'B']))


def _table():
    return pd.DataFrame({
        'Samples': ['D1', 'D2', 'D9'],
        'Age': ['9', '48', '60'],
        'Sex': ['male', 'female', 'female'],
        'Etiology of HF': ['NICM (familial)', 'No HF', 'No HF'],
        'BMI': ['16', '37.7', '24'],
    })


def test_read_table_csv_tsv_and_strip(tmp_path):
    t = _table()
    csv = tmp_path / 't.csv'
    t.to_csv(csv, index=False)
    tsv = tmp_path / 't.tsv'
    t.to_csv(tsv, sep='\t', index=False)
    for p in (csv, tsv):
        df = read_table(p)
        assert list(df.columns) == list(t.columns)
        assert df['Age'].tolist() == ['9', '48', '60']       # kept as text


def test_key_column_is_the_one_covering_sample_ids():
    assert match_column(_table(), ['D1', 'D2', 'D3']) == 'Samples'
    assert match_column(_table(), ['X', 'Y']) is None


def test_match_report_counts_both_directions():
    r = match_report(_table(), 'Samples', ['D1', 'D2', 'D3'])
    assert (r.n_samples, r.n_matched) == (3, 2)
    assert r.unmatched_samples == ['D3']
    assert r.unmatched_rows == ['D9']
    assert r.duplicate_keys == []


def test_match_report_flags_duplicate_keys():
    t = pd.concat([_table(), _table().iloc[[0]]])
    assert match_report(t, 'Samples', ['D1']).duplicate_keys == ['D1']


def test_suggested_columns_are_sex_and_age_only():
    assert suggested_columns(_table(), 'Samples') == ['Age', 'Sex']
    # already in obs -> not suggested again
    assert suggested_columns(_table(), 'Samples', existing=['sex']) == ['Age']


def test_column_summary():
    n, examples = column_summary(_table(), 'Sex')
    assert n == 2 and examples == 'male, female'


def test_safe_obs_name():
    assert safe_obs_name('Etiology of HF') == 'etiology_of_hf'
    assert safe_obs_name('Age', existing=['age']) == 'age_2'


def test_apply_writes_chosen_columns_per_cell_with_blanks_for_unmatched():
    a = _adata()
    written = apply_metadata(a, 'sample', _table(), 'Samples', ['Sex', 'Age'])
    assert written == {'sex': 'Sex', 'age': 'Age'}
    assert a.obs['sex'].astype(str).tolist() == ['male'] * 3 + ['female'] * 3 + [''] * 2
    assert a.obs['age'].astype(str).tolist() == ['9'] * 3 + ['48'] * 3 + [''] * 2
    assert 'bmi' not in a.obs.columns                       # not chosen, not written


def test_apply_honours_rename_and_first_duplicate_row():
    a = _adata()
    t = pd.concat([_table(), pd.DataFrame({'Samples': ['D1'], 'Age': ['99'], 'Sex': ['female'],
                                           'Etiology of HF': ['x'], 'BMI': ['1']})])
    written = apply_metadata(a, 'sample', t, 'Samples', ['Age'], rename={'Age': 'donor_age'})
    assert written == {'donor_age': 'Age'}
    assert a.obs.loc['c0', 'donor_age'] == '9'
