"""Sex inference per sample (kosmic.scrna.inspect.sex).

Female samples carry XIST and no Y-gene counts, male samples the
reverse; a sample with both is a mixture, and a sample with neither
gets no call. The recorded-sex comparison normalises the depositor's
labels and reports mismatches, missing records and inconsistent records.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from kosmic.scrna.inspect.sex import (
    FEMALE, MALE, UNCLEAR, compare_with_recorded, infer_sex,
    normalise_sex_value, recorded_sex_column, sex_genes, write_inferred_sex,
)

GENES = ['GAPDH', 'ACTB', 'XIST', 'DDX3Y', 'UTY', 'KDM5D', 'TP53']


def _adata(samples, xist, y, n_cells=50, seed=0, obs_extra=None):
    """One block of cells per sample; xist / y give the per-cell mean count
    of XIST and of each Y gene for that sample."""
    rng = np.random.default_rng(seed)
    rows, obs = [], []
    for s, xv, yv in zip(samples, xist, y):
        block = rng.poisson(5, size=(n_cells, len(GENES))).astype(np.float32)
        block[:, 2] = rng.poisson(xv, size=n_cells)
        block[:, 3:6] = rng.poisson(yv, size=(n_cells, 3))
        rows.append(block)
        obs += [s] * n_cells
    df = pd.DataFrame({'donor': obs})
    if obs_extra:
        for k, per_sample in obs_extra.items():
            df[k] = [per_sample[s] for s in obs]
    df.index = [f'c{i}' for i in range(len(df))]
    return ad.AnnData(X=np.vstack(rows), obs=df, var=pd.DataFrame(index=GENES))


def test_clean_female_and_male_are_called():
    a = _adata(['F1', 'M1'], xist=[3.0, 0.0], y=[0.0, 3.0])
    t = infer_sex(a, 'donor', species='human')
    assert t.loc['F1', 'sex_inferred'] == FEMALE
    assert t.loc['M1', 'sex_inferred'] == MALE
    assert (t['sex_flag'] == '').all()
    assert t.loc['F1', 'y_cpm'] == 0.0 and t.loc['M1', 'xist_cpm'] == 0.0
    assert t['n_cells'].tolist() == [50, 50]


def test_mixed_signal_is_flagged_but_still_called():
    """A female sample with a quarter of a male's Y signal (TWCM-11-104)."""
    a = _adata(['F_mixed'], xist=[3.0], y=[0.3])
    t = infer_sex(a, 'donor', species='human')
    assert t.loc['F_mixed', 'sex_inferred'] == FEMALE
    assert t.loc['F_mixed', 'sex_flag'] == 'mixed signal'


def test_no_signal_gets_no_call():
    a = _adata(['S'], xist=[0.0], y=[0.0])
    t = infer_sex(a, 'donor', species='human')
    assert t.loc['S', 'sex_inferred'] == UNCLEAR
    assert t.loc['S', 'sex_flag'] == 'no signal'


def test_genes_absent_gives_empty_table():
    a = _adata(['F1'], xist=[3.0], y=[0.0])
    a = a[:, ['GAPDH', 'ACTB', 'TP53']].copy()
    assert infer_sex(a, 'donor', species='human').empty


def test_mouse_gene_names():
    xist, y = sex_genes(['Gapdh', 'Xist', 'Ddx3y', 'Uty', 'Kdm5d', 'Eif2s3y'])
    assert xist == 'Xist' and y == ['Ddx3y', 'Uty', 'Kdm5d', 'Eif2s3y']


def test_normalise_sex_labels():
    assert normalise_sex_value('F') == FEMALE
    assert normalise_sex_value('Male') == MALE
    assert normalise_sex_value('unknown') is None
    assert normalise_sex_value(None) is None


def test_recorded_column_found_by_name():
    obs = pd.DataFrame({'donor': ['a'], 'Gender': ['F']})
    assert recorded_sex_column(obs) == 'Gender'
    assert recorded_sex_column(pd.DataFrame({'donor': ['a']})) is None


def test_compare_with_recorded_reports_mismatch_and_missing():
    a = _adata(['F1', 'M1', 'M2'], xist=[3.0, 0.0, 0.0], y=[0.0, 3.0, 3.0],
               obs_extra={'sex': {'F1': 'female', 'M1': 'F', 'M2': 'unknown'}})
    t = compare_with_recorded(infer_sex(a, 'donor', species='human'), a.obs, 'donor')
    assert t.loc['F1', 'sex_check'] == ''
    assert t.loc['M1', 'sex_check'] == 'mismatch' and t.loc['M1', 'sex_recorded'] == FEMALE
    assert t.loc['M2', 'sex_check'] == 'not recorded'


def test_compare_without_any_recorded_column():
    a = _adata(['F1'], xist=[3.0], y=[0.0])
    t = compare_with_recorded(infer_sex(a, 'donor', species='human'), a.obs, 'donor')
    assert t.loc['F1', 'sex_check'] == 'not recorded'


def test_inconsistent_record_within_a_sample():
    a = _adata(['F1'], xist=[3.0], y=[0.0])
    a.obs['sex'] = ['female'] * 25 + ['male'] * 25
    t = compare_with_recorded(infer_sex(a, 'donor', species='human'), a.obs, 'donor')
    assert t.loc['F1', 'sex_check'] == 'inconsistent record'


def test_write_inferred_sex_is_per_cell():
    a = _adata(['F1', 'M1'], xist=[3.0, 0.0], y=[0.0, 3.0])
    write_inferred_sex(a, 'donor', infer_sex(a, 'donor', species='human'))
    assert a.obs['sex_inferred'].astype(str).tolist() == [FEMALE] * 50 + [MALE] * 50


@pytest.mark.parametrize('sparse', [True, False])
def test_sparse_and_dense_agree(sparse):
    import scipy.sparse as sp
    a = _adata(['F1', 'M1'], xist=[3.0, 0.0], y=[0.0, 3.0])
    if sparse:
        a.X = sp.csr_matrix(a.X)
    t = infer_sex(a, 'donor', species='human')
    assert t['sex_inferred'].tolist() == [FEMALE, MALE]
