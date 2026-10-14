"""Cluster-level quality flags (kosmic.scrna.annotate.cluster_qc).

Each rule is exercised on its own; the thresholds themselves were set on
the five DCM studies, where the flags pick out exactly the doublet and
low-quality clusters and nothing else.
"""
import anndata as ad
import numpy as np
import pandas as pd

from kosmic.scrna.annotate.cluster_qc import (
    doublet_score_column, flag_clusters, prediction_agreement,
)
from kosmic.scrna.annotate.labels import cluster_label_summary


def _study(n_per=(100, 100, 100)):
    rows = []
    for c, n in enumerate(n_per):
        rows += [str(c)] * n
    n = len(rows)
    obs = pd.DataFrame({
        'leiden': pd.Categorical(rows),
        'cell_type': pd.Categorical(['CM'] * n_per[0] + ['Fib'] * n_per[1] + ['Fib'] * n_per[2]),
        'n_genes_by_counts': np.r_[np.full(n_per[0], 2000), np.full(n_per[1], 1000), np.full(n_per[2], 1000)],
    }, index=[f'c{i}' for i in range(n)])
    a = ad.AnnData(obs=obs)
    a.uns['rank_genes_groups'] = {'names': pd.DataFrame({
        '0': ['TTN', 'RYR2', 'MYH7'] + [f'G{i}' for i in range(7)],
        '1': ['COL1A1', 'DCN'] + [f'G{i}' for i in range(8)],
        '2': ['COL1A1'] + [f'G{i}' for i in range(9)],
    }).to_records(index=False)}
    return a


def _details(**per_cluster):
    d = {}
    for c, (score, margin) in per_cluster.items():
        d[c] = {'best_score': score, 'margin': margin, 'runner_up': 'Other'}
    return d


def test_clean_study_has_no_flags():
    a = _study()
    f = flag_clusters(a, details=_details(**{'0': (0.95, 0.5), '1': (0.96, 0.4), '2': (0.9, 0.3)}))
    assert list(f) == ['0', '1', '2'] and not any(x.flagged for x in f.values())


def test_weak_score_and_close_runner_up():
    a = _study()
    f = flag_clusters(a, details=_details(**{'0': (0.80, 0.5), '1': (0.96, 0.06), '2': (0.9, 0.3)}))
    assert f['0'].text().startswith('Weak reference match (correlation 0.80)')
    assert 'Close second type: Other (correlation margin 0.06)' in f['1'].text()
    assert not f['2'].flagged


def test_stress_markers():
    a = _study()
    names = pd.DataFrame(a.uns['rank_genes_groups']['names'])
    names['2'] = ['MT-CO1', 'MT-ND4', 'RPL13', 'RPS6', 'MT-CYB', 'G1', 'G2', 'G3', 'G4', 'G5']
    a.uns['rank_genes_groups'] = {'names': names.to_records(index=False)}
    f = flag_clusters(a)
    assert f['2'].text() == ('Stressed or damaged: 5 of its top 10 marker genes are '
                            'mitochondrial or ribosomal') and not f['1'].flagged


def test_doublet_score_and_gene_count_support():
    a = _study()
    a.obs['doublet_score'] = np.r_[np.full(100, 0.05), np.full(100, 0.05), np.full(100, 0.30)]
    a.obs.loc[a.obs['leiden'] == '2', 'n_genes_by_counts'] = 1500     # 1.5x the other Fib
    assert doublet_score_column(a.obs) == 'doublet_score'
    f = flag_clusters(a)
    assert 'Possible doublets: doublet score 6.0x the dataset median (doublet_score)' in f['2'].text()
    assert '1.5x the genes per cell of other Fib clusters' in f['2'].text()
    assert not f['0'].flagged                                          # 2000 genes but CM: fine


def test_gene_count_alone_never_flags():
    a = _study()
    a.obs.loc[a.obs['leiden'] == '2', 'n_genes_by_counts'] = 3000
    assert not flag_clusters(a)['2'].flagged


def test_author_labels_unlabelled_and_mixed():
    a = _study()
    author = ['cardiomyocyte'] * 100 + ['fibroblast'] * 55 + ['pericyte'] * 45 + ['unknown'] * 80 + ['fibroblast'] * 20
    a.obs['cell_type_author'] = pd.Categorical(author)
    f = flag_clusters(a)
    assert not f['0'].flagged
    assert f['1'].text() == "Authors' labels split: the most common covers 55%"
    assert f['2'].text() == 'Authors left most of its cells unlabelled'


def test_unclustered_cells_are_ignored():
    a = _study()
    a.obs['leiden'] = pd.Categorical(list(a.obs['leiden'].astype(str))[:-50] + [np.nan] * 50)
    f = flag_clusters(a, details=_details(**{'0': (0.95, 0.5), '1': (0.96, 0.4), '2': (0.9, 0.3)}))
    assert list(f) == ['0', '1', '2'] and not any(x.flagged for x in f.values())


# --- agreement with each cell's own prediction ------------------------------

def _voted(preds, label):
    """One cluster whose cells all carry 'label', as after majority voting."""
    n = len(preds)
    return pd.DataFrame({'leiden': pd.Categorical(['0'] * n),
                         'cell_type': [label] * n,
                         'predicted_labels': preds},
                        index=[f'c{i}' for i in range(n)])


def test_agreement_is_read_from_the_predictions_not_the_voted_label():
    # The T/NK case: one CD4 subtype wins the vote at 30%, but CD8 cells
    # outnumber CD4 cells once their subtypes are counted together.
    preds = (['CD4+T_naive'] * 30 + ['CD8+T_trans'] * 25 + ['CD8+T_em'] * 20
             + ['NK_CD16hi'] * 25)
    a = prediction_agreement(_voted(preds, 'CD4+ T Cell'), 'leiden')['0']
    assert a['agreement'] == 0.30
    assert (a['runner_up'], a['runner_up_share']) == ('CD8+ T Cell', 0.45)
    assert a['confidence'] == 'Below threshold'


def test_fine_labels_are_compared_as_they_are():
    preds = ['vCM1'] * 90 + ['vCM2'] * 10
    a = prediction_agreement(_voted(preds, 'vCM1'), 'leiden')['0']
    assert a['agreement'] == 0.9 and a['confidence'] == 'High'


def test_a_hand_named_cluster_is_left_out():
    preds = ['CD4+T_naive'] * 50 + ['NK_CD16hi'] * 50
    assert prediction_agreement(_voted(preds, 'T/NK Cell'), 'leiden') == {}


def test_no_predictions_no_agreement():
    obs = _voted(['vCM1'] * 10, 'vCM1').drop(columns='predicted_labels')
    assert prediction_agreement(obs, 'leiden') == {}


def test_a_mixed_cluster_is_flagged_on_agreement_not_on_score():
    # A CellTypist probability of 0.7 is not a weak reference correlation.
    preds = ['Mast'] * 60 + ['B'] * 40
    obs = _voted(preds, 'Mast Cell')
    a = ad.AnnData(obs=obs)
    details = {'0': {'best_score': 0.7,
                     **prediction_agreement(obs, 'leiden')['0']}}
    text = flag_clusters(a, details=details)['0'].text()
    assert text == 'Mixed: 60% of cells predicted as its label (next: B Cell 40%)'

    clean = _voted(['Mast'] * 95 + ['B'] * 5, 'Mast Cell')
    details = {'0': {'best_score': 0.7,
                     **prediction_agreement(clean, 'leiden')['0']}}
    assert not flag_clusters(ad.AnnData(obs=clean), details=details)['0'].flagged


def test_the_record_keeps_the_agreement():
    preds = ['Mast'] * 60 + ['B'] * 40
    obs = _voted(preds, 'Mast Cell')
    details = prediction_agreement(obs, 'leiden')
    assert cluster_label_summary(obs, 'leiden', details)['0']['agreement'] == 0.6


def test_each_cell_is_compared_with_its_own_label():
    # Two populations merged into one cluster at a lower resolution, each
    # correctly labelled: the cluster is not a mix of wrong labels.
    obs = pd.DataFrame({
        'leiden': pd.Categorical(['0'] * 100),
        'cell_type': ['Lymphatic Endothelial Cell'] * 60 + ['Adipocyte'] * 40,
        'predicted_labels': ['EC8_ln'] * 55 + ['Meso'] * 5 + ['Adip1'] * 40,
    }, index=[f'c{i}' for i in range(100)])
    a = prediction_agreement(obs, 'leiden')['0']
    assert a['agreement'] == 0.95
    assert (a['runner_up'], a['runner_up_share']) == ('Mesothelial Cell', 0.05)
