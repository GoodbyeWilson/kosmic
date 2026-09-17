"""Cluster-level quality flags (kosmic.scrna.annotate.cluster_qc).

Each rule is exercised on its own; the thresholds themselves were set on
the five DCM studies, where the flags pick out exactly the doublet and
low-quality clusters and nothing else.
"""
import anndata as ad
import numpy as np
import pandas as pd

from kosmic.scrna.annotate.cluster_qc import doublet_score_column, flag_clusters


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
    assert f['0'].text().startswith('weak match (0.80)')
    assert 'runner-up close (Other, margin 0.06)' in f['1'].text()
    assert not f['2'].flagged


def test_stress_markers():
    a = _study()
    names = pd.DataFrame(a.uns['rank_genes_groups']['names'])
    names['2'] = ['MT-CO1', 'MT-ND4', 'RPL13', 'RPS6', 'MT-CYB', 'G1', 'G2', 'G3', 'G4', 'G5']
    a.uns['rank_genes_groups'] = {'names': names.to_records(index=False)}
    f = flag_clusters(a)
    assert f['2'].text() == '5 of top 10 markers are MT/ribosomal' and not f['1'].flagged


def test_doublet_score_and_gene_count_support():
    a = _study()
    a.obs['doublet_score'] = np.r_[np.full(100, 0.05), np.full(100, 0.05), np.full(100, 0.30)]
    a.obs.loc[a.obs['leiden'] == '2', 'n_genes_by_counts'] = 1500     # 1.5x the other Fib
    assert doublet_score_column(a.obs) == 'doublet_score'
    f = flag_clusters(a)
    assert 'doublet score 6.0x study median (doublet_score)' in f['2'].text()
    assert '1.5x the genes of other Fib' in f['2'].text()
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
    assert f['1'].text() == "authors' labels mixed (top 55%)"
    assert f['2'].text() == 'authors left most cells unlabelled'


def test_unclustered_cells_are_ignored():
    a = _study()
    a.obs['leiden'] = pd.Categorical(list(a.obs['leiden'].astype(str))[:-50] + [np.nan] * 50)
    f = flag_clusters(a, details=_details(**{'0': (0.95, 0.5), '1': (0.96, 0.4), '2': (0.9, 0.3)}))
    assert list(f) == ['0', '1', '2'] and not any(x.flagged for x in f.values())
