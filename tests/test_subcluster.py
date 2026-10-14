"""Subclustering one population of a study (ADR-008).

The study has three clusters; cluster '1' holds two planted fibroblast
states, B1 and B2, that the study-wide clustering merged.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from kosmic.scrna.cluster import subcluster
from kosmic.scrna.annotate.labels import set_cell_type

SETTINGS = {'n_hvg': 150, 'n_pcs': 15, 'n_neighbors': 15, 'resolution': 0.3}


def _study():
    rng = np.random.default_rng(0)
    g = 200
    base = rng.gamma(0.6, 1.0, g)
    kinds = {'A': 200, 'B1': 150, 'B2': 150, 'C': 150}
    rows, kind = [], []
    for k, n in kinds.items():
        p = base.copy()
        p[{'A': slice(0, 20), 'B1': slice(20, 40), 'B2': slice(20, 40),
           'C': slice(60, 80)}[k]] *= 12
        if k == 'B1':
            p[40:50] *= 10
        if k == 'B2':
            p[50:60] *= 10
        rows.append(rng.multinomial(1500, p / p.sum(), size=n))
        kind += [k] * n
    n = len(kind)
    cluster = {'A': '0', 'B1': '1', 'B2': '1', 'C': '2'}
    label = {'A': 'Cardiomyocyte', 'B1': 'Fibroblast', 'B2': 'Fibroblast',
             'C': 'Endothelial'}
    obs = pd.DataFrame({
        'kind': kind,
        'leiden': pd.Categorical([cluster[k] for k in kind], categories=['0', '1', '2']),
        'cell_type': pd.Categorical([label[k] for k in kind]),
        'sample': [f's{i % 3}' for i in range(n)],
    }, index=[f'c{i}' for i in range(n)])
    return ad.AnnData(X=sp.csr_matrix(np.vstack(rows).astype(np.float32)), obs=obs,
                      var=pd.DataFrame(index=[f'G{j}' for j in range(g)]))


@pytest.fixture(scope='module')
def focused():
    """The study and its cluster '1', embedded and clustered on its own."""
    a = _study()
    focus = subcluster.select_population(a, 'leiden', ['1'])
    subcluster.embed_and_cluster(focus, SETTINGS)
    return a, focus


def test_a_population_gets_its_own_clusters(focused):
    a, focus = focused
    sub = focus.adata
    purity = pd.crosstab(sub.obs['leiden'], sub.obs['kind']).max(axis=1).sum() / sub.n_obs
    assert sub.obs['leiden'].nunique() >= 2 and purity > 0.95
    assert (focus.parent == '1').all() and focus.prefix == '1'
    assert focus.settings['n_pcs'] == 15 and focus.settings['resolution'] == 0.3
    assert 'X_umap' in sub.obsm and 'rank_genes_groups' in sub.uns
    assert np.allclose(a.X.data, np.round(a.X.data))     # the study keeps its counts


def test_apply_writes_ids_and_only_the_labels_that_changed(focused):
    study, focus = focused
    a = study.copy()
    sub = focus.adata
    first = sub.obs['leiden'].cat.categories[0]
    renamed = (sub.obs['leiden'] == first).to_numpy()
    set_cell_type(sub.obs, renamed, 'Fib_state')
    edits = [{'cluster': first, 'from': 'Fibroblast', 'to': 'Fib_state',
              'n_cells': int(renamed.sum())}]
    p = subcluster.apply_focus(a, focus, manual_edits=edits)

    ids = a.obs['leiden_sub'].astype(str)
    assert (ids[a.obs['leiden'] != '1'] == a.obs['leiden'][a.obs['leiden'] != '1'].astype(str)).all()
    assert set(ids[a.obs['leiden'] == '1']) == {f'1.{c}' for c in sub.obs['leiden'].cat.categories}
    assert (a.obs.loc[sub.obs_names[renamed], 'cell_type'] == 'Fib_state').all()
    assert (a.obs.loc[sub.obs_names[~renamed], 'cell_type'] == 'Fibroblast').all()
    others = a.obs['leiden'] != '1'
    assert a.obs.loc[others, 'cell_type'].astype(str).equals(
        study.obs.loc[others, 'cell_type'].astype(str))

    assert p['column'] == 'leiden' and p['groups'] == ['1'] and p['n_cells'] == 300
    assert p['n_cells_relabelled'] == int(renamed.sum())
    assert p['annotation'] == {'method': 'manual'}
    assert p['manual_edits'][0]['cluster'] == f'1.{first}'
    assert p['subclusters'][f'1.{first}'] == {'label': 'Fib_state', 'n_cells': int(renamed.sum())}

    assert subcluster.clear_subclusters(a) and 'leiden_sub' not in a.obs
    set_cell_type(sub.obs, renamed, 'Fibroblast')           # leave the fixture as found


def test_groups_follow_the_column_order_and_record_their_cells(focused):
    study, _ = focused
    a = study.copy()
    focus = subcluster.select_population(a, 'leiden', ['2', '1'])
    assert focus.groups == ['1', '2'] and focus.prefix == '1+2'
    subcluster.embed_and_cluster(focus, SETTINGS)
    p = subcluster.apply_focus(a, focus)
    assert p['n_cells_by_group'] == {'1': 300, '2': 150}
    assert all(sum(s['from'].values()) == s['n_cells'] for s in p['subclusters'].values())
    assert all(i.startswith('1+2.') for i in a.obs.loc[a.obs['leiden'] != '0', 'leiden_sub'])


def test_excluded_and_too_few_cells(focused):
    study, _ = focused
    a = study.copy()
    a.obs['_role'] = 'disease'
    a.obs.loc[a.obs_names[200:205], '_role'] = 'exclude'      # five fibroblasts
    focus = subcluster.select_population(a, 'cell_type', ['Fibroblast'])
    assert focus.adata.n_obs == 295 and focus.n_excluded == 5
    with pytest.raises(ValueError, match='at least'):
        subcluster.select_population(a.copy()[:240], 'leiden', ['1'])


def test_reclustering_keeps_labels_and_drops_what_named_the_old_subclusters(focused):
    _, focus = focused
    sub = focus.adata
    sub.uns['manual_relabels'] = {'cluster': ['0']}
    sub.uns['cluster_annotation_details'] = {'0': {}}
    labels = sub.obs['cell_type'].copy()
    n = subcluster.recluster(focus, 1.0)
    assert n == sub.obs['leiden'].nunique() and focus.settings['resolution'] == 1.0
    assert 'manual_relabels' not in sub.uns and 'cluster_annotation_details' not in sub.uns
    assert sub.obs['cell_type'].equals(labels)
    subcluster.recluster(focus, SETTINGS['resolution'])


def test_defaults_come_from_the_record_then_the_study():
    a = _study()
    a.uns['neighbors'] = {'params': {'n_neighbors': 20, 'n_pcs': 25}}
    record = {'stages': [
        {'stage': 'embedding', 'params': {'n_hvg': 3000, 'harmony': True,
                                          'harmony_batch_key': 'sample'}},
        {'stage': 'cluster', 'params': {'resolution': 0.8}}]}
    d = subcluster.focus_defaults(a, record)
    assert (d['n_hvg'], d['n_pcs'], d['n_neighbors'], d['resolution']) == (3000, 25, 20, 0.8)
    assert d['harmony'] and d['harmony_batch_key'] == 'sample'
    assert subcluster.focus_defaults(a)['resolution'] == 0.5


def test_methods_list_every_subclustering(tmp_path):
    """Each subclustering refines a different population, so the methods
    document keeps them all rather than only the last."""
    from kosmic import provenance
    for groups in (['1'], ['2']):
        provenance.record_stage(tmp_path, 's', 'subcluster', {'groups': groups})
    provenance.record_stage(tmp_path, 's', 'cluster', {'resolution': 0.5})
    provenance.record_stage(tmp_path, 's', 'cluster', {'resolution': 1.0})
    text = provenance.render_methods(provenance.load(tmp_path))
    assert text.count('Subclustering') == 2 and text.count('Clustering') == 1


def test_apply_writes_the_fine_labels_focus_mode_gave(focused):
    """A CellTypist pass in focus mode writes 'cell_type_fine' on the subset;
    without copying it the study kept the parent cluster's winning label."""
    study, focus = focused
    a = study.copy()
    a.obs['cell_type_fine'] = pd.Categorical(['FB1'] * a.n_obs)
    sub = focus.adata
    half = sub.n_obs // 2
    sub.obs['cell_type_fine'] = pd.Categorical(
        ['FB4_activated'] * half + ['FB1'] * (sub.n_obs - half))
    subcluster.apply_focus(a, focus)
    got = a.obs.loc[sub.obs_names, 'cell_type_fine'].astype(str).to_numpy()
    assert (got == sub.obs['cell_type_fine'].astype(str).to_numpy()).all()
    others = ~a.obs_names.isin(sub.obs_names)
    assert (a.obs.loc[others, 'cell_type_fine'] == 'FB1').all()
    del sub.obs['cell_type_fine']                       # leave the fixture as found
