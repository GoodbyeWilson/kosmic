"""Sentinels for the per-study provenance sidecar.

Covers the record/load round-trip, the structural fingerprint token (what it
does and does not react to), and the three staleness states -- including the
'none' state that guarantees an absent sidecar never flags an old file.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

from kosmic import provenance


def _adata(n=40, g=10, roles=None):
    X = csr_matrix(np.random.default_rng(0).poisson(3, size=(n, g)).astype('float32'))
    obs = pd.DataFrame({'sample': [f's{i % 4}' for i in range(n)],
                        'condition': ['ctrl'] * (n // 2) + ['dis'] * (n - n // 2)})
    if roles is not None:
        obs['_role'] = roles
    obs.index = [f'c{i}' for i in range(n)]
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f'G{i}' for i in range(g)]))


def test_record_load_roundtrip(tmp_path):
    a = _adata()
    fp = provenance.compute_fingerprint(a, 'sample', 'condition')
    provenance.record_stage(tmp_path, 'S1', 'qc', {'min_genes': 200}, fingerprint=fp)

    rec = provenance.load(tmp_path)
    assert rec['study'] == 'S1'
    assert rec['stages'][0]['stage'] == 'qc'
    assert rec['stages'][0]['params']['min_genes'] == 200
    assert rec['stages'][0]['produced_token'] == fp['token']
    assert rec['fingerprint']['token'] == fp['token']


def test_staleness_none_match_stale(tmp_path):
    a = _adata(roles=['disease'] * 20 + ['control'] * 20)
    # No sidecar yet -> 'none': an absent record must never flag an old file.
    assert provenance.staleness(tmp_path, a, 'sample', 'condition') == 'none'

    fp = provenance.compute_fingerprint(a, 'sample', 'condition')
    provenance.record_stage(tmp_path, 'S1', 'roles', {}, fingerprint=fp)
    assert provenance.staleness(tmp_path, a, 'sample', 'condition') == 'match'

    # Re-processed to fewer cells -> the data changed -> stale.
    a2 = _adata(n=30, roles=['disease'] * 15 + ['control'] * 15)
    assert provenance.staleness(tmp_path, a2, 'sample', 'condition') == 'stale'


def test_token_reflects_de_relevant_structure_only():
    a = _adata(roles=['disease'] * 20 + ['control'] * 20)
    base = provenance.compute_fingerprint(a, 'sample', 'condition')['token']

    # A decontX layer is a new count source DE can consume -> new token.
    a.layers['decontX_counts'] = a.X.copy()
    assert provenance.compute_fingerprint(a, 'sample', 'condition')['token'] != base

    # A different role split -> new token.
    a.obs['_role'] = ['disease'] * 10 + ['control'] * 30
    tok_layer = provenance.compute_fingerprint(a, 'sample', 'condition')['token']
    a.obs['_role'] = ['disease'] * 20 + ['control'] * 20
    # Back to the original split (with the extra layer) -> distinct from role change.
    assert provenance.compute_fingerprint(a, 'sample', 'condition')['token'] != tok_layer

    # A cluster column DE ignores must NOT change the token.
    a.obs['leiden'] = ['0'] * 40
    with_cluster = provenance.compute_fingerprint(a, 'sample', 'condition')['token']
    a2 = _adata(roles=['disease'] * 20 + ['control'] * 20)
    a2.layers['decontX_counts'] = a2.X.copy()
    assert provenance.compute_fingerprint(a2, 'sample', 'condition')['token'] == with_cluster


def test_render_and_write_methods(tmp_path):
    a = _adata(roles=['disease'] * 20 + ['control'] * 20)
    fp = provenance.compute_fingerprint(a, 'sample', 'condition')
    provenance.record_stage(tmp_path, 'S1', 'qc',
                            {'min_genes': 200, 'max_mt': 10.0}, fingerprint=fp)
    provenance.record_stage(tmp_path, 'S1', 'setup',
                            {'role_map': {'A': 'disease', 'B': 'control'}},
                            fingerprint=fp)

    text = provenance.render_methods(provenance.load(tmp_path))
    assert 'Study: S1' in text
    assert 'Quality control' in text and 'min_genes: 200' in text
    assert 'Sample / role setup' in text

    # methods.md is written next to the sidecar and matches the render.
    md = tmp_path / 'methods.md'
    assert md.is_file()
    assert md.read_text(encoding='utf-8') == text

    # An absent record renders empty, so callers fall back to their own text.
    assert provenance.render_methods(None) == ''


def test_subset_lineage_points_to_source(tmp_path):
    """A cell-type subset records a 'source' pointer, and load_lineage walks it
    back to the parent's full processing (oldest ancestor first)."""
    from kosmic.paths import processed_data_dir

    parent_dir = processed_data_dir(tmp_path / 'FullStudy')
    a = _adata(roles=['disease'] * 20 + ['control'] * 20)
    fp = provenance.compute_fingerprint(a, 'sample', 'condition')
    provenance.record_stage(parent_dir, 'FullStudy', 'load',
                            {'n_cells': 40}, fingerprint=fp)
    provenance.record_stage(parent_dir, 'FullStudy', 'qc',
                            {'min_genes': 200}, fingerprint=fp)
    parent_h5ad = parent_dir / 'FullStudy.h5ad'

    sub_dir = processed_data_dir(tmp_path / 'EndoSubset')
    provenance.record_stage(sub_dir, 'EndoSubset', 'subset',
                            {'cell_type_column': 'cell_type',
                             'cell_types': ['Endothelial'],
                             'n_cells_before': 40, 'n_cells_after': 12},
                            source=str(parent_h5ad))

    chain = provenance.load_lineage(sub_dir)
    assert [rec.get('study') for _d, rec in chain] == ['FullStudy', 'EndoSubset']

    # The readable render names the parent it was derived from.
    text = provenance.render_methods(provenance.load(sub_dir))
    assert 'Derived from: FullStudy.h5ad' in text
    assert 'Cell-type subset' in text


def test_combined_cross_study_methods(tmp_path):
    """The meta view aggregates each study's provenance plus the pooling step."""
    from kosmic.paths import processed_data_dir, meta_output_dir
    from kosmic.meta_analysis.io import build_combined_methods

    for study in ('StudyA', 'StudyB'):
        a = _adata(roles=['disease'] * 20 + ['control'] * 20)
        fp = provenance.compute_fingerprint(a, 'sample', 'condition')
        pdir = processed_data_dir(tmp_path / study)
        provenance.record_stage(pdir, study, 'qc', {'min_genes': 200}, fingerprint=fp)
        provenance.record_stage(pdir, study, 'gene_de', {'method': 'deseq2'},
                                consumed_token=fp['token'])

    # Meta pooling step lives in the project-level meta sidecar.
    provenance.record_stage(meta_output_dir(tmp_path), 'Proj', 'meta_gene',
                            {'pooling_methods': ['reml'], 'n_studies': 2})

    text = build_combined_methods(tmp_path, ['StudyA', 'StudyB'])
    assert 'Study: StudyA' in text and 'Study: StudyB' in text
    assert 'Gene differential expression' in text
    assert 'Meta-analysis' in text and 'pooling_methods: reml' in text

    # No provenance anywhere -> empty, so the viewer can fall back to a notice.
    assert build_combined_methods(tmp_path / 'nope') == ''


def test_combined_methods_flags_changed_study(tmp_path):
    """A study reprocessed after the meta ran is flagged stale in the view."""
    from kosmic.paths import processed_data_dir, meta_output_dir
    from kosmic.meta_analysis.io import build_combined_methods

    a = _adata(roles=['disease'] * 20 + ['control'] * 20)
    fp = provenance.compute_fingerprint(a, 'sample', 'condition')
    pdir = processed_data_dir(tmp_path / 'StudyA')
    provenance.record_stage(pdir, 'StudyA', 'qc', {}, fingerprint=fp)
    # Meta ran against this study token.
    provenance.record_stage(meta_output_dir(tmp_path), 'Proj', 'meta_gene',
                            {'study_tokens': {'StudyA': fp['token']}})
    assert 'changed since the last meta-analysis' not in build_combined_methods(tmp_path)

    # Re-QC to fewer cells -> new token -> the view flags StudyA.
    a2 = _adata(n=30, roles=['disease'] * 15 + ['control'] * 15)
    fp2 = provenance.compute_fingerprint(a2, 'sample', 'condition')
    provenance.record_stage(pdir, 'StudyA', 'qc', {}, fingerprint=fp2)
    text = build_combined_methods(tmp_path)
    assert 'changed since the last meta-analysis' in text and 'StudyA' in text


def test_every_recorded_stage_has_a_readable_title():
    """A stage with no title renders as its raw key in methods.md.

    The GUI is the only thing that names stages, so this drifts silently:
    'gene_names', 'filter_dataset' and 'sample_conditions' were all being
    recorded with no entry here.
    """
    from kosmic.provenance import _STAGE_TITLES

    recorded_by_gui = {
        'load', 'qc', 'normalize', 'cluster', 'decontx', 'subset',
        'embedding', 'filter_dataset', 'gene_names', 'setup', 'roles',
        'sample_conditions', 'gene_de', 'pathway_de',
        'meta_gene', 'meta_pathway',
    }
    missing = sorted(recorded_by_gui - set(_STAGE_TITLES))
    assert not missing, f"stages with no readable title: {missing}"


def test_a_filter_stage_renders_its_counts(tmp_path):
    from kosmic import provenance

    provenance.record_stage(
        tmp_path, 'S1', 'filter_dataset',
        params={'column': 'cell_type', 'kept_values': ['EC'],
                'n_cells_before': 12000, 'n_cells_after': 900})
    text = provenance.render_methods(provenance.load(tmp_path))
    assert 'Pre-analysis cell filter' in text
    assert 'filter_dataset' not in text          # the key never leaks through
