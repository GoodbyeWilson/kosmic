"""Tests for kosmic.de.batch -- per-cell-type batch DE.

Covers the parts that decide what gets run and where it lands: slug
safety against the meta-analysis discovery regex, donor-count
eligibility (which must ignore 'exclude' cells and mirror the
min_cells drop create_pseudobulk applies), the pseudobulk CSV schema,
and the loop's behaviour when one cell type fails.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")

from kosmic.de import batch as B  # noqa: E402
from kosmic.meta_analysis.io import discover_de_results  # noqa: E402
from kosmic.paths import de_stats_dir  # noqa: E402


def _adata(n_genes=40, cells_per_donor=30, seed=0):
    """Two cell types, 3 disease + 3 control donors, plus an excluded arm."""
    rng = np.random.default_rng(seed)
    rows = []
    for role, donors in (('disease', ['D1', 'D2', 'D3']),
                         ('control', ['C1', 'C2', 'C3']),
                         ('exclude', ['X1', 'X2'])):
        for donor in donors:
            for cell_type in ('Endothelial', 'Fibroblast'):
                for _ in range(cells_per_donor):
                    rows.append((donor, role, cell_type))
    obs = pd.DataFrame(rows, columns=['donor', '_role', 'cell_type'])
    obs['condition'] = obs['_role'].map(
        {'disease': 'DCM', 'control': 'Donor', 'exclude': 'DoxCM'})
    obs.index = [f'cell{i}' for i in range(len(obs))]
    X = rng.poisson(6, size=(len(obs), n_genes)).astype(np.float32)
    var = pd.DataFrame(index=[f'GENE{i}' for i in range(n_genes)])
    return anndata.AnnData(X=X, obs=obs, var=var)


# --- slugify -------------------------------------------------------------

@pytest.mark.parametrize('name,expected', [
    ('Endothelial', 'Endothelial'),
    ('Endothelial cell', 'Endothelial_cell'),
    ('Endothelial cell (vascular)', 'Endothelial_cell_vascular'),
    ('T/NK cells', 'T_NK_cells'),
    ('  spaced  ', 'spaced'),
    ('!!!', 'unnamed'),
])
def test_slugify_shapes(name, expected):
    assert B.slugify_cell_type(name) == expected


def test_slug_survives_the_discovery_regex():
    """A 'DE' in the cell type must not truncate the accession.

    'discover_de_results' splits on the first '_DE_', so a slug that
    contained one would hand the Meta workspace a half accession.
    """
    from kosmic.meta_analysis.io import _DE_RESULTS_PATTERN

    for name in ('DE cells', 'my DE type', 'type DE', 'Endothelial'):
        accession = f"GSE1_{B.slugify_cell_type(name)}"
        match = _DE_RESULTS_PATTERN.match(f"{accession}_DE_deseq2.csv")
        assert match is not None
        assert match.group(1) == accession


def test_method_label():
    assert B.method_label('ttest', True) == 'welch_cpm_eb'
    assert B.method_label('ttest', False) == 'welch_cpm'
    assert B.method_label('deseq2') == 'deseq2'
    assert B.method_label('unknown') == 'unknown'


# --- planning ------------------------------------------------------------

def test_plan_lists_both_cell_types_as_eligible():
    plans = B.plan_cell_types(_adata(), 'cell_type', 'donor', min_cells=10)
    assert [p.cell_type for p in plans] == ['Endothelial', 'Fibroblast']
    assert all(p.eligible for p in plans)
    assert all(p.n_disease_samples == 3 and p.n_control_samples == 3
               for p in plans)


def test_plan_ignores_excluded_cells():
    """The excluded arm must not inflate cell or donor counts."""
    adata = _adata(cells_per_donor=30)
    plans = {p.cell_type: p for p in
             B.plan_cell_types(adata, 'cell_type', 'donor', min_cells=10)}
    # 6 included donors x 30 cells, not 8 x 30.
    assert plans['Endothelial'].n_cells == 180
    assert plans['Endothelial'].n_samples == 6


def test_plan_marks_thin_cell_types_ineligible():
    adata = _adata(cells_per_donor=30)
    # A third type present in only one donor per arm.
    obs = adata.obs.copy()
    obs.loc[obs.index[:5], 'cell_type'] = 'Rare'
    adata.obs = obs
    plans = {p.cell_type: p for p in
             B.plan_cell_types(adata, 'cell_type', 'donor', min_cells=10)}
    assert plans['Rare'].eligible is False
    assert 'donors per arm' in plans['Rare'].reason


def test_plan_min_cells_drops_thin_donors():
    """A donor below min_cells is dropped, matching create_pseudobulk."""
    adata = _adata(cells_per_donor=30)
    plans = B.plan_cell_types(adata, 'cell_type', 'donor', min_cells=31)
    assert all(p.eligible is False for p in plans)
    assert all(p.n_disease_samples == 0 for p in plans)


def test_plan_honours_an_explicit_cell_type_list():
    plans = B.plan_cell_types(
        _adata(), 'cell_type', 'donor', cell_types=['Endothelial', 'Ghost'],
        min_cells=10)
    assert [p.cell_type for p in plans] == ['Endothelial', 'Ghost']
    assert plans[1].eligible is False
    assert plans[1].reason == 'no cells with this label'


def test_plan_requires_the_role_column():
    adata = _adata()
    del adata.obs['_role']
    with pytest.raises(ValueError, match='_role'):
        B.plan_cell_types(adata, 'cell_type', 'donor')


def test_plan_rejects_missing_columns():
    adata = _adata()
    with pytest.raises(ValueError, match='Cell-type column'):
        B.plan_cell_types(adata, 'nope', 'donor')
    with pytest.raises(ValueError, match='Sample column'):
        B.plan_cell_types(adata, 'cell_type', 'nope')


# --- pseudobulk frame ----------------------------------------------------

def test_pseudobulk_frame_carries_role():
    sample_df = pd.DataFrame({
        'sample': ['D1', 'C1'],
        'condition': ['DCM', 'Donor'],
        'n_cells': [10, 12],
        'role': ['disease', 'control'],
    })
    df = B.pseudobulk_frame(np.array([[1.0, 2.0], [3.0, 4.0]]),
                            sample_df, ['A', 'B'])
    assert list(df.columns) == ['condition', 'n_cells', 'role', 'A', 'B']
    assert list(df.index) == ['D1', 'C1']
    assert df['role'].tolist() == ['disease', 'control']


def test_pseudobulk_frame_without_role():
    sample_df = pd.DataFrame({
        'sample': ['D1'], 'condition': ['DCM'], 'n_cells': [10]})
    df = B.pseudobulk_frame(np.array([[1.0]]), sample_df, ['A'])
    assert list(df.columns) == ['condition', 'n_cells', 'A']


# --- the loop ------------------------------------------------------------

def test_run_writes_discoverable_per_cell_type_results(tmp_path):
    study = tmp_path / 'GSE1'
    runs = B.run_de_by_cell_type(
        _adata(), 'cell_type', 'donor', 'condition', study, 'GSE1',
        de_method='ttest', min_cells=10, detection_min_pct=0.0)

    assert [r.status for r in runs] == ['ok', 'ok']
    assert {r.accession for r in runs} == {'GSE1_Endothelial', 'GSE1_Fibroblast'}
    for r in runs:
        assert r.de_path.exists()
        assert r.pseudobulk_path.exists()
        assert r.n_samples == 6

    found = {e['accession'] for e in discover_de_results(tmp_path)
             if e['file_type'] == 'de_results'}
    assert found == {'GSE1_Endothelial', 'GSE1_Fibroblast'}


def test_run_written_pseudobulk_has_a_role_column(tmp_path):
    runs = B.run_de_by_cell_type(
        _adata(), 'cell_type', 'donor', 'condition', tmp_path / 'GSE1',
        'GSE1', cell_types=['Endothelial'], de_method='ttest',
        min_cells=10, detection_min_pct=0.0)
    pb = pd.read_csv(runs[0].pseudobulk_path, index_col=0)
    assert 'role' in pb.columns
    assert set(pb['role']) == {'disease', 'control'}
    # The excluded arm never reaches the pseudobulk.
    assert len(pb) == 6


def test_run_skips_ineligible_without_touching_disk(tmp_path):
    adata = _adata()
    obs = adata.obs.copy()
    obs.loc[obs.index[:5], 'cell_type'] = 'Rare'
    adata.obs = obs
    runs = {r.cell_type: r for r in B.run_de_by_cell_type(
        adata, 'cell_type', 'donor', 'condition', tmp_path / 'GSE1', 'GSE1',
        de_method='ttest', min_cells=10, detection_min_pct=0.0)}
    assert runs['Rare'].status == 'skipped'
    assert runs['Rare'].de_path is None
    written = {p.name for p in de_stats_dir(tmp_path / 'GSE1').iterdir()}
    assert not any('Rare' in name for name in written)


def test_run_continues_after_one_cell_type_fails(tmp_path, monkeypatch):
    calls = {'n': 0}

    from kosmic.de import de_analysis

    original = de_analysis.run_de_pipeline

    def flaky(adata, *args, **kwargs):
        calls['n'] += 1
        if calls['n'] == 1:
            raise RuntimeError('boom')
        return original(adata, *args, **kwargs)

    monkeypatch.setattr(de_analysis, 'run_de_pipeline', flaky)

    runs = B.run_de_by_cell_type(
        _adata(), 'cell_type', 'donor', 'condition', tmp_path / 'GSE1',
        'GSE1', de_method='ttest', min_cells=10, detection_min_pct=0.0)
    assert [r.status for r in runs] == ['failed', 'ok']
    assert 'boom' in runs[0].message


def test_run_reports_progress_per_cell_type(tmp_path):
    seen = []
    B.run_de_by_cell_type(
        _adata(), 'cell_type', 'donor', 'condition', tmp_path / 'GSE1',
        'GSE1', de_method='ttest', min_cells=10, detection_min_pct=0.0,
        progress_callback=seen.append)
    assert any('[1/2] Endothelial' in m for m in seen)
    assert any('[2/2] Fibroblast' in m for m in seen)


def test_run_can_be_cancelled(tmp_path):
    state = {'stop': False}

    def cancelled():
        return state['stop']

    def on_progress(_msg):
        state['stop'] = True

    runs = B.run_de_by_cell_type(
        _adata(), 'cell_type', 'donor', 'condition', tmp_path / 'GSE1',
        'GSE1', de_method='ttest', min_cells=10, detection_min_pct=0.0,
        progress_callback=on_progress, cancelled=cancelled)
    assert runs[-1].status == 'skipped'
    assert 'cancelled' in runs[-1].message


def test_run_write_false_computes_without_writing(tmp_path):
    runs = B.run_de_by_cell_type(
        _adata(), 'cell_type', 'donor', 'condition', tmp_path / 'GSE1',
        'GSE1', de_method='ttest', min_cells=10, detection_min_pct=0.0,
        write=False)
    assert all(r.status == 'ok' and r.de_path is None for r in runs)
    assert not (tmp_path / 'GSE1').exists()


# --- summaries -----------------------------------------------------------

def test_batch_summary_shape(tmp_path):
    runs = B.run_de_by_cell_type(
        _adata(), 'cell_type', 'donor', 'condition', tmp_path / 'GSE1',
        'GSE1', de_method='ttest', min_cells=10, detection_min_pct=0.0,
        write=False)
    summary = B.batch_summary(runs)
    assert len(summary) == 2
    assert {'cell_type', 'accession', 'status', 'n_significant'} <= set(summary.columns)


def test_common_tested_genes_intersects():
    a = pd.DataFrame({'names': ['A', 'B', 'C']})
    b = pd.DataFrame({'names': ['B', 'C', 'D']})
    assert B.common_tested_genes([a, b]) == ['B', 'C']


def test_common_tested_genes_ignores_empty():
    a = pd.DataFrame({'names': ['A', 'B']})
    assert B.common_tested_genes([a, pd.DataFrame(), None]) == ['A', 'B']
    assert B.common_tested_genes([]) == []


def test_common_tested_genes_accepts_runs():
    runs = [
        B.CellTypeRun('X', 'X', 'S_X', 'ok',
                      de_results=pd.DataFrame({'names': ['A', 'B']})),
        B.CellTypeRun('Y', 'Y', 'S_Y', 'ok',
                      de_results=pd.DataFrame({'names': ['B', 'C']})),
    ]
    assert B.common_tested_genes(runs) == ['B']


def test_write_false_needs_no_output_dir():
    """A dry run has nowhere to write and should not need a path."""
    runs = B.run_de_by_cell_type(
        _adata(), 'cell_type', 'donor', 'condition', None, 'GSE1',
        de_method='ttest', min_cells=10, detection_min_pct=0.0, write=False)
    assert [r.status for r in runs] == ['ok', 'ok']
    assert all(r.de_path is None and r.significant_path is None for r in runs)
    assert all(r.de_results is not None for r in runs)


def test_write_true_still_requires_output_dir():
    with pytest.raises(ValueError, match='output_dir'):
        B.run_de_by_cell_type(
            _adata(), 'cell_type', 'donor', 'condition', None, 'GSE1',
            de_method='ttest', min_cells=10, detection_min_pct=0.0, write=True)


# --- hypothesis mode must reach the batch runner -------------------------

def test_batch_honours_hypothesis_mode(tmp_path):
    """fdr_genes narrows the BH denominator; the fit stays genome-wide.

    The batch runner used to hardcode the pathway sets away and run
    discovery whatever mode was selected, so a hypothesis run came back
    genome-wide with an empty 'pathways' column.
    """
    adata = _adata(n_genes=60)
    panel = {'PanelA': ['GENE0', 'GENE1', 'GENE2', 'GENE3', 'GENE4']}
    committed = set(panel['PanelA'])

    runs = B.run_de_by_cell_type(
        adata, 'cell_type', 'donor', 'condition', tmp_path / 'GSE1', 'GSE1',
        cell_types=['Endothelial'], de_method='ttest', min_cells=10,
        detection_min_pct=0.0,
        pathway_gene_sets=panel, fdr_genes=committed)

    de = runs[0].de_results
    assert set(de['names']) <= committed          # results subset to the panel
    assert (de['pathways'] == 'PanelA').all()     # and annotated with it


def test_batch_without_gene_sets_is_genome_wide(tmp_path):
    runs = B.run_de_by_cell_type(
        _adata(n_genes=60), 'cell_type', 'donor', 'condition',
        tmp_path / 'GSE1', 'GSE1', cell_types=['Endothelial'],
        de_method='ttest', min_cells=10, detection_min_pct=0.0)
    de = runs[0].de_results
    assert len(de) > 5
    assert (de['pathways'].astype(str) == '').all()
