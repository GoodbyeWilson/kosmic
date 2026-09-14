"""Sentinel tests for DE-side enrichment + pathway scoring features.

- ORA (Fisher's exact + BH FDR over a gene-set library)
- GO elim enrichment (run-and-shape sentinel; correctness against
  goatools reference is implicit via shared dependency).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

from kosmic.de.ora import run_ora
from kosmic.de.de_analysis import compute_pathway_scores, run_pseudobulk_de


# ----------------------------------------------------------------------
# ORA: a study set hitting one library entry exclusively must produce
# a small p for that entry.
# ----------------------------------------------------------------------
def test_ora_finds_synthetic_enrichment():
    """50 study genes drawn entirely from PathA (50 members) — PathA must
    come back with a tiny p. Paths with no study-set hits are filtered
    out by the min_genes guard, so we don't see PathB/PathC in the result."""
    library = {
        'PathA': [f'A{i:02d}' for i in range(50)],
        'PathB': [f'B{i:02d}' for i in range(50)],
        'PathC': [f'C{i:02d}' for i in range(50)],
    }
    background = (library['PathA'] + library['PathB']
                  + library['PathC'] + [f'X{i:03d}' for i in range(500)])
    study = library['PathA']

    result = run_ora(study, background, library, min_genes=5)

    pa = result[result['Term'] == 'PathA'].iloc[0]
    assert pa['P_value'] < 1e-20
    assert pa['Gene_Count'] == 50
    assert 'FDR' in result.columns


# ----------------------------------------------------------------------
# GO elim enrichment: shape sentinel using bundled reference data.
# Correctness against goatools is implicit -- we delegate to that
# library; this test just proves the pipeline wires up and produces
# defensible output for a planted gene set.
# ----------------------------------------------------------------------
def test_go_elim_enrichment_recovers_glycolysis_term():
    """A pure glycolysis gene set must surface a glycolysis GO term as
    significant. This validates the elim algorithm + bundled gene2go +
    obo data are all wired up correctly."""
    from kosmic.de.go_enrichment import run_enrichment

    study = [
        'HK1', 'HK2', 'GPI', 'PFKL', 'PFKM', 'ALDOA', 'TPI1',
        'GAPDH', 'PGK1', 'PGAM1', 'ENO1', 'PKM', 'LDHA',
    ]
    background = study + [
        'TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS', 'PTEN',
        'AKT1', 'MTOR', 'BAX', 'BCL2', 'CDK4', 'RB1',
    ]

    result = run_enrichment(study, background, ontology='BP',
                             method='elim', alpha=0.05, min_genes=5)

    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    assert {'GO_ID', 'Term', 'P_value'} <= set(result.columns)
    # The top-ranked term for a pure glycolysis input must be glycolysis.
    top = result.sort_values('P_value').iloc[0]
    assert 'glycolysis' in top['Term'].lower()
    assert top['P_value'] < 0.01


# ----------------------------------------------------------------------
# Pathway scoring: counts_layer routes a decontaminated layer into the
# per-donor score, and the VST normalisation runs on the same inputs.
# ----------------------------------------------------------------------
def _toy_pathway_adata():
    rng = np.random.default_rng(0)
    n, ng = 60, 60
    raw = rng.negative_binomial(15, 0.3, size=(n, ng)).astype('float32')
    genes = [f'GENE{i}' for i in range(ng)]
    obs = pd.DataFrame({'sample': [f'S{i // 5}' for i in range(n)],
                        'condition': ['ctrl'] * 30 + ['dis'] * 30})
    obs.index = [f'c{i}' for i in range(n)]
    a = AnnData(X=csr_matrix(raw), obs=obs, var=pd.DataFrame(index=genes))
    a.layers['counts'] = csr_matrix(raw)
    a.layers['decontX_counts'] = csr_matrix(np.rint(raw * 0.6))  # distinct values
    return a, genes


def test_pathway_detection_and_coverage_gates():
    """Genes below the detection floor drop out of a pathway's score, and a
    pathway is dropped when too few of its annotated genes stay detected."""
    rng = np.random.default_rng(0)
    n_per, ng = 8, 10
    donors = [f"D{i}" for i in range(3)] + [f"C{i}" for i in range(3)]
    roles = ["disease"] * 3 + ["control"] * 3

    X = np.zeros((len(donors) * n_per, ng), dtype="float32")
    obs_sample, obs_role = [], []
    for di, (d, r) in enumerate(zip(donors, roles)):
        X[di * n_per:(di + 1) * n_per, :5] = rng.integers(3, 9, size=(n_per, 5))
        # G5-G9 stay zero -> undetectable in this cell type.
        obs_sample += [d] * n_per
        obs_role += [r] * n_per
    obs = pd.DataFrame({"sample": obs_sample, "condition": obs_role, "_role": obs_role})
    obs.index = [f"c{i}" for i in range(len(obs))]
    genes = [f"G{i}" for i in range(ng)]
    a = AnnData(X=csr_matrix(X), obs=obs, var=pd.DataFrame(index=genes))

    pw = {"PW_good": ["G0", "G1", "G2", "G3"],   # all detected
          "PW_bad":  ["G5", "G6", "G7", "G0"]}    # only G0 detected -> coverage 1/4

    _, _, names = compute_pathway_scores(
        a, pw, "sample", "condition", min_cells=1,
        gene_detection_pct=0.05, min_coverage=0.30, min_genes=3)
    assert names == ["PW_good"]

    # Disabling the detection filter lets PW_bad's undetectable genes count again.
    _, _, names_off = compute_pathway_scores(
        a, pw, "sample", "condition", min_cells=1,
        gene_detection_pct=0.0, min_coverage=0.30, min_genes=3)
    assert set(names_off) == {"PW_good", "PW_bad"}

    # The coverage report explains the scored set: every pathway a row, with a
    # status and coverage that mirror the gate compute_pathway_scores applies.
    from kosmic.de.de_analysis import pathway_coverage_report
    rep = pathway_coverage_report(
        a, pw, gene_detection_pct=0.05, min_coverage=0.30, min_genes=3)
    assert set(rep['pathway']) == {"PW_good", "PW_bad"}
    good = rep[rep['pathway'] == "PW_good"].iloc[0]
    bad = rep[rep['pathway'] == "PW_bad"].iloc[0]
    assert good['status'] == 'scored' and good['detected'] == 4
    assert bad['status'].startswith('dropped') and bad['detected'] == 1


def test_pathway_scores_counts_layer_routes_and_vst_runs():
    a, genes = _toy_pathway_adata()
    pw = {'PW1': genes[:20], 'PW2': genes[20:40]}

    raw_sm, _, names = compute_pathway_scores(
        a, pw, 'sample', 'condition', min_cells=1, normalization='deseq2')
    dec_sm, _, _ = compute_pathway_scores(
        a, pw, 'sample', 'condition', min_cells=1, normalization='deseq2',
        counts_layer='decontX_counts')

    assert names == ['PW1', 'PW2']
    # The decontaminated layer is a different count source -> different score.
    assert not np.allclose(raw_sm, dec_sm)

    vst_sm, _, vnames = compute_pathway_scores(
        a, pw, 'sample', 'condition', min_cells=1, normalization='vst')
    assert vnames == names
    assert vst_sm.shape == raw_sm.shape
    assert np.all(np.isfinite(vst_sm)) and vst_sm.min() > 0  # log2-like scale


def test_pseudobulk_expression_matrix_matches_method():
    """The plot/heatmap source is a donor x gene frame of finite, depth-
    corrected values, normalised to match the DE method (VST vs log-CPM)."""
    from kosmic.de.de_analysis import pseudobulk_expression_matrix
    a, genes = _toy_pathway_adata()

    vst = pseudobulk_expression_matrix(a, 'sample', 'condition',
                                       normalization='vst', min_cells=1)
    cpm = pseudobulk_expression_matrix(a, 'sample', 'condition',
                                       normalization='cpm', min_cells=1)

    for df in (vst, cpm):
        assert not df.empty
        assert list(df.columns[:2]) == ['sample', 'condition']
        assert set(genes).issubset(set(df.columns))
        assert df['sample'].nunique() == len(df)  # one row per donor
        assert np.all(np.isfinite(df[genes].to_numpy(dtype=float)))
    # The two normalisations are genuinely different quantities.
    assert not np.allclose(vst[genes].to_numpy(float), cpm[genes].to_numpy(float))


def test_expression_matrix_zero_floor():
    """log2 norm-counts (deseq2) and log2 CPM put a zero-count gene at 0; VST
    does not (its non-zero floor is why the bars use the former, not VST)."""
    from kosmic.de.de_analysis import pseudobulk_expression_matrix
    n_per, ng = 8, 6
    rng = np.random.default_rng(0)
    X = np.zeros((6 * n_per, ng), dtype='float32')
    X[:, :5] = rng.integers(20, 80, size=(6 * n_per, 5))  # G5 stays all-zero
    obs = pd.DataFrame({
        'sample': [f'd{i}' for i in range(6) for _ in range(n_per)],
        'condition': ['ctrl'] * (3 * n_per) + ['dis'] * (3 * n_per)})
    obs.index = [f'c{i}' for i in range(6 * n_per)]
    a = AnnData(X=csr_matrix(X), obs=obs,
                var=pd.DataFrame(index=[f'G{i}' for i in range(ng)]))

    common = dict(sample_col='sample', condition_col='condition', min_cells=1)
    cpm = pseudobulk_expression_matrix(a, normalization='cpm', **common)
    dsq = pseudobulk_expression_matrix(a, normalization='deseq2', **common)
    vst = pseudobulk_expression_matrix(a, normalization='vst', **common)

    assert np.allclose(cpm['G5'].to_numpy(float), 0.0)
    assert np.allclose(dsq['G5'].to_numpy(float), 0.0)
    assert (vst['G5'].to_numpy(float) > 0).all()  # VST floor is non-zero


def test_pathway_gene_standardization():
    """The gene-standardised VST scorer z-scores each gene across donors, so
    the module score is a different, mean-zero (equal-weight) quantity -- not
    the log2 VST score."""
    a, genes = _toy_pathway_adata()
    pw = {'PW1': genes[:20]}

    plain, _, _ = compute_pathway_scores(
        a, pw, 'sample', 'condition', min_cells=1, normalization='vst')
    std, _, _ = compute_pathway_scores(
        a, pw, 'sample', 'condition', min_cells=1, normalization='vst',
        standardize_genes=True)

    assert std.shape == plain.shape
    assert not np.allclose(plain, std)
    # Per-gene z-scoring makes the per-donor module score mean-zero.
    assert abs(float(std[:, 0].mean())) < 1e-6


# ----------------------------------------------------------------------
# SE floor: the config knob bounds the reported SE without touching the
# effect, so inverse-variance meta-analysis weights are not silently capped.
# ----------------------------------------------------------------------
def test_pseudobulk_de_se_floor_bounds_reported_se():
    # Near-identical donor scores -> tiny natural SE, well below any floor.
    scores = np.array([[5.00], [5.001], [5.002],
                       [5.50], [5.501], [5.502]])
    sample_df = pd.DataFrame({
        'sample': [f'd{i}' for i in range(6)],
        'condition': ['ctrl'] * 3 + ['dis'] * 3,
        'role': ['control'] * 3 + ['disease'] * 3,
    })
    hi = run_pseudobulk_de(scores, sample_df, ['PW'], pre_transformed=True,
                           moderate=False, se_floor=0.5, min_expressing_samples=1)
    lo = run_pseudobulk_de(scores, sample_df, ['PW'], pre_transformed=True,
                           moderate=False, se_floor=0.0, min_expressing_samples=1)

    assert hi['se'].iloc[0] >= 0.5 - 1e-9   # floored up
    assert lo['se'].iloc[0] < 0.5           # not floored
    # The effect (log2FC) is unchanged by the SE floor.
    assert np.isclose(hi['logfoldchanges'].iloc[0], lo['logfoldchanges'].iloc[0])


# ----------------------------------------------------------------------
# Spillover diagnostic: honours counts_layer, so it can be re-run on the
# decontaminated layer as a before/after validation.
# ----------------------------------------------------------------------
def test_pathway_spillover_counts_layer_routes():
    from kosmic.de.spillover import pathway_spillover

    rng = np.random.default_rng(1)
    n, ng = 120, 40  # 12 samples x 10 cells (>= DE min_cells)
    raw = rng.negative_binomial(15, 0.3, size=(n, ng)).astype('float32')
    genes = [f'GENE{i}' for i in range(ng)]
    obs = pd.DataFrame({
        'sample': [f'S{i // 10}' for i in range(n)],
        'condition': ['ctrl'] * 60 + ['dis'] * 60,
        'pct_counts_test': rng.random(n) * 10,
    })
    obs.index = [f'c{i}' for i in range(n)]
    a = AnnData(X=csr_matrix(raw), obs=obs, var=pd.DataFrame(index=genes))
    a.layers['counts'] = csr_matrix(raw)
    a.layers['decontX_counts'] = csr_matrix(np.rint(raw * 0.6))

    common = dict(sample_col='sample', condition_col='condition',
                  disease_label='dis', control_label='ctrl',
                  contam_key='pct_counts_test', method='deseq2')
    raw_r = pathway_spillover(a, 'PW', genes[:20], **common)
    dec_r = pathway_spillover(a, 'PW', genes[:20], counts_layer='decontX_counts', **common)

    # The decontaminated layer feeds a different per-donor score into the check.
    assert not np.allclose(raw_r['score'], dec_r['score'])


# ----------------------------------------------------------------------
# Pre-DE detection filter: keep genes detected in >= X% cells in >=1
# condition (OR), so on/off genes survive but the noise floor is dropped.
# ----------------------------------------------------------------------
def test_detection_prefilter():
    from kosmic.de.de_analysis import genes_passing_detection

    n = 20  # 10 disease, 10 control
    X = np.zeros((n, 4))
    X[:10, 0] = 5           # GENE0: all disease, 0 control -> on/off, kept via OR
    X[0, 1] = 3   # GENE1: 10% disease, 10% control
    X[10, 1] = 3
    # GENE2: detected in no cells -> dropped
    X[:, 3] = 2             # GENE3: all cells -> kept
    obs = pd.DataFrame({'_role': ['disease'] * 10 + ['control'] * 10})
    obs.index = [f'c{i}' for i in range(n)]
    a = AnnData(X=csr_matrix(X), obs=obs,
                var=pd.DataFrame(index=['GENE0', 'GENE1', 'GENE2', 'GENE3']))
    genes = ['GENE0', 'GENE1', 'GENE2', 'GENE3']

    assert set(genes_passing_detection(a, genes, 0.05)) == {'GENE0', 'GENE1', 'GENE3'}
    # GENE1 (10%) drops at a 20% floor; the on/off GENE0 still survives.
    assert set(genes_passing_detection(a, genes, 0.20)) == {'GENE0', 'GENE3'}
    # Disabled (0) and missing roles both return everything.
    assert genes_passing_detection(a, genes, 0) == genes
    a2 = a.copy()
    del a2.obs['_role']
    assert genes_passing_detection(a2, genes, 0.5) == genes
