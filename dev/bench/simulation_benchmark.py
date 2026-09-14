"""
Simulation-based Meta-Analysis Benchmark
=========================================
Generates synthetic scRNA-seq datasets with known ground truth DE genes,
runs DESeq2 per study, pools with all fast vectorised methods, and
evaluates performance (sensitivity, FDR, AUC-ROC, etc.).

Extracted from ``scripts/compare_top_methods_fast.py`` so both
the CLI script and the GUI dialog can reuse the same logic.
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from kosmic import DEFAULT_FDR


# ======================================================================
# Configuration
# ======================================================================

@dataclass
class SimulationConfig:
    """All tuneable parameters for the benchmark simulation."""

    n_genes: int = 15_000
    n_true_de: int = 500
    n_datasets: int = 6
    samples_per_group: tuple = (3, 8)
    cells_per_sample: tuple = (80, 400)
    nb_dispersion_range: tuple = (0.05, 0.8)
    dropout_rate: float = 0.05
    between_study_tau: float = 0.15
    n_perms: int = 1_000
    run_cc: bool = False
    n_cc_perms: int = 100
    rng_seed: int = 42

    # Spurious signal: fake DE injected into spurious_datasets/n_datasets studies.
    # These genes are NOT in the ground truth; a robust method should reject them.
    spurious_datasets: int = 0
    n_spurious_genes: int = 500
    spurious_lfc_range: tuple = (0.5, 1.5)

    @property
    def n_null(self) -> int:
        return self.n_genes - self.n_true_de


# ======================================================================
# Ground truth generation
# ======================================================================

def generate_ground_truth(config: SimulationConfig):
    """Create ground-truth log2FC and DE-gene indices.

    Returns
    -------
    true_log2fc : ndarray (n_genes,)
    de_indices  : ndarray (n_true_de,) sorted
    is_true_de  : ndarray bool (n_genes,)
    gene_names  : list of str
    spurious_info : dict or None
        If spurious_datasets > 0: {spurious_indices, spurious_lfc,
        is_spurious}.  These genes are NOT in is_true_de -- they are
        null genes that will get fake signal in a subset of studies.
    """
    rng = np.random.default_rng(config.rng_seed)
    true_log2fc = np.zeros(config.n_genes)
    de_indices = rng.choice(config.n_genes, config.n_true_de, replace=False)
    de_indices.sort()

    for idx in de_indices:
        r = rng.random()
        if r < 0.30:
            magnitude = rng.uniform(0.8, 2.0)
        elif r < 0.70:
            magnitude = rng.uniform(0.3, 0.8)
        else:
            magnitude = rng.uniform(0.1, 0.3)
        direction = rng.choice([-1, 1])
        true_log2fc[idx] = direction * magnitude

    is_true_de = np.zeros(config.n_genes, dtype=bool)
    is_true_de[de_indices] = True

    gene_names = [f'GENE_{i:05d}' for i in range(config.n_genes)]

    # Spurious signal setup
    spurious_info = None
    if config.spurious_datasets > 0 and config.n_spurious_genes > 0:
        # Spurious genes are picked from the null set (never true DE).
        null_indices = np.where(~is_true_de)[0]
        n_pick = min(config.n_spurious_genes, len(null_indices))
        spurious_indices = rng.choice(null_indices, n_pick, replace=False)
        spurious_indices.sort()

        spurious_lfc = np.zeros(config.n_genes)
        for idx in spurious_indices:
            magnitude = rng.uniform(*config.spurious_lfc_range)
            direction = rng.choice([-1, 1])
            spurious_lfc[idx] = direction * magnitude

        is_spurious = np.zeros(config.n_genes, dtype=bool)
        is_spurious[spurious_indices] = True

        spurious_info = {
            'spurious_indices': spurious_indices,
            'spurious_lfc': spurious_lfc,
            'is_spurious': is_spurious,
        }

    return true_log2fc, de_indices, is_true_de, gene_names, spurious_info


# ======================================================================
# Dataset simulation
# ======================================================================

def _simulate_single_study(study_idx, config, de_indices, true_log2fc,
                           gene_names, extra_lfc=None):
    """Simulate one NB-based scRNA-seq study with known effects.

    Parameters
    ----------
    extra_lfc : ndarray (n_genes,) or None
        Additional log2FC added on top of the true DE effects (e.g. spurious signal).
    """
    import anndata as ad

    rng = np.random.default_rng(config.rng_seed + study_idx + 1)
    n_genes = config.n_genes

    study_log2fc = np.zeros(n_genes)
    for idx in de_indices:
        study_log2fc[idx] = true_log2fc[idx] + rng.normal(0, config.between_study_tau)

    if extra_lfc is not None:
        study_log2fc = study_log2fc + extra_lfc

    log_baseline = rng.normal(-4.0, 2.0, n_genes)
    baseline_mu = np.clip(np.exp(log_baseline), 0.001, 200)
    gene_dispersion = np.exp(rng.uniform(
        np.log(config.nb_dispersion_range[0]),
        np.log(config.nb_dispersion_range[1]),
        n_genes,
    ))

    n_per_group = rng.integers(*config.samples_per_group, endpoint=True)
    obs_rows = []
    X_blocks = []

    for group in ['control', 'disease']:
        for s in range(n_per_group):
            n_cells = rng.integers(*config.cells_per_sample)
            mu = baseline_mu.copy()
            if group == 'disease':
                mu = mu * (2.0 ** study_log2fc)
            r = gene_dispersion
            p = np.clip(r / (r + mu), 1e-10, 1 - 1e-10)
            cells = rng.negative_binomial(r, p, size=(n_cells, n_genes)).astype(np.float32)
            dropout_mask = rng.random((n_cells, n_genes)) < config.dropout_rate
            cells[dropout_mask] = 0
            X_blocks.append(cells)
            sample_id = f'study{study_idx}_{group}_{s}'
            for _ in range(n_cells):
                obs_rows.append({
                    'sample': sample_id,
                    'condition': 'DCM' if group == 'disease' else 'NF',
                })

    X = np.vstack(X_blocks)
    obs = pd.DataFrame(obs_rows)
    obs.index = [f's{study_idx}_cell_{i}' for i in range(len(obs))]
    var = pd.DataFrame(index=gene_names)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs['sample'] = adata.obs['sample'].astype('category')
    adata.obs['condition'] = adata.obs['condition'].astype('category')
    return adata


def generate_simulated_datasets(config, true_log2fc, de_indices, gene_names,
                                progress_cb=None, spurious_info=None):
    """Generate *config.n_datasets* synthetic AnnData objects.

    Parameters
    ----------
    progress_cb : callable(step, current, total, detail), optional
    spurious_info : dict or None
        From generate_ground_truth().  If provided and
        config.spurious_datasets > 0, the first *spurious_datasets*
        studies get the spurious LFC injected.
    """
    datasets = []
    for i in range(config.n_datasets):
        extra_lfc = None
        is_spurious_study = (
            spurious_info is not None
            and i < config.spurious_datasets
        )
        if is_spurious_study:
            extra_lfc = spurious_info['spurious_lfc']

        adata = _simulate_single_study(
            i, config, de_indices, true_log2fc, gene_names,
            extra_lfc=extra_lfc)
        datasets.append(adata)

        tag = ' [SPURIOUS]' if is_spurious_study else ''
        if progress_cb:
            progress_cb('simulate', i + 1, config.n_datasets,
                        f'{adata.n_obs:,} cells, '
                        f'{adata.obs["sample"].nunique()} samples{tag}')
    return datasets


# ======================================================================
# DESeq2 per study
# ======================================================================

def run_de_on_simulated(datasets, config, progress_cb=None):
    """Run DESeq2 on each simulated AnnData, return MA-ready DataFrames.

    Returns
    -------
    ma_inputs : list of DataFrame
        Each with columns needed by the per-method pool functions.
    """
    from kosmic.de.de_analysis import run_de_pipeline
    from kosmic.numerical import bh_fdr

    ma_inputs = []
    for i, adata in enumerate(datasets):
        # Synthesise _role from the simulator's 'DCM' / 'NF' literals.
        if '_role' not in adata.obs.columns:
            cond = adata.obs['condition'].astype(str)
            adata.obs['_role'] = pd.Categorical(
                np.where(cond == 'DCM', 'disease',
                         np.where(cond == 'NF', 'control', 'exclude')),
                categories=['control', 'disease', 'exclude'])
        de_results, _, _, _ = run_de_pipeline(
            adata, sample_col='sample', condition_col='condition',
            pathway_gene_sets={}, min_cells=3,
            de_method='deseq2', full_genome=True, moderate=False,
            progress_callback=lambda msg: None,
        )
        ma_df = de_results[['names', 'logfoldchanges', 'pvals', 'se',
                            'disease_samples', 'control_samples']].copy()
        ma_df['dataset'] = f'study_{i}'
        valid = ma_df['pvals'].notna() & (ma_df['pvals'] > 0)
        ma_df['pvals_adj'] = 1.0
        if valid.sum() > 0:
            ma_df.loc[valid, 'pvals_adj'] = bh_fdr(ma_df.loc[valid, 'pvals'])
        ma_inputs.append(ma_df)
        if progress_cb:
            progress_cb('de_analysis', i + 1, len(datasets),
                        f'{len(de_results)} genes')
    return ma_inputs


# ======================================================================
# Fast pooling methods
# ======================================================================

def get_fast_pooling_methods():
    """Return ``{name: callable(datasets) -> DataFrame}`` for the analytical poolers.

    Methods are closed-form (z-test / t-test / Irwin-Hall / Beta). For
    empirical calibration use ``cc_permutation.consensus_cc_permutation``.

    Note: ``min_studies=1`` so the benchmark evaluates all genes.
    """
    from kosmic.meta_analysis.pooling.dl import dl_fast
    from kosmic.meta_analysis.pooling.reml import reml_fast
    from kosmic.meta_analysis.pooling.sumrank import sumrank_fast
    from kosmic.meta_analysis.pooling.gwop import gwop_fast
    from kosmic.meta_analysis.pooling.fisher import fisher_fast
    from kosmic.meta_analysis.pooling.stouffer import stouffer_fast

    ms = 1  # benchmark evaluates all genes

    # Analytical pooling methods, mirroring the UI method set.
    methods = {
        'DL_Analytical':               lambda ds: dl_fast(ds, min_studies=ms),
        'REML_Analytical':             lambda ds: reml_fast(ds, min_studies=ms),
        'DL_HKSJ_Analytical':          lambda ds: dl_fast(ds, min_studies=ms, hksj=True),
        'REML_HKSJ_Analytical':        lambda ds: reml_fast(ds, min_studies=ms, hksj=True),
        'DL_Top50_Analytical':         lambda ds: dl_fast(ds, min_studies=ms, proportion_top=0.5),
        'REML_Top50_Analytical':       lambda ds: reml_fast(ds, min_studies=ms, proportion_top=0.5),
        'SumRank_Analytical':          lambda ds: sumrank_fast(ds, min_studies=ms),
        'SumRank_Top50_Analytical':    lambda ds: sumrank_fast(ds, min_studies=ms, proportion_top=0.5),
        'gwOP_Analytical':             lambda ds: gwop_fast(ds, min_studies=ms),
        'Fisher_Analytical':           lambda ds: fisher_fast(ds, min_studies=ms),
        'Stouffer_Analytical':         lambda ds: stouffer_fast(ds, min_studies=ms),
    }

    return methods


FAST_METHOD_NAMES = [
    'DL_Analytical', 'REML_Analytical',
    'DL_HKSJ_Analytical', 'REML_HKSJ_Analytical',
    'DL_Top50_Analytical', 'REML_Top50_Analytical',
    'SumRank_Analytical', 'SumRank_Top50_Analytical',
    'gwOP_Analytical', 'Fisher_Analytical', 'Stouffer_Analytical',
]


# ======================================================================
# Evaluation
# ======================================================================

def evaluate(meta_df, true_log2fc, is_true_de, gene_names,
             is_spurious=None):
    """Compute all performance metrics for one pooled result.

    Parameters
    ----------
    is_spurious : ndarray bool (n_genes,), optional
        Genes that had spurious signal injected. Adds spurious metrics if given.

    Returns
    -------
    dict or None
        ``None`` if no genes have a finite p-value.
    """
    from sklearn.metrics import roc_auc_score, average_precision_score
    from kosmic.numerical import bh_fdr, neg_log10

    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    n = len(gene_names)
    pooled_lfc = np.full(n, np.nan)
    pooled_pval = np.full(n, np.nan)

    pval_col = 'pvals_pooled' if 'pvals_pooled' in meta_df.columns else 'pval'
    lfc_col = 'logfoldchanges' if 'logfoldchanges' in meta_df.columns else 'mean_logfc'

    for _, row in meta_df.iterrows():
        gene = row['names']
        if gene in gene_to_idx:
            idx = gene_to_idx[gene]
            if lfc_col in row.index:
                pooled_lfc[idx] = row[lfc_col]
            if pval_col in row.index:
                pooled_pval[idx] = row[pval_col]

    has_result = np.isfinite(pooled_pval)
    n_evaluated = has_result.sum()
    if n_evaluated == 0:
        return None

    valid = has_result & (pooled_pval > 0)
    padj = np.ones(n)
    if valid.sum() > 0:
        padj[valid] = bh_fdr(pooled_pval[valid])

    significant = has_result & (padj < DEFAULT_FDR)
    n_sig = significant.sum()
    tp = (significant & is_true_de).sum()
    fp = (significant & ~is_true_de).sum()

    sensitivity = tp / max(is_true_de[has_result].sum(), 1)
    empirical_fdr = fp / max(n_sig, 1)
    precision = tp / max(tp + fp, 1)

    labels = is_true_de[has_result].astype(int)
    scores = neg_log10(pooled_pval[has_result])

    try:
        auc_roc = roc_auc_score(labels, scores)
    except ValueError:
        auc_roc = np.nan
    try:
        auc_pr = average_precision_score(labels, scores)
    except ValueError:
        auc_pr = np.nan

    de_has_result = is_true_de & has_result & np.isfinite(pooled_lfc)
    rho = np.nan
    if de_has_result.sum() > 5:
        rho, _ = spearmanr(true_log2fc[de_has_result], pooled_lfc[de_has_result])

    sign_concordance = np.nan
    if de_has_result.sum() > 0:
        sign_match = np.sign(true_log2fc[de_has_result]) == np.sign(pooled_lfc[de_has_result])
        sign_concordance = sign_match.mean()

    top_k_values = {}
    ranked_indices = np.argsort(pooled_pval[has_result])
    for k in [50, 100, 200, 500]:
        if k > n_evaluated:
            continue
        top_k_genes = np.where(has_result)[0][ranked_indices[:k]]
        top_k_values[f'top{k}_precision'] = is_true_de[top_k_genes].sum() / k

    # FPR: fraction of null genes called significant.
    n_null_tested = (~is_true_de & has_result).sum()
    fpr = fp / max(n_null_tested, 1)

    spurious_metrics = {}
    if is_spurious is not None and is_spurious.any():
        n_spurious_total = is_spurious.sum()
        spurious_called = (significant & is_spurious).sum()
        spurious_metrics['n_spurious_called'] = int(spurious_called)
        spurious_metrics['n_spurious_total'] = int(n_spurious_total)
        spurious_metrics['spurious_call_rate'] = (
            spurious_called / max(n_spurious_total, 1))

    return {
        'n_significant': int(n_sig),
        'n_false_positives': int(fp),
        'false_positive_rate': fpr,
        'sensitivity': sensitivity,
        'empirical_fdr': empirical_fdr,
        'precision': precision,
        'auc_roc': auc_roc,
        'auc_pr': auc_pr,
        'spearman_rho': rho,
        'sign_concordance': sign_concordance,
        **top_k_values,
        **spurious_metrics,
    }


# ======================================================================
# Per-study significance (baseline)
# ======================================================================

def compute_per_study_significance(ma_inputs, is_true_de, gene_names,
                                   fdr_threshold=DEFAULT_FDR):
    """Compute per-study and union significance stats.

    Returns
    -------
    dict
        ``per_study``: list with ``n_sig``, ``n_tp``, ``n_fp`` per study.
        ``union_sig`` / ``union_tp`` / ``union_fp``: counts across studies.
    """
    from kosmic.numerical import bh_fdr

    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    n_genes = len(gene_names)

    per_study = []
    sig_in_any = np.zeros(n_genes, dtype=bool)

    for i, df in enumerate(ma_inputs):
        pvals = df['pvals'].values.copy()
        padj = bh_fdr(pvals)

        study_sig = np.zeros(n_genes, dtype=bool)
        for j, row in df.iterrows():
            gene = row['names']
            if gene in gene_to_idx and padj[j] < fdr_threshold:
                study_sig[gene_to_idx[gene]] = True

        sig_in_any |= study_sig
        n_sig = study_sig.sum()
        n_tp = (study_sig & is_true_de).sum()
        n_fp = (study_sig & ~is_true_de).sum()
        per_study.append({
            'study': i, 'n_sig': int(n_sig),
            'n_tp': int(n_tp), 'n_fp': int(n_fp),
        })

    return {
        'per_study': per_study,
        'union_sig': int(sig_in_any.sum()),
        'union_tp': int((sig_in_any & is_true_de).sum()),
        'union_fp': int((sig_in_any & ~is_true_de).sum()),
    }


# ======================================================================
# Orchestrator
# ======================================================================

def _create_pseudobulk_from_adata(adata, gene_names):
    """Aggregate cells per sample into a sum-counts pseudobulk matrix.

    Mirrors ``create_pseudobulk(..., aggregate='sum')`` for benchmark use.

    Returns
    -------
    expression : ndarray (n_samples, n_genes)
    conditions : ndarray of str (n_samples,)
    gene_names : list of str
    samples : list of str
    """
    samples = adata.obs['sample'].cat.categories.tolist()
    n_genes = adata.shape[1]
    expression = np.zeros((len(samples), n_genes))
    conditions = []

    for i, sample in enumerate(samples):
        mask = adata.obs['sample'] == sample
        X = adata[mask].X
        if hasattr(X, 'toarray'):
            X = X.toarray()
        expression[i] = np.asarray(X).sum(axis=0).ravel()
        conditions.append(adata.obs.loc[mask, 'condition'].iloc[0])

    return expression, np.array(conditions), gene_names, samples


def run_welch_de_on_simulated(datasets, config, gene_names, eb_moderate=False,
                              progress_cb=None):
    """Run Welch t-test on pseudobulk from simulated AnnData.

    Returns
    -------
    ma_inputs : list of DataFrame
    pseudobulk_data : list of dict (for CC permutation)
    """
    from kosmic.meta_analysis.cc_permutation import fast_ttest_de
    from kosmic.numerical import bh_fdr

    ma_inputs = []
    pseudobulk_data = []

    for i, adata in enumerate(datasets):
        expression, conditions, gnames, samples = _create_pseudobulk_from_adata(
            adata, gene_names)

        # Engines take masks, not labels.
        is_disease = (conditions == 'DCM')
        is_control = (conditions == 'NF')
        two_group = is_disease | is_control
        logfcs, pvals, ses = fast_ttest_de(
            expression, is_disease, two_group,
            eb_moderate=eb_moderate)

        ma_df = pd.DataFrame({
            'names': gnames,
            'logfoldchanges': logfcs,
            'pvals': pvals,
            'se': ses,
            'dataset': f'study_{i}',
        })
        ma_df['pvals_adj'] = bh_fdr(pvals)
        ma_inputs.append(ma_df)

        # Synthesise per-sample roles for downstream cc_permutation.
        roles = np.where(
            conditions == 'DCM', 'disease',
            np.where(conditions == 'NF', 'control', 'exclude'),
        )
        pseudobulk_data.append({
            'expression': expression,
            'conditions': conditions,
            'gene_names': gnames,
            'role': roles,
            'dataset_name': f'study_{i}',
        })

        if progress_cb:
            label = 'welch_eb_de' if eb_moderate else 'welch_de'
            progress_cb(label, i + 1, len(datasets),
                        f'{len(samples)} samples')

    return ma_inputs, pseudobulk_data


def get_cc_perm_methods(pseudobulk_data, ma_inputs_welch, ma_inputs_welch_eb,
                        n_perms=100):
    """Return CC-perm benchmark methods for each pooling x DE combination.

    Parameters
    ----------
    pseudobulk_data : list of dict
        Pseudobulk data for CC permutation.
    ma_inputs_welch : list of DataFrame
        Welch t-test DE results (plain).
    ma_inputs_welch_eb : list of DataFrame
        Welch t-test (EB) DE results.
    n_perms : int
        Number of CC permutations.
    """
    from kosmic.meta_analysis.pooling.dl import dl_fast
    from kosmic.meta_analysis.pooling.reml import reml_fast
    from kosmic.meta_analysis.pooling.sumrank import sumrank_fast
    from kosmic.meta_analysis.pooling.gwop import gwop_fast
    from kosmic.meta_analysis.cc_permutation import consensus_cc_permutation

    ms = 1  # min_studies

    def _normalize(df):
        if 'pval' in df.columns and 'pvals_pooled' not in df.columns:
            df['pvals_pooled'] = df['pval']
        if 'mean_logfc' in df.columns and 'logfoldchanges' not in df.columns:
            df['logfoldchanges'] = df['mean_logfc']
        return df

    # Benchmark only methods present in both the UI and CC perm dispatch.
    # Each entry: (pooling_name, pool_fn, dispatch_key)
    pool_specs = [
        ('DL',            lambda ds: dl_fast(ds, min_studies=ms),               'dl'),
        ('REML',          lambda ds: reml_fast(ds, min_studies=ms),             'reml'),
        ('DL_HKSJ',       lambda ds: dl_fast(ds, min_studies=ms, hksj=True),    'dl_hksj'),
        ('REML_HKSJ',     lambda ds: reml_fast(ds, min_studies=ms, hksj=True),  'reml_hksj'),
        ('DL_Top50',      lambda ds: dl_fast(ds, min_studies=ms, proportion_top=0.5),
                          'dl_top50'),
        ('REML_Top50',    lambda ds: reml_fast(ds, min_studies=ms, proportion_top=0.5),
                          'reml_top50'),
        ('SumRank',       lambda ds: sumrank_fast(ds, min_studies=ms),          'sumrank'),
        ('SumRank_Top50', lambda ds: sumrank_fast(ds, min_studies=ms, proportion_top=0.5),
                          'sumrank_top50'),
        ('gwOP',          lambda ds: gwop_fast(ds, min_studies=ms),             'gwop'),
    ]

    methods = {}

    # Cache consensus_cc_permutation per DE config so all methods sharing
    # the same ma_inputs + de_method run it once.
    for de_label, ma_inputs, de_method in (
        ('Welch',   ma_inputs_welch,    'welch_cpm'),
        ('WelchEB', ma_inputs_welch_eb, 'welch_cpm_eb'),
    ):
        cache = {}
        def _make_runner(name, pool_fn, key, cache=cache,
                         ma=ma_inputs, dem=de_method):
            def _run(ds):
                if not cache:
                    observed_by_key = {
                        k: _normalize(pf(ma))
                        for (_, pf, k) in pool_specs
                    }
                    cache.update(consensus_cc_permutation(
                        pseudobulk_data,
                        method_meta_dfs=observed_by_key,
                        n_perms=n_perms, seed=42,
                        de_method=dem, min_studies=ms,
                    ))
                return cache.get(key, _normalize(pool_fn(ma)))
            return _run

        for pooling_name, pool_fn, key in pool_specs:
            methods[f'{pooling_name}_CC_{de_label}'] = _make_runner(
                f'{pooling_name}_CC_{de_label}', pool_fn, key)

    return methods


def run_simulation_benchmark(config: SimulationConfig,
                             progress_cb: Optional[Callable] = None):
    """Full benchmark pipeline: simulate -> DE -> pool -> evaluate.

    Parameters
    ----------
    config : SimulationConfig
    progress_cb : callable(step, current, total, detail), optional

    Returns
    -------
    dict
        ``config``, ``true_log2fc``, ``is_true_de``, ``gene_names``,
        ``ma_results`` (per-method DataFrame), ``metrics`` (per-method dict),
        ``results_df`` (DataFrame of all metrics).
    """
    # 1. Ground truth.
    true_log2fc, de_indices, is_true_de, gene_names, spurious_info = \
        generate_ground_truth(config)

    # 2. Simulate datasets.
    datasets = generate_simulated_datasets(
        config, true_log2fc, de_indices, gene_names, progress_cb,
        spurious_info=spurious_info)

    # 3. DESeq2 (analytical / gene-perm methods).
    ma_inputs = run_de_on_simulated(datasets, config, progress_cb)

    # 3b. Welch t-test variants (CC perm methods).
    ma_inputs_welch = None
    ma_inputs_welch_eb = None
    pseudobulk_data = None
    if config.run_cc:
        ma_inputs_welch, pseudobulk_data = run_welch_de_on_simulated(
            datasets, config, gene_names, eb_moderate=False,
            progress_cb=progress_cb)
        ma_inputs_welch_eb, _ = run_welch_de_on_simulated(
            datasets, config, gene_names, eb_moderate=True,
            progress_cb=progress_cb)

    # 3c. Per-study significance baseline.
    study_stats = compute_per_study_significance(
        ma_inputs, is_true_de, gene_names)

    # 4. Pooling (analytical methods on DESeq2 inputs).
    methods = get_fast_pooling_methods()

    # 4b. CC perm methods (on Welch inputs).
    if config.run_cc and pseudobulk_data is not None:
        cc_methods = get_cc_perm_methods(
            pseudobulk_data, ma_inputs_welch, ma_inputs_welch_eb,
            n_perms=config.n_cc_perms)
        methods.update(cc_methods)
    ma_results = {}
    method_names = list(methods.keys())
    total = len(method_names)

    for i, name in enumerate(method_names):
        if progress_cb:
            progress_cb('pooling', i, total, name)
        try:
            meta_df = methods[name](ma_inputs)
            # Normalise column names.
            if 'pval' in meta_df.columns and 'pvals_pooled' not in meta_df.columns:
                meta_df['pvals_pooled'] = meta_df['pval']
            if 'mean_logfc' in meta_df.columns and 'logfoldchanges' not in meta_df.columns:
                meta_df['logfoldchanges'] = meta_df['mean_logfc']
            ma_results[name] = meta_df
        except Exception:
            ma_results[name] = pd.DataFrame()  # failed -> empty frame
        if progress_cb:
            progress_cb('pooling', i + 1, total, name)

    # 5. Evaluate.
    is_spurious = spurious_info['is_spurious'] if spurious_info else None
    all_metrics = {}
    for i, (name, meta_df) in enumerate(ma_results.items()):
        if meta_df.empty:
            continue
        metrics = evaluate(meta_df, true_log2fc, is_true_de, gene_names,
                           is_spurious=is_spurious)
        if metrics:
            metrics['method'] = name
            all_metrics[name] = metrics
        if progress_cb:
            progress_cb('evaluate', i + 1, len(ma_results), name)

    results_df = pd.DataFrame(list(all_metrics.values())) if all_metrics else pd.DataFrame()

    return {
        'config': config,
        'true_log2fc': true_log2fc,
        'is_true_de': is_true_de,
        'gene_names': gene_names,
        'spurious_info': spurious_info,
        'study_stats': study_stats,
        'ma_results': ma_results,
        'metrics': all_metrics,
        'results_df': results_df,
    }


# ======================================================================
# Spurious sweep
# ======================================================================

@dataclass
class SweepConfig:
    """Parameters for the multifactor spurious sweep."""

    n_genes: int = 15_000
    n_true_de: int = 0          # 0 = pure calibration
    n_spurious_genes: int = 500
    n_perms: int = 1_000
    rng_seed: int = 42
    samples_per_group: tuple = (3, 8)
    cells_per_sample: tuple = (80, 400)
    nb_dispersion_range: tuple = (0.05, 0.8)
    dropout_rate: float = 0.05
    between_study_tau: float = 0.15

    # Sweep axes (k_spurious is automatically 1..N//2 for each N).
    n_datasets_values: tuple = (4, 6, 8, 10)
    magnitude_levels: tuple = (
        (0.1, 0.3),    # subtle
        (0.3, 0.8),    # moderate
        (0.8, 1.5),    # strong
        (1.5, 3.0),    # extreme
    )


def run_multifactor_sweep(sweep_cfg: SweepConfig,
                          progress_cb: Optional[Callable] = None):
    """Multifactor sweep over k_contaminated x signal_magnitude x n_datasets.

    For each (N, magnitude) pair: simulate N clean + N contaminated studies
    once, run DESeq2, then for each k in 1..N//2 pick k contaminated +
    (N-k) clean inputs and pool.

    Parameters
    ----------
    sweep_cfg : SweepConfig
    progress_cb : callable(step, current, total, detail), optional

    Returns
    -------
    dict
        ``sweep_cfg`` and ``sweep_df`` (columns: method, k_spurious,
        n_datasets, magnitude_label, spurious_call_rate,
        n_spurious_called, ...).
    """
    mag_labels = []
    for lo, hi in sweep_cfg.magnitude_levels:
        mag_labels.append(f'{lo}-{hi}')

    # Count total simulation + DESeq2 runs for progress.
    total_sim = 0
    total_de_runs = 0
    for n_ds in sweep_cfg.n_datasets_values:
        total_sim += 2 * n_ds * len(sweep_cfg.magnitude_levels)
        total_de_runs += 2 * n_ds * len(sweep_cfg.magnitude_levels)
    sim_done = 0

    # Count total pooling combos for progress.
    total_pooling = 0
    for n_ds in sweep_cfg.n_datasets_values:
        n_k = n_ds // 2
        total_pooling += n_k * len(sweep_cfg.magnitude_levels) * 9  # 9 methods

    de_done = 0
    pooling_done = 0
    sweep_rows = []

    methods = get_fast_pooling_methods()
    method_names = list(methods.keys())

    for n_ds in sweep_cfg.n_datasets_values:
        base_cfg = SimulationConfig(
            n_genes=sweep_cfg.n_genes,
            n_true_de=sweep_cfg.n_true_de,
            n_datasets=n_ds,
            samples_per_group=sweep_cfg.samples_per_group,
            cells_per_sample=sweep_cfg.cells_per_sample,
            nb_dispersion_range=sweep_cfg.nb_dispersion_range,
            dropout_rate=sweep_cfg.dropout_rate,
            between_study_tau=sweep_cfg.between_study_tau,
            n_perms=sweep_cfg.n_perms,
            rng_seed=sweep_cfg.rng_seed,
            spurious_datasets=n_ds,
            n_spurious_genes=sweep_cfg.n_spurious_genes,
        )

        for mag_idx, (mag_lo, mag_hi) in enumerate(sweep_cfg.magnitude_levels):
            mag_label = mag_labels[mag_idx]

            # Generate ground truth with this magnitude
            mag_cfg = SimulationConfig(
                n_genes=sweep_cfg.n_genes,
                n_true_de=sweep_cfg.n_true_de,
                n_datasets=n_ds,
                samples_per_group=sweep_cfg.samples_per_group,
                cells_per_sample=sweep_cfg.cells_per_sample,
                nb_dispersion_range=sweep_cfg.nb_dispersion_range,
                dropout_rate=sweep_cfg.dropout_rate,
                between_study_tau=sweep_cfg.between_study_tau,
                n_perms=sweep_cfg.n_perms,
                rng_seed=sweep_cfg.rng_seed,
                spurious_datasets=n_ds,
                n_spurious_genes=sweep_cfg.n_spurious_genes,
                spurious_lfc_range=(mag_lo, mag_hi),
            )
            true_log2fc, de_indices, is_true_de, gene_names, spurious_info = \
                generate_ground_truth(mag_cfg)
            is_spurious = spurious_info['is_spurious']

            # Simulate clean + contaminated.
            clean_datasets = []
            for i in range(n_ds):
                adata = _simulate_single_study(
                    i, base_cfg, de_indices, true_log2fc, gene_names,
                    extra_lfc=None)
                clean_datasets.append(adata)
                sim_done += 1
                if progress_cb:
                    progress_cb('simulate', sim_done, total_sim,
                                f'N={n_ds}, mag={mag_label}, clean {i}')

            contaminated_datasets = []
            for i in range(n_ds):
                adata = _simulate_single_study(
                    i, base_cfg, de_indices, true_log2fc, gene_names,
                    extra_lfc=spurious_info['spurious_lfc'])
                contaminated_datasets.append(adata)
                sim_done += 1
                if progress_cb:
                    progress_cb('simulate', sim_done, total_sim,
                                f'N={n_ds}, mag={mag_label}, contam {i}')

            # DESeq2 on all studies, clean + contaminated.
            all_ds = clean_datasets + contaminated_datasets

            def _de_progress(step, current, total, detail,
                             _base=de_done, _n=n_ds, _mag=mag_label):
                if progress_cb:
                    progress_cb('de_analysis', _base + current,
                                total_de_runs,
                                f'N={_n}, mag={_mag}: DESeq2 {current}/{total}')

            ma_all = run_de_on_simulated(all_ds, base_cfg,
                                         progress_cb=_de_progress)
            de_done += 2 * n_ds
            ma_clean = ma_all[:n_ds]
            ma_contaminated = ma_all[n_ds:]

            # Sweep k for this (N, magnitude).
            k_max = n_ds // 2
            for k in range(1, k_max + 1):
                mixed = list(ma_contaminated[:k]) + list(ma_clean[k:])

                for method_name in method_names:
                    if progress_cb:
                        progress_cb('pooling', pooling_done, total_pooling,
                                    f'N={n_ds}, mag={mag_label}, '
                                    f'k={k}, {method_name}')
                    try:
                        meta_df = methods[method_name](mixed)
                        if 'pval' in meta_df.columns and 'pvals_pooled' not in meta_df.columns:
                            meta_df['pvals_pooled'] = meta_df['pval']
                        if 'mean_logfc' in meta_df.columns and 'logfoldchanges' not in meta_df.columns:
                            meta_df['logfoldchanges'] = meta_df['mean_logfc']

                        metrics = evaluate(meta_df, true_log2fc, is_true_de,
                                           gene_names, is_spurious=is_spurious)
                        if metrics:
                            metrics['method'] = method_name
                            metrics['k_spurious'] = k
                            metrics['n_datasets'] = n_ds
                            metrics['magnitude'] = mag_label
                            sweep_rows.append(metrics)
                    except Exception:
                        pass
                    pooling_done += 1

    if progress_cb:
        progress_cb('evaluate', 1, 1, 'Done')

    sweep_df = pd.DataFrame(sweep_rows) if sweep_rows else pd.DataFrame()

    return {
        'sweep_cfg': sweep_cfg,
        'sweep_df': sweep_df,
    }


# ======================================================================
# Save
# ======================================================================

def save_results(results: dict, output_dir: str) -> str:
    """Write results CSV to *output_dir*, return the file path."""
    os.makedirs(output_dir, exist_ok=True)
    timestamp = time.strftime('%Y%m%d_%H%M%S')
    csv_path = os.path.join(output_dir, f'benchmark_{timestamp}.csv')
    results['results_df'].to_csv(csv_path, index=False)
    return csv_path
