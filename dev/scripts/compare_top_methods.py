"""
Quick benchmark: top 4 meta-analysis pooling methods
=====================================================
Runs only DL Random, DL REML, DL Top50 Perm, and DL GeneGate Top50.
Uses the same simulation framework as compare_gene_methods.py but
skips the slower/weaker methods.

Usage:
  conda activate kirk
  python scripts/compare_top_methods.py
"""

import sys
import os
import time
import warnings
from kosmic import DEFAULT_FDR

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # repo root: dev/<this dir>/<file>
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMBA_NUM_THREADS'] = '1'
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
from scipy.stats import spearmanr

# ============================================================================
# 1. Simulation parameters
# ============================================================================

N_GENES = 15000
N_TRUE_DE = 500
N_NULL = N_GENES - N_TRUE_DE
N_DATASETS = 6
SAMPLES_PER_GROUP = (3, 8)
CELLS_PER_SAMPLE = (80, 400)
NB_DISPERSION_RANGE = (0.05, 0.8)
DROPOUT_RATE = 0.05
BETWEEN_STUDY_TAU = 0.15
RNG_SEED = 42

GENE_NAMES = [f'GENE_{i:05d}' for i in range(N_GENES)]

print("=" * 70)
print("  Top Methods Benchmark (DL, REML, Top50 Perm, GeneGate)")
print("=" * 70)
print(f"  {N_GENES} genes ({N_TRUE_DE} truly DE, {N_NULL} null)")
print(f"  {N_DATASETS} studies, {SAMPLES_PER_GROUP[0]}-{SAMPLES_PER_GROUP[1]} samples/group")
print()

# ============================================================================
# 2. Generate true effects
# ============================================================================

rng = np.random.default_rng(RNG_SEED)

true_log2fc = np.zeros(N_GENES)
de_indices = rng.choice(N_GENES, N_TRUE_DE, replace=False)
de_indices.sort()

for i, idx in enumerate(de_indices):
    r = rng.random()
    if r < 0.30:
        magnitude = rng.uniform(0.8, 2.0)
    elif r < 0.70:
        magnitude = rng.uniform(0.3, 0.8)
    else:
        magnitude = rng.uniform(0.1, 0.3)
    direction = rng.choice([-1, 1])
    true_log2fc[idx] = direction * magnitude

is_true_de = np.zeros(N_GENES, dtype=bool)
is_true_de[de_indices] = True

print(f"True DE genes: {N_TRUE_DE}")
print(f"  Strong (|lfc|>0.8): {np.sum(np.abs(true_log2fc[de_indices]) > 0.8)}")
print(f"  Moderate (0.3-0.8): {np.sum((np.abs(true_log2fc[de_indices]) >= 0.3) & (np.abs(true_log2fc[de_indices]) <= 0.8))}")
print(f"  Weak (<0.3):        {np.sum(np.abs(true_log2fc[de_indices]) < 0.3)}")
print()

# ============================================================================
# 3. Simulate datasets
# ============================================================================

def simulate_study(study_idx, rng):
    study_log2fc = np.zeros(N_GENES)
    for idx in de_indices:
        study_log2fc[idx] = true_log2fc[idx] + rng.normal(0, BETWEEN_STUDY_TAU)

    log_baseline = rng.normal(-4.0, 2.0, N_GENES)
    baseline_mu = np.clip(np.exp(log_baseline), 0.001, 200)
    gene_dispersion = np.exp(rng.uniform(
        np.log(NB_DISPERSION_RANGE[0]), np.log(NB_DISPERSION_RANGE[1]), N_GENES))

    n_per_group = rng.integers(*SAMPLES_PER_GROUP, endpoint=True)
    obs_rows = []
    X_blocks = []

    for group in ['control', 'disease']:
        for s in range(n_per_group):
            n_cells = rng.integers(*CELLS_PER_SAMPLE)
            mu = baseline_mu.copy()
            if group == 'disease':
                mu = mu * (2.0 ** study_log2fc)
            r = gene_dispersion
            p = np.clip(r / (r + mu), 1e-10, 1 - 1e-10)
            cells = rng.negative_binomial(r, p, size=(n_cells, N_GENES)).astype(np.float32)
            dropout_mask = rng.random((n_cells, N_GENES)) < DROPOUT_RATE
            cells[dropout_mask] = 0
            X_blocks.append(cells)
            sample_id = f'study{study_idx}_{group}_{s}'
            for _ in range(n_cells):
                obs_rows.append({'sample': sample_id, 'condition': 'DCM' if group == 'disease' else 'NF'})

    X = np.vstack(X_blocks)
    obs = pd.DataFrame(obs_rows)
    obs.index = [f's{study_idx}_cell_{i}' for i in range(len(obs))]
    var = pd.DataFrame(index=GENE_NAMES)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs['sample'] = adata.obs['sample'].astype('category')
    adata.obs['condition'] = adata.obs['condition'].astype('category')
    return adata


print("Generating datasets...", flush=True)
t0 = time.time()
datasets_adata = []
for i in range(N_DATASETS):
    study_rng = np.random.default_rng(RNG_SEED + i + 1)
    adata = simulate_study(i, study_rng)
    datasets_adata.append(adata)
    print(f"  Study {i}: {adata.n_obs:,} cells, {adata.obs['sample'].nunique()} samples", flush=True)
print(f"  Done in {time.time() - t0:.1f}s\n", flush=True)

# ============================================================================
# 4. Run DESeq2 on each study
# ============================================================================

from kosmic.de.de_analysis import run_de_pipeline

print("Running DESeq2 on each study...", flush=True)
de_results_list = []
for i, adata in enumerate(datasets_adata):
    t0 = time.time()
    try:
        if '_role' not in adata.obs.columns:
            import numpy as _np
            import pandas as _pd
            cond = adata.obs['condition'].astype(str)
            adata.obs['_role'] = _pd.Categorical(
                _np.where(cond == 'DCM', 'disease',
                          _np.where(cond == 'NF', 'control', 'exclude')),
                categories=['control', 'disease', 'exclude'])
        de_results, _, _, _ = run_de_pipeline(
            adata, sample_col='sample', condition_col='condition',
            pathway_gene_sets={}, min_cells=3,
            de_method='deseq2', full_genome=True, moderate=False,
            progress_callback=lambda msg: None,
        )
        de_results_list.append(de_results)
        print(f"  Study {i}: {len(de_results)} genes, {time.time()-t0:.1f}s", flush=True)
    except Exception as e:
        print(f"  Study {i}: FAILED -- {e}", flush=True)

print(flush=True)

# ============================================================================
# 5. Convert to MA input format
# ============================================================================

ma_inputs = []
for i, de_df in enumerate(de_results_list):
    ma_df = de_df[['names', 'logfoldchanges', 'pvals', 'se',
                    'disease_samples', 'control_samples']].copy()
    ma_df['dataset'] = f'study_{i}'

    from kosmic.numerical import bh_fdr
    ma_df['pvals_adj'] = bh_fdr(ma_df['pvals'])
    ma_inputs.append(ma_df)

print(f"Prepared {len(ma_inputs)} studies for meta-analysis\n", flush=True)

# ============================================================================
# 6. Run the 4 pooling methods
# ============================================================================

from kosmic.meta_analysis.pooling.dl import dl_fast
from kosmic.meta_analysis.pooling.reml import reml_fast as reml_random_effects

_perm_t0 = None

def _perm_progress(current, total):
    global _perm_t0
    if current == 1:
        _perm_t0 = time.time()
    elapsed = time.time() - _perm_t0 if _perm_t0 else 0
    rate = current / elapsed if elapsed > 0 else 0
    eta = (total - current) / rate if rate > 0 else 0
    pct = current / total * 100
    bar_len = 30
    filled = int(bar_len * current / total)
    bar = '#' * filled + '-' * (bar_len - filled)
    print(f"\r      [{bar}] {current}/{total} ({pct:.0f}%) "
          f"{rate:.1f} perm/s, ETA {eta:.0f}s   ", end='', flush=True)
    if current == total:
        print(f"\r      [{bar}] {total}/{total} (100%) "
              f"done in {elapsed:.1f}s               ", flush=True)

def _gene_progress(current, total):
    global _perm_t0
    if current == 1:
        _perm_t0 = time.time()
    elapsed = time.time() - _perm_t0 if _perm_t0 else 0
    rate = current / elapsed if elapsed > 0 else 0
    eta = (total - current) / rate if rate > 0 else 0
    pct = current / total * 100
    bar_len = 30
    filled = int(bar_len * current / total)
    bar = '#' * filled + '-' * (bar_len - filled)
    print(f"\r      [{bar}] {current}/{total} genes ({pct:.0f}%) "
          f"ETA {eta:.0f}s   ", end='', flush=True)
    if current == total:
        print(f"\r      [{bar}] {total}/{total} genes (100%) "
              f"done in {elapsed:.1f}s               ", flush=True)

POOLING_METHODS = {
    'DL_random': lambda ds: dl_fast(ds, min_studies=1, progress_callback=_gene_progress),
    'DL_REML': lambda ds: reml_random_effects(ds, min_studies=1, progress_callback=_gene_progress),
}

print("Running meta-analysis pooling...", flush=True)
ma_results = {}

for pool_name, pool_fn in POOLING_METHODS.items():
    print(f"  {pool_name}...", flush=True)
    t0 = time.time()
    try:
        meta_df = pool_fn(ma_inputs)

        # Normalise columns
        if 'pval' in meta_df.columns and 'pvals_pooled' not in meta_df.columns:
            meta_df['pvals_pooled'] = meta_df['pval']
        if 'mean_logfc' in meta_df.columns and 'logfoldchanges' not in meta_df.columns:
            meta_df['logfoldchanges'] = meta_df['mean_logfc']

        ma_results[pool_name] = meta_df
        print(f"    -> {len(meta_df)} genes, {time.time()-t0:.1f}s", flush=True)
    except Exception as e:
        print(f"    -> FAILED: {e}", flush=True)
        import traceback
        traceback.print_exc()

print(flush=True)

# ============================================================================
# 7. Evaluate
# ============================================================================

from sklearn.metrics import roc_auc_score, average_precision_score


def evaluate(meta_df, true_log2fc, is_true_de, gene_names):
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    pooled_lfc = np.full(len(gene_names), np.nan)
    pooled_pval = np.full(len(gene_names), np.nan)

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

    from kosmic.numerical import bh_fdr, neg_log10
    valid = has_result & (pooled_pval > 0)
    padj = np.ones(len(gene_names))
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
    if de_has_result.sum() > 5:
        rho, _ = spearmanr(true_log2fc[de_has_result], pooled_lfc[de_has_result])
    else:
        rho = np.nan

    if de_has_result.sum() > 0:
        sign_match = np.sign(true_log2fc[de_has_result]) == np.sign(pooled_lfc[de_has_result])
        sign_concordance = sign_match.mean()
    else:
        sign_concordance = np.nan

    top_k_values = {}
    ranked_indices = np.argsort(pooled_pval[has_result])
    for k in [50, 100, 200, 500]:
        if k > n_evaluated:
            continue
        top_k_genes = np.where(has_result)[0][ranked_indices[:k]]
        top_k_values[f'top{k}_precision'] = is_true_de[top_k_genes].sum() / k

    return {
        'n_significant': int(n_sig),
        'sensitivity': sensitivity,
        'empirical_fdr': empirical_fdr,
        'precision': precision,
        'auc_roc': auc_roc,
        'auc_pr': auc_pr,
        'spearman_rho': rho,
        'sign_concordance': sign_concordance,
        **top_k_values,
    }


# ============================================================================
# 8. Results table
# ============================================================================

print("=" * 70)
print("  RESULTS")
print("=" * 70)
print()

print(f"{'Pooling':<22} {'Sens':>6} {'eFDR':>6} "
      f"{'Prec':>6} {'AUC-ROC':>8} {'AUC-PR':>8} "
      f"{'Spearman':>9} {'Sign%':>6} "
      f"{'Top50':>6} {'Top100':>7} {'Top200':>7} {'Top500':>7}")
print("-" * 120)

for pool_name, meta_df in ma_results.items():
    metrics = evaluate(meta_df, true_log2fc, is_true_de, GENE_NAMES)
    if metrics is None:
        print(f"  {pool_name}: NO RESULTS")
        continue

    top50 = f"{metrics.get('top50_precision', float('nan')):.3f}"
    top100 = f"{metrics.get('top100_precision', float('nan')):.3f}"
    top200 = f"{metrics.get('top200_precision', float('nan')):.3f}"
    top500 = f"{metrics.get('top500_precision', float('nan')):.3f}"

    print(f"{pool_name:<22} "
          f"{metrics['sensitivity']:>6.3f} {metrics['empirical_fdr']:>6.3f} "
          f"{metrics['precision']:>6.3f} {metrics['auc_roc']:>8.4f} {metrics['auc_pr']:>8.4f} "
          f"{metrics['spearman_rho']:>9.4f} {metrics['sign_concordance']:>6.3f} "
          f"{top50:>6} {top100:>7} {top200:>7} {top500:>7}")

print()
print("Done.")
