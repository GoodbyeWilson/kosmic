"""
Simulation: Compare gene-level meta-analysis methods
======================================================
Realistic synthetic scRNA-seq with known truly-DE genes.
Runs through KOSMIC's actual DE pipeline (welch_ttest, welch_ttest_eb, deseq2),
then pools with each meta-analysis method and measures performance.

Combinations tested:
  DE methods:  welch_ttest_eb (EB-moderated), deseq2
  Pooling:     DerSimonian-Laird (random), Hartung-Knapp,
               SumRank (Irwin-Hall), SumRank (Permutation, 1000 perms),
               rOP (rth Ordered P-value)

Metrics:
  - Sensitivity (TPR) at FDR < 0.05
  - Empirical FDR (actual false discovery proportion among "significant" genes)
  - AUC-ROC and AUC-PR
  - Spearman correlation of pooled log2FC with true log2FC
  - Rank recovery (proportion of true DE genes in top-K)

Usage:
  conda activate kirk
  python scripts/compare_gene_methods.py
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
N_TRUE_DE = 500          # number of truly DE genes
N_NULL = N_GENES - N_TRUE_DE
N_DATASETS = 6           # number of studies
SAMPLES_PER_GROUP = (3, 8)  # range of pseudobulk samples per condition
CELLS_PER_SAMPLE = (80, 400)
NB_DISPERSION_RANGE = (0.05, 0.8)
DROPOUT_RATE = 0.05
BETWEEN_STUDY_TAU = 0.15  # heterogeneity in true effect across studies

# True effect sizes for DE genes — mix of strong, moderate, weak, both directions
# Drawn from a distribution to be realistic
RNG_SEED = 42

# Gene names
GENE_NAMES = [f'GENE_{i:05d}' for i in range(N_GENES)]

print("=" * 70)
print("  Gene-Level Meta-Analysis Method Benchmark")
print("=" * 70)
print(f"  {N_GENES} genes ({N_TRUE_DE} truly DE, {N_NULL} null)")
print(f"  {N_DATASETS} studies, {SAMPLES_PER_GROUP[0]}-{SAMPLES_PER_GROUP[1]} samples/group")
print(f"  {CELLS_PER_SAMPLE[0]}-{CELLS_PER_SAMPLE[1]} cells/sample")
print(f"  Between-study tau: {BETWEEN_STUDY_TAU}")
print()


# ============================================================================
# 2. Generate true effects
# ============================================================================

rng = np.random.default_rng(RNG_SEED)

# True log2FC for DE genes: mixture of effect sizes
# 30% strong (|lfc| 0.8-2.0), 40% moderate (0.3-0.8), 30% weak (0.1-0.3)
true_log2fc = np.zeros(N_GENES)

de_indices = rng.choice(N_GENES, N_TRUE_DE, replace=False)
de_indices.sort()

for i, idx in enumerate(de_indices):
    r = rng.random()
    if r < 0.30:
        # Strong effect
        magnitude = rng.uniform(0.8, 2.0)
    elif r < 0.70:
        # Moderate effect
        magnitude = rng.uniform(0.3, 0.8)
    else:
        # Weak effect
        magnitude = rng.uniform(0.1, 0.3)

    # Random direction (up or down)
    direction = rng.choice([-1, 1])
    true_log2fc[idx] = direction * magnitude

is_true_de = np.zeros(N_GENES, dtype=bool)
is_true_de[de_indices] = True

print(f"True DE genes: {N_TRUE_DE}")
print(f"  Strong (|lfc|>0.8): {np.sum(np.abs(true_log2fc[de_indices]) > 0.8)}")
print(f"  Moderate (0.3-0.8): {np.sum((np.abs(true_log2fc[de_indices]) >= 0.3) & (np.abs(true_log2fc[de_indices]) <= 0.8))}")
print(f"  Weak (<0.3):        {np.sum(np.abs(true_log2fc[de_indices]) < 0.3)}")
print(f"  Up-regulated:       {np.sum(true_log2fc[de_indices] > 0)}")
print(f"  Down-regulated:     {np.sum(true_log2fc[de_indices] < 0)}")
print()


# ============================================================================
# 3. Simulate datasets (cell-level counts → AnnData)
# ============================================================================

def simulate_study(study_idx, rng):
    """Generate one study's AnnData with known DE effects + between-study noise."""
    # This study's per-gene effect (true effect + heterogeneity)
    study_log2fc = np.zeros(N_GENES)
    for idx in de_indices:
        study_log2fc[idx] = true_log2fc[idx] + rng.normal(0, BETWEEN_STUDY_TAU)

    # Baseline expression (realistic scRNA-seq)
    log_baseline = rng.normal(-4.0, 2.0, N_GENES)
    baseline_mu = np.clip(np.exp(log_baseline), 0.001, 200)

    gene_dispersion = np.exp(rng.uniform(
        np.log(NB_DISPERSION_RANGE[0]),
        np.log(NB_DISPERSION_RANGE[1]),
        N_GENES,
    ))

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

            # Dropout
            dropout_mask = rng.random((n_cells, N_GENES)) < DROPOUT_RATE
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
    n_cells = adata.n_obs
    n_samples = adata.obs['sample'].nunique()
    print(f"  Study {i}: {n_cells:,} cells, {n_samples} samples", flush=True)
print(f"  Done in {time.time() - t0:.1f}s\n", flush=True)


# ============================================================================
# 4. Run DE analysis with KOSMIC's actual pipeline
# ============================================================================

from kosmic.de.de_analysis import run_de_pipeline

DE_METHODS = ['deseq2']
# DE_METHODS = ['welch_ttest_eb', 'deseq2']  # uncomment to include Welch t-test


def run_de_for_study(adata, method, moderate=False):
    """Run one DE method on one study using KOSMIC's pipeline."""
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
            adata,
            sample_col='sample',
            condition_col='condition',
            pathway_gene_sets={},
            min_cells=3,
            de_method=method,
            full_genome=True,
            moderate=moderate,
            progress_callback=lambda msg: None,
        )
        return de_results
    except Exception as e:
        print(f"    DE failed ({method}): {e}")
        return None


print("Running DE analysis on each study...", flush=True)
# de_results_by_method[method_name] = list of DataFrames, one per study
de_results_by_method = {m: [] for m in DE_METHODS}

for i, adata in enumerate(datasets_adata):
    t0 = time.time()
    print(f"  Study {i}:", end='', flush=True)

    if 'deseq2' in DE_METHODS:
        t1 = time.time()
        de = run_de_for_study(adata, 'deseq2')
        if de is not None:
            de_results_by_method['deseq2'].append(de)
            print(f" deseq2({len(de)} genes, {time.time()-t1:.1f}s)", end='', flush=True)
        else:
            print(f" deseq2(FAILED)", end='', flush=True)

    if 'welch_ttest_eb' in DE_METHODS:
        t1 = time.time()
        de_eb = run_de_for_study(adata, 'welch_ttest', moderate=True)
        if de_eb is not None:
            de_results_by_method['welch_ttest_eb'].append(de_eb)
            print(f" eb({len(de_eb)} genes, {time.time()-t1:.1f}s)", end='', flush=True)

    print(f"  total: {time.time()-t0:.1f}s", flush=True)

print(flush=True)


# ============================================================================
# 5. Convert DE results to meta-analysis input format
# ============================================================================

def de_to_ma_input(de_df, study_name):
    """Convert KOSMIC DE results DataFrame to meta-analysis input format."""
    ma_df = de_df[['names', 'logfoldchanges', 'pvals', 'se',
                    'disease_samples', 'control_samples']].copy()
    ma_df['dataset'] = study_name

    # Add pvals_adj for volcano
    from kosmic.numerical import bh_fdr
    ma_df['pvals_adj'] = bh_fdr(ma_df['pvals'])

    return ma_df


# ============================================================================
# 6. Run all meta-analysis pooling methods
# ============================================================================

from kosmic.meta_analysis.pooling.dl import dl_fast
from kosmic.meta_analysis.pooling.reml import reml_fast as reml_random_effects_fast
from kosmic.meta_analysis.pooling.sumrank import sumrank_fast as sumrank_irwin_hall_fast
from kosmic.meta_analysis.pooling.fisher import fisher_fast
from kosmic.meta_analysis.pooling.rop import rop_fast
from kosmic.meta_analysis.pooling.wop import wop_fast as wop_analytical

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
    'DL_random': lambda datasets: dl_fast(datasets, min_studies=1, progress_callback=_gene_progress),
    'Hartung_Knapp': lambda datasets: dl_fast(datasets, min_studies=1, hksj=True, progress_callback=_gene_progress),
    'SumRank_IH': lambda datasets: sumrank_irwin_hall_fast(datasets, min_studies=1, progress_callback=_gene_progress),
    'rOP': lambda datasets: rop_fast(datasets, min_studies=1, progress_callback=_gene_progress),
    'DL_REML_random': lambda datasets: reml_random_effects_fast(datasets, min_studies=1, progress_callback=_gene_progress),
    'wOP': lambda datasets: wop_analytical(datasets, min_studies=1, progress_callback=_gene_progress),
    'Fisher_Top50': lambda datasets: fisher_fast(
        datasets, min_studies=1, progress_callback=_gene_progress),
}


print("Running meta-analysis pooling...", flush=True)
# results[(de_method, pooling_method)] = meta_df
ma_results = {}

for de_method in DE_METHODS:
    de_list = de_results_by_method[de_method]
    if len(de_list) < 2:
        print(f"  Skipping {de_method} — only {len(de_list)} studies succeeded")
        continue

    # Convert to MA input format
    ma_inputs = []
    for i, de_df in enumerate(de_list):
        ma_df = de_to_ma_input(de_df, f'study_{i}')
        ma_inputs.append(ma_df)

    print(f"  {de_method} ({len(ma_inputs)} studies):", end='', flush=True)

    for pool_name, pool_fn in POOLING_METHODS.items():
        print(f"    {pool_name}...", end=' ', flush=True)
        t0 = time.time()
        print("preparing...", end=' ', flush=True)
        try:
            meta_df = pool_fn(ma_inputs)

            # Normalise column names
            if 'pval' in meta_df.columns and 'pvals_pooled' not in meta_df.columns:
                meta_df['pvals_pooled'] = meta_df['pval']
            if 'mean_logfc' in meta_df.columns and 'logfoldchanges' not in meta_df.columns:
                meta_df['logfoldchanges'] = meta_df['mean_logfc']

            ma_results[(de_method, pool_name)] = meta_df
            print(f"    -> {pool_name}: {len(meta_df)} genes, {time.time()-t0:.1f}s", flush=True)
        except Exception as e:
            print(f"    -> {pool_name}: FAILED — {e}", flush=True)

    print(flush=True)

print(flush=True)


# ============================================================================
# 7. Evaluate performance
# ============================================================================

from sklearn.metrics import roc_auc_score, average_precision_score


def evaluate(meta_df, true_log2fc, is_true_de, gene_names):
    """Compute performance metrics for one MA result."""
    # Map results back to gene indices
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

    # Only evaluate genes that appear in the MA results
    has_result = np.isfinite(pooled_pval)
    n_evaluated = has_result.sum()

    if n_evaluated == 0:
        return None

    # BH FDR correction on pooled p-values
    from kosmic.numerical import bh_fdr, neg_log10
    valid = has_result & (pooled_pval > 0)
    padj = np.ones(len(gene_names))
    if valid.sum() > 0:
        padj[valid] = bh_fdr(pooled_pval[valid])

    # Significant at FDR < 0.05
    significant = has_result & (padj < DEFAULT_FDR)
    n_sig = significant.sum()

    # True/false positives among significant
    tp = (significant & is_true_de).sum()
    fp = (significant & ~is_true_de).sum()
    fn = (is_true_de & has_result & ~significant).sum()

    sensitivity = tp / max(is_true_de[has_result].sum(), 1)
    empirical_fdr = fp / max(n_sig, 1)
    precision = tp / max(tp + fp, 1)

    # AUC-ROC and AUC-PR (using -log10(pval) as score)
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

    # Spearman correlation of pooled log2FC with true log2FC (DE genes only)
    de_has_result = is_true_de & has_result & np.isfinite(pooled_lfc)
    if de_has_result.sum() > 5:
        rho, _ = spearmanr(true_log2fc[de_has_result], pooled_lfc[de_has_result])
    else:
        rho = np.nan

    # Sign concordance (DE genes only)
    if de_has_result.sum() > 0:
        sign_match = np.sign(true_log2fc[de_has_result]) == np.sign(pooled_lfc[de_has_result])
        sign_concordance = sign_match.mean()
    else:
        sign_concordance = np.nan

    # Top-K recovery (what proportion of top-K by p-value are truly DE?)
    ranked_indices = np.argsort(pooled_pval[has_result])
    top_k_values = {}
    for k in [50, 100, 200, 500]:
        if k > n_evaluated:
            continue
        top_k_genes = np.where(has_result)[0][ranked_indices[:k]]
        top_k_tp = is_true_de[top_k_genes].sum()
        top_k_values[f'top{k}_precision'] = top_k_tp / k

    return {
        'n_genes_evaluated': int(n_evaluated),
        'n_significant': int(n_sig),
        'true_positives': int(tp),
        'false_positives': int(fp),
        'false_negatives': int(fn),
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

all_results = []

for (de_method, pool_method), meta_df in sorted(ma_results.items()):
    metrics = evaluate(meta_df, true_log2fc, is_true_de, GENE_NAMES)
    if metrics is None:
        print(f"  {de_method} + {pool_method}: NO RESULTS")
        continue

    metrics['de_method'] = de_method
    metrics['pool_method'] = pool_method
    all_results.append(metrics)

# Print summary table
if all_results:
    results_df = pd.DataFrame(all_results)

    # Reorder columns
    col_order = ['de_method', 'pool_method', 'n_genes_evaluated', 'n_significant',
                 'sensitivity', 'empirical_fdr', 'precision', 'auc_roc', 'auc_pr',
                 'spearman_rho', 'sign_concordance',
                 'true_positives', 'false_positives', 'false_negatives']
    top_k_cols = [c for c in results_df.columns if c.startswith('top')]
    col_order += sorted(top_k_cols)
    col_order = [c for c in col_order if c in results_df.columns]
    results_df = results_df[col_order]

    # Print formatted table
    print(f"{'DE Method':<18} {'Pooling':<18} {'Sens':>6} {'eFDR':>6} "
          f"{'Prec':>6} {'AUC-ROC':>8} {'AUC-PR':>8} "
          f"{'Spearman':>9} {'Sign%':>6} "
          f"{'Top50':>6} {'Top100':>7} {'Top200':>7} {'Top500':>7}")
    print("-" * 140)

    for _, row in results_df.iterrows():
        top50 = f"{row.get('top50_precision', float('nan')):.3f}" if 'top50_precision' in row else '  N/A'
        top100 = f"{row.get('top100_precision', float('nan')):.3f}" if 'top100_precision' in row else '  N/A'
        top200 = f"{row.get('top200_precision', float('nan')):.3f}" if 'top200_precision' in row else '  N/A'
        top500 = f"{row.get('top500_precision', float('nan')):.3f}" if 'top500_precision' in row else '  N/A'

        print(f"{row['de_method']:<18} {row['pool_method']:<18} "
              f"{row['sensitivity']:>6.3f} {row['empirical_fdr']:>6.3f} "
              f"{row['precision']:>6.3f} {row['auc_roc']:>8.4f} {row['auc_pr']:>8.4f} "
              f"{row['spearman_rho']:>9.4f} {row['sign_concordance']:>6.3f} "
              f"{top50:>6} {top100:>7} {top200:>7} {top500:>7}")

    print()

    # Save to CSV
    out_path = os.path.join(os.path.dirname(__file__), 'gene_method_benchmark_results.csv')
    results_df.to_csv(out_path, index=False)
    print(f"Results saved to: {out_path}")

    # Best method summary
    print()
    print("BEST METHODS:")
    for metric in ['sensitivity', 'auc_roc', 'auc_pr', 'spearman_rho']:
        best_idx = results_df[metric].idxmax()
        best = results_df.loc[best_idx]
        print(f"  Best {metric}: {best['de_method']} + {best['pool_method']} = {best[metric]:.4f}")

    best_fdr_idx = results_df['empirical_fdr'].idxmin()
    best_fdr = results_df.loc[best_fdr_idx]
    print(f"  Best FDR control: {best_fdr['de_method']} + {best_fdr['pool_method']} = {best_fdr['empirical_fdr']:.4f}")
else:
    print("No results to report!")

print("\nDone.")
