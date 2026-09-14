"""
FAST benchmark: top 4 meta-analysis pooling methods (vectorized + parallel)
=============================================================================
Uses the per-method pool functions (dl.py / reml.py / sumrank.py / rop.py / ...).
Same simulation, same metrics, ~100x faster pooling.

Usage:
  conda activate kirk
  python scripts/compare_top_methods_fast.py
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
from scipy.stats import spearmanr

# ============================================================================
# Simulation parameters
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


def simulate_study(study_idx, rng, de_indices, true_log2fc):
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


def evaluate(meta_df, true_log2fc, is_true_de, gene_names):
    from sklearn.metrics import roc_auc_score, average_precision_score
    from kosmic.numerical import bh_fdr, neg_log10

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


_progress_t0 = None

def _progress(current, total):
    global _progress_t0
    if current == 1:
        _progress_t0 = time.time()
    elapsed = time.time() - _progress_t0 if _progress_t0 else 0
    rate = current / elapsed if elapsed > 0 else 0
    eta = (total - current) / rate if rate > 0 else 0
    pct = current / total * 100
    bar_len = 30
    filled = int(bar_len * current / total)
    bar = '#' * filled + '-' * (bar_len - filled)
    print(f"\r      [{bar}] {current}/{total} ({pct:.0f}%) "
          f"{rate:.1f}/s, ETA {eta:.0f}s   ", end='', flush=True)
    if current == total:
        print(f"\r      [{bar}] {total}/{total} (100%) "
              f"done in {elapsed:.1f}s               ", flush=True)


def main():
    print("=" * 70)
    print("  Top Methods FAST Benchmark (vectorized + parallel)")
    print("=" * 70)
    print(f"  {N_GENES} genes ({N_TRUE_DE} truly DE, {N_NULL} null)")
    print(f"  {N_DATASETS} studies, {SAMPLES_PER_GROUP[0]}-{SAMPLES_PER_GROUP[1]} samples/group")
    print()

    # Generate true effects
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

    # Generate datasets
    print("Generating datasets...", flush=True)
    t0 = time.time()
    datasets_adata = []
    for i in range(N_DATASETS):
        study_rng = np.random.default_rng(RNG_SEED + i + 1)
        adata = simulate_study(i, study_rng, de_indices, true_log2fc)
        datasets_adata.append(adata)
        print(f"  Study {i}: {adata.n_obs:,} cells, {adata.obs['sample'].nunique()} samples", flush=True)
    print(f"  Done in {time.time() - t0:.1f}s\n", flush=True)

    # Run DESeq2
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

    # Convert to MA input format
    from kosmic.numerical import bh_fdr

    ma_inputs = []
    for i, de_df in enumerate(de_results_list):
        ma_df = de_df[['names', 'logfoldchanges', 'pvals', 'se',
                        'disease_samples', 'control_samples']].copy()
        ma_df['dataset'] = f'study_{i}'
        ma_df['pvals_adj'] = bh_fdr(ma_df['pvals'])
        ma_inputs.append(ma_df)

    print(f"Prepared {len(ma_inputs)} studies for meta-analysis\n", flush=True)

    # Run FAST pooling methods
    from kosmic.meta_analysis.pooling.dl import dl_fast
    from kosmic.meta_analysis.pooling.reml import reml_fast
    from kosmic.meta_analysis.pooling.sumrank import sumrank_fast as sumrank_irwin_hall_fast

    pooling_methods = {
        'DL_random_FAST': lambda ds: dl_fast(
            ds, min_studies=1, progress_callback=_progress),
        'DL_REML_FAST': lambda ds: reml_fast(
            ds, min_studies=1, progress_callback=_progress),
        'DL_HKSJ_FAST': lambda ds: dl_fast(
            ds, min_studies=1, hksj=True, progress_callback=_progress),
        'REML_HKSJ_FAST': lambda ds: reml_fast(
            ds, min_studies=1, hksj=True, progress_callback=_progress),
        'DL_Top50_FAST': lambda ds: dl_fast(
            ds, min_studies=1, proportion_top=0.5, progress_callback=_progress),
        'SumRank_IH_FAST': lambda ds: sumrank_irwin_hall_fast(
            ds, min_studies=1, progress_callback=_progress),
    }

    print("Running FAST meta-analysis pooling...", flush=True)
    ma_results = {}

    for pool_name, pool_fn in pooling_methods.items():
        print(f"  {pool_name}...", flush=True)
        t0 = time.time()
        try:
            meta_df = pool_fn(ma_inputs)

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

    # ====================================================================
    # Case-control permutation (gold standard) -- optional, slow
    # ====================================================================

    RUN_CC = '--cc' in sys.argv
    N_CC_PERMS = 100  # override with --cc-perms N
    for i, arg in enumerate(sys.argv):
        if arg == '--cc-perms' and i + 1 < len(sys.argv):
            N_CC_PERMS = int(sys.argv[i + 1])

    if not RUN_CC:
        print("Skipping case-control permutation (add --cc flag to enable)\n", flush=True)
    else:
        print(f"Running case-control permutation ({N_CC_PERMS} perms)...", flush=True)
        print(f"  Each perm: shuffle labels + DESeq2 on {N_DATASETS} studies", flush=True)

        from scipy.stats import ecdf as scipy_ecdf

        # Observed statistics
        obs_meta = dl_fast(ma_inputs, min_studies=1)
        obs_sr = sumrank_irwin_hall_fast(ma_inputs, min_studies=1)

        obs_z_map = {}
        for _, row in obs_meta.iterrows():
            se = max(row['se'], 1e-10)
            obs_z_map[row['names']] = abs(row['logfoldchanges'] / se)

        obs_rs_map = {}
        if not obs_sr.empty:
            obs_rs_map = dict(zip(obs_sr['names'], obs_sr['rank_sum']))

        # Run CC permutations
        t0 = time.time()
        k_to_null_dl = {}
        k_to_null_sr = {}

        for perm_i in range(N_CC_PERMS):
            rng_cc = np.random.default_rng(RNG_SEED + 10000 + perm_i)
            perm_de_list = []

            for adata in datasets_adata:
                adata_perm = adata.copy()
                samples = adata_perm.obs['sample'].values
                unique_samples = np.unique(samples)
                sample_conds = {s: adata_perm.obs.loc[samples == s, 'condition'].values[0]
                                for s in unique_samples}
                sample_names = list(sample_conds.keys())
                shuffled = rng_cc.permutation([sample_conds[s] for s in sample_names])
                new_map = dict(zip(sample_names, shuffled))
                adata_perm.obs['condition'] = pd.Categorical(
                    [new_map[s] for s in samples])

                try:
                    cond_perm = adata_perm.obs['condition'].astype(str)
                    adata_perm.obs['_role'] = pd.Categorical(
                        np.where(cond_perm == 'DCM', 'disease',
                                 np.where(cond_perm == 'NF', 'control', 'exclude')),
                        categories=['control', 'disease', 'exclude'])
                    de_res, _, _, _ = run_de_pipeline(
                        adata_perm, sample_col='sample', condition_col='condition',
                        pathway_gene_sets={}, min_cells=3,
                        de_method='deseq2', full_genome=True, moderate=False,
                        progress_callback=lambda msg: None)
                    ma_df = de_res[['names', 'logfoldchanges', 'pvals', 'se',
                                    'disease_samples', 'control_samples']].copy()
                    ma_df['dataset'] = f'perm_{perm_i}'
                    perm_de_list.append(ma_df)
                except Exception:
                    pass

            if len(perm_de_list) >= 2:
                # DL null
                pm = dl_fast(perm_de_list, min_studies=1)
                if not pm.empty:
                    z_abs = np.abs(pm['logfoldchanges'].values / pm['se'].clip(1e-10).values)
                    k_arr = pm['n_studies'].values.astype(int)
                    for k in np.unique(k_arr):
                        k_to_null_dl.setdefault(int(k), []).extend(z_abs[k_arr == k].tolist())

                # SumRank null
                psr = sumrank_irwin_hall_fast(perm_de_list, min_studies=1)
                if not psr.empty:
                    rs = psr['rank_sum'].values
                    k_arr_sr = psr['n_studies'].values.astype(int)
                    for k in np.unique(k_arr_sr):
                        k_to_null_sr.setdefault(int(k), []).extend(rs[k_arr_sr == k].tolist())

            elapsed = time.time() - t0
            rate = (perm_i + 1) / elapsed if elapsed > 0 else 0
            eta = (N_CC_PERMS - perm_i - 1) / rate if rate > 0 else 0
            pct = (perm_i + 1) / N_CC_PERMS * 100
            filled = int(30 * (perm_i + 1) / N_CC_PERMS)
            bar = '#' * filled + '-' * (30 - filled)
            print(f"\r  [{bar}] {perm_i+1}/{N_CC_PERMS} ({pct:.0f}%) "
                  f"{rate:.2f}/s, ETA {eta:.0f}s   ", end='', flush=True)

        print(f"\r  [{'#' * 30}] {N_CC_PERMS}/{N_CC_PERMS} (100%) "
              f"done in {time.time() - t0:.1f}s               ", flush=True)

        # Calibrate DL
        cc_dl = obs_meta.copy()
        cc_dl['pvals_analytical'] = cc_dl['pvals_pooled'].copy()
        cal_dl = np.ones(len(cc_dl))
        for k, nv in k_to_null_dl.items():
            if not nv:
                continue
            null_cdf = scipy_ecdf(np.array(nv)).cdf
            k_mask = cc_dl['n_studies'].values == k
            for j in np.where(k_mask)[0]:
                gene = cc_dl.iloc[j]['names']
                cal_dl[j] = max(float(1.0 - null_cdf.evaluate(obs_z_map.get(gene, 0.0))),
                                1.0 / (len(nv) + 1))
        cc_dl['pvals_pooled'] = cal_dl
        ma_results['CC_Perm_DL'] = cc_dl
        print(f"  -> CC_Perm_DL: {len(cc_dl)} genes", flush=True)

        # Calibrate SumRank
        if not obs_sr.empty and k_to_null_sr:
            cal_sr = np.ones(len(obs_sr))
            for k, nv in k_to_null_sr.items():
                if not nv:
                    continue
                null_cdf = scipy_ecdf(np.array(nv)).cdf
                k_mask = obs_sr['n_studies'].values == k
                for j in np.where(k_mask)[0]:
                    gene = obs_sr.iloc[j]['names']
                    cal_sr[j] = max(float(null_cdf.evaluate(obs_rs_map.get(gene, 1.0))),
                                    1.0 / (len(nv) + 1))
            obs_sr['pval'] = cal_sr
            ma_results['CC_Perm_SumRank'] = obs_sr
            print(f"  -> CC_Perm_SumRank: {len(obs_sr)} genes", flush=True)

        print(f"  Total: {time.time()-t0:.1f}s\n", flush=True)

    # Results
    print("=" * 70)
    print("  RESULTS")
    print("=" * 70)
    print()

    print(f"{'Pooling':<26} {'Sens':>6} {'eFDR':>6} "
          f"{'Prec':>6} {'AUC-ROC':>8} {'AUC-PR':>8} "
          f"{'Spearman':>9} {'Sign%':>6} "
          f"{'Top50':>6} {'Top100':>7} {'Top200':>7} {'Top500':>7}")
    print("-" * 130)

    for pool_name, meta_df in ma_results.items():
        metrics = evaluate(meta_df, true_log2fc, is_true_de, GENE_NAMES)
        if metrics is None:
            print(f"  {pool_name}: NO RESULTS")
            continue

        top50 = f"{metrics.get('top50_precision', float('nan')):.3f}"
        top100 = f"{metrics.get('top100_precision', float('nan')):.3f}"
        top200 = f"{metrics.get('top200_precision', float('nan')):.3f}"
        top500 = f"{metrics.get('top500_precision', float('nan')):.3f}"

        print(f"{pool_name:<26} "
              f"{metrics['sensitivity']:>6.3f} {metrics['empirical_fdr']:>6.3f} "
              f"{metrics['precision']:>6.3f} {metrics['auc_roc']:>8.4f} {metrics['auc_pr']:>8.4f} "
              f"{metrics['spearman_rho']:>9.4f} {metrics['sign_concordance']:>6.3f} "
              f"{top50:>6} {top100:>7} {top200:>7} {top500:>7}")

    # Save results
    all_metrics = []
    for pool_name, meta_df in ma_results.items():
        metrics = evaluate(meta_df, true_log2fc, is_true_de, GENE_NAMES)
        if metrics:
            metrics['method'] = pool_name
            all_metrics.append(metrics)

    if all_metrics:
        results_df = pd.DataFrame(all_metrics)
        out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
        os.makedirs(out_dir, exist_ok=True)

        timestamp = time.strftime('%Y%m%d_%H%M%S')
        csv_path = os.path.join(out_dir, f'benchmark_{timestamp}.csv')
        results_df.to_csv(csv_path, index=False)
        print(f"Results saved to {csv_path}")

        # Also save per-method pooled results for inspection
        for pool_name, meta_df in ma_results.items():
            method_path = os.path.join(out_dir, f'{pool_name}_{timestamp}.csv')
            cols = [c for c in ['names', 'logfoldchanges', 'se', 'pvals_pooled',
                                'pval', 'n_studies', 'direction', 'pathways']
                    if c in meta_df.columns]
            meta_df[cols].to_csv(method_path, index=False)

        print(f"Per-method results saved to {out_dir}/")

    print()
    print("Done.")


if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()
