"""
Benchmark: DL pooling x 3 calibrations, DESeq2 DE.

Generates synthetic data, runs DESeq2, pools with DL + analytical/gene-perm/CC-perm.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # repo root: dev/<this dir>/<file>
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import time

from dev.bench.simulation_benchmark import (
    SimulationConfig, generate_ground_truth, generate_simulated_datasets,
    run_de_on_simulated, run_welch_de_on_simulated,
    evaluate, compute_per_study_significance,
)
from kosmic.meta_analysis.pooling.dl import dl_fast
from kosmic.meta_analysis.cc_permutation import consensus_cc_permutation


def progress(step, current, total, detail=''):
    print(f"  [{step}] {current}/{total} {detail}")


def normalize(df):
    if 'pval' in df.columns and 'pvals_pooled' not in df.columns:
        df['pvals_pooled'] = df['pval']
    if 'mean_logfc' in df.columns and 'logfoldchanges' not in df.columns:
        df['logfoldchanges'] = df['mean_logfc']
    return df


def main():
    config = SimulationConfig(
        n_genes=15_000,
        n_true_de=500,
        n_datasets=6,
        samples_per_group=(3, 8),
        cells_per_sample=(80, 400),
        n_perms=1000,
        rng_seed=42,
    )

    print(f"=== Benchmark: {config.n_genes} genes, "
          f"{config.n_true_de} true DE, {config.n_datasets} datasets ===\n")

    # 1. Ground truth
    print("Generating ground truth...")
    true_log2fc, de_indices, is_true_de, gene_names, _ = generate_ground_truth(config)

    # 2. Simulate datasets
    print("Simulating datasets...")
    datasets = generate_simulated_datasets(
        config, true_log2fc, de_indices, gene_names, progress)

    # 3a. DESeq2 DE
    print("\nRunning DESeq2...")
    ma_deseq = run_de_on_simulated(datasets, config, progress)

    # 3b. Welch DE + pseudobulk
    print("\nRunning Welch t-test...")
    ma_welch, pseudobulk_data = run_welch_de_on_simulated(
        datasets, config, gene_names, eb_moderate=False, progress_cb=progress)

    # Per-study baselines
    study_deseq = compute_per_study_significance(ma_deseq, is_true_de, gene_names)
    study_welch = compute_per_study_significance(ma_welch, is_true_de, gene_names)
    print(f"\nPer-study baseline (DESeq2): union {study_deseq['union_sig']} sig "
          f"({study_deseq['union_tp']} TP, {study_deseq['union_fp']} FP)")
    print(f"Per-study baseline (Welch):  union {study_welch['union_sig']} sig "
          f"({study_welch['union_tp']} TP, {study_welch['union_fp']} FP)")

    # 4. Run methods -- all use matched DE for observed and null
    ms = 1

    from kosmic.meta_analysis.pooling.sumrank import sumrank_fast

    dl_fn = lambda ds, min_studies=ms: dl_fast(ds, min_studies=min_studies)
    sr_fn = lambda ds, min_studies=ms: sumrank_fast(
        ds, min_studies=min_studies)

    methods = {}

    # --- DESeq2 observed, analytical ---
    methods['DL_Analytical_DESeq2'] = ('deseq', lambda ds: dl_fn(ds))
    methods['SumRank_Analytical_DESeq2'] = ('deseq', lambda ds: sr_fn(ds))

    # --- Welch observed, analytical ---
    methods['DL_Analytical_Welch'] = ('welch', lambda ds: dl_fn(ds))
    methods['SumRank_Analytical_Welch'] = ('welch', lambda ds: sr_fn(ds))

    # --- Welch observed, CC perm (matched: Welch for both) ---
    # Runs DE once per permutation and pools DL + SumRank from the
    # shared result -- one call, both methods.
    def _run_dl_and_sumrank_cc():
        observed_dl = normalize(dl_fn(ma_welch))
        observed_sr = normalize(sr_fn(ma_welch))
        return consensus_cc_permutation(
            pseudobulk_data,
            method_meta_dfs={'dl': observed_dl, 'sumrank': observed_sr},
            n_perms=config.n_perms, seed=42,
            de_method='welch_cpm',
            min_studies=1)
    _cc_cache = {}
    def _get_cc():
        if not _cc_cache:
            _cc_cache.update(_run_dl_and_sumrank_cc())
        return _cc_cache
    methods['DL_CC_Welch']      = ('welch', lambda ds: _get_cc()['dl'])
    methods['SumRank_CC_Welch'] = ('welch', lambda ds: _get_cc()['sumrank'])

    results = {}
    for name, (input_type, fn) in methods.items():
        ma = ma_deseq if input_type.startswith('deseq') else ma_welch
        print(f"\nRunning {name}...")
        t0 = time.time()
        try:
            meta_df = normalize(fn(ma))
            elapsed = time.time() - t0
            metrics = evaluate(meta_df, true_log2fc, is_true_de, gene_names)
            if metrics:
                metrics['method'] = name
                metrics['time_s'] = round(elapsed, 1)
                results[name] = metrics
                print(f"  {elapsed:.1f}s | sig={metrics['n_significant']} | "
                      f"FP={metrics['n_false_positives']} | "
                      f"FDR={metrics['empirical_fdr']:.3f} | "
                      f"sens={metrics['sensitivity']:.3f} | "
                      f"AUC={metrics['auc_roc']:.3f}")
        except Exception as e:
            import traceback
            print(f"  FAILED: {e}")
            traceback.print_exc()

    # Summary table
    if results:
        print("\n" + "=" * 80)
        print("SUMMARY")
        print("=" * 80)
        df = pd.DataFrame(results.values())
        cols = ['method', 'n_significant', 'n_false_positives', 'false_positive_rate',
                'empirical_fdr', 'sensitivity', 'auc_roc', 'auc_pr',
                'top50_precision', 'top100_precision', 'top200_precision',
                'time_s']
        cols = [c for c in cols if c in df.columns]
        print(df[cols].to_string(index=False))


if __name__ == '__main__':
    main()
