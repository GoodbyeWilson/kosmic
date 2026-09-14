"""
Parameter Sweep: k (datasets) x tau (heterogeneity)
====================================================
Runs a comprehensive benchmark sweep across different numbers of studies
and heterogeneity levels with multiple replicates per combination.

For the KOSMIC paper: shows how each pooling method behaves across realistic
conditions rather than a single fixed scenario.

Usage:
  conda activate kirk
  python scripts/sweep_k_tau.py               # full sweep (125 runs)
  python scripts/sweep_k_tau.py --quick        # small test (8 runs)
"""

import sys
import os
import time
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # repo root: dev/<this dir>/<file>
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMBA_NUM_THREADS'] = '1'
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

from dev.bench.simulation_benchmark import (
    SimulationConfig,
    generate_ground_truth,
    generate_simulated_datasets,
    run_de_on_simulated,
    get_fast_pooling_methods,
    evaluate,
)

# ======================================================================
# Sweep parameters
# ======================================================================

K_VALUES = [3, 4, 6, 8, 12]
TAU_VALUES = [0.0, 0.05, 0.15, 0.30, 0.50]
N_REPLICATES = 5

# Fixed simulation parameters (same as compare_top_methods_fast.py)
N_GENES = 15_000
N_TRUE_DE = 500
N_PERMS = 1_000


def run_single(k, tau, seed, methods):
    """Simulate k datasets with given tau, run DESeq2, pool, evaluate.

    Parameters
    ----------
    k : int
        Number of datasets.
    tau : float
        Between-study heterogeneity.
    seed : int
        RNG seed for reproducibility.
    methods : dict
        {name: callable(datasets) -> DataFrame} from get_fast_pooling_methods.

    Returns
    -------
    list of dict
        One dict per method with metrics + (k, tau, seed, method) columns.
    """
    config = SimulationConfig(
        n_genes=N_GENES,
        n_true_de=N_TRUE_DE,
        n_datasets=k,
        samples_per_group=(3, 8),
        cells_per_sample=(80, 400),
        nb_dispersion_range=(0.05, 0.8),
        dropout_rate=0.05,
        between_study_tau=tau,
        n_perms=N_PERMS,
        rng_seed=seed,
    )

    # Ground truth
    true_log2fc, de_indices, is_true_de, gene_names, _ = generate_ground_truth(config)

    # Simulate datasets
    datasets = generate_simulated_datasets(config, true_log2fc, de_indices, gene_names)

    # DESeq2
    ma_inputs = run_de_on_simulated(datasets, config)

    # Pool + evaluate
    rows = []
    for method_name, method_fn in methods.items():
        try:
            meta_df = method_fn(ma_inputs)
            # Normalise column names
            if 'pval' in meta_df.columns and 'pvals_pooled' not in meta_df.columns:
                meta_df['pvals_pooled'] = meta_df['pval']
            if 'mean_logfc' in meta_df.columns and 'logfoldchanges' not in meta_df.columns:
                meta_df['logfoldchanges'] = meta_df['mean_logfc']

            metrics = evaluate(meta_df, true_log2fc, is_true_de, gene_names)
            if metrics:
                metrics['k'] = k
                metrics['tau'] = tau
                metrics['seed'] = seed
                metrics['method'] = method_name
                rows.append(metrics)
        except Exception as exc:
            print(f"      FAILED {method_name}: {exc}")

    return rows


def main():
    quick = '--quick' in sys.argv

    if quick:
        k_values = [4, 6]
        tau_values = [0.0, 0.15]
        n_replicates = 2
    else:
        k_values = K_VALUES
        tau_values = TAU_VALUES
        n_replicates = N_REPLICATES

    total_combos = len(k_values) * len(tau_values) * n_replicates
    total_runs_label = f"{len(k_values)} k x {len(tau_values)} tau x {n_replicates} seeds = {total_combos}"

    print("=" * 70)
    print("  Parameter Sweep: k x tau")
    print("=" * 70)
    print(f"  k values:  {k_values}")
    print(f"  tau values: {tau_values}")
    print(f"  Replicates: {n_replicates}")
    print(f"  Total: {total_runs_label}")
    print(f"  Genes: {N_GENES} ({N_TRUE_DE} DE)")
    print(f"  Perms: {N_PERMS}")
    print()

    methods = get_fast_pooling_methods(n_perms=N_PERMS)
    print(f"  Methods ({len(methods)}): {', '.join(methods.keys())}")
    print()

    # Output directory
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(out_dir, exist_ok=True)
    timestamp = time.strftime('%Y%m%d_%H%M%S')
    csv_path = os.path.join(out_dir, f'sweep_k_tau_{timestamp}.csv')

    all_rows = []
    combo_idx = 0
    t0_total = time.time()

    for k in k_values:
        for tau in tau_values:
            for rep in range(n_replicates):
                combo_idx += 1
                seed = 42 + rep * 1000  # spread seeds

                # Progress
                elapsed = time.time() - t0_total
                rate = combo_idx / elapsed if elapsed > 1 else 0
                remaining = total_combos - combo_idx
                eta = remaining / rate if rate > 0 else 0

                print(f"  [{combo_idx}/{total_combos}] "
                      f"k={k}, tau={tau:.2f}, seed={seed}  "
                      f"(ETA {eta/60:.1f}min)", flush=True)

                t0 = time.time()
                rows = run_single(k, tau, seed, methods)
                dt = time.time() - t0

                n_methods_ok = len(rows)
                print(f"    -> {n_methods_ok} methods, {dt:.1f}s", flush=True)

                all_rows.extend(rows)

                # Save intermediate results (crash recovery)
                if all_rows:
                    df = pd.DataFrame(all_rows)
                    df.to_csv(csv_path, index=False)

    # Final summary
    total_time = time.time() - t0_total
    print()
    print("=" * 70)
    print(f"  DONE -- {total_combos} combinations in {total_time/60:.1f} min")
    print(f"  Results: {csv_path}")
    print("=" * 70)

    if not all_rows:
        print("  No results collected!")
        return

    df = pd.DataFrame(all_rows)
    df.to_csv(csv_path, index=False)

    # Print summary table: mean +/- std over seeds, grouped by (k, tau)
    print()
    print("  Summary (mean over seeds):")
    print()

    metric_cols = ['sensitivity', 'empirical_fdr', 'precision', 'auc_roc', 'auc_pr']
    available_metrics = [c for c in metric_cols if c in df.columns]

    for method in df['method'].unique():
        print(f"\n  --- {method} ---")
        mdf = df[df['method'] == method]
        summary = mdf.groupby(['k', 'tau'])[available_metrics].agg(['mean', 'std'])
        # Flatten column names
        summary.columns = [f'{m}_{s}' for m, s in summary.columns]
        print(f"  {'k':>3} {'tau':>5}  "
              f"{'Sens':>12} {'eFDR':>12} {'Prec':>12} "
              f"{'AUC-ROC':>12} {'AUC-PR':>12}")
        print(f"  {'---':>3} {'-----':>5}  "
              f"{'------------':>12} {'------------':>12} {'------------':>12} "
              f"{'------------':>12} {'------------':>12}")
        for (k_val, tau_val), row in summary.iterrows():
            parts = [f"  {k_val:>3} {tau_val:>5.2f}  "]
            for m in available_metrics:
                mean_val = row.get(f'{m}_mean', float('nan'))
                std_val = row.get(f'{m}_std', float('nan'))
                if np.isnan(std_val) or std_val == 0:
                    parts.append(f"{mean_val:>6.3f}      ")
                else:
                    parts.append(f"{mean_val:>6.3f}+/-{std_val:.3f}")
            print(''.join(parts))

    print(f"\n  Full results saved to: {csv_path}")


if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()
