"""Compare per-study DE methods on synthetic data with known truth.

Generates one synthetic study with planted DE genes, runs DESeq2 and
Welch t-test (raw + EB-moderated), and reports sensitivity / empirical
FDR / AUC for each. Answers: 'on KOSMIC's pseudobulk DE pipeline,
which method gives the best power-FDR tradeoff?'

Usage:

    python scripts/compare_de_methods.py
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # repo root: dev/<this dir>/<file>
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

from dev.bench.simulation_benchmark import (
    SimulationConfig, generate_ground_truth, generate_simulated_datasets,
    run_de_on_simulated, run_welch_de_on_simulated, evaluate,
)


def main():
    config = SimulationConfig(
        n_genes=5_000,
        n_true_de=200,
        n_datasets=1,
        samples_per_group=(8, 8),
        rng_seed=42,
    )

    print(f"Generating one study with {config.n_true_de} planted DE genes "
          f"({config.n_genes} total)...")
    true_log2fc, de_indices, gene_names = generate_ground_truth(config)
    is_true_de = np.zeros(len(gene_names), dtype=bool)
    is_true_de[de_indices] = True

    datasets = generate_simulated_datasets(
        config, true_log2fc, de_indices, gene_names)

    print('Running DESeq2...')
    deseq2_df = run_de_on_simulated(datasets, config)[0]

    print('Running Welch (raw)...')
    welch_raw_df = run_welch_de_on_simulated(
        datasets, config, gene_names, eb_moderate=False)[0][0]

    print('Running Welch (EB-moderated)...')
    welch_eb_df = run_welch_de_on_simulated(
        datasets, config, gene_names, eb_moderate=True)[0][0]

    rows = []
    for name, df in [
        ('DESeq2', deseq2_df),
        ('Welch (raw)', welch_raw_df),
        ('Welch + EB', welch_eb_df),
    ]:
        m = evaluate(df, true_log2fc, is_true_de, gene_names)
        rows.append({
            'method':         name,
            'n_significant':  m['n_significant'],
            'sensitivity':    round(m['sensitivity'], 3),
            'empirical_fdr':  round(m['empirical_fdr'], 3),
            'precision':      round(m['precision'], 3),
            'auc_roc':        round(m['auc_roc'], 3),
            'auc_pr':         round(m['auc_pr'], 3),
        })

    print()
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == '__main__':
    main()
