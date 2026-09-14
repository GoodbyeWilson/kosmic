"""
Compare simulated data properties against a real scRNA-seq dataset.
Checks sparsity, mean-variance, library sizes, expression distributions.

Uses sparse-aware operations to avoid OOM on large datasets.

Usage:
  python simulations/validate_simulation.py
"""

import numpy as np
import scipy.sparse as sp
import scanpy as sc
import sys

sys.path.insert(0, '.')
from dev.scripts.compare_pathway_methods import simulate_dataset, N_GENES

REAL_H5AD = "GSE183852/processed_data/GSE183852_DCM_Nuclei_normalized.h5ad"


def stats_block(name, X):
    """Print key stats for a (cells x genes) matrix. Sparse-aware."""
    is_sparse = sp.issparse(X)
    n_cells, n_genes = X.shape
    total_entries = n_cells * n_genes

    print(f"\n{'=' * 60}")
    print(f"  {name}")
    print(f"{'=' * 60}")
    print(f"  Shape:           {n_cells:,} cells x {n_genes:,} genes")
    print(f"  Storage:         {'sparse' if is_sparse else 'dense'}")

    # Sparsity
    if is_sparse:
        n_nonzero = X.nnz
    else:
        n_nonzero = np.count_nonzero(X)
    sparsity = (1 - n_nonzero / total_entries) * 100
    print(f"  Sparsity:        {sparsity:.1f}% zeros")

    # Library sizes (sum per cell) — works on sparse
    lib_sizes = np.asarray(X.sum(axis=1)).ravel()
    print(f"  Library size:    median={np.median(lib_sizes):.1f}  "
          f"mean={np.mean(lib_sizes):.1f}  "
          f"SD={np.std(lib_sizes):.1f}")
    print(f"                   range=[{lib_sizes.min():.1f}, {lib_sizes.max():.1f}]")

    # Gene means and variances — sparse-safe
    gene_means = np.asarray(X.mean(axis=0)).ravel()

    if is_sparse:
        # var = E[X^2] - E[X]^2, sparse-friendly
        X_sq = X.copy()
        X_sq.data **= 2
        gene_vars = np.asarray(X_sq.mean(axis=0)).ravel() - gene_means ** 2
    else:
        gene_vars = X.var(axis=0)

    print(f"  Gene means:      median={np.median(gene_means):.4f}  "
          f"mean={np.mean(gene_means):.4f}")
    print(f"                   range=[{gene_means.min():.5f}, {gene_means.max():.1f}]")
    print(f"  Gene variance:   median={np.median(gene_vars):.4f}  "
          f"mean={np.mean(gene_vars):.4f}")

    # Mean-variance relationship (binned)
    print(f"\n  Mean-variance relationship (binned):")
    print(f"  {'Mean range':<25} {'Avg Mean':>10} {'Avg Var':>10} {'Var/Mean':>10} {'n genes':>8}")
    print(f"  {'-' * 65}")
    bins = [0, 0.01, 0.1, 0.5, 1, 5, 50, 500, 10000]
    for i in range(len(bins) - 1):
        mask = (gene_means >= bins[i]) & (gene_means < bins[i + 1])
        n_in_bin = int(mask.sum())
        if n_in_bin > 0:
            avg_m = gene_means[mask].mean()
            avg_v = gene_vars[mask].mean()
            ratio = avg_v / avg_m if avg_m > 0 else 0
            print(f"  [{bins[i]:>7.2f}, {bins[i + 1]:>7.1f})  "
                  f"{avg_m:>10.4f} {avg_v:>10.4f} {ratio:>10.2f} {n_in_bin:>8}")

    # Non-zero expression distribution (sample to avoid OOM)
    print(f"\n  Non-zero expression distribution:")
    if is_sparse:
        nonzero_vals = X.data
    else:
        nonzero_vals = X[X > 0]

    if len(nonzero_vals) > 0:
        # Sample if too many
        if len(nonzero_vals) > 5_000_000:
            rng = np.random.default_rng(0)
            sample = rng.choice(nonzero_vals, 5_000_000, replace=False)
        else:
            sample = nonzero_vals

        percentiles = np.percentile(sample, [25, 50, 75, 90, 95, 99])
        print(f"  Non-zero count:  {len(nonzero_vals):,} ({len(nonzero_vals) / total_entries * 100:.1f}% of matrix)")
        print(f"  Percentiles:     25th={percentiles[0]:.3f}  50th={percentiles[1]:.3f}  "
              f"75th={percentiles[2]:.3f}")
        print(f"                   90th={percentiles[3]:.2f}  95th={percentiles[4]:.2f}  "
              f"99th={percentiles[5]:.2f}")
        print(f"  Max value:       {nonzero_vals.max():.2f}")


def main():
    # --- Real data ---
    print("Loading real dataset...")
    adata = sc.read_h5ad(REAL_H5AD)

    stats_block("REAL DATA: GSE183852 (normalised)", adata.X)

    if adata.raw is not None:
        stats_block("REAL DATA: GSE183852 (raw counts)", adata.raw.X)

    # --- Simulated data ---
    print("\nGenerating simulated dataset (20,000 genes)...")
    rng = np.random.default_rng(42)
    ds = simulate_dataset(0, rng)
    sim_X = ds['cell_level']
    stats_block("SIMULATED DATA (raw counts)", sim_X)

    # Pseudobulk summary
    pb_sum = ds['pseudobulk_sum']
    n_samples = pb_sum.shape[0]
    conds = ds['conditions']
    n_ctrl = sum(1 for c in conds if c == 'control')
    n_dis = sum(1 for c in conds if c == 'disease')
    print(f"\n  Pseudobulk: {n_samples} samples ({n_ctrl} control, {n_dis} disease)")
    print(f"  Cells per sample: {ds['cells_per_sample']}")
    lib = pb_sum.sum(axis=1)
    print(f"  Sum-agg lib sizes: median={np.median(lib):.0f}  "
          f"range=[{lib.min():.0f}, {lib.max():.0f}]")


if __name__ == '__main__':
    main()
