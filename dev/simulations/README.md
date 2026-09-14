# Simulations

Benchmarking scripts for KOSMIC's meta-analysis pooling methods.

## Scripts

### `compare_gene_methods.py`
Full benchmark of all 13 pooling methods (slow methods included).
Uses simulated scRNA-seq data with known ground truth.

```bash
conda activate kirk
python scripts/compare_gene_methods.py
```

### `compare_top_methods.py`
Reduced benchmark of 4 top methods using the original (non-vectorized) code.
Useful for validating that fast versions produce the same results.

```bash
python scripts/compare_top_methods.py
```

### `compare_top_methods_fast.py`
Fast benchmark using vectorized numpy + parallel permutations.
Runs 8 methods in seconds. Optionally runs case-control permutation
(gold standard, takes hours).

```bash
# Fast methods only (seconds/minutes)
python scripts/compare_top_methods_fast.py

# Include case-control permutation (hours)
python scripts/compare_top_methods_fast.py --cc
```

## Documentation

- **`METHODS.md`** -- detailed writeup of each pooling method, permutation
  strategies (gene-label vs case-control), benchmark results, novel
  contributions, and limitations. Read this first.

## Key files

| File | Purpose |
|------|---------|
| `src/meta_analysis/pooling.py` | Per-gene helpers (hartung_knapp, wop_analytical, etc.) |
| `src/meta_analysis/pooling_core.py` | Shared matrix cores used by every per-method file |
| `src/meta_analysis/dl.py` / `reml.py` / `sumrank.py` / `rop.py` / `fisher.py` / `stouffer.py` | Per-method entry points |
| `src/meta_analysis/pooling_benchmark.py` | In-app benchmark dialog support |
| `src/meta_analysis/random_effects.py` | Core DL + SE estimation |

## Simulation parameters

- 15,000 genes, 500 truly DE (3.3%)
- 6 studies, 3-8 samples/condition, 80-400 cells/sample
- Negative binomial counts, 5% dropout
- Between-study heterogeneity tau = 0.15
- DE via KOSMIC's actual DESeq2 pipeline
