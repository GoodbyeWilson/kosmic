# dev/ — internal tooling

Nothing here is part of the application. It is excluded from the
installable package (`pyproject.toml`), not covered by the lint gate,
and may be rough. It exists so the two of us can benchmark methods,
validate against simulated data, and rebuild reference files.

Run everything from the repo root with the venv active. Each script
finds the repo root as three directories up from itself, so moving a
file within `dev/` breaks that — adjust the `sys.path.insert` line.

## bench/

Benchmark harnesses that other scripts import.

- `simulation_benchmark.py` — generates simulated multi-study datasets
  with planted effects; used by the `compare_*` and `sweep_*` scripts.
- `gwas_overlap.py` — nearest-gene GWAS overlap. Tested
  (`tests/test_gwas_overlap.py`) but **not yet wired into the app**; the
  Validation step's GWAS tab is a placeholder pointing here.
- `proteomics_validation.py` — Olink / proteomics cross-checks.

## scripts/

Method comparisons and reference builders.

- `compare_de_methods.py`, `compare_gene_methods.py`,
  `compare_pathway_methods.py`, `compare_top_methods*.py` — head-to-head
  runs of the pooling and DE methods on simulated data. Write CSVs
  alongside themselves or into `results/`.
- `benchmark_dl_welch.py`, `sweep_k_tau.py` — DL vs Welch, and how the
  methods behave as study count k and heterogeneity τ² vary.
- `build_reference.py`, `build_all_references.py` — rebuild the bundled
  reference files under `kosmic/reference/`.
- `generate_gene_positions.py` — regenerates the GRCh37 gene-position
  table used by the GWAS overlap.
- `capture_help_screenshots.py` — regenerates every screenshot the
  help pages embed, from a live window on a real project. Run after
  any UI change; never hand-capture.
- `smoke_drive.py` — clicks through every workspace and step on a
  real project and reports any exception with the step it happened
  on. The last thing to run before a release.

## simulations/

Validation and one-off data gathering.

- `validate_simulation.py`, `validate_deseq2_fast.py` — check the
  fast DESeq2 path and the simulation framework against ground truth.
  `METHODS.md` describes the simulation design;
  `gene_method_benchmark_results.csv` is the last full run.
- `find_hf_geo_datasets.py` — query GEO for heart-failure snRNA-seq
  datasets.
- `query_covid_pbmc_datasets.py`, `merge_covid_datasets.py` — assemble
  the COVID-19 PBMC cohort that SumRank was published on, to reproduce
  that benchmark against our pooling methods.
- `heart_umap.py` — renders the heart-shaped UMAP used as the logo.

### Large data (gitignored, local only)

- `simulations/results/` — outputs of the scripts above; regenerate by
  re-running them.
- `simulations/Heart Logo/` — the raw and processed data behind the
  logo UMAP (~1.3 GB). Not needed for anything else. If you need it,
  ask rather than regenerate.
