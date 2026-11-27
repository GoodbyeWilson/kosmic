"""Exact against vectorised case-control permutation, through the shipped code (ADR-019).

Runs `consensus_cc_permutation` twice on the same studies and seed -- once with
the exact engine (PyDESeq2 re-run on every shuffle) and once with the vectorised
engine -- so both use the same stratified shuffles, then compares:

  * the null values of each permutation (how many differ, and by how much);
  * the null tail fractions at p < 0.05, 0.01 and 0.001;
  * the genes called at BH FDR 5% from the calibrated p-values.

The permutation p-value of a gene depends on the whole pooled null, so the
tail fractions and the gene counts are the comparison that matters; single null
values differ at float precision (PyDESeq2 and the array code agree to about
1e-7 relative) and, for a small fraction of genes, more.

Read-only on the project; writes only to <outdir>.

usage (repo root, venv active):
    set KOSMIC_DCM_PROJECT=<project folder holding the study folders>
    python dev/validation/cc_exact_vs_vec.py <cell type folder name> <n_perms> <outdir>
    e.g. ... Ventricular_Cardiomyocyte 14 out_cm

The exact engine takes minutes per permutation for a large cell type; start with
about 14 permutations.
"""
import os
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
warnings.filterwarnings("ignore")

from kosmic import META_IF_FILTER_STAT  # noqa: E402
from kosmic.meta_analysis.cc_permutation import consensus_cc_permutation  # noqa: E402
from kosmic.meta_analysis.independent_filter import (  # noqa: E402
    apply_filter, filter_datasets_to_universe)
from kosmic.meta_analysis.io import load_de_results, load_pseudobulk_set  # noqa: E402
from kosmic.meta_analysis.pooling_dispatch import get_analytical_pool_fn  # noqa: E402
from kosmic.numerical import bh_fdr  # noqa: E402

ROOT = Path(os.environ.get("KOSMIC_DCM_PROJECT", "D:/!!_KOSMICDATA2/New_QC_DecontX"))
STUDIES = ["Chaffin", "Guo", "Koenig", "Reichart", "Youness"]
METHODS = ["dl", "reml", "sumrank", "sumrank_top50", "gwop", "fisher", "stouffer"]
MIN_STUDIES = 3


def main(cell_type, n_perms, outdir):
    folder = {s: ROOT / s / "results" / "de_analysis" / "statistics" / "decontX_counts" / cell_type
              for s in STUDIES}
    entries = [{"name": f"{s}_{cell_type}",
                "df": load_de_results(str(folder[s] / f"{s}_{cell_type}_DE_deseq2.csv"),
                                      dataset_name=f"{s}_{cell_type}")} for s in STUDIES]
    universe = apply_filter(entries, threshold=0.05, stat=META_IF_FILTER_STAT,
                            min_studies=MIN_STUDIES)
    dfs = [e["df"] for e in filter_datasets_to_universe(entries, universe)]
    observed = {k: get_analytical_pool_fn(k, MIN_STUDIES)(dfs, min_studies=MIN_STUDIES)
                for k in METHODS}
    pb = load_pseudobulk_set([folder[s] / f"{s}_{cell_type}_pseudobulk.csv" for s in STUDIES],
                             labels=[e["name"] for e in entries])
    out = Path(outdir)
    result = {}
    for engine in ("deseq2_vec", "deseq2_exact"):
        result[engine] = consensus_cc_permutation(
            pb, observed, n_perms=n_perms, seed=42, de_method=engine, min_studies=MIN_STUDIES,
            status_cb=print, resume_dir=out / engine, n_workers=max(1, (os.cpu_count() or 2) - 2))
    print(f"\n{'method':14s} {'FDR5% exact':>11s} {'vec':>6s} | tail p<.05 exact/vec | "
          f"p<.001 exact/vec | null values differing by >0.01")
    for m in METHODS:
        calls = {e: int((bh_fdr(result[e][m]['pvals_pooled'].to_numpy(float)) < 0.05).sum())
                 for e in result}
        nulls = {}
        for e in result:
            chunks = []
            for i in range(n_perms):
                with np.load(out / e / f"perm_{i:05d}.npz") as z:
                    chunks.append(z[m].astype(float))
            nulls[e] = np.concatenate(chunks)
        a, b = nulls["deseq2_exact"], nulls["deseq2_vec"]
        tail = lambda v, p: float((v >= -np.log10(p)).mean())  # noqa: E731
        print(f"{m:14s} {calls['deseq2_exact']:11d} {calls['deseq2_vec']:6d} | "
              f"{tail(a, .05):.5f}/{tail(b, .05):.5f}    | {tail(a, .001):.6f}/{tail(b, .001):.6f} | "
              f"{float((np.abs(a - b) > 0.01).mean()):.5f}")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
