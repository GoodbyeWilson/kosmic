"""Exact against fast case-control permutation, on real data (ADR-018).

Compares three calibrations of the pooled meta-analysis p-value on the DCM
endothelial cells, five studies:

  exact     - labels shuffled within each study, the study's own DE re-run
              on every shuffle (design ~sex_inferred + condition, dispersions
              re-estimated), then pooled with KOSMIC's matrix poolers.
  fast CC   - the shipped fast engine on the SAME shuffles, then pooled.
  analytic  - the pooler's own p-value, judged on the exact null.

The numbers in docs/ADR-018 and in the Statistical approach help page come
from this script (200 permutations; genes at BH FDR 5%: REML exact 4,427 /
fast 4,287, top-50% SumRank 305 / 193, Stouffer 6,789 / 6,462).

Read-only on the project; writes only to <outdir>. Two script-local patches:
single-core PyDESeq2 fits and apeGLM shrinkage skipped (nothing pools the
shrunk columns). About 170 s per permutation per core.

usage (repo root, venv active; set KOSMIC_DCM_PROJECT to the project folder):
    python dev/validation/cc_exact_vs_fast.py run <n_perms> <workers> <outdir>
    python dev/validation/cc_exact_vs_fast.py boot <outdir>      # paired bootstrap
    python dev/validation/cc_exact_vs_fast.py analyse <outdir>   # plain summary
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import pickle
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
warnings.filterwarnings("ignore")

ROOT = os.environ.get("KOSMIC_DCM_PROJECT", "D:/!!_KOSMIC DATA/New_QC_DecontX")
STUDIES = ["Chaffin", "Guo", "Koenig", "Reichart", "Youness"]
CT = "Endothelial_Cell"
METHODS = ["dl", "reml", "sumrank", "sumrank_top50", "gwop", "fisher", "stouffer"]
SEED = 20261119


def patch_pydeseq2():
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    o1, o2 = DeseqDataSet.__init__, DeseqStats.__init__

    def i1(self, *a, **k):
        k["n_cpus"] = 1
        o1(self, *a, **k)

    def i2(self, *a, **k):
        k["n_cpus"] = 1
        o2(self, *a, **k)
    DeseqDataSet.__init__, DeseqStats.__init__ = i1, i2
    DeseqStats.lfc_shrink = lambda self, *a, **k: None   # unshrunk columns are what is pooled


def proper_de(st, role_perm):
    """KOSMIC's per-study DE (design ~sex_inferred + condition) with permuted roles."""
    from kosmic.de.de_analysis import _run_deseq2
    n = len(st["genes"])
    lfc, se, p = (np.full(n, np.nan) for _ in range(3))
    sdf = st["sample_df"].copy()
    sdf["role"] = role_perm
    try:
        df = _run_deseq2(st["expr"], sdf, st["genes"], covariates=["sex_inferred"])
    except Exception:
        return lfc, se, p, False
    if df is None or len(df) == 0:
        return lfc, se, p, False
    pos = np.array([st["gene_pos"][g] for g in df["names"]])
    lfc[pos] = df["logfoldchanges"].to_numpy(float)
    se[pos] = df["se"].to_numpy(float)
    p[pos] = df["pvals"].to_numpy(float)
    return lfc, se, p, True


def matrices(per_study, ds_gene_indices, n_genes):
    n_ds = len(per_study)
    eff = np.full((n_genes, n_ds), np.nan)
    var = np.full((n_genes, n_ds), np.nan)
    pv = np.full((n_genes, n_ds), np.nan)
    for i, (l, p, s) in enumerate(per_study):
        gi = ds_gene_indices[i]
        ok = np.isfinite(l) & np.isfinite(s) & (s > 0) & np.isfinite(p)
        eff[gi[ok], i] = l[ok]
        var[gi[ok], i] = s[ok] ** 2
        pv[gi[ok], i] = p[ok]
    valid = np.isfinite(eff) & np.isfinite(var) & (var > 0)
    return eff, var, pv, valid


def pool(eff, var, pv, valid, min_studies):
    from kosmic.meta_analysis.cc_permutation import _MATRIX_POOLERS
    gene_mask = valid.sum(axis=1) >= min_studies
    out = {}
    for m in METHODS:
        fn, kind = _MATRIX_POOLERS[m]
        if kind == "weighted_rank":
            pp = fn(pv, eff, var, valid)
        elif kind == "rank":
            pp = fn(pv, eff, valid)
        elif kind == "pvals":
            pp = fn(eff, var, valid, pval_matrix=pv)
        else:
            pp = fn(eff, var, valid)
        x = np.full(len(pp), np.nan)
        ok = gene_mask & np.isfinite(pp) & (pp > 0)
        x[ok] = -np.log10(pp[ok])
        out[m] = x.astype(np.float32)
    return out


def build_state(outdir):
    import anndata as ad
    from kosmic import MIN_STUDIES
    from kosmic.meta_analysis.cc_permutation import load_pseudobulk
    from kosmic.meta_analysis.deseq2_fast import precompute_deseq2
    patch_pydeseq2()
    studies, log = [], []
    for s in STUDIES:
        base = f"{ROOT}/{s}/results/de_analysis/statistics/decontX_counts/{CT}"
        expr, conds, genes, samples, n_cells, roles = load_pseudobulk(f"{base}/{s}_{CT}_pseudobulk.csv")
        genes = list(genes)
        a = ad.read_h5ad(f"{ROOT}/{s}/processed_data/{s}.h5ad", backed="r")
        sx = (a.obs[["sample", "sex_inferred"]].drop_duplicates("sample")
              .set_index("sample")["sex_inferred"].astype(str))
        a.file.close()
        sex = [sx.get(str(x), "unknown") for x in samples]
        role = np.asarray(roles).astype(str)
        sdf = pd.DataFrame({"sample": [str(x) for x in samples], "role": role, "sex_inferred": sex})
        is_dis, two = role == "disease", np.isin(role, ["disease", "control"])
        t0 = time.time()
        pc = precompute_deseq2(expr, is_dis, two)
        st = dict(name=s, expr=np.asarray(expr, float), genes=genes,
                  gene_pos={g: i for i, g in enumerate(genes)}, sample_df=sdf,
                  role=role, is_dis=is_dis, two=two, precomp=pc)
        # real saved DE (observed labels, sex-adjusted) = observed per-study input to the pooling
        d = pd.read_csv(f"{base}/{s}_{CT}_DE_deseq2.csv", usecols=["names", "logfoldchanges", "se", "pvals"])
        l, se, p = (np.full(len(genes), np.nan) for _ in range(3))
        pos = np.array([st["gene_pos"].get(g, -1) for g in d["names"]])
        m = pos >= 0
        l[pos[m]], se[pos[m]], p[pos[m]] = (d["logfoldchanges"].to_numpy(float)[m],
                                            d["se"].to_numpy(float)[m], d["pvals"].to_numpy(float)[m])
        st["obs"] = (l, se, p)
        # harness check: our proper_de on the OBSERVED roles must reproduce the saved csv
        pl, pse, pp, ok = proper_de(st, role)
        both = np.isfinite(l) & np.isfinite(pl) & np.isfinite(p) & np.isfinite(pp)
        log.append(f"{s}: {len(role)} donors ({is_dis.sum()} dis/{(role=='control').sum()} ctl/"
                   f"{(role=='exclude').sum()} excl), sex levels {sorted(set(sex))}, genes {len(genes)}; "
                   f"harness vs saved csv: lfc r={np.corrcoef(l[both], pl[both])[0,1]:.6f}, "
                   f"-log10p r={np.corrcoef(-np.log10(p[both]), -np.log10(pp[both]))[0,1]:.6f} "
                   f"({both.sum()} genes) [{time.time()-t0:.0f}s]")
        print(log[-1], flush=True)
        studies.append(st)
    gene_to_idx = {}
    for st in studies:
        for g in st["genes"]:
            gene_to_idx.setdefault(g, len(gene_to_idx))
    for st in studies:
        st["gi"] = np.array([gene_to_idx[g] for g in st["genes"]])
    state = dict(studies=studies, n_genes=len(gene_to_idx), min_studies=MIN_STUDIES,
                 gene_names=list(gene_to_idx))
    # observed pooled -log10 p from the real per-study DE
    eff, var, pv, valid = matrices([(st["obs"][0], st["obs"][2], st["obs"][1]) for st in studies],
                                   [st["gi"] for st in studies], state["n_genes"])
    obs = pool(eff, var, pv, valid, MIN_STUDIES)
    np.savez(Path(outdir) / "observed.npz", genes=np.array(state["gene_names"]),
             n_valid=valid.sum(1), **obs)
    return state


_S = None


def _init(state):
    global _S
    patch_pydeseq2()
    _S = state


def one_perm(perm_i):
    S = _S
    studies = S["studies"]
    ps_proper, ps_fast, fails = [], [], 0
    t_proper = t_fast = 0.0
    from kosmic.meta_analysis.deseq2_fast import fast_deseq2_perm
    for si, st in enumerate(studies):
        rng = np.random.default_rng([SEED, perm_i, si])
        idx = rng.permutation(len(st["role"]))
        t0 = time.time()
        l, se, p, ok = proper_de(st, st["role"][idx])
        t_proper += time.time() - t0
        fails += (not ok)
        ps_proper.append((l, p, se))
        t0 = time.time()
        lf, pf, sf = fast_deseq2_perm(st["precomp"], st["is_dis"][idx], st["two"][idx])
        t_fast += time.time() - t0
        ps_fast.append((lf, pf, sf))
    gis = [st["gi"] for st in studies]
    out = {}
    for tag, per in (("proper", ps_proper), ("fast", ps_fast)):
        eff, var, pv, valid = matrices(per, gis, S["n_genes"])
        for m, x in pool(eff, var, pv, valid, S["min_studies"]).items():
            out[f"{tag}_{m}"] = x
    return perm_i, out, t_proper, t_fast, fails


def run(n_perms, workers, outdir):
    from concurrent.futures import ProcessPoolExecutor, as_completed
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    state = build_state(outdir)
    todo = [i for i in range(n_perms) if not (outdir / f"perm_{i:04d}.npz").exists()]
    print(f"{len(todo)} permutations to run on {workers} workers", flush=True)
    t0, done = time.time(), 0
    with ProcessPoolExecutor(max_workers=workers, initializer=_init, initargs=(state,)) as ex:
        futs = [ex.submit(one_perm, i) for i in todo]
        for f in as_completed(futs):
            i, out, tp, tf, fails = f.result()
            np.savez_compressed(outdir / f"perm_{i:04d}.npz", t_proper=tp, t_fast=tf, fails=fails, **out)
            done += 1
            print(f"perm {i}: proper {tp:.0f}s fast {tf:.2f}s fails {fails} | {done}/{len(todo)} "
                  f"elapsed {time.time()-t0:.0f}s", flush=True)
    print("RUN DONE", flush=True)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def analyse(outdir):
    from scipy.stats import ks_2samp
    outdir = Path(outdir)
    obs = np.load(outdir / "observed.npz", allow_pickle=True)
    files = sorted(outdir.glob("perm_*.npz"))
    perms = [np.load(f) for f in files]
    print(f"{len(perms)} permutations | proper DE {np.mean([p['t_proper'] for p in perms]):.0f} s "
          f"vs fast {np.mean([p['t_fast'] for p in perms]):.2f} s per permutation (5 studies) | "
          f"failed study fits: {int(sum(p['fails'] for p in perms))}")
    for m in METHODS:
        nul = {t: np.concatenate([p[f"{t}_{m}"][np.isfinite(p[f"{t}_{m}"])] for p in perms]).astype(float)
               for t in ("proper", "fast")}
        o = obs[m].astype(float)
        ok = np.isfinite(o)
        print(f"\n=== {m} === null values: proper {nul['proper'].size}, fast {nul['fast'].size}; "
              f"observed genes {ok.sum()}")
        print("  analytical p under PROPER null (should equal the nominal level if calibrated):")
        for thr in (0.05, 0.01, 0.001, 1e-4):
            x = -np.log10(thr)
            print(f"    P(p < {thr:g}): proper {np.mean(nul['proper'] > x):.5f} | "
                  f"fast {np.mean(nul['fast'] > x):.5f}")
        print(f"  KS(proper vs fast null): D={ks_2samp(nul['proper'], nul['fast']).statistic:.4f}")
        sp, sf = np.sort(nul["proper"]), np.sort(nul["fast"])

        def cal(sorted_null, x):
            n = sorted_null.size
            return np.maximum((n - np.searchsorted(sorted_null, x, side="left")) / n, 1.0 / (n + 1))
        cp, cf = cal(sp, o[ok]), cal(sf, o[ok])
        pa = 10.0 ** (-o[ok])
        qa, qp, qf = bh(pa), bh(cp), bh(cf)
        sa, spp, sff = qa < 0.05, qp < 0.05, qf < 0.05
        print(f"  observed genes with BH FDR<0.05: analytical {sa.sum()} | fast-CC {sff.sum()} | "
              f"proper-CC {spp.sum()} | fast&proper {(sff & spp).sum()}")
        print("  analytical p -> calibrated p (proper | fast) at the observed genes' quantiles:")
        for target in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
            j = np.argmin(np.abs(np.log10(pa) - np.log10(target)))
            print(f"    analytical {pa[j]:.1e} -> proper {cp[j]:.1e} | fast {cf[j]:.1e}")


def analyse_boot(outdir, n_boot=1000):
    """Paired bootstrap over permutations: the same resampled permutation set is used for the
    proper and the fast null, so differences between them carry a proper interval."""
    outdir = Path(outdir)
    obs = np.load(outdir / "observed.npz", allow_pickle=True)
    perms = [np.load(f) for f in sorted(outdir.glob("perm_*.npz"))]
    P = len(perms)
    print(f"{P} permutations | proper DE {np.mean([p['t_proper'] for p in perms]):.0f} s vs fast "
          f"{np.mean([p['t_fast'] for p in perms]):.2f} s per permutation (5 studies) | "
          f"failed study fits: {int(sum(p['fails'] for p in perms))} | bootstrap draws {n_boot}")
    step = 0.02
    edges = np.append(np.arange(0, 16 + step / 2, step), np.inf)
    nb = len(edges) - 1
    rng = np.random.default_rng(7)
    boots = [rng.integers(0, P, P) for _ in range(n_boot)]
    thr = [0.05, 0.01, 0.001]
    ci = lambda a: (np.percentile(a, 2.5), np.percentile(a, 97.5))

    for m in METHODS:
        o = obs[m].astype(float)
        ok = np.isfinite(o)
        o = o[ok]
        pa = 10.0 ** (-o)
        H = {}
        for t in ("proper", "fast"):
            H[t] = np.array([np.histogram(p[f"{t}_{m}"][np.isfinite(p[f"{t}_{m}"])], bins=edges)[0]
                             for p in perms], float)
        jb = np.minimum((o / step).astype(int), nb - 1)

        def surv(h):  # P(null -log10p >= lower edge of bin)
            tot = h.sum()
            s = np.cumsum(h[::-1])[::-1] / tot
            return s, tot

        def calp(h):
            s, tot = surv(h)
            return np.maximum(s[jb], 1.0 / (tot + 1))

        def tailfrac(h, x):
            s, _ = surv(h)
            return s[min(int(np.ceil(x / step)), nb - 1)]

        def stats(idx):
            hp, hf = H["proper"][idx].sum(0), H["fast"][idx].sum(0)
            cp, cf = calp(hp), calp(hf)
            return (np.array([tailfrac(hp, -np.log10(t)) for t in thr]),
                    np.array([tailfrac(hf, -np.log10(t)) for t in thr]),
                    int((bh(cp) < 0.05).sum()), int((bh(cf) < 0.05).sum()),
                    np.log10(cf / cp)[[np.argmin(np.abs(np.log10(pa) - np.log10(a)))
                                       for a in (1e-3, 1e-4, 1e-5)]])
        point = stats(np.arange(P))
        B = [stats(b) for b in boots]
        n_an = int((bh(pa) < 0.05).sum())
        print(f"\n=== {m} ===  observed genes {ok.sum()} | null values proper {int(H['proper'].sum())}")
        print("  null P(analytical p < x)   [95% CI]            fast - proper  [95% CI]")
        for k, t in enumerate(thr):
            bp = np.array([b[0][k] for b in B]); bf = np.array([b[1][k] for b in B])
            d = bf - bp
            print(f"    x={t:<6g} proper {point[0][k]:.5f} [{ci(bp)[0]:.5f}, {ci(bp)[1]:.5f}]   "
                  f"fast {point[1][k]:.5f}   diff {point[1][k]-point[0][k]:+.5f} [{ci(d)[0]:+.5f}, {ci(d)[1]:+.5f}]")
        bpn = np.array([b[2] for b in B]); bfn = np.array([b[3] for b in B])
        print(f"  genes at BH FDR<0.05: analytical {n_an} | proper-CC {point[2]} [{ci(bpn)[0]:.0f}, {ci(bpn)[1]:.0f}] | "
              f"fast-CC {point[3]} [{ci(bfn)[0]:.0f}, {ci(bfn)[1]:.0f}] | fast-proper {point[3]-point[2]:+d} "
              f"[{ci(bfn-bpn)[0]:+.0f}, {ci(bfn-bpn)[1]:+.0f}]")
        for k, a in enumerate((1e-3, 1e-4, 1e-5)):
            bl = np.array([b[4][k] for b in B])
            print(f"  analytical p~{a:g}: log10(fast-calibrated / proper-calibrated) = {point[4][k]:+.2f} "
                  f"[{ci(bl)[0]:+.2f}, {ci(bl)[1]:+.2f}]")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
    elif sys.argv[1] == "boot":
        analyse_boot(sys.argv[2])
    else:
        analyse(sys.argv[2])
