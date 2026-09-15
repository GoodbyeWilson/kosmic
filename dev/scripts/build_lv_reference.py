"""Build a KOSMIC centroid reference from a large integrated cardiac atlas, streaming.

The bundled builder ('kosmic.scrna.annotate.reference.build_reference_from_h5ad')
loads the whole h5ad; the HeartMap (42 GB) and Gao (96 GB) objects do not fit in
memory, so this script walks the CSR count matrix in row chunks and accumulates
what the reference needs: per-gene mean/variance of CP10k-normalised expression
(for scanpy 'seurat'-flavour HVG selection) and per-cell-type sums of
log1p(CP10k) expression (for the centroids). The output JSON has exactly the
format 'kosmic/scrna/annotate/reference.py' reads, plus a 'description' field
recording every filter and count so the build is reproducible.

Gene symbols are harmonised to HGNC with KOSMIC's own lookup before HVG
selection, summing columns that map to one approved symbol, so the reference
speaks the same symbols as a harmonised KOSMIC study.

Usage, from the repo root:

    python dev/scripts/build_lv_reference.py heartmap "D:/dcm source/SCP3689 (HeartMap)/HeartMap_V1.0_raw_counts.h5ad"
    python dev/scripts/build_lv_reference.py gao "D:/dcm source/gao/GSE290367_integrated_RNA_adata.h5ad"

The first argument names the atlas preset (which fixes the matrix slot,
label column, row filters and label mapping); the second is the h5ad. The
output goes to kosmic/reference/atlases/ by default, where the annotation
step discovers it; '--out' overrides. '--max-rows N' walks only the first
N rows, for a smoke test. The two presets, and the run on 2026-09-15:

  heartmap  Broad SCP3689 HeartMap_V1.0_raw_counts.h5ad (42 GB), X = CellBender
            raw counts, labels 'original_cell_type'. Filters: LV only; AtrialCM
            and Doublet dropped; Koenig paediatric donors dropped (a no-op on
            the released object). 1.21 M nuclei, 19 min.
  gao       GEO GSE290367_integrated_RNA_adata.h5ad (96 GB, decompressed),
            layers/counts, labels 'final_cell_type'. Filters: region LV;
            postnatal; not paediatric HF; age >= 18 where numeric. 2.01 M
            nuclei, 7 min.

Both keep every disease state. Memory is small (per-gene and per-type
sums only); the run is disk-bound.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", ".."))  # repo root: dev/scripts/<file>
from kosmic.scrna.annotate.reference import REFERENCES_DIR  # noqa: E402
from kosmic.scrna.inspect.gene_names import load_hgnc_lookup  # noqa: E402
N_TOP_GENES = 2000
CHUNK_ROWS = 20_000
MIN_CELLS_PER_TYPE = 10

# Koenig paediatric donors (DATASETS.md §4); matched as substrings of donor_id.
KOENIG_PAEDIATRIC = ("TWCM-11-93", "TWCM-13-181", "TWCM-13-198", "TWCM-13-235")

HEARTMAP_LABELS = {
    "VentricularCM": "Cardiomyocyte",
    "Fibroblast": "Fibroblast",
    "Endothelial": "Endothelial",
    "Pericyte": "Pericyte",
    "Macrophage": "Myeloid",
    "VSMC": "vSMC",
    "Lymphocyte": "Lymphoid",
    "Endocardial": "Endocardial",
    "Neuronal": "Neuronal",
    "Adipocyte": "Adipocyte",
    "LymphaticEndothelial": "LEC",
    "Mast": "Mast",
    "Epicardial": "Epicardial",
    # AtrialCM and Doublet are dropped.
}

# Per-atlas presets: which matrix and label column to read, and the output
# name. The input path is a command-line argument.
ATLASES = {
    "heartmap": dict(
        matrix="X",
        label_col="original_cell_type",
        name="HeartMap LV broad (Datar 2026)",
        source="doi:10.1038/s44161-026-00831-5; Broad SCP3689 HeartMap_V1.0_raw_counts.h5ad",
        out="heartmap_lv_broad.json.gz",
    ),
    "gao": dict(
        matrix="layers/counts",
        label_col="final_cell_type",
        name="Gao/Wu LV broad (Gao 2026)",
        source="doi:10.1186/s13059-026-04061-7; GEO GSE290367_integrated_RNA_adata.h5ad",
        out="gao_lv_broad.json.gz",
    ),
}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_categorical(obs, col):
    g = obs[col]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]], dtype=object)
        codes = g["codes"][:]
        out = np.full(codes.shape[0], "", dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        return out
    vals = g[:]
    return np.array([v.decode() if isinstance(v, bytes) else v for v in vals], dtype=object)


def read_numeric(obs, col):
    g = obs[col]
    if isinstance(g, h5py.Group):
        return None
    try:
        return np.asarray(g[:], dtype=float)
    except Exception:
        return None


def build_row_filter(atlas, obs):
    """Return (keep_mask, labels, filter_report)."""
    n = obs[obs.attrs.get("_index", "_index")].shape[0]
    keep = np.ones(n, dtype=bool)
    rep = {"n_total": int(n)}
    if atlas == "heartmap":
        organ = read_categorical(obs, "organ__ontology_label")
        m = organ == "heart left ventricle"
        rep["n_after_LV"] = int(m.sum())
        keep &= m
        raw = read_categorical(obs, "original_cell_type")
        labels = np.array([HEARTMAP_LABELS.get(x, "") for x in raw], dtype=object)
        m = labels != ""
        rep["n_after_dropping_AtrialCM_Doublet"] = int((keep & m).sum())
        keep &= m
        donor = read_categorical(obs, "donor_id")
        m = np.array([not any(p in d for p in KOENIG_PAEDIATRIC) for d in donor])
        rep["n_after_dropping_Koenig_paediatric"] = int((keep & m).sum())
        keep &= m
        rep["filters"] = [
            "organ__ontology_label == 'heart left ventricle'",
            "original_cell_type not in {AtrialCM, Doublet}",
            f"donor_id not containing {list(KOENIG_PAEDIATRIC)} (Koenig paediatric)",
            "all disease states retained",
        ]
        rep["label_map"] = HEARTMAP_LABELS
    else:
        region = read_categorical(obs, "region")
        m = region == "LV"
        rep["n_after_LV"] = int(m.sum())
        keep &= m
        age_status = read_categorical(obs, "age_status")
        m = age_status == "postnatal"
        rep["n_after_postnatal"] = int((keep & m).sum())
        keep &= m
        disease = read_categorical(obs, "disease")
        m = disease != "pediatric HF"
        rep["n_after_dropping_pediatric_HF"] = int((keep & m).sum())
        keep &= m
        age = read_numeric(obs, "age")
        if age is not None:
            m = ~(np.isfinite(age) & (age < 18))
            rep["n_after_age_ge_18"] = int((keep & m).sum())
            keep &= m
        labels = read_categorical(obs, "final_cell_type")
        rep["filters"] = [
            "region == 'LV'",
            "age_status == 'postnatal'",
            "disease != 'pediatric HF'",
            "age >= 18 where age is numeric",
            "all disease states retained",
        ]
    rep["n_kept"] = int(keep.sum())
    return keep, labels, rep


def harmonisation_matrix(var_names, lookup):
    """Sparse (n_old x n_new) matrix that renames and sums colliding symbols."""
    new_names = [lookup.get(g, g) for g in var_names]
    uniq = pd.unique(pd.Series(new_names))
    col = pd.Series(np.arange(len(uniq)), index=uniq)
    j = col[new_names].to_numpy()
    G = sp.csr_matrix((np.ones(len(var_names)), (np.arange(len(var_names)), j)),
                      shape=(len(var_names), len(uniq)))
    n_renamed = sum(1 for a, b in zip(var_names, new_names) if a != b)
    n_merged = len(var_names) - len(uniq)
    return G, list(uniq), n_renamed, n_merged


def select_hvg_seurat(mean_cp10k, var_cp10k, n_top, n_bins=20):
    """scanpy 'seurat' flavour on normalised (not logged) data."""
    mean = mean_cp10k.copy()
    var = var_cp10k.copy()
    mean[mean == 0] = 1e-12
    dispersion = var / mean
    dispersion[dispersion == 0] = np.nan
    dispersion = np.log(dispersion)
    mean = np.log1p(mean)
    df = pd.DataFrame({"mean": mean, "disp": dispersion})
    df["bin"] = pd.cut(df["mean"], bins=n_bins)
    grp = df.groupby("bin", observed=True)["disp"]
    mu = grp.transform("mean")
    sd = grp.transform("std")
    single = grp.transform("size") == 1
    sd = sd.where(~single, mu)  # scanpy: bins with one gene get std := mean
    df["disp_norm"] = (df["disp"] - mu) / sd
    df["disp_norm"] = df["disp_norm"].fillna(-np.inf)
    order = np.argsort(-df["disp_norm"].to_numpy(), kind="stable")
    return np.sort(order[:n_top])


def main(atlas, path, out=None, max_rows=None):
    cfg = ATLASES[atlas]
    path = str(path)
    log(f"opening {path}")
    lookup, _ = load_hgnc_lookup()
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        var = f["var"]
        var_names = [v.decode() if isinstance(v, bytes) else v
                     for v in var[var.attrs.get("_index", "_index")][:]]
        keep, labels, rep = build_row_filter(atlas, obs)
        log(f"rows kept: {rep['n_kept']:,} of {rep['n_total']:,}")

        G, new_genes, n_renamed, n_merged = harmonisation_matrix(var_names, lookup)
        n_genes = len(new_genes)
        log(f"genes: {len(var_names):,} -> {n_genes:,} after HGNC ({n_renamed:,} renamed, {n_merged:,} merged)")

        types = sorted({t for t in labels[keep]})
        type_idx = {t: i for i, t in enumerate(types)}
        lab_codes = np.array([type_idx.get(t, -1) for t in labels], dtype=np.int64)

        M = f[cfg["matrix"]]
        assert M.attrs.get("encoding-type", "csr_matrix") == "csr_matrix", "expected CSR"
        indptr = M["indptr"][:]
        n_rows = indptr.shape[0] - 1
        if max_rows:
            n_rows = min(n_rows, int(max_rows))
        data_ds, ind_ds = M["data"], M["indices"]

        sum_cp10k = np.zeros(n_genes)
        sumsq_cp10k = np.zeros(n_genes)
        type_sum_log = np.zeros((len(types), n_genes))
        type_n = np.zeros(len(types), dtype=np.int64)
        n_seen = 0
        t0 = time.time()

        for a in range(0, n_rows, CHUNK_ROWS):
            b = min(a + CHUNK_ROWS, n_rows)
            rows_keep = keep[a:b]
            if not rows_keep.any():
                continue
            s, e = indptr[a], indptr[b]
            X = sp.csr_matrix((data_ds[s:e].astype(np.float64), ind_ds[s:e], indptr[a:b + 1] - s),
                              shape=(b - a, len(var_names)))
            X = X[rows_keep]
            codes = lab_codes[a:b][rows_keep]
            X = X @ G                                  # HGNC rename + merge
            libsize = np.asarray(X.sum(axis=1)).ravel()
            libsize[libsize == 0] = 1.0
            X = sp.diags(1e4 / libsize) @ X            # CP10k
            sum_cp10k += np.asarray(X.sum(axis=0)).ravel()
            sumsq_cp10k += np.asarray(X.multiply(X).sum(axis=0)).ravel()
            X.data = np.log1p(X.data)                  # log1p(CP10k)
            for ti in np.unique(codes):
                m = codes == ti
                type_sum_log[ti] += np.asarray(X[m].sum(axis=0)).ravel()
                type_n[ti] += int(m.sum())
            n_seen += int(rows_keep.sum())
            if (a // CHUNK_ROWS) % 10 == 0:
                rate = n_seen / max(time.time() - t0, 1)
                log(f"  {n_seen:,}/{rep['n_kept']:,} nuclei  ({rate:,.0f}/s)")

    mean = sum_cp10k / n_seen
    var_ = (sumsq_cp10k - n_seen * mean ** 2) / max(n_seen - 1, 1)
    var_ = np.clip(var_, 0, None)
    hvg_idx = select_hvg_seurat(mean, var_, N_TOP_GENES)
    genes = [new_genes[i] for i in hvg_idx]

    valid = [i for i, t in enumerate(types) if type_n[i] >= MIN_CELLS_PER_TYPE]
    cell_types = [types[i] for i in valid]
    centroids = [(type_sum_log[i][hvg_idx] / type_n[i]).tolist() for i in valid]
    cell_counts = [int(type_n[i]) for i in valid]

    rep.update({
        "built": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input_file": path,
        "matrix": cfg["matrix"],
        "genes_in_object": len(var_names),
        "genes_after_hgnc": n_genes,
        "hgnc_renamed": n_renamed,
        "hgnc_merged": n_merged,
        "normalisation": "CP10k then log1p; centroids are mean log1p(CP10k) per type over all kept nuclei",
        "hvg": f"scanpy 'seurat' flavour on CP10k, top {N_TOP_GENES}, 20 mean bins",
        "script": "dev/scripts/build_lv_reference.py",
        "max_rows": int(max_rows) if max_rows else None,
        "cells_per_type": dict(zip(cell_types, cell_counts)),
    })
    ref = {
        "name": cfg["name"],
        "species": "human",
        "source": cfg["source"],
        "description": json.dumps(rep),
        "n_cell_types": len(cell_types),
        "n_genes": len(genes),
        "cell_types": cell_types,
        "genes": genes,
        "centroids": centroids,
        "cell_counts": cell_counts,
    }
    out = Path(out) if out else REFERENCES_DIR / cfg["out"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt") as fh:
        json.dump(ref, fh)
    log(f"wrote {out}  ({out.stat().st_size/1e6:.1f} MB)")
    for t, n in zip(cell_types, cell_counts):
        log(f"   {t:16} {n:>10,}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("atlas", choices=sorted(ATLASES), help="atlas preset")
    ap.add_argument("input", help="the atlas h5ad (CSR count matrix)")
    ap.add_argument("--out", help="output .json.gz (default: kosmic/reference/atlases/<preset>)")
    ap.add_argument("--max-rows", type=int, help="walk only the first N rows (smoke test)")
    args = ap.parse_args()
    if not Path(args.input).is_file():
        ap.error(f"not a file: {args.input}")
    main(args.atlas, args.input, out=args.out, max_rows=args.max_rows)
