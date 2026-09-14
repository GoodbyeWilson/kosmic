#!/usr/bin/env python
"""
Build lightweight reference files from popV HuggingFace h5ad downloads.

Downloads minified_ref_adata.h5ad files from HuggingFace popV repos,
subsamples, computes per-cell-type mean expression centroids, and saves
as compact JSON.gz files (~1-5 MB each).

Usage:
    python -u scripts/build_all_references.py

Downloads are ~15 GB for Tabula Sapiens All Cells + smaller for organ-specific.
Each is processed then deleted. Final output is <50 MB total.
"""

import sys
import os
import time
import gc

sys.stdout.reconfigure(line_buffering=True)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))  # repo root

import json
import gzip
import numpy as np
import pandas as pd
from pathlib import Path
import urllib.request
import tempfile

OUTPUT_DIR = Path(__file__).resolve().parents[2] / 'data' / 'references'  # repo root
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

HF_BASE = 'https://huggingface.co'

# What to download and build
REFERENCES = [
    # Full multi-organ atlases
    {
        'name': 'Tabula Sapiens',
        'hf_repo': 'popV/tabula_sapiens_All_Cells',
        'species': 'human',
        'source': 'Tabula Sapiens Consortium, Science 2022',
        'filename': 'tabula_sapiens.json.gz',
    },
    {
        'name': 'Tabula Muris',
        'hf_repo': 'popV/tabula_muris_All',
        'species': 'mouse',
        'source': 'Tabula Muris Consortium, Nature 2018',
        'filename': 'tabula_muris.json.gz',
    },
    # Organ-specific (from Tabula Sapiens subsets)
    {
        'name': 'Heart Cell Atlas',
        'hf_repo': 'popV/tabula_sapiens_Heart',
        'species': 'human',
        'source': 'Tabula Sapiens - Heart',
        'filename': 'heart_atlas.json.gz',
    },
    {
        'name': 'Human Lung Cell Atlas',
        'hf_repo': 'popV/tabula_sapiens_Lung',
        'species': 'human',
        'source': 'Tabula Sapiens - Lung',
        'filename': 'human_lung.json.gz',
    },
    {
        'name': 'PBMC / Blood',
        'hf_repo': 'popV/tabula_sapiens_Blood',
        'species': 'human',
        'source': 'Tabula Sapiens - Blood',
        'filename': 'azimuth_pbmc.json.gz',
    },
    {
        'name': 'Human Brain',
        'hf_repo': 'popV/tabula_sapiens_Neural',
        'species': 'human',
        'source': 'Tabula Sapiens - Neural',
        'filename': 'human_brain.json.gz',
    },
    {
        'name': 'Human Gut',
        'hf_repo': 'popV/tabula_sapiens_Large_Intestine',
        'species': 'human',
        'source': 'Tabula Sapiens - Large Intestine',
        'filename': 'human_gut.json.gz',
    },
    {
        'name': 'Kidney Cell Atlas',
        'hf_repo': 'popV/tabula_sapiens_Kidney',
        'species': 'human',
        'source': 'Tabula Sapiens - Kidney',
        'filename': 'kidney_atlas.json.gz',
    },
    {
        'name': 'Human Pancreas',
        'hf_repo': 'popV/tabula_sapiens_Pancreas',
        'species': 'human',
        'source': 'Tabula Sapiens - Pancreas',
        'filename': 'human_pancreas.json.gz',
    },
    # Mouse organ-specific
    {
        'name': 'Mouse Heart',
        'hf_repo': 'popV/tabula_muris_Heart_10x',
        'species': 'mouse',
        'source': 'Tabula Muris - Heart (10x)',
        'filename': 'mouse_heart.json.gz',
    },
    {
        'name': 'Mouse Brain',
        'hf_repo': 'popV/tabula_muris_Brain_non-myeloid_cells',
        'species': 'mouse',
        'source': 'Tabula Muris - Brain',
        'filename': 'mouse_brain.json.gz',
    },
]


def download_with_progress(url, dest_path):
    """Download a file with a text progress bar."""
    req = urllib.request.Request(url)
    req.add_header('User-Agent', 'Kirk-scRNA-pipeline/1.0')
    response = urllib.request.urlopen(req)
    total_size = int(response.headers.get('Content-Length', 0))
    total_mb = total_size / (1024 * 1024)

    block_size = 1024 * 1024  # 1 MB
    downloaded = 0
    last_print = 0

    with open(dest_path, 'wb') as f:
        while True:
            chunk = response.read(block_size)
            if not chunk:
                break
            f.write(chunk)
            downloaded += len(chunk)
            dl_mb = downloaded / (1024 * 1024)

            if dl_mb - last_print >= 50 or downloaded == total_size:
                if total_size > 0:
                    pct = downloaded / total_size * 100
                    print(f"      {dl_mb:.0f}/{total_mb:.0f} MB ({pct:.0f}%)")
                else:
                    print(f"      {dl_mb:.0f} MB")
                last_print = dl_mb

    return downloaded


def build_reference(h5ad_path, name, species, source, output_filename,
                    n_top_genes=2000, min_cells=10, max_cells_per_type=500):
    """Build a lightweight reference from a local h5ad file."""
    import scanpy as sc

    output_path = OUTPUT_DIR / output_filename

    print(f"    Loading h5ad...")
    t0 = time.time()
    adata = sc.read_h5ad(h5ad_path, backed='r')
    print(f"    Loaded: {adata.shape[0]:,} x {adata.shape[1]:,} ({time.time()-t0:.0f}s)")

    # Find cell type column
    cell_type_col = None
    for col in ['cell_type', 'celltype', 'CellType', 'cell_ontology_class', 'author_cell_type']:
        if col in adata.obs.columns:
            cell_type_col = col
            break

    if cell_type_col is None:
        print(f"    ERROR: No cell type column. Available: {list(adata.obs.columns)[:15]}")
        adata.file.close()
        return None

    print(f"    Cell type column: '{cell_type_col}'")

    type_counts = adata.obs[cell_type_col].value_counts()
    valid_types = type_counts[type_counts >= min_cells].index.tolist()
    cell_types = sorted(valid_types)
    print(f"    {len(cell_types)} cell types")

    # Subsample
    print(f"    Subsampling ({max_cells_per_type}/type max)...")
    rng = np.random.default_rng(42)
    sample_indices = []
    for ct in cell_types:
        ct_idx = adata.obs.index[adata.obs[cell_type_col] == ct].tolist()
        if len(ct_idx) > max_cells_per_type:
            ct_idx = rng.choice(ct_idx, size=max_cells_per_type, replace=False).tolist()
        sample_indices.extend(ct_idx)

    print(f"    Loading {len(sample_indices):,} cells into memory...")
    t0 = time.time()
    adata_sub = adata[sample_indices].to_memory()
    adata.file.close()
    del adata
    gc.collect()
    print(f"    In memory: {adata_sub.shape} ({time.time()-t0:.0f}s)")

    # Convert Ensembl IDs to gene symbols if feature_name column exists
    if 'feature_name' in adata_sub.var.columns:
        symbols = adata_sub.var['feature_name'].values
        n_valid = (symbols != '').sum() if hasattr(symbols, '__len__') else 0
        if n_valid > 0:
            # Use gene symbols as var_names, handle duplicates
            adata_sub.var_names = symbols.astype(str)
            adata_sub.var_names_make_unique()
            print(f"    Converted to gene symbols (feature_name column)")
        else:
            print(f"    WARNING: feature_name column exists but is empty")
    else:
        # Check if var_names are already gene symbols (not ENSG/ENSMUSG)
        sample_name = str(adata_sub.var_names[0])
        if sample_name.startswith('ENS'):
            print(f"    WARNING: var_names are Ensembl IDs but no feature_name column found")

    # Normalize if needed
    if hasattr(adata_sub.X, 'data') and len(adata_sub.X.data) > 0:
        max_val = float(adata_sub.X.data.max())
    elif hasattr(adata_sub.X, 'max'):
        max_val = float(adata_sub.X.max())
    else:
        max_val = 100

    if max_val > 20:
        sc.pp.normalize_total(adata_sub, target_sum=1e4)
        sc.pp.log1p(adata_sub)
        print(f"    Normalized (raw max was {max_val:.0f})")
    else:
        print(f"    Already normalized (max {max_val:.2f})")

    # HVGs
    print(f"    Finding {n_top_genes} HVGs...")
    try:
        sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_top_genes, flavor='seurat_v3')
    except Exception:
        try:
            sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_top_genes)
        except Exception:
            gene_means = np.asarray(adata_sub.X.mean(axis=0)).flatten()
            top_idx = np.argsort(gene_means)[-n_top_genes:]
            adata_sub.var['highly_variable'] = False
            adata_sub.var.iloc[top_idx, adata_sub.var.columns.get_loc('highly_variable')] = True

    hvg_mask = adata_sub.var['highly_variable']
    genes = list(adata_sub.var_names[hvg_mask])
    print(f"    {len(genes)} HVGs")

    # Centroids
    print(f"    Computing centroids...")
    centroids = []
    cell_counts = []
    for ct in cell_types:
        mask = adata_sub.obs[cell_type_col] == ct
        cell_counts.append(int(type_counts[ct]))
        X_ct = adata_sub[mask][:, hvg_mask].X
        if hasattr(X_ct, 'toarray'):
            X_ct = X_ct.toarray()
        centroids.append(np.asarray(X_ct.mean(axis=0)).flatten().tolist())

    ref = {
        'name': name, 'species': species, 'source': source,
        'n_cell_types': len(cell_types), 'n_genes': len(genes),
        'cell_types': cell_types, 'genes': genes,
        'centroids': centroids, 'cell_counts': cell_counts,
    }

    with gzip.open(output_path, 'wt') as f:
        json.dump(ref, f)

    size_mb = output_path.stat().st_size / (1024 * 1024)
    print(f"    DONE: {output_filename} ({size_mb:.1f} MB, {len(cell_types)} types, {len(genes)} genes)")

    del adata_sub
    gc.collect()
    return ref


def main():
    print(f"Building {len(REFERENCES)} reference atlases from HuggingFace popV\n")

    for i, ref_info in enumerate(REFERENCES, 1):
        name = ref_info['name']
        output_path = OUTPUT_DIR / ref_info['filename']

        print(f"[{i}/{len(REFERENCES)}] {name}")

        if output_path.exists():
            size_mb = output_path.stat().st_size / (1024 * 1024)
            print(f"  SKIP: already exists ({size_mb:.1f} MB)\n")
            continue

        # Download
        url = f"{HF_BASE}/{ref_info['hf_repo']}/resolve/main/minified_ref_adata.h5ad?download=true"
        print(f"  Downloading from {ref_info['hf_repo']}...")

        tmp_dir = tempfile.mkdtemp()
        tmp_path = os.path.join(tmp_dir, 'atlas.h5ad')

        try:
            t0 = time.time()
            download_with_progress(url, tmp_path)
            file_mb = os.path.getsize(tmp_path) / (1024 * 1024)
            print(f"    Downloaded: {file_mb:.0f} MB ({time.time()-t0:.0f}s)")

            build_reference(
                tmp_path, name=name, species=ref_info['species'],
                source=ref_info['source'], output_filename=ref_info['filename'],
            )
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
            if os.path.exists(tmp_dir):
                os.rmdir(tmp_dir)

        print()

    # Summary
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    total_size = 0
    for f in sorted(OUTPUT_DIR.glob('*.json.gz')):
        size_mb = f.stat().st_size / (1024 * 1024)
        total_size += size_mb
        with gzip.open(f, 'rt') as fh:
            meta = json.load(fh)
        print(f"  {f.name:30s} {size_mb:5.1f} MB  {meta['n_cell_types']:3d} types  {meta['n_genes']:5d} genes  ({meta['species']})")
    print(f"  {'TOTAL':30s} {total_size:5.1f} MB")


if __name__ == '__main__':
    main()
