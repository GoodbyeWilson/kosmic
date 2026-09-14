# Format Converters for scRNA-seq Data
# Pure-Python functions for converting between data formats:
# - Dense matrix.txt.gz → h5ad
# - CSV/TSV count matrix → h5ad
# - MTX (Market Matrix) + genes/barcodes → h5ad
# - 10X folder (matrix.mtx + features.tsv + barcodes.tsv) → AnnData
# - Assemble h5ad from R-exported intermediate files
# - Format detection heuristics

import gzip
import os
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from scipy import sparse


# Format detection

def detect_file_format(filename: str) -> str:
    """Detect scRNA-seq file format from a filename. Returns a human-readable string."""
    f = filename.lower()
    if f.endswith('.rds') or f.endswith('.rds.gz') or f.endswith('.robj.gz') or f.endswith('.robj'):
        return 'Seurat RDS'
    elif f.endswith('.h5ad'):
        return 'h5ad (ready)'
    elif f.endswith('.h5seurat'):
        return 'h5Seurat'
    elif 'raw.tar' in f or f.endswith('.tar'):
        return '10X RAW'
    elif f.endswith('.h5'):
        return 'HDF5'
    elif f.endswith('.mtx.gz') or f.endswith('.mtx'):
        return 'MTX'
    elif f.endswith(('.csv.gz', '.csv', '.tsv.gz', '.tsv', '.txt.gz', '.txt')):
        if any(x in f for x in ['metadata', 'meta.', 'cellinfo', 'cluster', 'annotation']):
            return 'Metadata'
        elif any(x in f for x in ['count', 'matrix', 'expression']):
            return 'Count Matrix'
        else:
            return 'Other'
    else:
        return 'Other'


def looks_like_barcode(names: list) -> bool:
    """Check if a list of names looks like cell barcodes (long strings with dashes/numbers)."""
    return any('-' in str(n) and len(str(n)) > 10 for n in names)


def detect_matrix_orientation(df: pd.DataFrame) -> str:
    """Detect whether 'df' has rows=cells or rows=genes.

    Returns 'cells_x_genes' or 'genes_x_cells' (defaults to the
    latter when ambiguous, since GEO files usually ship that way).
    """
    sample_row_names = list(df.index[:10])
    sample_col_names = list(df.columns[:10])

    rows_are_cells = looks_like_barcode(sample_row_names)
    cols_are_cells = looks_like_barcode(sample_col_names)

    if rows_are_cells and not cols_are_cells:
        return 'cells_x_genes'
    elif cols_are_cells and not rows_are_cells:
        return 'genes_x_cells'
    else:
        return 'genes_x_cells'


# Dense matrix → h5ad

def convert_matrix_to_h5ad(matrix_file, sample_name=None, condition=None,
                           output_file=None, progress_callback=None):
    """Convert a dense matrix.txt(.gz) file to h5ad.

    Parameters
    ----------
    matrix_file : str
        Path to matrix.txt or matrix.txt.gz.
    sample_name : str, optional
        Cell-metadata sample name (auto-extracted from filename if None).
    condition : str, optional
        Condition label added to cell metadata.
    output_file : str, optional
        Output h5ad path. 'None' skips writing; '''' auto-generates from input.
    progress_callback : callable(msg), optional

    Returns
    -------
    anndata.AnnData
    """
    import scanpy as sc

    def _status(msg):
        if progress_callback:
            progress_callback(msg)

    if str(matrix_file).endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'

    _status(f"Reading {Path(matrix_file).name}...")

    with opener(str(matrix_file), mode) as f:
        header = f.readline().rstrip('\n\r').split('\t')
        cell_barcodes = header[1:]

        _status(f"Found {len(cell_barcodes)} cells...")

        genes = []
        expr_data = []

        for line_num, line in enumerate(f):
            parts = line.strip().split('\t')
            genes.append(parts[0])
            values = [float(x) if x else 0 for x in parts[1:]]
            expr_data.append(values)

            if (line_num + 1) % 5000 == 0:
                _status(f"Processed {line_num + 1} genes...")

    _status(f"Building matrix ({len(genes)} genes x {len(cell_barcodes)} cells)...")

    expr_matrix = np.array(expr_data, dtype=np.float32)
    expr_matrix = expr_matrix.T  # cells x genes
    expr_sparse = sparse.csr_matrix(expr_matrix)

    adata = sc.AnnData(X=expr_sparse)
    adata.obs_names = cell_barcodes
    adata.var_names = genes
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    if sample_name is None:
        sample_name = _extract_sample_name(matrix_file)
    adata.obs['sample'] = sample_name

    if condition:
        adata.obs['condition'] = condition

    if output_file is not None:
        if output_file == '':
            base_name = os.path.basename(str(matrix_file))
            base_name = base_name.replace('.txt.gz', '').replace('.txt', '')
            output_file = f"{base_name}.h5ad"
        _status(f"Writing {Path(output_file).name}...")
        adata.write_h5ad(output_file)

    return adata


# CSV/TSV → h5ad

def csv_to_h5ad(csv_path, output_path=None):
    """Convert a CSV/TSV count matrix to h5ad.

    Auto-detects delimiter, matrix orientation, and sample info from barcode patterns.

    Parameters
    ----------
    csv_path : str or Path
        CSV / TSV (gzip-compressed accepted).
    output_path : str or Path, optional
        Defaults to '<input_name>.h5ad' in the same directory.

    Returns
    -------
    anndata.AnnData
    """
    import anndata

    csv_path = Path(csv_path)

    # Peek at header to detect delimiter
    if str(csv_path).endswith('.gz'):
        with gzip.open(csv_path, 'rt') as f:
            header = f.readline().strip()
    else:
        with open(csv_path, 'r') as f:
            header = f.readline().strip()

    sep = '\t' if '\t' in header else ','

    df = pd.read_csv(csv_path, sep=sep, index_col=0)

    # Detect orientation
    orientation = detect_matrix_orientation(df)

    if orientation == 'cells_x_genes':
        cell_names = df.index.tolist()
        gene_names = df.columns.tolist()
        X = df.values
    else:
        cell_names = df.columns.tolist()
        gene_names = df.index.tolist()
        X = df.values.T

    X_sparse = sparse.csr_matrix(X.astype(np.float32))

    adata = anndata.AnnData(X=X_sparse)
    adata.obs_names = cell_names
    adata.var_names = gene_names

    # Try to extract sample info from barcode patterns
    if '_' in cell_names[0]:
        samples = [c.rsplit('_', 1)[0] if '_' in c else 'unknown' for c in cell_names]
        adata.obs['sample'] = samples
    elif cell_names[0].count('-') >= 1:
        suffixes = [c.split('-')[-1] if '-' in c else '1' for c in cell_names]
        if len(set(suffixes)) > 1 and len(set(suffixes)) < 50:
            adata.obs['sample'] = [f"sample_{s}" for s in suffixes]

    if output_path is None:
        csv_name = csv_path.stem
        if csv_name.endswith('.csv'):
            csv_name = csv_name[:-4]
        output_path = csv_path.parent / f"{csv_name}.h5ad"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)

    return adata


# MTX → h5ad

def mtx_to_h5ad(mtx_path, output_path=None, genes_path=None, barcodes_path=None,
                metadata_path=None, clusters_path=None):
    """Convert MTX (Matrix Market) + gene/barcode files to h5ad.

    Auto-detects genes and barcodes files in the same directory when not given.

    Parameters
    ----------
    mtx_path : str or Path
        '.mtx' or '.mtx.gz'.
    output_path : str or Path, optional
    genes_path, barcodes_path : str or Path, optional
        Paths to features / barcodes TSV files; auto-detected if omitted.
    metadata_path : str or Path, optional
        Cell metadata CSV/TSV.
    clusters_path : str or Path, optional
        Cluster / embedding coordinates CSV/TSV.

    Returns
    -------
    anndata.AnnData
    """
    import scanpy as sc
    from scipy.io import mmread

    mtx_path = Path(mtx_path)
    mtx_dir = mtx_path.parent

    # Auto-detect features file.
    if genes_path is None:
        for pattern in ['genes.tsv*', 'features.tsv*', '*genes*.txt*', '*features*.txt*']:
            matches = list(mtx_dir.glob(pattern))
            if matches:
                genes_path = matches[0]
                break

    # Auto-detect barcodes file.
    if barcodes_path is None:
        for pattern in ['barcodes.tsv*', '*barcodes*.txt*', '*cells*.txt*']:
            matches = list(mtx_dir.glob(pattern))
            if matches:
                barcodes_path = matches[0]
                break

    # Load the matrix.
    X = mmread(str(mtx_path))
    X = X.T.tocsr().astype(np.float32)  # genes x cells -> cells x genes

    # Load gene names.
    genes = None
    if genes_path and Path(genes_path).exists():
        genes_df = pd.read_csv(genes_path, sep='\t', header=None)
        if genes_df.shape[1] >= 2:
            genes = genes_df.iloc[:, 1].astype(str).tolist()
        else:
            genes = genes_df.iloc[:, 0].astype(str).tolist()

    # Load barcodes.
    barcodes = None
    if barcodes_path and Path(barcodes_path).exists():
        barcodes_df = pd.read_csv(barcodes_path, sep='\t', header=None)
        barcodes = barcodes_df.iloc[:, 0].astype(str).tolist()

    adata = sc.AnnData(X=X)

    if genes and len(genes) == adata.n_vars:
        adata.var_names = genes
    else:
        adata.var_names = [f"Gene_{i}" for i in range(adata.n_vars)]

    if barcodes and len(barcodes) == adata.n_obs:
        adata.obs_names = barcodes
    else:
        adata.obs_names = [f"Cell_{i}" for i in range(adata.n_obs)]

    # Optional cell metadata.
    if metadata_path and Path(metadata_path).exists():
        metadata_path = Path(metadata_path)
        sep = '\t' if str(metadata_path).endswith('.tsv') else ','
        try:
            meta = pd.read_csv(metadata_path, sep=sep, index_col=0)
        except (pd.errors.ParserError, UnicodeDecodeError):
            meta = pd.read_csv(metadata_path, sep='\t', index_col=0)

        common_cells = adata.obs_names.intersection(meta.index)
        if len(common_cells) > 0:
            meta = meta.loc[adata.obs_names.intersection(meta.index)]
            for col in meta.columns:
                if col in adata.obs_names:
                    continue
                adata.obs[col] = meta[col].reindex(adata.obs_names).values
        elif len(meta) == adata.n_obs:
            for col in meta.columns:
                adata.obs[col] = meta[col].values

    # Optional cluster / embedding coordinates.
    if clusters_path and Path(clusters_path).exists():
        clusters_path = Path(clusters_path)
        sep = '\t' if str(clusters_path).endswith('.tsv') else ','
        try:
            coords = pd.read_csv(clusters_path, sep=sep, index_col=0)
        except (pd.errors.ParserError, UnicodeDecodeError):
            coords = pd.read_csv(clusters_path, sep='\t', index_col=0)

        for embed_name in ['umap', 'tsne', 'UMAP', 'tSNE', 'TSNE']:
            cols = [c for c in coords.columns if embed_name.lower() in c.lower()]
            if len(cols) >= 2:
                embed_key = f"X_{embed_name.lower()}"
                coord_vals = coords[cols[:2]].values
                if len(coord_vals) == adata.n_obs:
                    adata.obsm[embed_key] = coord_vals.astype(np.float32)

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(output_path)

    return adata


# 10X folder → AnnData

def load_10x_folder(folder):
    """Load 10X Genomics data from a Cell Ranger output folder.

    Handles both gzipped (v3+) and uncompressed (v2) formats, including
    Cell Ranger v2's genome-subfolder nesting.

    Parameters
    ----------
    folder : str or Path
        Folder containing 'matrix.mtx[.gz]', 'features/genes.tsv[.gz]',
        'barcodes.tsv[.gz]'.

    Returns
    -------
    anndata.AnnData
    """
    import scanpy as sc
    import anndata
    from scipy.io import mmread

    folder = Path(folder)

    has_gzipped = (folder / "matrix.mtx.gz").exists()
    has_uncompressed = (folder / "matrix.mtx").exists()

    # Cell Ranger v2 nests files under a genome subfolder (e.g. mm10/, hg19/).
    if not has_gzipped and not has_uncompressed:
        subdirs = [d for d in folder.iterdir() if d.is_dir() and not d.name.startswith('.')]
        if len(subdirs) == 1:
            folder = subdirs[0]
            has_gzipped = (folder / "matrix.mtx.gz").exists()
            has_uncompressed = (folder / "matrix.mtx").exists()

    if has_gzipped:
        return sc.read_10x_mtx(str(folder), var_names='gene_symbols', cache=False)

    if has_uncompressed:
        matrix_file = folder / "matrix.mtx"

        if (folder / "features.tsv").exists():
            features_file = folder / "features.tsv"
        elif (folder / "genes.tsv").exists():
            features_file = folder / "genes.tsv"
        else:
            raise FileNotFoundError(f"No features.tsv or genes.tsv in {folder}")

        if (folder / "barcodes.tsv").exists():
            barcodes_file = folder / "barcodes.tsv"
        else:
            raise FileNotFoundError(f"No barcodes.tsv in {folder}")

        X = mmread(str(matrix_file)).T.tocsr().astype(np.float32)

        barcodes = pd.read_csv(barcodes_file, header=None, sep='\t')[0].tolist()

        features_df = pd.read_csv(features_file, header=None, sep='\t')
        if features_df.shape[1] >= 2:
            gene_symbols = features_df[1].tolist()
        else:
            gene_symbols = features_df[0].tolist()

        adata = anndata.AnnData(X=X)
        adata.obs_names = barcodes
        adata.var_names = gene_symbols
        adata.var_names_make_unique()

        return adata

    raise FileNotFoundError(f"No matrix.mtx or matrix.mtx.gz found in {folder}")


# R-export intermediate files → h5ad

def normalise_gene_index(adata):
    """Put gene symbols on ``var_names`` and Ensembl IDs in ``var['ensembl_id']``.

    Depositors index on whichever identifier suits them: CellxGene requires
    Ensembl IDs on the index with symbols in ``feature_name``; CellRanger and
    Seurat exports index on symbols, sometimes carrying Ensembl IDs alongside
    in ``gene_ids``. Combining studies joins on ``var_names``, so a mismatch
    is silent and total -- an Ensembl-indexed study shares *no* genes with a
    symbol-indexed one, and an inner join returns nothing.

    Symbols win as the index because everything downstream matches on them:
    reference annotation, marker genes, pathway sets. The Ensembl ID is kept
    as a column where the file supplies one (Seurat objects do not), since it
    is the stabler identifier for joins across reference versions.

    Modified in place; returns ``adata`` for chaining.
    """
    import pandas as pd

    def _looks_ensembl(values) -> bool:
        head = [str(v) for v in list(values)[:20]]
        return bool(head) and sum(v.startswith("ENSG") for v in head) > len(head) / 2

    if _looks_ensembl(adata.var_names):
        symbol_col = next((c for c in ("feature_name", "gene_symbol", "symbol")
                           if c in adata.var.columns), None)
        if symbol_col is not None:
            if "ensembl_id" not in adata.var.columns:
                adata.var["ensembl_id"] = list(adata.var_names)
            adata.var_names = pd.Index(
                adata.var[symbol_col].astype(str).values)
            adata.var_names_make_unique()
    elif "ensembl_id" not in adata.var.columns:
        # Symbol-indexed already: surface any Ensembl column under one name.
        ens_col = next((c for c in ("gene_ids", "gene_id", "feature_id")
                        if c in adata.var.columns
                        and _looks_ensembl(adata.var[c])), None)
        if ens_col is not None:
            adata.var["ensembl_id"] = adata.var[ens_col].astype(str).values
    return adata


def assemble_h5ad_from_r_export(intermediate_dir, prefix, h5ad_path):
    """Assemble an h5ad from R-exported (Seurat-conversion) intermediate files.

    Expected files in 'intermediate_dir':

    - '{prefix}_expression.mtx' -- sparse matrix (genes x cells)
    - '{prefix}_genes.txt'      -- gene names
    - '{prefix}_cells.txt'      -- cell barcodes
    - '{prefix}_metadata.csv'   -- cell metadata
    - '{prefix}_variable_features.txt'  (optional) HVG list
    - '{prefix}_reduction_*.csv' (optional) UMAP / tSNE / PCA embeddings

    Parameters
    ----------
    intermediate_dir : str or Path
    prefix : str
        File prefix (typically the sample name).
    h5ad_path : str or Path
        Output path.

    Returns
    -------
    anndata.AnnData
    """
    import scanpy as sc
    from scipy.io import mmread

    intermediate_dir = Path(intermediate_dir)
    h5ad_path = Path(h5ad_path)

    metadata = pd.read_csv(intermediate_dir / f"{prefix}_metadata.csv", index_col=0, low_memory=False)

    with open(intermediate_dir / f"{prefix}_genes.txt") as f:
        genes = [line.strip() for line in f]

    with open(intermediate_dir / f"{prefix}_cells.txt") as f:
        cells = [line.strip() for line in f]

    X = mmread(intermediate_dir / f"{prefix}_expression.mtx").T
    X = X.tocsr().astype(np.float32)

    adata = sc.AnnData(X=X)
    adata.obs_names = cells
    adata.var_names = genes

    for col in metadata.columns:
        adata.obs[col] = metadata[col].values

    # Standardise common column names.
    if 'Condition' in adata.obs.columns:
        adata.obs['condition'] = adata.obs['Condition']
    if 'Names' in adata.obs.columns:
        adata.obs['cell_type'] = adata.obs['Names']
    if 'orig.ident' in adata.obs.columns:
        adata.obs['sample'] = adata.obs['orig.ident']

    # Optional HVG list.
    var_features_file = intermediate_dir / f"{prefix}_variable_features.txt"
    if var_features_file.exists():
        try:
            with open(var_features_file) as f:
                var_features = [line.strip() for line in f if line.strip()]
            adata.var['highly_variable'] = adata.var_names.isin(var_features)
        except OSError:
            pass

    # Optional UMAP / t-SNE / PCA embeddings.
    reduction_files = list(intermediate_dir.glob(f"{prefix}_reduction_*.csv"))
    for red_file in reduction_files:
        try:
            red_name = red_file.stem.replace(f"{prefix}_reduction_", "")
            embeddings = pd.read_csv(red_file, index_col=0)

            obsm_key = f"X_{red_name}"
            if red_name.lower() == "umap":
                obsm_key = "X_umap"
            elif red_name.lower() == "tsne":
                obsm_key = "X_tsne"
            elif red_name.lower() == "pca":
                obsm_key = "X_pca"

            embeddings = embeddings.loc[cells]
            adata.obsm[obsm_key] = embeddings.values.astype(np.float32)
        except (OSError, KeyError, ValueError):
            pass

    h5ad_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(h5ad_path)

    return adata


# AnnData summary builder (shared by multiple workers)

def build_adata_summary(adata, file_path=None) -> dict:
    """Build a summary dict for a loaded AnnData object.

    Parameters
    ----------
    adata : anndata.AnnData
    file_path : str or Path, optional
        Used for 'file_size_mb'.

    Returns
    -------
    dict
        Keys: n_cells, n_genes, file_size_mb, obs_columns, var_columns,
        obsm_keys, obs_summary, is_normalized, has_raw, max_expression, sparsity.
    """
    summary = {
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'file_size_mb': 0.0,
        'obs_columns': list(adata.obs.columns),
        'var_columns': list(adata.var.columns),
        'obsm_keys': list(adata.obsm.keys()) if hasattr(adata, 'obsm') else [],
        'obs_summary': {},
        'is_normalized': False,
        'has_raw': hasattr(adata, 'raw') and adata.raw is not None,
    }

    if file_path and Path(file_path).exists():
        summary['file_size_mb'] = Path(file_path).stat().st_size / (1024 * 1024)

    # Heuristic: max < 20 suggests already log-normalised.
    try:
        if hasattr(adata.X, 'max'):
            max_val = adata.X.max()
        elif hasattr(adata.X, 'data') and len(adata.X.data) > 0:
            max_val = adata.X.data.max()
        else:
            max_val = 0
        summary['max_expression'] = float(max_val)
        summary['is_normalized'] = max_val < 20
    except Exception:
        summary['max_expression'] = 0
        summary['is_normalized'] = False

    if hasattr(adata.X, 'nnz'):
        total = adata.n_obs * adata.n_vars
        summary['sparsity'] = 1 - (adata.X.nnz / total) if total > 0 else 0
    else:
        summary['sparsity'] = None

    # Per-obs-column summary (dtype, n_unique, sample values).
    for col in list(adata.obs.columns):
        col_data = adata.obs[col]
        col_summary = {
            'dtype': str(col_data.dtype),
            'n_unique': int(col_data.nunique()),
            'sample_values': [str(v) for v in col_data.dropna().unique()[:5]],
        }
        if col_data.nunique() <= 50:
            col_summary['value_counts'] = {
                str(k): int(v) for k, v in col_data.value_counts().head(10).items()
            }
        else:
            col_summary['value_counts'] = {}
        summary['obs_summary'][col] = col_summary

    return summary


# h5ad loader with summary

def load_h5ad_with_summary(h5ad_path):
    """Load an h5ad and return '(AnnData, summary)'.

    Parameters
    ----------
    h5ad_path : str or Path

    Returns
    -------
    tuple of (anndata.AnnData, dict)
    """
    import scanpy as sc

    h5ad_path = Path(h5ad_path)
    adata = sc.read_h5ad(h5ad_path)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    # If any var_names are Ensembl IDs, swap to gene symbols from feature_name
    # column. Stash the original Ensembl IDs into 'ensembl_id' first so the
    # stable identifier survives the swap (provenance for reproducible joins).
    import pandas as pd
    from collections import Counter
    ensg_in_var = any(v.startswith('ENSG') for v in adata.var_names)
    if ensg_in_var and 'feature_name' in adata.var.columns:
        if 'ensembl_id' not in adata.var.columns:
            adata.var['ensembl_id'] = list(adata.var_names)
        adata.var_names = adata.var['feature_name'].values
        adata.var_names_make_unique()

    # Also rename raw slot var_names if any are still Ensembl IDs (used by CellTypist)
    if adata.raw is not None:
        ensg_in_raw = any(v.startswith('ENSG') for v in adata.raw.var_names)
        if ensg_in_raw and 'feature_name' in adata.raw.var.columns:
            if 'ensembl_id' not in adata.raw.var.columns:
                adata.raw._var['ensembl_id'] = list(adata.raw.var_names)
            base_names = list(adata.raw.var['feature_name'].values)
            counts = Counter()
            unique_names = []
            for name in base_names:
                counts[name] += 1
                unique_names.append(f"{name}-{counts[name]-1}" if counts[name] > 1 else name)
            adata.raw._var.index = pd.Index(unique_names)

    summary = build_adata_summary(adata, h5ad_path)
    return adata, summary


# Prefixed 10X file detection and loading

def find_10x_files(folder) -> List[dict]:
    """Find 10X files in a folder (standard and GSM-prefixed formats).

    Standard:  'matrix.mtx.gz', 'features.tsv.gz', 'barcodes.tsv.gz'.
    Prefixed:  '{prefix}_matrix.mtx.gz' etc. (or '.'-separated).

    Returns
    -------
    list of dict
        Each entry has 'folder', 'prefix', 'sample', and (for
        prefixed sets) 'matrix', 'features', 'barcodes'.
    """
    import re

    folder = Path(folder)

    # Standard layout (Cell Ranger v3+).
    standard_matrix = folder / "matrix.mtx.gz"
    standard_features = folder / "features.tsv.gz"
    standard_genes = folder / "genes.tsv.gz"
    standard_barcodes = folder / "barcodes.tsv.gz"

    if standard_matrix.exists() and standard_barcodes.exists():
        if standard_features.exists() or standard_genes.exists():
            return [{'folder': folder, 'prefix': None, 'sample': folder.name}]

    # GSM-prefixed layout.
    mtx_files = (
        list(folder.glob("*_matrix.mtx.gz")) + list(folder.glob("*_matrix.mtx")) +
        list(folder.glob("*.matrix.mtx.gz")) + list(folder.glob("*.matrix.mtx"))
    )

    samples = []
    for mtx_file in mtx_files:
        prefix_match = re.match(r'^(.+?)[._]matrix\.mtx', mtx_file.name)
        if prefix_match:
            sample_prefix = prefix_match.group(1)

            # Determine separator
            if (folder / f"{sample_prefix}.features.tsv.gz").exists() or \
               (folder / f"{sample_prefix}.barcodes.tsv.gz").exists():
                sep = "."
            else:
                sep = "_"

            features_file = folder / f"{sample_prefix}{sep}features.tsv.gz"
            genes_file = folder / f"{sample_prefix}{sep}genes.tsv.gz"
            barcodes_file = folder / f"{sample_prefix}{sep}barcodes.tsv.gz"

            if barcodes_file.exists() and (features_file.exists() or genes_file.exists()):
                samples.append({
                    'folder': folder,
                    'prefix': sample_prefix,
                    'sample': sample_prefix,
                    'matrix': mtx_file,
                    'features': features_file if features_file.exists() else genes_file,
                    'barcodes': barcodes_file,
                })

    return samples


def load_prefixed_10x(sample_info: dict):
    """Load 10X data from prefixed files (e.g. 'GSM123_matrix.mtx.gz').

    Parameters
    ----------
    sample_info : dict
        Must contain 'matrix', 'features', 'barcodes' Path objects.

    Returns
    -------
    anndata.AnnData
    """
    import gzip
    from anndata import AnnData
    from scipy.io import mmread

    matrix_path = Path(sample_info['matrix'])
    features_path = Path(sample_info['features'])
    barcodes_path = Path(sample_info['barcodes'])

    # Read matrix
    if str(matrix_path).endswith('.gz'):
        with gzip.open(matrix_path, 'rb') as f:
            matrix = mmread(f).T.tocsr()
    else:
        matrix = mmread(matrix_path).T.tocsr()

    # Read barcodes
    if str(barcodes_path).endswith('.gz'):
        with gzip.open(barcodes_path, 'rt') as f:
            barcodes = [line.strip().split('\t')[0] for line in f]
    else:
        with open(barcodes_path, 'r') as f:
            barcodes = [line.strip().split('\t')[0] for line in f]

    # Read features/genes
    if str(features_path).endswith('.gz'):
        with gzip.open(features_path, 'rt') as f:
            features = [line.strip().split('\t') for line in f]
    else:
        with open(features_path, 'r') as f:
            features = [line.strip().split('\t') for line in f]

    if len(features[0]) >= 2:
        gene_names = [f[1] if len(f) > 1 and f[1] else f[0] for f in features]
        gene_ids = [f[0] for f in features]
    else:
        gene_names = [f[0] for f in features]
        gene_ids = gene_names

    adata = AnnData(X=matrix)
    adata.obs_names = barcodes
    adata.var_names = gene_names
    adata.var['gene_ids'] = gene_ids
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    return adata


def load_10x_multi_sample(
    folder,
    output_path,
    progress_callback=None,
):
    """Load multiple 10X samples (standard or prefixed) and combine into one AnnData.

    Parameters
    ----------
    folder : str or Path
        Top-level folder containing 10X files or per-sample subfolders.
    output_path : str or Path
    progress_callback : callable(msg), optional

    Returns
    -------
    tuple of (anndata.AnnData, dict)
        Combined AnnData and summary dict.
    """
    import scanpy as sc

    folder = Path(folder)
    output_path = Path(output_path)

    # Search both the top folder and one level of subdirectories.
    samples = find_10x_files(folder)
    for subdir in folder.iterdir():
        if subdir.is_dir():
            samples.extend(find_10x_files(subdir))

    if not samples:
        raise FileNotFoundError(
            "No 10X files found. Expected matrix.mtx.gz, features.tsv.gz, "
            "barcodes.tsv.gz (or prefixed versions)."
        )

    if progress_callback:
        progress_callback(f"Found {len(samples)} sample(s)")

    adata_list = []
    sample_names = []

    for sample_info in samples:
        sample_name = sample_info['sample']
        if progress_callback:
            progress_callback(f"Loading {sample_name}...")

        try:
            if sample_info.get('prefix'):
                adata = load_prefixed_10x(sample_info)
            else:
                adata = sc.read_10x_mtx(
                    str(sample_info['folder']),
                    var_names='gene_symbols',
                    cache=False,
                )
            adata.obs['sample'] = sample_name
            adata_list.append(adata)
            sample_names.append(sample_name)
        except (OSError, ValueError, KeyError):
            continue

    if not adata_list:
        raise RuntimeError("Failed to load any samples.")

    if len(adata_list) > 1:
        adata_combined = sc.concat(adata_list, label='sample', keys=sample_names)
        adata_combined.obs_names_make_unique()
    else:
        adata_combined = adata_list[0]

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata_combined.write_h5ad(output_path)

    summary = build_adata_summary(adata_combined, output_path)
    summary['n_samples'] = len(sample_names)
    summary['samples'] = sample_names

    return adata_combined, summary


# TSV/CSV expression matrix → AnnData (chunked + featureCounts support)

def _extract_sample_name(file_path):
    """Extract a meaningful sample name from a filename."""
    p = Path(file_path)
    file_name = p.stem
    for ext in ['.tsv', '.csv', '.txt', '.gz']:
        if file_name.endswith(ext):
            file_name = file_name[:-len(ext)]

    generic_names = ['matrix', 'gene', 'counts', 'count', 'data', 'expression', 'raw', 'csv', 'tsv']

    sample_name = file_name
    if '_' in file_name:
        parts = file_name.split('_')
        # Skip GSM prefix, take first non-generic part.
        for part in parts:
            if part and not part.startswith('GSM') and part.lower() not in generic_names:
                sample_name = part
                break

    if sample_name.lower() in generic_names or sample_name == file_name:
        parent_name = p.parent.name
        if parent_name and parent_name.lower() not in ['processed', 'data', 'raw', 'output', 'results', 'matrices', 'raw_data', 'processed_data']:
            sample_name = parent_name

    return sample_name


def tsv_to_h5ad(tsv_path, progress_callback=None):
    """Convert a TSV/CSV expression matrix to AnnData.

    Auto-detects separator, featureCounts format, matrix orientation,
    and switches to chunked loading for files >100 MB. Does not write
    to disk.

    Parameters
    ----------
    tsv_path : str or Path
        TSV / CSV (gzip-compressed accepted).
    progress_callback : callable(msg), optional

    Returns
    -------
    anndata.AnnData
    """
    import scanpy as sc

    tsv_path = Path(tsv_path)

    def _status(msg):
        if progress_callback:
            progress_callback(msg)

    # Read first line to detect format.
    compression = 'gzip' if str(tsv_path).endswith('.gz') else None
    if compression == 'gzip':
        import gzip as _gzip
        with _gzip.open(tsv_path, 'rt') as f:
            first_line = f.readline()
    else:
        with open(tsv_path, 'r') as f:
            first_line = f.readline()

    skip_rows = 1 if first_line.startswith('#') else 0

    # Detect separator.
    file_str = str(tsv_path).lower()
    if '.csv' in file_str:
        sep = ','
    elif '.tsv' in file_str:
        sep = '\t'
    else:
        tabs = first_line.count('\t')
        commas = first_line.count(',')
        sep = '\t' if tabs > commas else ','

    _status(f"Detected separator: {'comma' if sep == ',' else 'tab'}")

    # >100 MB switches to chunked loading.
    file_size_mb = tsv_path.stat().st_size / (1024 * 1024)
    use_chunked = file_size_mb > 100

    if use_chunked:
        _status(f"Large file ({file_size_mb:.0f} MB) — using chunked loading...")
        X, gene_names, cell_names = _tsv_chunked_load(
            tsv_path, sep, skip_rows, compression, progress_callback
        )
    else:
        df = pd.read_csv(tsv_path, sep=sep, skiprows=skip_rows, compression=compression)
        _status(f"Raw shape: {df.shape[0]} rows x {df.shape[1]} cols")

        featurecounts_cols = ['Chr', 'Start', 'End', 'Strand', 'Length']
        is_featurecounts = all(col in df.columns for col in featurecounts_cols)

        if is_featurecounts:
            _status("Detected featureCounts format")
            X, gene_names, cell_names = _parse_featurecounts(df)
        else:
            df = df.set_index(df.columns[0])
            X, gene_names, cell_names = _orient_matrix(df, progress_callback)

    # Create AnnData
    adata = sc.AnnData(X=X)
    adata.obs_names = cell_names
    adata.var_names = gene_names
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    # Set sample name
    sample_name = _extract_sample_name(tsv_path)
    adata.obs['sample'] = sample_name
    adata.obs['file'] = tsv_path.name

    _status(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes (sample: {sample_name})")
    return adata


def _tsv_chunked_load(tsv_path, sep, skip_rows, compression, progress_callback=None):
    """Load a large TSV in chunks, building a sparse matrix incrementally."""
    def _status(msg):
        if progress_callback:
            progress_callback(msg)

    # Read header
    header_df = pd.read_csv(tsv_path, sep=sep, skiprows=skip_rows,
                            compression=compression, nrows=0)
    columns = header_df.columns.tolist()

    # Count rows
    _status("Counting rows...")
    if compression == 'gzip':
        import gzip as _gzip
        with _gzip.open(tsv_path, 'rt') as f:
            total_rows = sum(1 for _ in f) - 1 - skip_rows
    else:
        with open(tsv_path, 'r') as f:
            total_rows = sum(1 for _ in f) - 1 - skip_rows
    _status(f"Total rows: {total_rows:,}")

    # Build sparse matrix from chunks
    chunk_size = 5000
    chunks_data = []
    chunks_indices = []
    chunks_indptr = [0]
    gene_names = []
    current_ptr = 0

    chunk_iter = pd.read_csv(tsv_path, sep=sep, skiprows=skip_rows,
                             compression=compression, chunksize=chunk_size,
                             index_col=0, low_memory=True)

    rows_processed = 0
    for chunk in chunk_iter:
        rows_processed += len(chunk)
        _status(f"Loading: {rows_processed:,}/{total_rows:,} rows")

        gene_names.extend(chunk.index.astype(str).tolist())

        for row_idx in range(len(chunk)):
            row_data = chunk.iloc[row_idx].values.astype(np.float32)
            non_zero_mask = row_data != 0
            non_zero_vals = row_data[non_zero_mask]
            non_zero_cols = np.where(non_zero_mask)[0]
            chunks_data.extend(non_zero_vals)
            chunks_indices.extend(non_zero_cols)
            current_ptr += len(non_zero_vals)
            chunks_indptr.append(current_ptr)

        del chunk

    _status("Building sparse matrix...")
    n_genes = len(gene_names)
    n_cells = len(columns) - 1
    cell_names = columns[1:]

    X = sparse.csr_matrix(
        (np.array(chunks_data, dtype=np.float32),
         np.array(chunks_indices, dtype=np.int32),
         np.array(chunks_indptr, dtype=np.int32)),
        shape=(n_genes, n_cells)
    )

    _status("Transposing to cells x genes...")
    X = X.T.tocsr()

    return X, gene_names, cell_names


def _parse_featurecounts(df):
    """Parse a featureCounts formatted DataFrame."""
    if 'gene_name' in df.columns:
        gene_names = df['gene_name'].astype(str).tolist()
    else:
        gene_names = df.iloc[:, 0].astype(str).tolist()

    # Find where count data starts
    count_start_idx = 0
    for idx, col in enumerate(df.columns):
        if col in ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'gene_name']:
            count_start_idx = idx + 1
        else:
            break

    cell_cols = df.columns[count_start_idx:].tolist()
    cell_names = []
    for col in cell_cols:
        if '/' in col:
            col = col.split('/')[-1]
        if '.bam' in col:
            col = col.replace('.curated.minimap2.bam', '').replace('.bam', '')
        cell_names.append(col)

    count_data = df.iloc[:, count_start_idx:].values.T.astype(np.float32)
    X = sparse.csr_matrix(count_data)

    return X, gene_names, cell_names


def _orient_matrix(df, progress_callback=None):
    """Determine orientation of a standard expression matrix and return (X, genes, cells)."""
    sample_row_names = df.index[:10].astype(str).tolist()
    sample_col_names = df.columns[:10].astype(str).tolist()

    def _looks_like_barcode(name):
        acgt_count = sum(1 for c in name.upper() if c in 'ACGT')
        return len(name) >= 10 and acgt_count / len(name) > 0.5

    def _looks_like_gene(name):
        if not name:
            return False
        # Strip trailing barcode suffix like -1, -2 before checking
        core = name.rsplit('-', 1)[0] if '-' in name else name
        return ('Rik' in core or 'Gm' in core or
                core[0].isdigit() or core.isupper() or
                (any(c.isdigit() for c in core) and len(core) < 15))

    rows_are_genes = False

    if df.shape[0] > df.shape[1] * 3:
        rows_are_genes = True

    # Barcode detection takes priority over gene heuristic
    row_barcodes = any(_looks_like_barcode(n) for n in sample_row_names)
    col_barcodes = any(_looks_like_barcode(n) for n in sample_col_names)

    if col_barcodes:
        rows_are_genes = True
    elif row_barcodes:
        rows_are_genes = False
    elif any(_looks_like_gene(n) for n in sample_row_names if len(n) > 0):
        rows_are_genes = True

    if rows_are_genes:
        X = sparse.csr_matrix(df.values.T.astype(np.float32))
        gene_names = df.index.astype(str).tolist()
        cell_names = df.columns.astype(str).tolist()
    else:
        X = sparse.csr_matrix(df.values.astype(np.float32))
        gene_names = df.columns.astype(str).tolist()
        cell_names = df.index.astype(str).tolist()

    return X, gene_names, cell_names


# 10X HDF5 (.h5) → AnnData (with CellBender extras)

def h5_to_h5ad(h5_path):
    """
    Convert a 10X HDF5 (.h5) file to AnnData.

    Reads via scanpy.read_10x_h5() and extracts CellBender extras
    (latent embeddings, cell probability, RT efficiency, ambient expression)
    if present.

    Parameters
    ----------
    h5_path : str or Path
        Path to the .h5 file.

    Returns
    -------
    anndata.AnnData
        Does NOT write to disk or merge.
    """
    import scanpy as sc
    import h5py

    h5_path = Path(h5_path)
    adata = sc.read_10x_h5(str(h5_path))

    # Extract sample name
    sample_name = _extract_sample_name(h5_path)
    adata.obs['sample'] = sample_name
    adata.obs['source_file'] = h5_path.name
    adata.var_names_make_unique()

    # Try to extract CellBender extras
    try:
        with h5py.File(str(h5_path), 'r') as f:
            if 'matrix' in f:
                m = f['matrix']
                if 'latent_gene_encoding' in m:
                    latent = m['latent_gene_encoding'][:]
                    adata.obsm['X_cellbender_latent'] = latent.astype(np.float32)
                if 'latent_cell_probability' in m:
                    adata.obs['cellbender_cell_prob'] = m['latent_cell_probability'][:].astype(np.float32)
                if 'latent_RT_efficiency' in m:
                    adata.obs['cellbender_RT_efficiency'] = m['latent_RT_efficiency'][:].astype(np.float32)
                if 'ambient_expression' in m:
                    adata.var['ambient_expression'] = m['ambient_expression'][:].astype(np.float32)
    except Exception:
        pass  # CellBender extras are optional

    return adata
