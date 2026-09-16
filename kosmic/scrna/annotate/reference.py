# Reference-Based Cell Type Annotation
# Lightweight label transfer using pre-computed reference centroids.
#
# Reference format (JSON):
# {
#     "name": "Human Lung Cell Atlas",
#     "species": "human",
#     "source": "doi:10.1038/s41591-023-02327-2",
#     "n_cell_types": 61,
#     "n_genes": 2000,
#     "cell_types": ["AT1", "AT2", ...],
#     "genes": ["GENE1", "GENE2", ...],
#     "centroids": [[mean_expr_type1_gene1, ...], ...],  # n_types x n_genes
#     "cell_counts": [1234, 5678, ...],  # cells per type in source atlas
# }

import json
import gzip
from kosmic.scrna.inspect.roles import cluster_labels
import numpy as np
import pandas as pd
from importlib.resources import files
from pathlib import Path


REFERENCES_DIR = Path(str(files("kosmic.reference.atlases")))

# Built-in reference catalog. Empty on purpose: every '.json.gz' in
# kosmic/reference/atlases/ -- the two shipped LV references included -- is
# discovered by 'list_available_references' and describes itself (name,
# source, and a 'description' recording the build). Build another with
# dev/scripts/build_reference.py or build_lv_reference.py.
BUILTIN_REFERENCES = {}


def list_available_references():
    """Return list of available reference names (built-in + any custom in kosmic/reference/atlases/)."""
    available = {}
    for name, info in BUILTIN_REFERENCES.items():
        ref_path = REFERENCES_DIR / info['file']
        available[name] = {
            **info,
            'installed': ref_path.exists(),
            'builtin': True,
        }

    # Also scan for custom .json.gz files
    if REFERENCES_DIR.exists():
        for f in REFERENCES_DIR.glob('*.json.gz'):
            if f.name not in {v['file'] for v in BUILTIN_REFERENCES.values()}:
                try:
                    meta = _load_reference_metadata(f)
                    available[meta.get('name', f.stem)] = {
                        'file': f.name,
                        'species': meta.get('species', 'unknown'),
                        'description': meta.get('description', ''),
                        'source': meta.get('source', ''),
                        'installed': True,
                        'builtin': False,
                    }
                except Exception:
                    pass
        for f in REFERENCES_DIR.glob('*.json'):
            if f.name not in {v['file'] for v in BUILTIN_REFERENCES.values()}:
                try:
                    meta = _load_reference_metadata(f)
                    available[meta.get('name', f.stem)] = {
                        'file': f.name,
                        'species': meta.get('species', 'unknown'),
                        'description': meta.get('description', ''),
                        'source': meta.get('source', ''),
                        'installed': True,
                        'builtin': False,
                    }
                except Exception:
                    pass

    return available


def _load_reference_metadata(path):
    """Load just the metadata fields from a reference file (not centroids)."""
    path = Path(path)
    if path.suffix == '.gz':
        with gzip.open(path, 'rt') as f:
            # Read just enough to get metadata — for large files, parse incrementally
            data = json.load(f)
    else:
        with open(path) as f:
            data = json.load(f)
    return {k: v for k, v in data.items() if k not in ('centroids', 'genes')}


def load_reference(name_or_path):
    """
    Load a reference file by name (from catalog) or file path.

    Returns
    -------
    dict with keys: name, species, cell_types, genes, centroids (numpy array),
    cell_counts, source, n_cell_types, n_genes
    """
    path = Path(name_or_path)

    # A name rather than a path: resolve through whatever is discovered in
    # kosmic/reference/atlases/ (the same list the Annotate dropdown
    # shows), then the built-in catalogue.
    if not path.exists():
        available = list_available_references()
        if name_or_path in available:
            path = REFERENCES_DIR / available[name_or_path]['file']
        elif name_or_path in BUILTIN_REFERENCES:
            path = REFERENCES_DIR / BUILTIN_REFERENCES[name_or_path]['file']

    if not path.exists():
        raise FileNotFoundError(
            f"Reference not found: {name_or_path}. Available: "
            + (", ".join(sorted(list_available_references())) or "none"))

    if path.suffix == '.gz':
        with gzip.open(path, 'rt') as f:
            data = json.load(f)
    else:
        with open(path) as f:
            data = json.load(f)

    data['centroids'] = np.array(data['centroids'], dtype=np.float32)
    return data


def annotate_by_reference(adata, reference, cluster_col='leiden',
                          min_correlation=0.3, progress_callback=None):
    """
    Annotate cells by correlating expression with reference centroids.

    For each cluster, computes mean expression over shared genes, then
    correlates with each reference cell type centroid. The best-correlating
    type is assigned if above min_correlation.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have been clustered. Uses raw expression if available.
    reference : dict
        Loaded reference (from load_reference).
    cluster_col : str
        Cluster column name.
    min_correlation : float
        Minimum Pearson correlation to assign a type (else 'Unknown').
    progress_callback : callable, optional
        Called with status strings.

    Returns
    -------
    tuple of (cluster_types, annotation_details)
        cluster_types: {cluster_id: cell_type_name}
        annotation_details: {cluster_id: {correlation, runner_up, ...}}
    """
    from kosmic.scrna.inspect.detection import detect_species

    ref_genes = reference['genes']
    ref_centroids = reference['centroids']  # n_types x n_genes
    ref_cell_types = reference['cell_types']
    ref_species = reference.get('species', 'human')

    # Expression on the reference's scale: the centroids are mean
    # log1p(CP10k) per type, so the query is compared as log1p(CP10k)
    # too. (Earlier versions correlated mean raw counts from '.raw'
    # against the log-scale centroids.)
    from kosmic.scrna.counts import log_normalised
    expr_X, expr_var_names = log_normalised(adata)

    # Detect species of query data
    query_species = detect_species(expr_var_names)

    # Convert reference gene names to match query species if needed
    if ref_species == 'human' and query_species == 'mouse':
        from kosmic.scrna.inspect.detection import format_gene_for_species
        ref_genes_fmt = [format_gene_for_species(g, 'mouse') if g else g
                         for g in ref_genes]
    elif ref_species == 'mouse' and query_species == 'human':
        ref_genes_fmt = [g.upper() for g in ref_genes]
    else:
        ref_genes_fmt = ref_genes

    # Find shared genes
    expr_gene_set = set(expr_var_names)
    shared_mask = [g in expr_gene_set for g in ref_genes_fmt]
    shared_ref_indices = [i for i, m in enumerate(shared_mask) if m]
    shared_genes = [ref_genes_fmt[i] for i in shared_ref_indices]

    n_shared = len(shared_genes)
    if progress_callback:
        progress_callback(f"Found {n_shared}/{len(ref_genes)} shared genes with reference")

    if n_shared < 20:
        raise ValueError(
            f"Too few shared genes ({n_shared}) between query and reference. "
            f"Need at least 20. Check species compatibility."
        )

    # Subset reference centroids to shared genes
    ref_centroids_shared = ref_centroids[:, shared_ref_indices]

    # Expression matrix for the shared genes
    col_idx = [expr_var_names.index(g) for g in shared_genes]
    X = expr_X[:, col_idx]
    expr_df = pd.DataFrame(
        X.toarray() if hasattr(X, 'toarray') else np.asarray(X),
        index=adata.obs_names,
        columns=shared_genes,
    )

    # Detect cluster column
    if cluster_col not in adata.obs.columns:
        for alt in ['leiden', 'clusters', 'seurat_clusters', 'louvain']:
            if alt in adata.obs.columns:
                cluster_col = alt
                break
        else:
            raise ValueError("No cluster column found in adata.obs")

    clusters = cluster_labels(adata, cluster_col)

    cluster_types = {}
    annotation_details = {}

    for i, cluster in enumerate(clusters):
        if progress_callback and i % 5 == 0:
            progress_callback(f"Annotating cluster {i+1}/{len(clusters)}...")

        mask = adata.obs[cluster_col] == cluster
        n_cells = int(mask.sum())

        # Compute mean expression for this cluster
        cluster_mean = expr_df.loc[mask].mean(axis=0).values

        # Correlate with each reference centroid
        correlations = {}
        for j, ct in enumerate(ref_cell_types):
            ref_vec = ref_centroids_shared[j]
            # Pearson correlation
            r = _pearson_correlation(cluster_mean, ref_vec)
            correlations[ct] = r

        # Sort by correlation
        sorted_corr = sorted(correlations.items(), key=lambda x: x[1], reverse=True)
        best_type, best_corr = sorted_corr[0]
        runner_up, runner_up_corr = sorted_corr[1] if len(sorted_corr) > 1 else ('None', 0)

        if best_corr < min_correlation:
            assigned_type = 'Unknown'
            confidence = 'Low correlation'
        elif best_corr - runner_up_corr < 0.05:
            assigned_type = best_type
            confidence = 'Ambiguous'
        else:
            assigned_type = best_type
            confidence = 'High'

        cluster_types[str(cluster)] = assigned_type
        annotation_details[str(cluster)] = {
            'top_scores': sorted_corr[:5],
            'confidence': confidence,
            'n_cells': n_cells,
            'best_score': best_corr,
            'runner_up': runner_up,
            'runner_up_score': runner_up_corr,
            'margin': best_corr - runner_up_corr,
            'n_shared_genes': n_shared,
        }

    return cluster_types, annotation_details


def _pearson_correlation(x, y):
    """Pearson correlation between two vectors. Returns 0 if degenerate."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    xm = x - x.mean()
    ym = y - y.mean()
    denom = np.sqrt((xm ** 2).sum() * (ym ** 2).sum())
    if denom < 1e-12:
        return 0.0
    return float(np.dot(xm, ym) / denom)


def run_reference_annotation(adata, reference_name_or_path, cluster_col='leiden',
                              min_correlation=0.3, progress_callback=None):
    """
    End-to-end reference-based annotation.

    Loads the reference, runs label transfer, and writes results
    to adata.obs['cell_type'] and adata.obs['cell_type_score'].

    Parameters
    ----------
    adata : anndata.AnnData
        Must have cluster labels in obs.
    reference_name_or_path : str
        Built-in reference name or path to custom reference JSON.
    cluster_col : str
        Cluster column name.
    min_correlation : float
        Minimum Pearson correlation for assignment.
    progress_callback : callable, optional
        Called with status strings.

    Returns
    -------
    tuple of (adata, cluster_types, annotation_details)
    """
    if progress_callback:
        progress_callback(f"Loading reference: {reference_name_or_path}...")

    reference = load_reference(reference_name_or_path)

    if progress_callback:
        n_types = reference.get('n_cell_types', len(reference['cell_types']))
        progress_callback(
            f"Reference loaded: {reference['name']} "
            f"({n_types} cell types, {len(reference['genes'])} genes)"
        )

    cluster_types, annotation_details = annotate_by_reference(
        adata, reference,
        cluster_col=cluster_col,
        min_correlation=min_correlation,
        progress_callback=progress_callback,
    )

    # Map to cells
    if cluster_types and cluster_col in adata.obs.columns:
        adata.obs['cell_type'] = adata.obs[cluster_col].astype(str).map(cluster_types)
        adata.obs['cell_type'] = adata.obs['cell_type'].fillna('Unknown')

        # Store correlation as score
        cluster_scores = {}
        for cid, details in annotation_details.items():
            cluster_scores[cid] = details['best_score']
        adata.obs['cell_type_score'] = (
            adata.obs[cluster_col].astype(str).map(cluster_scores).fillna(0.0)
        )

    if progress_callback:
        n_assigned = sum(1 for v in cluster_types.values() if v != 'Unknown')
        progress_callback(
            f"Reference annotation complete: {n_assigned}/{len(cluster_types)} "
            f"clusters assigned ({len(set(cluster_types.values()) - {'Unknown'})} types)"
        )

    return adata, cluster_types, annotation_details


def build_reference_from_h5ad(h5ad_path, output_path, cell_type_col='cell_type',
                               name=None, species='human', source='',
                               n_top_genes=2000, min_cells_per_type=10):
    """
    Build a lightweight reference JSON from an annotated h5ad atlas.

    Parameters
    ----------
    h5ad_path : str
        Path to annotated h5ad file with cell_type labels.
    output_path : str
        Output path for .json.gz reference file.
    cell_type_col : str
        Column in obs containing cell type labels.
    name : str
        Reference name.
    species : str
        'human' or 'mouse'.
    source : str
        Citation / DOI.
    n_top_genes : int
        Number of highly variable genes to include.
    min_cells_per_type : int
        Minimum cells per type (types with fewer cells are excluded).

    Returns
    -------
    dict with reference metadata
    """
    import scanpy as sc

    adata = sc.read_h5ad(h5ad_path)

    if cell_type_col not in adata.obs.columns:
        raise ValueError(f"Column '{cell_type_col}' not found in obs. "
                         f"Available: {list(adata.obs.columns)}")

    # Filter small types
    type_counts = adata.obs[cell_type_col].value_counts()
    valid_types = type_counts[type_counts >= min_cells_per_type].index.tolist()
    adata = adata[adata.obs[cell_type_col].isin(valid_types)].copy()

    # Normalize if needed
    if hasattr(adata.X, 'data') and len(adata.X.data) > 0:
        max_val = adata.X.data.max()
    elif hasattr(adata.X, 'max'):
        max_val = adata.X.max()
    else:
        max_val = 100
    if max_val > 20:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # Find HVGs
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='seurat_v3'
                                 if max_val > 20 else 'seurat')
    hvg_mask = adata.var['highly_variable']
    genes = list(adata.var_names[hvg_mask])

    # Compute centroids (mean expression per cell type)
    cell_types = sorted(valid_types)
    centroids = []
    cell_counts = []

    for ct in cell_types:
        mask = adata.obs[cell_type_col] == ct
        n = int(mask.sum())
        cell_counts.append(n)

        X_ct = adata[mask][:, hvg_mask].X
        if hasattr(X_ct, 'toarray'):
            X_ct = X_ct.toarray()
        mean_expr = X_ct.mean(axis=0).tolist()
        if hasattr(mean_expr, 'tolist'):
            mean_expr = mean_expr.tolist()
        # Handle numpy matrix (returns nested list)
        if isinstance(mean_expr, list) and len(mean_expr) == 1:
            mean_expr = mean_expr[0]
        centroids.append(mean_expr)

    ref = {
        'name': name or Path(h5ad_path).stem,
        'species': species,
        'source': source,
        'n_cell_types': len(cell_types),
        'n_genes': len(genes),
        'cell_types': cell_types,
        'genes': genes,
        'centroids': centroids,
        'cell_counts': cell_counts,
    }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, 'wt') as f:
        json.dump(ref, f)

    return {k: v for k, v in ref.items() if k != 'centroids'}
