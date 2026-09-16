# Cell Type Annotation
# Marker-based cell type scoring and cluster annotation.

import re
from kosmic.scrna.inspect.roles import cluster_labels
import numpy as np

from kosmic.scrna.inspect.detection import detect_species, format_gene_for_species


# CellTypist broad label collapsing

# Known prefix → broad name mappings for CellTypist models.
# Longest-prefix-first matching ensures e.g. "CD4+T" matches before "CD4".
_CELLTYPIST_BROAD_MAP = {
    'aCM': 'Atrial Cardiomyocyte',
    'vCM': 'Ventricular Cardiomyocyte',
    'EC7_endocardial': 'Endocardial Cell',
    'EC8_ln': 'Lymphatic Endothelial Cell',
    'EC': 'Endothelial Cell',
    'FB': 'Fibroblast',
    'PC': 'Pericyte',
    'NC': 'Neural Crest',
    'SMC': 'Smooth Muscle Cell',
    'Adip': 'Adipocyte',
    'CD4+T': 'CD4+ T Cell',
    'CD8+T': 'CD8+ T Cell',
    'CD14+Mo': 'CD14+ Monocyte',
    'CD16+Mo': 'CD16+ Monocyte',
    'NK': 'NK Cell',
    'AVN': 'AV Node Cell',
    'SAN': 'SA Node Cell',
    'LYVE1': 'LYVE1+ Macrophage',
    'B_': 'B Cell',
    'B': 'B Cell',
    'Purkinje': 'Purkinje Cell',
    'Neut': 'Neutrophil',
    'Mast': 'Mast Cell',
    'Meso': 'Mesothelial Cell',
    'MoMP': 'Monocyte-Macrophage',
    'DC': 'Dendritic Cell',
    'ILC': 'Innate Lymphoid Cell',
    'gdT': 'Gamma-delta T Cell',
    'MAIT': 'MAIT Cell',
    'T/NK': 'T/NK Cell',
}

_SORTED_PREFIXES = sorted(_CELLTYPIST_BROAD_MAP.keys(), key=len, reverse=True)


def collapse_celltypist_label(label):
    """Collapse a CellTypist fine-grained label to a broad category.

    Uses known prefix mappings first, then falls back to stripping
    trailing digits and underscore-suffixes.

    Non-string input is returned unchanged. Cells left out of the
    embedding carry no label, and ``Categorical.map`` additionally calls
    the mapper with ``np.nan`` to work out what NA should become -- so
    this is reached with a float even when every real label is a string.
    """
    if not isinstance(label, str):
        return label

    for prefix in _SORTED_PREFIXES:
        if label.startswith(prefix):
            return _CELLTYPIST_BROAD_MAP[prefix]

    # Fallback: strip trailing digit+suffix pattern (e.g. "Type3_sub" → "Type")
    broad = re.sub(r'\d+(_\w+)?$', '', label)
    if broad and broad != label:
        return broad

    return label


# Default marker genes for common cell types
DEFAULT_MARKERS = {
    'Endothelial': ['PECAM1', 'CDH5', 'VWF', 'KDR', 'ENG'],
    'SMC': ['ACTA2', 'MYH11', 'TAGLN', 'CNN1', 'SMTN'],
    'Pericyte': ['PDGFRB', 'RGS5', 'ABCC9', 'KCNJ8', 'DES'],
    'Fibroblast': ['DCN', 'LUM', 'COL1A1', 'COL1A2', 'PDGFRA'],
    'Macrophage': ['CD68', 'CD163', 'CSF1R', 'MARCO', 'MSR1'],
    'T_cell': ['CD3D', 'CD3E', 'CD3G', 'CD2', 'IL7R'],
    'B_cell': ['CD19', 'MS4A1', 'CD79A', 'CD79B', 'PAX5'],
    'Neuron': ['RBFOX3', 'MAP2', 'TUBB3', 'NEFL', 'SYP'],
    'Astrocyte': ['GFAP', 'AQP4', 'S100B', 'ALDH1L1', 'SLC1A3'],
    'Oligodendrocyte': ['MBP', 'PLP1', 'MOG', 'MAG', 'OLIG2'],
    'Microglia': ['CX3CR1', 'P2RY12', 'TMEM119', 'AIF1', 'ITGAM'],
    'Cardiomyocyte': ['TNNT2', 'MYH7', 'MYH6', 'ACTC1', 'RYR2'],
    'Epithelial': ['EPCAM', 'KRT18', 'KRT19', 'CDH1', 'MUC1'],
}


def score_marker_genes(adata, markers=None, skip_gene_conversion=False,
                       progress_callback=None, max_rank=1500, min_markers=3):
    """
    Score cells for marker gene sets using UCell (rank-based) scoring.

    For each cell, genes are ranked by expression (descending).  The score
    for a marker set is 1 − (U / U_max), where U is the Mann-Whitney U
    statistic comparing marker ranks to the maximum possible ranks.
    Scores range from 0 (markers not enriched) to 1 (all markers at the
    very top of the ranking).

    This method is agnostic to the number of markers — a 3-marker set and
    a 30-marker set produce directly comparable scores.

    Reference: Andreatta & Carmona, Comput Struct Biotechnol J (2021).

    Parameters
    ----------
    adata : anndata.AnnData
        Input data with expression values.
    markers : dict, optional
        {cell_type: [gene1, gene2, ...]}. Defaults to DEFAULT_MARKERS.
    skip_gene_conversion : bool
        If True, use marker names as-is without species conversion.
    progress_callback : callable, optional
        Called with status strings.
    max_rank : int
        Only consider the top N genes per cell for ranking.  Genes ranked
        below this are all tied at max_rank+1.  Speeds up computation and
        reduces noise from lowly-expressed genes.  Default 1500.
    min_markers : int
        Minimum markers found in data to score a cell type.  Types with
        fewer are skipped.  Default 3.

    Returns
    -------
    tuple of (adata, scored_types)
        scored_types: dict of {cell_type: n_markers_found}
    """
    if markers is None:
        markers = DEFAULT_MARKERS

    species = detect_species(list(adata.var_names))

    # Clear stale score columns from previous runs
    old_score_cols = [c for c in adata.obs.columns if c.endswith('_score')]
    if old_score_cols:
        adata.obs.drop(columns=old_score_cols, inplace=True)

    # Log-normalised expression on the full gene list
    from kosmic.scrna.counts import log_normalised
    X, gene_names = log_normalised(adata)

    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    n_cells, n_genes = X.shape

    if progress_callback:
        progress_callback("Resolving marker genes...")

    # Resolve all marker gene sets
    resolved_markers = {}
    scored_types = {}
    for cell_type, marker_list in markers.items():
        if skip_gene_conversion:
            markers_fmt = marker_list
        else:
            markers_fmt = [format_gene_for_species(m, species) for m in marker_list]

        found_idx = [gene_to_idx[m] for m in markers_fmt if m in gene_to_idx]
        scored_types[cell_type] = len(found_idx)
        if len(found_idx) >= min_markers:
            resolved_markers[cell_type] = np.array(found_idx)

    # Report skipped types
    skipped = [f"{ct} ({n})" for ct, n in scored_types.items()
               if 0 < n < min_markers]
    if skipped and progress_callback:
        progress_callback(f"Skipped {len(skipped)} types with <{min_markers} markers: "
                          f"{', '.join(skipped[:5])}")

    if not resolved_markers:
        return adata, scored_types

    n_resolved = len(resolved_markers)
    if progress_callback:
        progress_callback(f"Scoring {n_resolved} cell types (UCell)...")

    # Process in chunks to manage memory
    chunk_size = 5000
    all_scores = {ct: np.empty(n_cells, dtype=np.float64)
                  for ct in resolved_markers}

    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        if progress_callback:
            progress_callback(
                f"Ranking genes ({start:,}–{end:,} of {n_cells:,} cells)...",
                start / n_cells,
            )

        # Get dense chunk
        X_chunk = X[start:end]
        if hasattr(X_chunk, 'toarray'):
            X_chunk = X_chunk.toarray()
        X_chunk = np.asarray(X_chunk, dtype=np.float64)

        # Rank genes per cell (descending — rank 1 = highest expression)
        # Use argpartition for speed: only sort the top max_rank genes
        chunk_cells = X_chunk.shape[0]
        ranks = np.full((chunk_cells, n_genes), max_rank + 1, dtype=np.float64)

        for i in range(chunk_cells):
            row = X_chunk[i]
            if max_rank < n_genes:
                # Partial sort: find top max_rank genes
                top_idx = np.argpartition(row, -max_rank)[-max_rank:]
                # Full sort only within top genes
                sorted_within = top_idx[np.argsort(row[top_idx])[::-1]]
                ranks[i, sorted_within] = np.arange(1, max_rank + 1, dtype=np.float64)
            else:
                order = np.argsort(row)[::-1]
                ranks[i, order] = np.arange(1, n_genes + 1, dtype=np.float64)

        # Score each cell type using UCell formula
        for cell_type, marker_idx in resolved_markers.items():
            m = len(marker_idx)
            marker_ranks = ranks[:, marker_idx]  # (chunk_cells, m)
            # U statistic: sum of ranks minus the minimum possible sum
            rank_sum = marker_ranks.sum(axis=1)
            u = rank_sum - m * (m + 1) / 2
            u_max = m * max_rank  # worst case: all markers at max_rank
            # UCell score: 1 - U/U_max, clipped to [0, 1]
            scores = np.clip(1.0 - u / u_max, 0.0, 1.0)
            all_scores[cell_type][start:end] = scores

    # Store scores in adata
    for cell_type, scores in all_scores.items():
        adata.obs[f'{cell_type}_score'] = scores

    if progress_callback:
        progress_callback(f"Scored {n_resolved} cell types (UCell)", 1.0)

    return adata, scored_types


def assign_cell_types_per_cluster(adata, cluster_col='leiden',
                                   score_threshold=0.1, confidence_margin=0.05,
                                   force_assignment=False):
    """
    Assign cell types to clusters based on mean marker scores.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have score columns ({cell_type}_score) in obs.
    cluster_col : str
        Column in adata.obs with cluster labels.
    score_threshold : float
        Min mean score to assign a cell type (else 'Unknown').
    confidence_margin : float
        Min difference between top 2 scores for confident assignment.
    force_assignment : bool
        If True, assign even when below threshold.

    Returns
    -------
    tuple of (cluster_types, annotation_details)
        cluster_types: {cluster_id: cell_type_name}
        annotation_details: {cluster_id: {top_scores, confidence, n_cells, ...}}
    """
    # Find score columns (exclude cell_type_score which is an aggregate, not a type)
    score_cols = [c for c in adata.obs.columns
                  if c.endswith('_score') and c != 'cell_type_score']
    if not score_cols:
        return {}, {}

    # Detect cluster column
    if cluster_col not in adata.obs.columns:
        for alt in ['leiden', 'clusters', 'seurat_clusters', 'louvain']:
            if alt in adata.obs.columns:
                cluster_col = alt
                break
        else:
            return {}, {}

    cluster_types = {}
    annotation_details = {}

    for cluster in cluster_labels(adata, cluster_col):
        mask = adata.obs[cluster_col] == cluster
        cluster_cells = adata.obs.loc[mask]
        n_cells = int(mask.sum())

        # Compute mean scores for each cell type
        mean_scores = {}
        for col in score_cols:
            ct_name = col.replace('_score', '')
            mean_scores[ct_name] = cluster_cells[col].mean()

        if not mean_scores:
            cluster_types[str(cluster)] = 'Unknown'
            continue

        # Sort by score
        sorted_scores = sorted(mean_scores.items(), key=lambda x: x[1], reverse=True)
        best_type, best_score = sorted_scores[0]
        runner_up, runner_up_score = sorted_scores[1] if len(sorted_scores) > 1 else ('None', 0)
        margin = best_score - runner_up_score

        # Assignment logic
        if best_score < score_threshold and not force_assignment:
            assigned_type = 'Unknown'
            confidence = 'Below threshold'
        elif best_score < score_threshold and force_assignment:
            assigned_type = best_type
            confidence = 'Below threshold (forced)'
        elif margin < confidence_margin:
            assigned_type = best_type
            confidence = 'Ambiguous'
        else:
            assigned_type = best_type
            confidence = 'High'

        cluster_types[str(cluster)] = assigned_type
        annotation_details[str(cluster)] = {
            'top_scores': sorted_scores[:3],
            'confidence': confidence,
            'n_cells': n_cells,
            'best_score': best_score,
            'runner_up': runner_up,
            'runner_up_score': runner_up_score,
            'margin': margin,
        }

    return cluster_types, annotation_details


def run_ora_per_cluster(adata, markers, cluster_col='leiden',
                        skip_gene_conversion=False, n_top_genes=200):
    """
    Run ORA per cluster using Fisher's exact test.

    For each cluster, identifies up-regulated genes (higher mean expression
    in the cluster than global mean), then tests for enrichment of each
    marker gene set using scipy's Fisher's exact test.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have cluster labels in obs.
    markers : dict
        {cell_type: [gene1, gene2, ...]}
    cluster_col : str
        Cluster column in adata.obs.
    skip_gene_conversion : bool
        If True, use marker names as-is without species conversion.
    n_top_genes : int
        Number of top up-regulated genes per cluster to use as the
        "expressed" set for Fisher's test.

    Returns
    -------
    dict
        {cluster_id: {cell_type: ora_pvalue, ...}}
    """
    from scipy.stats import fisher_exact

    species = detect_species(list(adata.var_names))

    # Detect cluster column
    if cluster_col not in adata.obs.columns:
        for alt in ['leiden', 'clusters', 'seurat_clusters', 'louvain']:
            if alt in adata.obs.columns:
                cluster_col = alt
                break
        else:
            return {}

    # Build resolved marker sets
    resolved_markers = {}
    for ct, gene_list in markers.items():
        if skip_gene_conversion:
            genes = set(gene_list)
        else:
            genes = set(format_gene_for_species(g, species) for g in gene_list)
        resolved_markers[ct] = genes

    # Log-normalised expression on the full gene list
    from kosmic.scrna.counts import log_normalised
    X, names = log_normalised(adata)
    var_names = np.array(names)

    n_genes = len(var_names)

    # Compute global mean expression per gene
    if hasattr(X, 'toarray'):
        global_mean = np.asarray(X.mean(axis=0)).flatten()
    else:
        global_mean = np.mean(X, axis=0).flatten()

    clusters = cluster_labels(adata, cluster_col)

    ora_results = {}
    for cluster in clusters:
        cluster_id = str(cluster)
        mask = (adata.obs[cluster_col] == cluster).values
        X_cluster = X[mask]

        # Compute cluster mean expression
        if hasattr(X_cluster, 'toarray'):
            cluster_mean = np.asarray(X_cluster.mean(axis=0)).flatten()
        else:
            cluster_mean = np.mean(X_cluster, axis=0).flatten()

        # Fold change: cluster mean / global mean (with pseudocount)
        fc = (cluster_mean + 1e-9) / (global_mean + 1e-9)

        # Select top up-regulated genes for this cluster
        top_idx = np.argsort(fc)[-n_top_genes:]
        top_genes = set(var_names[top_idx])

        # Fisher's exact test for each cell type marker set
        cluster_pvals = {}
        for ct, marker_genes in resolved_markers.items():
            # Genes in marker set that exist in the data
            markers_in_data = marker_genes & set(var_names)
            if not markers_in_data:
                continue

            # 2x2 contingency:
            #                  In top genes | Not in top genes
            # In marker set    |    a       |      b
            # Not in marker    |    c       |      d
            a = len(markers_in_data & top_genes)
            b = len(markers_in_data - top_genes)
            c = len(top_genes - markers_in_data)
            d = n_genes - a - b - c

            _, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
            cluster_pvals[ct] = pval

        ora_results[cluster_id] = cluster_pvals

    return ora_results


def annotate_clusters(adata, markers=None, cluster_col='leiden',
                      score_threshold=0.1, confidence_margin=0.05,
                      skip_gene_conversion=False, force_assignment=False):
    """
    End-to-end cluster annotation: score markers then assign cell types.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have cluster labels in obs.
    markers : dict, optional
        Marker gene dict. Defaults to DEFAULT_MARKERS.
    cluster_col : str
        Cluster column name.
    score_threshold, confidence_margin, force_assignment :
        Passed to assign_cell_types_per_cluster.
    skip_gene_conversion : bool
        Passed to score_marker_genes.

    Returns
    -------
    tuple of (adata, cluster_types, annotation_details)
        adata has 'cell_type' and 'cell_type_score' columns added to obs.
    """
    # Score
    adata, scored_types = score_marker_genes(
        adata, markers=markers, skip_gene_conversion=skip_gene_conversion
    )

    # Assign
    cluster_types, annotation_details = assign_cell_types_per_cluster(
        adata, cluster_col=cluster_col,
        score_threshold=score_threshold,
        confidence_margin=confidence_margin,
        force_assignment=force_assignment,
    )

    # Map cluster assignments to individual cells
    if cluster_types and cluster_col in adata.obs.columns:
        adata.obs['cell_type'] = adata.obs[cluster_col].astype(str).map(cluster_types)
        adata.obs['cell_type'] = adata.obs['cell_type'].fillna('Unknown')

        # Also store the best score per cell
        score_cols = [c for c in adata.obs.columns
                      if c.endswith('_score') and c != 'cell_type_score']
        if score_cols:
            adata.obs['cell_type_score'] = adata.obs[score_cols].max(axis=1)

    return adata, cluster_types, annotation_details


def run_celltypist(adata, model_name, majority_voting=True, progress_callback=None):
    """
    Run CellTypist ML-based cell type annotation.

    Handles model download, expression normalization from raw,
    majority voting with existing leiden clusters, and standard column names.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have been through clustering (leiden in obs) for majority voting.
    model_name : str
        CellTypist model name (e.g. 'Immune_All_Low.pkl').
    majority_voting : bool
        Whether to apply majority voting within leiden clusters.
    progress_callback : callable, optional
        Called with status strings.

    Returns
    -------
    anndata.AnnData
        adata with 'cell_type', 'cell_type_auto', 'cell_type_score',
        'predicted_labels', and optionally 'majority_voting' columns in obs.

    Raises
    ------
    ImportError
        If celltypist is not installed.
    """
    import celltypist
    from celltypist import models
    import scanpy as sc

    if progress_callback:
        progress_callback(f"Loading CellTypist model: {model_name}...")

    # Download model if not cached
    try:
        model = models.Model.load(model=model_name)
    except Exception:
        if progress_callback:
            progress_callback(f"Downloading model {model_name}...")
        models.download_models(model=model_name)
        model = models.Model.load(model=model_name)

    if progress_callback:
        progress_callback("Preparing expression data for CellTypist...")

    # CellTypist requires log1p normalised expression (target_sum=1e4).
    # Rebuild it from the counts rather than trusting X, which may have
    # been scaled by the cluster step.
    from kosmic.scrna.counts import counts_adata, has_counts_layer
    if has_counts_layer(adata):
        adata_ct = counts_adata(adata)
        sc.pp.normalize_total(adata_ct, target_sum=1e4)
        sc.pp.log1p(adata_ct)
    else:
        adata_ct = adata.copy()
        if hasattr(adata_ct.X, 'data') and len(adata_ct.X.data) > 0:
            max_val = adata_ct.X.data.max()
        elif hasattr(adata_ct.X, 'max'):
            max_val = adata_ct.X.max()
        else:
            max_val = 100
        if max_val > 20:
            sc.pp.normalize_total(adata_ct, target_sum=1e4)
            sc.pp.log1p(adata_ct)

    # Carry over obs (leiden clusters needed for majority_voting)
    adata_ct.obs = adata.obs.copy()

    # For very large datasets, CellTypist's internal scaling densifies the
    # sparse matrix causing OOM. Process in chunks that fit in memory,
    # then do majority voting on the full result.
    CHUNK_SIZE = 100_000  # cells per chunk (~3 GB for 3802 features × float64)
    chunked = adata_ct.n_obs > CHUNK_SIZE

    if chunked:
        if progress_callback:
            progress_callback(
                f"Annotating {adata_ct.n_obs:,} cells in chunks of {CHUNK_SIZE:,}..."
            )

        # Annotate in chunks (per-cell predictions, no majority voting yet)
        n_chunks = int(np.ceil(adata_ct.n_obs / CHUNK_SIZE))
        all_labels = []
        all_scores = []

        for chunk_i in range(n_chunks):
            start = chunk_i * CHUNK_SIZE
            end = min((chunk_i + 1) * CHUNK_SIZE, adata_ct.n_obs)
            if progress_callback:
                progress_callback(
                    f"Annotating chunk {chunk_i + 1}/{n_chunks} "
                    f"({start:,}-{end:,} of {adata_ct.n_obs:,} cells)..."
                )

            chunk = adata_ct[start:end].copy()
            chunk_pred = celltypist.annotate(chunk, model=model, majority_voting=False)
            chunk_result = chunk_pred.to_adata()

            all_labels.extend(chunk_result.obs['predicted_labels'].values)
            if 'conf_score' in chunk_result.obs.columns:
                all_scores.extend(chunk_result.obs['conf_score'].values)

            del chunk, chunk_pred, chunk_result

        # Build full predictions
        adata_ct.obs['predicted_labels'] = all_labels
        if all_scores:
            adata_ct.obs['conf_score'] = all_scores

        # Now do majority voting on the full set (just label assignment, no matrix ops)
        use_leiden = majority_voting and 'leiden' in adata_ct.obs.columns
        if use_leiden:
            if progress_callback:
                progress_callback("Running majority voting across clusters...")
            cluster_col = adata_ct.obs['leiden']
            labels = adata_ct.obs['predicted_labels']
            majority = {}
            # Only clusters with cells: a cell left out of the embedding has
            # no cluster, and an empty group has no majority to take.
            for cluster in cluster_labels(adata_ct, 'leiden'):
                votes = labels[cluster_col == cluster].value_counts()
                if len(votes):
                    majority[cluster] = votes.index[0]
            adata_ct.obs['majority_voting'] = cluster_col.map(majority)

        # Transfer to adata
        ct_cols = ['predicted_labels', 'majority_voting', 'conf_score']
        for col in ct_cols:
            if col in adata_ct.obs.columns:
                adata.obs[col] = adata_ct.obs[col].values
        del adata_ct
    else:
        if progress_callback:
            progress_callback(f"Annotating {adata_ct.n_obs:,} cells with {model_name}...")

        use_leiden = majority_voting and 'leiden' in adata_ct.obs.columns
        predictions = celltypist.annotate(
            adata_ct,
            model=model,
            majority_voting=majority_voting,
            over_clustering='leiden' if use_leiden else None,
        )
        del adata_ct

        adata_result = predictions.to_adata()
        del predictions

        # Transfer annotation columns back to original adata
        ct_cols = ['predicted_labels', 'over_clustering', 'majority_voting', 'conf_score']
        for col in ct_cols:
            if col in adata_result.obs.columns:
                adata.obs[col] = adata_result.obs[col].values
        del adata_result

    # Map to standard column names
    if majority_voting and 'majority_voting' in adata.obs.columns:
        adata.obs['cell_type'] = adata.obs['majority_voting']
    else:
        adata.obs['cell_type'] = adata.obs['predicted_labels']

    adata.obs['cell_type_auto'] = adata.obs['cell_type']

    if 'conf_score' in adata.obs.columns:
        adata.obs['cell_type_score'] = adata.obs['conf_score']
    else:
        adata.obs['cell_type_score'] = 0.0

    return adata
