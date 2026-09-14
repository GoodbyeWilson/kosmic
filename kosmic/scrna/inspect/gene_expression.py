# Gene Expression Analysis
# Single-gene expression stats and multi-gene dotplot statistics for scRNA-seq data.

import numpy as np
import pandas as pd


def analyze_gene_expression(adata, gene_name, sample_col, condition_col, use_raw=True):
    """Analyze expression of a single gene across conditions.

    Parameters
    ----------
    adata : anndata.AnnData
        Input data.
    gene_name : str
        Gene to analyze (case-insensitive lookup).
    sample_col : str
        Sample grouping column in obs.
    condition_col : str
        Condition column in obs.
    use_raw : bool
        If True and adata.raw is available, use raw expression values.

    Returns
    -------
    tuple of (per_cell_df, summary_df, resolved_gene)
        per_cell_df: DataFrame with cell_barcode, sample, condition, gene, expression.
        summary_df: DataFrame with per-condition summary stats.
        resolved_gene: str, the actual gene name found in the dataset.

    Raises
    ------
    ValueError
        If gene is not found in dataset.
    """
    # Determine data source
    if use_raw and adata.raw is not None:
        var_names = adata.raw.var_names
    else:
        var_names = adata.var_names

    # Case-insensitive gene lookup
    gene = _resolve_gene_name(gene_name, var_names)
    if gene is None:
        raise ValueError(f"Gene '{gene_name}' not found in dataset")

    # Extract expression
    if use_raw and adata.raw is not None:
        gene_data = adata.raw[:, gene].X
        if hasattr(gene_data, 'toarray'):
            gene_expr = gene_data.toarray().flatten()
        else:
            gene_expr = np.array(gene_data).flatten()
    else:
        gene_idx = adata.var_names.get_loc(gene)
        if hasattr(adata.X, 'toarray'):
            gene_expr = adata.X[:, gene_idx].toarray().flatten()
        else:
            gene_expr = adata.X[:, gene_idx].flatten()

    # Build per-cell DataFrame
    per_cell_df = pd.DataFrame({
        'cell_barcode': adata.obs_names,
        'sample': adata.obs[sample_col].values,
        'condition': adata.obs[condition_col].values,
        'gene': gene,
        'expression': gene_expr,
    })

    # Summary stats per condition
    conditions = adata.obs[condition_col].unique()
    summary_rows = []
    for condition in conditions:
        cond_data = per_cell_df[per_cell_df['condition'] == condition]
        n_cells = len(cond_data)
        mean_expr = cond_data['expression'].mean()
        std_expr = cond_data['expression'].std(ddof=1)
        sem = std_expr / np.sqrt(n_cells) if n_cells > 0 else 0
        pct_expressing = (cond_data['expression'] > 0).sum() / n_cells * 100 if n_cells > 0 else 0

        summary_rows.append({
            'condition': condition,
            'gene': gene,
            'n_cells': n_cells,
            'n_samples': cond_data['sample'].nunique(),
            'mean_expression': mean_expr,
            'sd': std_expr,
            'sem': sem,
            'pct_expressing': pct_expressing,
        })

    summary_df = pd.DataFrame(summary_rows)

    # Calculate log2FC if exactly 2 conditions
    if len(summary_rows) == 2:
        control_patterns = [
            'control', 'ctrl', 'healthy', 'normal', 'wt', 'wildtype',
            'wild-type', 'sham', 'baseline', 'nf', 'donor',
        ]
        cond1, cond2 = summary_rows[0], summary_rows[1]
        cond1_ctrl = any(p in str(cond1['condition']).lower() for p in control_patterns)
        cond2_ctrl = any(p in str(cond2['condition']).lower() for p in control_patterns)

        if cond2_ctrl and not cond1_ctrl:
            disease, control = cond1, cond2
        elif cond1_ctrl and not cond2_ctrl:
            disease, control = cond2, cond1
        else:
            disease, control = (cond1, cond2) if str(cond1['condition']) < str(cond2['condition']) else (cond2, cond1)

        pseudo = 0.001
        log2fc = np.log2((disease['mean_expression'] + pseudo) / (control['mean_expression'] + pseudo))
        summary_df['log2fc_vs_control'] = summary_df['condition'].apply(
            lambda x: log2fc if x == disease['condition'] else -log2fc
        )

    return per_cell_df, summary_df, gene


def compute_dotplot_stats(adata, gene_list, condition_col, conditions,
                          use_raw=True, expressing_only=False):
    """Compute per-gene, per-condition statistics for a dotplot.

    Parameters
    ----------
    adata : anndata.AnnData
        Input data.
    gene_list : list of str
        Genes to include (case-insensitive lookup).
    condition_col : str
        Condition column in obs.
    conditions : list of str
        Conditions to include.
    use_raw : bool
        If True and adata.raw available, use raw expression.
    expressing_only : bool
        If True, compute mean only over expressing cells (>0).

    Returns
    -------
    tuple of (stats_df, resolved_genes)
        stats_df: DataFrame with gene, condition, pct_expressing, mean_expression,
                  mean_all_cells, n_cells columns.
        resolved_genes: list of resolved gene names.
    """
    if use_raw and adata.raw is not None:
        var_names = adata.raw.var_names
    else:
        var_names = adata.var_names

    # Resolve gene names
    resolved_genes = []
    for gene in gene_list:
        gene = gene.strip()
        if not gene:
            continue
        resolved = _resolve_gene_name(gene, var_names)
        if resolved is not None:
            resolved_genes.append(resolved)

    if not resolved_genes:
        return pd.DataFrame(), []

    # Compute stats
    rows = []
    for gene in resolved_genes:
        if use_raw and adata.raw is not None:
            gene_data = adata.raw[:, gene].X
            if hasattr(gene_data, 'toarray'):
                expr = gene_data.toarray().flatten()
            else:
                expr = np.array(gene_data).flatten()
        else:
            gene_idx = var_names.get_loc(gene)
            if hasattr(adata.X, 'toarray'):
                expr = adata.X[:, gene_idx].toarray().flatten()
            else:
                expr = adata.X[:, gene_idx].flatten()

        cond_values = adata.obs[condition_col].values

        for cond in conditions:
            mask = cond_values == cond
            cond_expr = expr[mask]
            if len(cond_expr) == 0:
                continue

            n_cells = len(cond_expr)
            n_expressing = np.sum(cond_expr > 0)
            pct_expressing = (n_expressing / n_cells) * 100
            mean_all = np.mean(cond_expr)

            if expressing_only and n_expressing > 0:
                mean_for_color = np.mean(cond_expr[cond_expr > 0])
            else:
                mean_for_color = mean_all

            rows.append({
                'gene': gene,
                'condition': cond,
                'pct_expressing': pct_expressing,
                'mean_expression': mean_for_color,
                'mean_all_cells': mean_all,
                'n_cells': n_cells,
            })

    return pd.DataFrame(rows), resolved_genes


def _resolve_gene_name(gene_name, var_names):
    """Case-insensitive gene name resolution.

    Returns the matching gene name from var_names, or None if not found.
    """
    if gene_name in var_names:
        return gene_name

    for variant in [gene_name.upper(), gene_name.lower(), gene_name.capitalize()]:
        if variant in var_names:
            return variant

    matches = [g for g in var_names if g.lower() == gene_name.lower()]
    return matches[0] if matches else None
