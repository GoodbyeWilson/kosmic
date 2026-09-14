# Transcriptome Size Test
# Tests whether total per-cell expression differs between disease and control.
# A significant difference indicates that global normalization methods
# (CPM, DESeq2 size factors) may mask real biological signal.

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


def test_transcriptome_size_adata(adata, sample_col, condition_col,
                                  disease_label, control_label):
    """Test for transcriptome size difference from AnnData directly.

    Computes three per-cell metrics for each sample, then compares
    disease vs control samples with Welch's t-test:
    1. Mean total UMI per cell (transcriptional output)
    2. Mean genes detected per cell (transcriptional diversity)
    3. Mean UMI per gene per cell (expression depth)

    Parameters
    ----------
    adata : AnnData
        Raw or processed data.  Uses raw counts if available.
    sample_col, condition_col : str
        Column names in adata.obs.
    disease_label, control_label : str

    Returns
    -------
    result : dict with keys: n_disease, n_control, metrics (dict of
        metric_name -> {mean_disease, mean_control, ratio, t_stat, pval}),
        per_sample (list of dicts).
    """
    # Use raw counts if available
    if 'counts' in adata.layers:
        X = adata.layers['counts']
    elif adata.raw is not None:
        X = adata.raw.X
    else:
        X = adata.X

    samples = adata.obs[sample_col].unique()
    per_sample = []

    for sample_id in samples:
        mask = adata.obs[sample_col] == sample_id
        cells = X[mask.values] if hasattr(mask, 'values') else X[mask]
        if hasattr(cells, 'toarray'):
            cells = cells.toarray()
        cells = np.asarray(cells, dtype=float)

        total_per_cell = cells.sum(axis=1)
        genes_per_cell = (cells > 0).sum(axis=1)
        umi_per_gene = np.where(genes_per_cell > 0,
                                total_per_cell / genes_per_cell, 0)

        condition = adata.obs.loc[mask, condition_col].iloc[0]

        per_sample.append({
            'sample': str(sample_id),
            'condition': str(condition),
            'n_cells': int(mask.sum()),
            'mean_total_umi': float(total_per_cell.mean()),
            'mean_genes_detected': float(genes_per_cell.mean()),
            'mean_umi_per_gene': float(umi_per_gene.mean()),
        })

    sample_df = pd.DataFrame(per_sample)

    disease_mask = sample_df['condition'] == disease_label
    control_mask = sample_df['condition'] == control_label

    if disease_mask.sum() < 2 or control_mask.sum() < 2:
        return None

    metric_keys = [
        ('mean_total_umi', 'Transcripts per cell (UMI counts)'),
        ('mean_genes_detected', 'Genes detected per cell'),
        ('mean_umi_per_gene', 'UMI per gene per cell'),
    ]

    metrics = {}
    for col, label in metric_keys:
        d_vals = sample_df.loc[disease_mask, col].values
        c_vals = sample_df.loc[control_mask, col].values
        ratio = float(d_vals.mean() / max(c_vals.mean(), 1e-10))
        try:
            t_stat, pval = ttest_ind(d_vals, c_vals, equal_var=False)
        except ValueError:
            t_stat, pval = 0.0, 1.0  # n<2 or both groups zero-variance
        metrics[col] = {
            'label': label,
            'mean_disease': float(d_vals.mean()),
            'mean_control': float(c_vals.mean()),
            'ratio': ratio,
            't_stat': float(t_stat),
            'pval': float(pval),
        }

    return {
        'n_disease': int(disease_mask.sum()),
        'n_control': int(control_mask.sum()),
        'metrics': metrics,
        'per_sample': per_sample,
    }

