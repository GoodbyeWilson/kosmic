# Pathway scoring: one score per donor per pathway.
#
# Each donor's counts are summed, normalised with DESeq2's median-of-ratios
# size factors and log2(x + 1) transformed (DESeq2's normTransform), and the
# pathway's genes are averaged per donor (de_analysis.compute_pathway_scores).
# Disease and control donors are then compared with Welch's t-test
# (de_analysis.run_pseudobulk_de). The group difference of the score is the
# average of the genes' log2 differences.
#
# Earlier versions offered seven other scores (VST, z-scored VST, CPM, and
# four per-cell scores tested after Shapiro / Levene pre-tests). They are
# removed; their labels stay so the methods can describe runs made with them.

from pathlib import Path

import numpy as np
import pandas as pd


PER_DONOR_SCORES_FILE = 'pathway_scores_per_donor.csv'


def save_per_donor_scores(scores, path):
    """Persist the per-donor pathway score matrix ('names', 'matrix', 'sample_df')
    so figures can reload the exact values a scoring run produced instead of
    recomputing them. One row per donor: sample, condition, then a column per
    pathway. Returns the written path."""
    names = list(scores['names'])
    matrix = np.asarray(scores['matrix'], dtype=float)
    sdf = scores['sample_df']
    out = pd.DataFrame()
    for col in ('sample', 'condition'):
        if col in getattr(sdf, 'columns', []):
            out[col] = np.asarray(sdf[col])
    for j, pw in enumerate(names):
        out[pw] = matrix[:, j]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    return path


def load_per_donor_scores(path):
    """Reload per-donor pathway scores written by 'save_per_donor_scores' into the
    {'names', 'matrix', 'sample_df'} shape the plot expects, or None if the file
    is missing or lacks a condition column."""
    path = Path(path)
    if not path.exists():
        return None
    df = pd.read_csv(path)
    meta_cols = [c for c in ('sample', 'condition') if c in df.columns]
    names = [c for c in df.columns if c not in meta_cols]
    if not names or 'condition' not in df.columns:
        return None
    return {
        'names': names,
        'matrix': df[names].to_numpy(dtype=float),
        'sample_df': df[meta_cols].copy(),
    }


#: The scoring method, as recorded with each run.
SCORING_METHOD = "deseq2"

#: Plain description of the score and its test.
SCORING_DESCRIPTION = (
    "One score per donor for each pathway: the average of the pathway's "
    "genes after DESeq2 normalisation (median-of-ratios size factors, "
    "log2(x + 1)). Disease and control donors are compared with Welch's "
    "t-test. The effect is the difference in mean score, in log2 units.")

#: Labels of every method a run may have been made with, including those
#: earlier versions offered, so the methods can describe older runs.
SCORING_METHODS = {
    "deseq2": "Mean DESeq2-normalised log2 expression per donor, Welch's t-test",
    "deseq2_vst": "DESeq2 VST (per-donor, variance-stabilised)",
    "deseq2_vst_zscore": "DESeq2 VST, gene-standardised (per-donor, equal weight)",
    "pseudobulk_cpm": "Pseudobulk CPM-log2 (per-donor, full-genome normalised)",
    "mean": "Mean expression (per-cell, log-normalised)",
    "zscore": "Z-score normalised mean (per-cell)",
    "singscore": "singscore (rank-based, Foroutan 2018)",
    "score_genes": "Background-corrected mean (Tirosh et al. 2016)",
}


def per_donor_pathway_score(adata, pathway_name, genes, sample_col,
                            condition_col, method=SCORING_METHOD, counts_layer=None):
    """Per-donor score for one pathway: the score the pathway test uses
    (de_analysis.compute_pathway_scores with DESeq2 normalisation).

    Returns
    -------
    (scores, sample_df) : (numpy.ndarray, pandas.DataFrame)
        One score per donor; 'sample_df' has 'sample' and 'condition'
        row-aligned to 'scores'. Empty array + empty frame if the pathway
        fails coverage.
    """
    if method != SCORING_METHOD:
        raise ValueError(f"Pathway scoring method {method!r} is no longer offered.")
    from kosmic.de.de_analysis import compute_pathway_scores
    score_matrix, sample_df, names = compute_pathway_scores(
        adata, {pathway_name: list(genes)}, sample_col, condition_col,
        normalization='deseq2', counts_layer=counts_layer)
    if score_matrix.size == 0 or pathway_name not in names:
        return np.array([]), pd.DataFrame()
    return (np.asarray(score_matrix[:, names.index(pathway_name)], dtype=float),
            sample_df)
