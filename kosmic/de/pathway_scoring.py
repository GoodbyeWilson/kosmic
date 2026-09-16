# Pathway Scoring
# Multiple pathway scoring methods for scRNA-seq data with sample-level
# statistical testing.
#
# Methods:
#   - score_genes: background-corrected mean expression (scanpy/AddModuleScore style)
#   - mean: simple mean expression of pathway genes
#   - zscore: mean of z-scored gene expression
#   - singscore: rank-based (Foroutan et al. 2018)

from pathlib import Path

import numpy as np
import pandas as pd

from kosmic.scrna.inspect.detection import detect_species, format_gene_for_species
from kosmic import (
    DEFAULT_FDR,
    PATHWAY_MIN_GENES,
    PATHWAY_MIN_COVERAGE,
    PATHWAY_N_CTRL_GENES,
)


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


SCORING_METHODS = {
    "deseq2": "DESeq2 (per-donor log2FC, recommended)",
    "deseq2_vst": "DESeq2 VST (per-donor, variance-stabilised)",
    "deseq2_vst_zscore": "DESeq2 VST, gene-standardised (per-donor, equal weight)",
    "pseudobulk_cpm": "Pseudobulk CPM-log2 (per-donor, full-genome normalised)",
    "mean": "Mean expression (per-cell, log-normalised)",
    "zscore": "Z-score normalised mean (per-cell)",
    "singscore": "singscore (rank-based, Foroutan 2018)",
    "score_genes": "Background-corrected mean (Tirosh et al. 2016)",
}

SCORING_METHOD_HELP = {
    "deseq2": (
        "Is this pathway's activity different per donor, on a DESeq2 log2FC scale?\n"
        "Sums counts per donor, normalises with DESeq2 median-of-ratios size factors, "
        "log2-transforms, then averages the pathway's genes per donor and tests disease "
        "vs control across donors with a plain Welch t-test. The effect is the average "
        "per-gene log2 fold-change -- a genuine, cross-study-poolable log2FC. Collapsing "
        "to one value per donor absorbs inter-gene correlation (no VIF/rho) and powers "
        "the test by donors, not genes. Recommended."
    ),
    "deseq2_vst": (
        "Same donor-level pathway log2FC as DESeq2, but on the variance-stabilised\n"
        "scale. Sums counts per donor, applies the DESeq2 VST (fit on the full-genome\n"
        "pseudobulk), then averages the pathway's genes per donor (equal weight) and\n"
        "tests disease vs control with Welch's t-test. VST flattens the variance-mean\n"
        "trend so low-count genes do not inflate the mean's noise -- steadier per-donor\n"
        "scores and SEs for meta-analysis. Best for highly-expressed pathways, where VST\n"
        "is also a faithful log2FC."
    ),
    "deseq2_vst_zscore": (
        "Same per-donor VST score, but each gene is z-scored across donors before\n"
        "averaging, so every subunit contributes equally and a few high-abundance\n"
        "genes cannot dominate the module (a robustness check on the VST score). "
        "The effect it tests is a STANDARDISED score difference, not a log2FC -- do "
        "not report it as a fold change or pool it against log2FC studies. Use it to "
        "confirm the module is coherently shifted, not just driven by its loudest genes."
    ),
    "pseudobulk_cpm": (
        "Is the proportion of the transcriptome allocated to this pathway different?\n"
        "Sums counts per donor, CPM-normalises using ALL genes as library size, "
        "log2-transforms, then averages across pathway genes. Tests at the donor level "
        "with Welch's t-test. Best for meta-analysis compatibility."
    ),
    "mean": (
        "How much pathway mRNA does each cell make, relative to the rest?\n"
        "Mean of the pathway genes' log-normalised expression per cell, "
        "aggregated by donor. Library-size normalised (log1p CP10k), so it "
        "is not confounded by sequencing depth, but highly-expressed genes "
        "still dominate the mean."
    ),
    "zscore": (
        "Are pathway genes shifted toward the tails of their expression distributions?\n"
        "On log-normalised expression, each gene is z-scored across all cells "
        "first (removing magnitude differences between genes), then averaged "
        "within the pathway. Aggregated by donor. Sensitive to shifts but not "
        "dominated by highly-expressed genes."
    ),
    "singscore": (
        "Does this pathway rank higher or lower within each cell's transcriptome?\n"
        "Ranks all genes by expression within each cell (ties averaged) and scores "
        "the gene set by its mean rank, centred at the middle rank and scaled to "
        "[-0.5, +0.5] (singscore; Foroutan et al. 2018). +0.5 = the pathway's genes "
        "are the highest-expressed in the cell, -0.5 = the lowest. Rank-based, so "
        "robust to normalisation."
    ),
    "score_genes": (
        "Is this pathway more or less expressed than expression-matched control genes?\n"
        "Subtracts mean expression of randomly-selected control genes from the same "
        "expression bins (Tirosh et al. 2016 / scanpy score_genes). Sensitive to "
        "which control genes are selected and can overcorrect when there is global "
        "transcriptional change. Use with caution."
    ),
}

DEFAULT_METHOD = "deseq2"


def score_pathways(adata, pathway_gene_sets, method="score_genes",
                   min_coverage=PATHWAY_MIN_COVERAGE,
                   min_genes=PATHWAY_MIN_GENES, n_ctrl_genes=PATHWAY_N_CTRL_GENES):
    """Compute per-cell pathway scores using the specified method.

    Parameters
    ----------
    adata : anndata.AnnData
        Input data. Scored on the raw counts (see kosmic.scrna.counts).
    pathway_gene_sets : dict
        {pathway_name: [gene_list]}.
    method : str
        One of 'score_genes', 'mean', 'zscore', 'singscore'.
    min_coverage : float
        Minimum fraction of pathway genes that must be present (0-1).
        Defaults to 'PATHWAY_MIN_COVERAGE'.
    min_genes : int
        Minimum absolute number of pathway genes that must be present.
        Defaults to 'PATHWAY_MIN_GENES'. Both thresholds must be
        satisfied -- a 100%-coverage one-gene pathway still drops.
    n_ctrl_genes : int
        Number of control genes per bin for score_genes method.

    Returns
    -------
    tuple of (score_adata, available_pathways)
        score_adata: AnnData with '{pathway}_score' columns in obs.
        available_pathways: dict of {pathway_name: [available_genes]}.
    """
    # Score on the raw counts (layers['counts'], else .raw, else X).
    from kosmic.scrna.counts import counts_adata
    work = counts_adata(adata)
    work.obs = adata.obs.copy()

    # Per-cell module scores (score_genes / mean / zscore) expect
    # log-normalised values, but `adata.raw` holds RAW COUNTS in the
    # standard pipeline (snapshotted before normalisation so DE can read
    # counts). Log-normalise this throwaway copy if it still looks like
    # counts; `adata.raw` itself is never touched.
    _ensure_lognorm(work)

    species = detect_species(list(work.var_names))

    # Filter to pathways with sufficient coverage AND enough absolute genes
    available_pathways = {}
    for pathway_name, genes in pathway_gene_sets.items():
        genes_fmt = [format_gene_for_species(g, species) for g in genes]
        available = [g for g in genes_fmt if g in work.var_names]
        coverage = len(available) / len(genes) if genes else 0
        if len(available) >= min_genes and coverage >= min_coverage:
            available_pathways[pathway_name] = available

    if not available_pathways:
        return work, {}

    if method == "singscore":
        _score_singscore(work, available_pathways)
    elif method == "zscore":
        _score_zscore(work, available_pathways)
    elif method == "score_genes":
        _score_genes(work, available_pathways, n_ctrl_genes)
    else:  # "mean"
        _score_mean(work, available_pathways)

    return work, available_pathways


def _ensure_lognorm(work):
    """Log-normalise ``work`` in place if it still looks like raw counts.

    ``score_pathways`` sources the full gene set from the counts, which
    holds raw counts in the standard pipeline. Per-cell scoring expects
    log-normalised data, so detect counts (non-negative integers) and apply
    library-size + log1p; skip if already normalised (e.g. an imported raw
    slot that is float) to avoid double normalisation. ``work`` is a copy,
    so the study's counts are never mutated.
    """
    import scanpy as sc
    from scipy import sparse as sp

    X = work.X
    data = X.data if sp.issparse(X) else np.asarray(X).reshape(-1)
    if data.size == 0:
        return
    sample = data[:10000]
    looks_like_counts = bool(
        np.all(sample >= 0) and np.allclose(sample, np.rint(sample)))
    if looks_like_counts:
        sc.pp.normalize_total(work, target_sum=1e4)
        sc.pp.log1p(work)


def _score_mean(work, available_pathways):
    """Simple mean expression of pathway genes per cell."""
    for pathway_name, genes in available_pathways.items():
        expr = work[:, genes].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()
        work.obs[f'{pathway_name}_score'] = np.mean(expr.astype(np.float64), axis=1)


def _score_zscore(work, available_pathways):
    """Mean of z-scored gene expression per cell.

    Each gene is standardised across all cells first, then averaged
    within each pathway. This prevents highly-expressed genes from
    dominating the score.
    """
    # Collect all pathway genes
    all_genes = list(dict.fromkeys(g for gs in available_pathways.values() for g in gs))
    X = work[:, all_genes].X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    X = X.astype(np.float64)

    # Z-score each gene (column)
    means = np.mean(X, axis=0)
    stds = np.std(X, axis=0)
    stds[stds == 0] = 1.0
    X_z = (X - means) / stds

    gene_to_idx = {g: i for i, g in enumerate(all_genes)}

    for pathway_name, genes in available_pathways.items():
        indices = [gene_to_idx[g] for g in genes]
        work.obs[f'{pathway_name}_score'] = np.mean(X_z[:, indices], axis=1)


def _score_singscore(work, available_pathways):
    """singscore rank-based pathway score (Foroutan et al., BMC Bioinformatics 2018).

    Within each cell, all genes are ranked by expression (ties averaged).
    A gene set's score is its mean rank, centred at the middle rank and
    scaled to [-0.5, +0.5]:

        score = (mean_rank - (N+1)/2) / (N - n)

    with N genes total and n in the set. +0.5 means the set's genes are the
    highest-expressed in the cell, -0.5 the lowest, 0 the middle. The set is
    treated as undirected (an "up" signature). Being rank-based, the score is
    unaffected by the choice of normalisation.
    """
    from scipy.stats import rankdata

    X = work.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    N = X.shape[1]

    # Per-cell ranks over the full transcriptome, ties averaged (singscore).
    ranks = np.vstack([rankdata(X[i], method='average') for i in range(X.shape[0])])

    gene_to_idx = {g: i for i, g in enumerate(work.var_names)}
    mid = (N + 1) / 2.0
    for pathway_name, genes in available_pathways.items():
        idx = [gene_to_idx[g] for g in genes if g in gene_to_idx]
        n = len(idx)
        if n == 0 or N <= n:
            work.obs[f'{pathway_name}_score'] = 0.0
            continue
        mean_rank = ranks[:, idx].mean(axis=1)
        work.obs[f'{pathway_name}_score'] = (mean_rank - mid) / (N - n)


def _score_genes(work, available_pathways, n_ctrl_genes=PATHWAY_N_CTRL_GENES):
    """Background-corrected mean expression (AddModuleScore / score_genes).

    Delegates to scanpy's reference ``sc.tl.score_genes`` so the
    expression-matched control binning matches the community standard
    (equal-frequency rank bins) rather than a hand-rolled scheme. Runs on
    ``work.X`` -- log-normalised by ``score_pathways`` -- with
    ``use_raw=False`` so it never falls back to a raw slot.

    Reference: Tirosh et al. (2016), Science.
    """
    import scanpy as sc

    for pathway_name, genes in available_pathways.items():
        sc.tl.score_genes(
            work,
            gene_list=list(genes),
            ctrl_size=n_ctrl_genes,
            score_name=f'{pathway_name}_score',
            random_state=0,
            use_raw=False,
        )


def per_donor_pathway_score(adata, pathway_name, genes, sample_col,
                            condition_col, method, counts_layer=None):
    """Per-donor score for one pathway under the selected scoring method.

    Reuses each method's own per-donor computation so the returned score
    is the one that method actually tests on -- nothing is re-derived:

      - ``'pseudobulk_cpm'`` -> ``compute_pathway_scores`` (CPM-log2)
      - ``'deseq2'``         -> ``compute_pathway_scores(normalization='deseq2')``
        (median-of-ratios size-factor log2 -- the per-donor log2FC score)
      - ``'mean'`` / ``'zscore'`` / ``'singscore'`` / ``'score_genes'`` ->
        ``score_pathways`` per-cell scores averaged per donor (identical to
        the aggregation ``test_pathway_scores`` performs).

    Returns
    -------
    (scores, sample_df) : (numpy.ndarray, pandas.DataFrame)
        One score per donor; ``sample_df`` has ``'sample'`` and
        ``'condition'`` row-aligned to ``scores``. Empty array + empty
        frame if the pathway fails coverage.
    """
    gs = {pathway_name: list(genes)}

    _norm = {'deseq2': 'deseq2', 'deseq2_vst': 'vst',
             'deseq2_vst_zscore': 'vst', 'pseudobulk_cpm': 'cpm'}
    if method in _norm:
        from kosmic.de.de_analysis import compute_pathway_scores
        score_matrix, sample_df, names = compute_pathway_scores(
            adata, gs, sample_col, condition_col,
            normalization=_norm[method], counts_layer=counts_layer,
            standardize_genes=(method == 'deseq2_vst_zscore'))
        if score_matrix.size == 0 or pathway_name not in names:
            return np.array([]), pd.DataFrame()
        return (np.asarray(score_matrix[:, names.index(pathway_name)], dtype=float),
                sample_df)

    # Per-cell methods: score per cell, then average per donor.
    work, available = score_pathways(adata, gs, method=method)
    if pathway_name not in available:
        return np.array([]), pd.DataFrame()
    score_col = f'{pathway_name}_score'
    grp = work.obs.groupby(sample_col, observed=True)
    scores = grp[score_col].mean()
    conds = grp[condition_col].first().reindex(scores.index)
    sample_df = pd.DataFrame({
        'sample': scores.index.to_numpy(),
        'condition': conds.to_numpy(),
    })
    return scores.to_numpy(dtype=float), sample_df


def test_pathway_scores(score_adata, available_pathways, sample_col, condition_col,
                        control_label, disease_label):
    """Perform statistical tests on pathway scores at the sample level.

    For each pathway, aggregates scores per sample, then tests
    disease vs control using Shapiro -> Levene -> t/Welch/MWU.

    Parameters
    ----------
    score_adata : anndata.AnnData
        AnnData with '{pathway}_score' columns in obs.
    available_pathways : dict
        {pathway_name: [genes]}.
    sample_col : str
        Sample grouping column.
    condition_col : str
        Condition column.
    control_label, disease_label : str
        Condition labels.

    Returns
    -------
    pandas.DataFrame
        Columns: Pathway, Disease_Mean, Control_Mean, Effect_Size, P_Value,
        N_Disease, N_Control, Statistical_Test, P_Adjusted, Significant,
        Log10_P_Corrected.
    """
    from scipy import stats as sp_stats
    from scipy.stats import mannwhitneyu, shapiro, levene
    from kosmic.numerical import bh_fdr, neg_log10

    results = []
    for pathway_name in available_pathways:
        score_col = f'{pathway_name}_score'
        if score_col not in score_adata.obs.columns:
            continue

        # Aggregate by sample
        sample_agg = []
        for sample_id in score_adata.obs[sample_col].unique():
            mask = score_adata.obs[sample_col] == sample_id
            sample_data = score_adata.obs[mask]
            condition = sample_data[condition_col].iloc[0]
            sample_agg.append({
                'sample': sample_id,
                'condition': condition,
                'score': sample_data[score_col].mean(),
            })

        sample_df = pd.DataFrame(sample_agg)
        disease_scores = sample_df[sample_df['condition'] == disease_label]['score'].values
        control_scores = sample_df[sample_df['condition'] == control_label]['score'].values

        if len(disease_scores) < 2 or len(control_scores) < 2:
            continue

        # Statistical test selection
        test_used = "Mann-Whitney U"
        try:
            disease_normal = len(disease_scores) < 3 or shapiro(disease_scores)[1] > 0.05
            control_normal = len(control_scores) < 3 or shapiro(control_scores)[1] > 0.05

            if disease_normal and control_normal:
                _, levene_p = levene(disease_scores, control_scores)
                if levene_p > 0.05:
                    statistic, pvalue = sp_stats.ttest_ind(disease_scores, control_scores, equal_var=True)
                    test_used = "Student's t-test"
                else:
                    statistic, pvalue = sp_stats.ttest_ind(disease_scores, control_scores, equal_var=False)
                    test_used = "Welch's t-test"
            else:
                statistic, pvalue = mannwhitneyu(disease_scores, control_scores, alternative='two-sided')
        except (ValueError, RuntimeWarning):
            statistic, pvalue = mannwhitneyu(disease_scores, control_scores, alternative='two-sided')

        # Cohen's d
        pooled_std = np.sqrt(
            ((len(disease_scores) - 1) * np.var(disease_scores, ddof=1) +
             (len(control_scores) - 1) * np.var(control_scores, ddof=1)) /
            (len(disease_scores) + len(control_scores) - 2)
        )
        effect_size = (np.mean(disease_scores) - np.mean(control_scores)) / pooled_std if pooled_std > 0 else 0

        results.append({
            'Pathway': pathway_name,
            'Disease_Mean': np.mean(disease_scores),
            'Control_Mean': np.mean(control_scores),
            'Effect_Size': effect_size,
            'P_Value': pvalue,
            'N_Disease': len(disease_scores),
            'N_Control': len(control_scores),
            'Statistical_Test': test_used,
            'Disease_Scores': disease_scores.tolist(),
            'Control_Scores': control_scores.tolist(),
        })

    stats_df = pd.DataFrame(results)

    if len(stats_df) > 0:
        stats_df['P_Adjusted'] = bh_fdr(stats_df['P_Value'])
        stats_df['Significant'] = stats_df['P_Adjusted'] < DEFAULT_FDR
        stats_df['Log10_P_Corrected'] = neg_log10(stats_df['P_Adjusted'])

    return stats_df
