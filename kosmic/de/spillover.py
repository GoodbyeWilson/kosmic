# Ambient-spillover diagnostic for pseudobulk DE.
#
# Answers the reviewer question "is this cell type's DE just soup from the
# dominant cell type?" WITHOUT modifying the DE engine. Given a per-sample
# score (a pathway score or a single gene's pseudobulk expression), a
# per-sample contamination index (mean `pct_counts_<key>` per sample, from
# `kosmic.scrna.qc.signatures`), and disease/control labels, it reports how
# much of the disease effect is attributable to contamination.
#
# The method is a two-model comparison at the per-sample (per-donor) unit:
#   M0:  score ~ disease
#   M1:  score ~ disease + contamination
# The shrinkage of the disease coefficient from M0 to M1 is the fraction of
# the effect the contamination covariate can explain. This is deliberately
# conservative: because contamination is itself correlated with disease, M1
# can absorb genuine disease signal, so "survives adjustment" is a stringent
# test, not a lenient one.
#
# Pure logic (no Qt, no DE-engine coupling).
from __future__ import annotations

import numpy as np


def _ols(y, X):
    """Ordinary least squares. Returns (coefs, residual_df).

    `X` includes its own intercept column if wanted. Uses lstsq so it is
    robust to rank issues; residual df is n - rank(X).
    """
    coef, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    return coef, max(len(y) - rank, 1)


def spillover_attribution(score, contamination, is_disease):
    """Quantify how much of a disease effect on `score` is contamination.

    Parameters
    ----------
    score : array-like, shape (n_samples,)
        Per-sample score (pathway score or gene pseudobulk expression).
    contamination : array-like, shape (n_samples,)
        Per-sample contamination index (e.g. mean `pct_counts_cardiomyocyte`).
    is_disease : array-like of bool/int, shape (n_samples,)
        1/True for disease samples, 0/False for control.

    Returns
    -------
    dict
        raw_effect : disease coefficient with no adjustment (= mean diff).
        adjusted_effect : disease coefficient adjusting for contamination.
        pct_attributable : 100 * (raw - adjusted) / raw, clipped to [0, 100];
            NaN if the raw effect is ~0.
        raw_effect_sd, adjusted_effect_sd : same effects in SD-of-score units.
        adjusted_p : two-sided p-value for the adjusted disease coefficient.
        spearman_score_contam : Spearman rho between contamination and score.
        n_disease, n_control : sample counts.
    """
    from scipy import stats

    y = np.asarray(score, dtype=float).ravel()
    c = np.asarray(contamination, dtype=float).ravel()
    d = np.asarray(is_disease, dtype=float).ravel()

    n = len(y)
    if not (len(c) == len(d) == n):
        raise ValueError("score, contamination, is_disease must be same length")
    if n < 4:
        raise ValueError(f"need >=4 samples for the diagnostic, got {n}")

    ones = np.ones(n)

    # M0: score ~ 1 + disease
    X0 = np.column_stack([ones, d])
    coef0, _ = _ols(y, X0)
    raw_effect = float(coef0[1])

    # M1: score ~ 1 + disease + contamination
    X1 = np.column_stack([ones, d, c])
    coef1, df1 = _ols(y, X1)
    adjusted_effect = float(coef1[1])

    # p-value for adjusted disease coefficient.
    resid = y - X1 @ coef1
    sigma2 = float(resid @ resid) / df1
    try:
        xtx_inv = np.linalg.inv(X1.T @ X1)
        var_disease = sigma2 * xtx_inv[1, 1]
        # Guard against tiny negative values from near-collinear designs
        # (e.g. contamination almost perfectly predicting condition).
        se_disease = float(np.sqrt(var_disease)) if var_disease > 0 else 0.0
        t = adjusted_effect / se_disease if se_disease > 0 else np.inf
        adjusted_p = float(2 * stats.t.sf(abs(t), df1))
    except np.linalg.LinAlgError:
        adjusted_p = float('nan')

    if abs(raw_effect) < 1e-12:
        pct_attributable = float('nan')
    else:
        pct_attributable = float(
            np.clip((raw_effect - adjusted_effect) / raw_effect * 100.0, 0.0, 100.0))

    sd = float(np.std(y)) or 1.0
    rho = stats.spearmanr(c, y).statistic if n >= 3 else float('nan')

    return {
        'raw_effect': raw_effect,
        'adjusted_effect': adjusted_effect,
        'pct_attributable': pct_attributable,
        'raw_effect_sd': raw_effect / sd,
        'adjusted_effect_sd': adjusted_effect / sd,
        'adjusted_p': adjusted_p,
        'spearman_score_contam': float(rho),
        'n_disease': int(d.sum()),
        'n_control': int((1 - d).sum()),
    }


def pathway_spillover(adata, pathway_name, genes, sample_col, condition_col,
                      disease_label, control_label, contam_key,
                      method='pseudobulk_cpm', counts_layer=None):
    """Per-sample spillover diagnostic for one pathway.

    Scores the pathway per donor **using the selected scoring method**
    (via :func:`kosmic.de.pathway_scoring.per_donor_pathway_score`, so the
    diagnostic operates on the same score the Pathway Score plot shows),
    aggregates the contamination metric ``obs[contam_key]`` to a per-sample
    mean, restricts to the two conditions, and runs
    :func:`spillover_attribution`.

    Parameters
    ----------
    adata : anndata.AnnData
        Must carry ``obs[contam_key]`` (e.g. produced by the QC signature
        step) and ``obs[sample_col]`` / ``obs[condition_col]``.
    pathway_name : str
    genes : list of str
        Pathway member genes.
    sample_col, condition_col : str
    disease_label, control_label : str
    contam_key : str
        An ``obs`` column, e.g. ``'pct_counts_cardiomyocyte'``.
    method : str
        Scoring method (``'deseq2'``, ``'pseudobulk_cpm'``, ``'mean'``,
        ``'zscore'``, ``'singscore'``, ``'score_genes'``).
    counts_layer : str, optional
        Counts layer to score (e.g. ``'decontX_counts'``); scores the
        pathway on decontaminated counts so the diagnostic can be re-run
        as a before/after validation of decontamination.

    Returns
    -------
    dict
        The :func:`spillover_attribution` result, plus per-sample arrays
        ``score``, ``contamination``, ``is_disease``, ``samples`` and the
        echoed ``pathway`` / ``contam_key``. Raises ``ValueError`` with a
        clear message when a prerequisite is missing.
    """
    import numpy as np
    import pandas as pd
    from kosmic.de.pathway_scoring import per_donor_pathway_score

    if contam_key not in adata.obs.columns:
        raise ValueError(
            f"'{contam_key}' is not in obs -- score a contamination panel "
            "in the QC step first.")

    scores, sample_df = per_donor_pathway_score(
        adata, pathway_name, list(genes), sample_col, condition_col, method,
        counts_layer=counts_layer)
    if scores.size == 0:
        raise ValueError(
            f"pathway '{pathway_name}' could not be scored (too few member "
            "genes present, or coverage below threshold).")

    samples = sample_df['sample'].to_numpy()
    conditions = sample_df['condition'].astype(str).to_numpy()

    contam_by_sample = adata.obs.groupby(sample_col, observed=True)[contam_key].mean()
    contamination = pd.Series(samples).map(contam_by_sample).to_numpy(dtype=float)

    keep = (np.isin(conditions, [str(disease_label), str(control_label)])
            & ~np.isnan(contamination))
    scores, contamination = scores[keep], contamination[keep]
    is_disease = (conditions[keep] == str(disease_label)).astype(float)
    samples = samples[keep]

    res = spillover_attribution(scores, contamination, is_disease)
    res.update({
        'pathway': pathway_name,
        'contam_key': contam_key,
        'score': scores,
        'contamination': contamination,
        'is_disease': is_disease,
        'samples': samples,
    })
    return res


__all__ = ["spillover_attribution", "pathway_spillover"]
