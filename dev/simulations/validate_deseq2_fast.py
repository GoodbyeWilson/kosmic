"""Validate the fast CC-perm DESeq2 engine against full PyDESeq2.

Compares 'kosmic/meta_analysis/deseq2_fast.py' (cached-dispersion +
10-iter IRLS) against full PyDESeq2 (re-fit per call, see
'pydeseq2_reference_de' below) on simulated pseudobulk data.

Two checks:

1. Observed-statistic agreement on the true labels (logFC, SE, p).
2. Null-distribution agreement under K random label permutations.

Plus a convergence audit: how many genes fail to converge in the 3
IRLS iterations the perm engine uses, and whether their stats look
suspicious.

Run from the kirk env, repo root::

    python simulations/validate_deseq2_fast.py
"""
from __future__ import annotations

import os
import sys
import time
from kosmic import DEFAULT_FDR

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # repo root: dev/<this dir>/<file>

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, pearsonr, spearmanr

from kosmic.meta_analysis.deseq2_fast import (
    _fit_lfc_irls,
    _wald_test,
    fast_deseq2_perm,
    precompute_deseq2,
)


# ---------------------------------------------------------------------------
# Reference DESeq2 (full PyDESeq2)
# ---------------------------------------------------------------------------

def pydeseq2_reference_de(expression, is_disease, two_group_mask):
    """Full PyDESeq2 per call on pseudobulk sums. Slow -- re-fits
    dispersions every call -- but mathematically the canonical answer
    against which 'fast_deseq2_perm' is validated.

    Parameters
    ----------
    expression : ndarray (n_samples, n_genes)
        Sum counts per sample (pseudobulk).
    is_disease : ndarray of bool, shape (n_samples,)
        True iff this sample's role is 'disease'.
    two_group_mask : ndarray of bool, shape (n_samples,)
        True iff this sample is in disease union control.

    Returns
    -------
    logfcs : ndarray (n_genes,)
    pvals : ndarray (n_genes,)
    ses : ndarray (n_genes,)
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    n_genes = expression.shape[1]

    if two_group_mask.sum() < 2:
        return np.zeros(n_genes), np.ones(n_genes), np.full(n_genes, np.nan)

    expr = expression[two_group_mask]
    is_disease_sub = is_disease[two_group_mask]
    synth_conds = np.where(is_disease_sub, 'disease', 'control')
    sample_ids = [f"S{i}" for i in range(len(synth_conds))]

    counts_df = pd.DataFrame(
        expr.round().astype(int),
        index=sample_ids,
        columns=[f"G{i}" for i in range(n_genes)],
    )

    keep = ((counts_df > 0).sum(axis=0) >= 2) & (counts_df.sum(axis=0) > 0)
    if keep.sum() == 0:
        return np.zeros(n_genes), np.ones(n_genes), np.full(n_genes, np.nan)

    counts_filtered = counts_df.loc[:, keep]
    metadata_df = pd.DataFrame(
        {'condition': pd.Categorical(
            synth_conds, categories=['control', 'disease'])},
        index=sample_ids,
    )

    try:
        dds = DeseqDataSet(
            counts=counts_filtered, metadata=metadata_df,
            design="~condition", quiet=True,
        )
        dds.deseq2()
        ds = DeseqStats(
            dds, contrast=["condition", "disease", "control"],
            alpha=DEFAULT_FDR, quiet=True,
        )
        ds.summary()
        results = ds.results_df
    except (ValueError, RuntimeError, KeyError, AttributeError,
            np.linalg.LinAlgError):
        return np.zeros(n_genes), np.ones(n_genes), np.full(n_genes, np.nan)

    logfcs = np.zeros(n_genes)
    pvals = np.ones(n_genes)
    ses = np.full(n_genes, np.nan)
    for col_name in results.index:
        idx = int(col_name[1:])
        row = results.loc[col_name]
        lfc = row['log2FoldChange']
        se = row['lfcSE']
        p = row['pvalue']
        logfcs[idx] = lfc if np.isfinite(lfc) else 0.0
        ses[idx] = se if np.isfinite(se) else np.nan
        pvals[idx] = p if np.isfinite(p) else 1.0

    return logfcs, pvals, ses


# ---------------------------------------------------------------------------
# Pseudobulk simulator
# ---------------------------------------------------------------------------

def simulate_pseudobulk(n_genes: int = 1500,
                        n_per_group: int = 6,
                        n_de: int = 150,
                        lfc_magnitude: float = 1.0,
                        seed: int = 42) -> dict:
    """Generate one pseudobulk count matrix with known DE genes.

    NB-distributed counts; DESeq2-style dispersion-mean trend
    (alpha = 0.1 + 5/mu); library-size variation 0.5x-2x.
    """
    rng = np.random.default_rng(seed)

    log_baseline = rng.normal(4.0, 2.0, n_genes)
    baseline_mu = np.clip(np.exp(log_baseline), 1.0, 1e5)
    dispersions_true = 0.1 + 5.0 / baseline_mu
    dispersions_true = np.clip(dispersions_true, 0.01, 5.0)

    de_idx = rng.choice(n_genes, n_de, replace=False)
    de_idx.sort()
    direction = rng.choice([-1, 1], size=n_de)
    true_lfc = np.zeros(n_genes)
    true_lfc[de_idx] = direction * lfc_magnitude
    is_de = np.zeros(n_genes, dtype=bool)
    is_de[de_idx] = True

    n_samples = 2 * n_per_group
    lib_factor = rng.uniform(0.5, 2.0, n_samples)
    conditions = np.array(['control'] * n_per_group + ['disease'] * n_per_group)

    counts = np.zeros((n_samples, n_genes), dtype=np.int64)
    for s in range(n_samples):
        mu_s = baseline_mu * lib_factor[s]
        if conditions[s] == 'disease':
            mu_s = mu_s * (2.0 ** true_lfc)
        r = 1.0 / dispersions_true
        p = r / (r + mu_s)
        p = np.clip(p, 1e-10, 1 - 1e-10)
        counts[s] = rng.negative_binomial(r, p)

    return {
        'counts': counts.astype(np.float64),
        'conditions': conditions,
        'true_lfc': true_lfc,
        'is_de': is_de,
        'dispersions_true': dispersions_true,
    }


# ---------------------------------------------------------------------------
# Validation runs
# ---------------------------------------------------------------------------

def _masks_from_conditions(conds: np.ndarray) -> tuple:
    """Build (is_disease, two_group_mask) for the simulator's labels.

    Simulator uses literal 'disease' / 'control' strings; ADR-029 made
    engines mask-only. This helper resolves once at the script
    boundary, mirroring what the GUI does at the Inspect-tab boundary.
    """
    is_disease = (conds == 'disease')
    is_control = (conds == 'control')
    return is_disease, (is_disease | is_control)


def run_observed(sim: dict) -> dict:
    """Run both engines on the true labels and align gene indices."""
    counts = sim['counts']
    conds = sim['conditions']
    is_disease, two_group = _masks_from_conditions(conds)

    t0 = time.time()
    lfc_ref, p_ref, se_ref = pydeseq2_reference_de(
        counts, is_disease, two_group)
    t_ref = time.time() - t0

    t0 = time.time()
    pc = precompute_deseq2(counts, is_disease, two_group)
    t_pre = time.time() - t0

    t0 = time.time()
    lfc_fast, p_fast, se_fast = fast_deseq2_perm(
        pc, is_disease, two_group)
    t_fast = time.time() - t0

    return {
        'pc': pc,
        'lfc_ref': lfc_ref, 'p_ref': p_ref, 'se_ref': se_ref,
        'lfc_fast': lfc_fast, 'p_fast': p_fast, 'se_fast': se_fast,
        't_ref': t_ref, 't_pre': t_pre, 't_fast': t_fast,
    }


def run_perm_null(pc: dict, conditions: np.ndarray,
                  n_perms: int, seed: int = 0) -> dict:
    """Run N permutations through both engines on identical perm indices.

    Returns aligned matrices of (n_perms, n_genes) for stat comparison.
    """
    rng = np.random.default_rng(seed)
    n_genes = pc['n_genes_total']
    n_samples = len(conditions)

    is_disease_full, two_group_full = _masks_from_conditions(conditions)

    fast_p = np.zeros((n_perms, n_genes))
    fast_lfc = np.zeros((n_perms, n_genes))
    ref_p = np.zeros((n_perms, n_genes))
    ref_lfc = np.zeros((n_perms, n_genes))

    counts_full = np.zeros((n_samples, n_genes))
    counts_full[:, pc['gene_mask']] = pc['counts']

    t_fast = 0.0
    t_ref = 0.0

    for k in range(n_perms):
        perm_idx = rng.permutation(n_samples)
        is_disease_perm = is_disease_full[perm_idx]
        two_group_perm = two_group_full[perm_idx]

        t0 = time.time()
        lf_f, pp_f, _ = fast_deseq2_perm(pc, is_disease_perm, two_group_perm)
        t_fast += time.time() - t0

        t0 = time.time()
        lf_r, pp_r, _ = pydeseq2_reference_de(
            counts_full, is_disease_perm, two_group_perm)
        t_ref += time.time() - t0

        fast_lfc[k] = lf_f
        fast_p[k] = pp_f
        ref_lfc[k] = lf_r
        ref_p[k] = pp_r

    return {
        'fast_lfc': fast_lfc, 'fast_p': fast_p,
        'ref_lfc': ref_lfc, 'ref_p': ref_p,
        't_fast': t_fast, 't_ref': t_ref,
    }


def convergence_audit(pc: dict, conditions: np.ndarray) -> dict:
    """Compare 10-iter (production) vs 50-iter IRLS. Reports per-gene
    non-convergence and the magnitude of any β disagreement.

    Production also pre-filters genes with all-zero counts in either
    group; that mask is computed here too so the audit reflects what
    actually reaches IRLS."""
    counts = pc['counts']
    sf = pc['size_factors']
    disp = pc['dispersions']

    two_group = np.isin(conditions, ['disease', 'control'])
    is_disease = (conditions[two_group] == 'disease')
    counts_sub = counts[two_group]
    sf_sub = sf[two_group]

    ctrl_zero = counts_sub[~is_disease].sum(axis=0) == 0
    dis_zero = counts_sub[is_disease].sum(axis=0) == 0
    n_filtered = int((ctrl_zero | dis_zero).sum())

    _, lfc10, _, conv10 = _fit_lfc_irls(
        counts_sub, sf_sub, disp, is_disease, max_iter=10)
    _, lfc50, _, conv50 = _fit_lfc_irls(
        counts_sub, sf_sub, disp, is_disease, max_iter=50)

    # Restrict the convergence stat to genes that actually reach IRLS
    # in production (i.e. exclude the zero-in-one-group prefilter).
    fittable = ~(ctrl_zero | dis_zero)
    diff = np.abs(lfc10[fittable] - lfc50[fittable])
    return {
        'n_genes_kept': counts_sub.shape[1],
        'n_zero_in_one_group_filtered': n_filtered,
        'frac_converged_10': float(conv10[fittable].mean()) if fittable.any() else float('nan'),
        'frac_converged_50': float(conv50[fittable].mean()) if fittable.any() else float('nan'),
        'lfc_diff_max': float(diff.max()) if diff.size else 0.0,
        'lfc_diff_mean': float(diff.mean()) if diff.size else 0.0,
        'lfc_diff_p99': float(np.quantile(diff, 0.99)) if diff.size else 0.0,
    }


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def _agreement(name: str, ref: np.ndarray, fast: np.ndarray) -> dict:
    """Pearson + Spearman + max abs diff on a paired stat vector."""
    mask = np.isfinite(ref) & np.isfinite(fast)
    if mask.sum() < 3:
        return {'name': name, 'pearson': np.nan, 'spearman': np.nan,
                'max_abs': np.nan, 'mean_abs': np.nan}
    r_p = pearsonr(ref[mask], fast[mask]).statistic
    r_s = spearmanr(ref[mask], fast[mask]).statistic
    diff = np.abs(ref[mask] - fast[mask])
    return {'name': name, 'pearson': float(r_p), 'spearman': float(r_s),
            'max_abs': float(diff.max()), 'mean_abs': float(diff.mean())}


def report_observed(obs: dict, sim: dict) -> None:
    print('\n=== OBSERVED (true labels) ===')
    print(f'  PyDESeq2 reference:  {obs["t_ref"]:.2f}s')
    print(f'  precompute_deseq2:   {obs["t_pre"]:.2f}s')
    print(f'  fast_deseq2_perm:    {obs["t_fast"]:.3f}s')

    # Restrict to genes both engines kept (PyDESeq2 may drop some).
    valid = (np.isfinite(obs['p_ref']) & np.isfinite(obs['p_fast'])
             & (obs['p_ref'] > 0) & (obs['p_fast'] > 0))
    print(f'  comparable genes:    {valid.sum()} / {len(obs["p_ref"])}')

    rows = [
        _agreement('logFC',   obs['lfc_ref'][valid], obs['lfc_fast'][valid]),
        _agreement('SE',      obs['se_ref'][valid],  obs['se_fast'][valid]),
        _agreement('-log10p', -np.log10(obs['p_ref'][valid]),
                              -np.log10(obs['p_fast'][valid])),
    ]
    print('\n  Stat        Pearson   Spearman   max|diff|   mean|diff|')
    for r in rows:
        print(f'  {r["name"]:8s}    {r["pearson"]:.5f}   {r["spearman"]:.5f}    '
              f'{r["max_abs"]:.4f}      {r["mean_abs"]:.4f}')

    # DE recovery check
    is_de = sim['is_de']
    sig_ref = (obs['p_ref'] < DEFAULT_FDR) & is_de
    sig_fast = (obs['p_fast'] < DEFAULT_FDR) & is_de
    print(f'\n  True DE genes called sig (p<0.05):')
    print(f'    reference: {sig_ref.sum()} / {is_de.sum()}')
    print(f'    fast:      {sig_fast.sum()} / {is_de.sum()}')


def report_null(perm: dict, n_genes_total: int) -> None:
    print('\n=== NULL DISTRIBUTION (label permutations) ===')
    n_perms = perm['fast_lfc'].shape[0]
    print(f'  perms:               {n_perms}')
    print(f'  fast total:          {perm["t_fast"]:.2f}s '
          f'({perm["t_fast"] / n_perms * 1000:.1f} ms/perm)')
    print(f'  reference total:     {perm["t_ref"]:.2f}s '
          f'({perm["t_ref"] / n_perms * 1000:.1f} ms/perm)')
    print(f'  speedup:             {perm["t_ref"] / perm["t_fast"]:.0f}x')

    # Per-perm gene-wise correlation
    valid = (np.isfinite(perm['fast_p']) & np.isfinite(perm['ref_p'])
             & (perm['fast_p'] > 0) & (perm['ref_p'] > 0))

    perm_corrs_lfc = []
    perm_corrs_lp = []
    for k in range(n_perms):
        m = valid[k]
        if m.sum() < 10:
            continue
        perm_corrs_lfc.append(spearmanr(
            perm['fast_lfc'][k][m], perm['ref_lfc'][k][m]).statistic)
        perm_corrs_lp.append(spearmanr(
            -np.log10(perm['fast_p'][k][m]),
            -np.log10(perm['ref_p'][k][m])).statistic)

    perm_corrs_lfc = np.array(perm_corrs_lfc)
    perm_corrs_lp = np.array(perm_corrs_lp)
    print(f'\n  per-perm Spearman correlations (fast vs reference):')
    print(f'    logFC:    median={np.median(perm_corrs_lfc):.4f}  '
          f'min={perm_corrs_lfc.min():.4f}  '
          f'max={perm_corrs_lfc.max():.4f}')
    print(f'    -log10p:  median={np.median(perm_corrs_lp):.4f}  '
          f'min={perm_corrs_lp.min():.4f}  '
          f'max={perm_corrs_lp.max():.4f}')

    # Aggregate null p-value KS test (uniformity)
    all_fast_p = perm['fast_p'][valid]
    all_ref_p = perm['ref_p'][valid]
    print(f'\n  null p-value uniformity (KS vs U(0,1)):')
    fast_ks = ks_2samp(all_fast_p, np.random.uniform(0, 1, len(all_fast_p)))
    ref_ks = ks_2samp(all_ref_p, np.random.uniform(0, 1, len(all_ref_p)))
    print(f'    fast:       D={fast_ks.statistic:.4f}  p={fast_ks.pvalue:.3g}')
    print(f'    reference:  D={ref_ks.statistic:.4f}  p={ref_ks.pvalue:.3g}')

    # Two-sample KS: do fast and reference draw from the same null?
    cross_ks = ks_2samp(all_fast_p, all_ref_p)
    print(f'\n  fast vs reference null distributions (two-sample KS):')
    print(f'    D={cross_ks.statistic:.4f}  p={cross_ks.pvalue:.3g}')
    print(f'    (large p = the two engines produce indistinguishable nulls)')


def report_convergence(conv: dict) -> None:
    print('\n=== IRLS CONVERGENCE AUDIT (production: max_iter=10) ===')
    print(f'  genes after PyDESeq2 filter:                {conv["n_genes_kept"]}')
    print(f'  zero-in-one-group prefilter drops:          '
          f'{conv["n_zero_in_one_group_filtered"]}')
    print(f'  fittable genes converged in 10 iter:        '
          f'{conv["frac_converged_10"] * 100:.2f}%')
    print(f'  fittable genes converged in 50 iter:        '
          f'{conv["frac_converged_50"] * 100:.2f}%')
    print(f'  |lfc(10-iter) - lfc(50-iter)| on fittable genes:')
    print(f'    mean:  {conv["lfc_diff_mean"]:.6f}')
    print(f'    p99:   {conv["lfc_diff_p99"]:.6f}')
    print(f'    max:   {conv["lfc_diff_max"]:.6f}')
    print(f'    (genes still flagged non-converged by 10 iter are '
          'returned as NaN by production)')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(n_genes=1500, n_per_group=6, n_de=150, n_perms=20, seed=42):
    print('='*70)
    print(f'DESeq2-fast validation')
    print(f'  n_genes={n_genes}, n_per_group={n_per_group}, '
          f'n_de={n_de}, n_perms={n_perms}, seed={seed}')
    print('='*70)

    sim = simulate_pseudobulk(n_genes=n_genes, n_per_group=n_per_group,
                               n_de=n_de, seed=seed)

    obs = run_observed(sim)
    report_observed(obs, sim)

    conv = convergence_audit(obs['pc'], sim['conditions'])
    report_convergence(conv)

    perm = run_perm_null(obs['pc'], sim['conditions'],
                         n_perms=n_perms, seed=seed + 1)
    report_null(perm, n_genes)

    print('\n' + '='*70)
    print('VERDICT — summary of agreement metrics for the paper:')
    print('='*70)


if __name__ == '__main__':
    main()
