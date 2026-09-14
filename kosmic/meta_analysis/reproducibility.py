# Cross-dataset reproducibility analysis (pre-meta-analysis).
#
# Entry point is 'run_reproducibility'; everything else is internal.
# Three replication modes are computed in parallel:
#
#     strict         -- pvals_adj < threshold AND same direction
#     loose          -- pvals     < threshold AND same direction
#     direction_only -- same sign of log2FC, no p-value test
#
# All tables use the "testable" denominator convention: genes not
# measured in a comparison study don't count as failures.

from __future__ import annotations

from typing import List, Optional

import numpy as np
import pandas as pd
from kosmic import DEFAULT_FDR


REPLICATION_MODES = ('strict', 'loose', 'direction_only')

REPLICATION_LABELS = {
    'strict':         'Strict (FDR < threshold + same direction)',
    'loose':          'Loose (raw p < threshold + same direction)',
    'direction_only': 'Direction-only (same sign of log2FC)',
}

# direction_only -> 'direction' suffix for column names; the others pass through.
_MODE_SUFFIX = {
    'strict':         'strict',
    'loose':          'loose',
    'direction_only': 'direction',
}


def _significant_signed(df: pd.DataFrame,
                        mode: str,
                        threshold: float) -> dict:
    """Return '{gene: +1|-1}' for significantly-called genes (first occurrence wins)."""
    if df is None or len(df) == 0:
        return {}

    lfc = df['logfoldchanges'].to_numpy(dtype=float)
    if mode == 'direction_only':
        sig_mask = np.ones(len(df), dtype=bool)
    elif mode == 'loose':
        if 'pvals' not in df.columns:
            return {}
        p = df['pvals'].to_numpy(dtype=float)
        sig_mask = np.isfinite(p) & (p < threshold)
    else:  # strict
        if 'pvals_adj' not in df.columns:
            return {}
        padj = df['pvals_adj'].to_numpy(dtype=float)
        sig_mask = np.isfinite(padj) & (padj < threshold)

    keep = sig_mask & np.isfinite(lfc) & (lfc != 0)
    if not keep.any():
        return {}

    sub = pd.DataFrame({
        'name': df['names'].astype(str).to_numpy()[keep],
        'sign': np.sign(lfc[keep]).astype(int),
    })
    sub = sub[sub['name'] != ''].drop_duplicates(subset='name', keep='first')
    return dict(zip(sub['name'], sub['sign']))


def _all_mode_sig_maps(datasets: List[pd.DataFrame],
                       threshold: float) -> dict:
    """'{mode: [per-study {gene: +1|-1}]}' for all three modes."""
    return {
        mode: [_significant_signed(d, mode, threshold) for d in datasets]
        for mode in REPLICATION_MODES
    }


def _compute_pairwise_agreement(datasets: List[pd.DataFrame],
                                labels: List[str],
                                mode: str = 'strict',
                                threshold: float = DEFAULT_FDR) -> pd.DataFrame:
    """K x K row-wise agreement matrix.

    Cell (i, j) = fraction of study i's sig genes that are also sig in
    j with the same direction, over genes measured in j. NaN when study
    i has no sig genes.
    """
    if mode not in REPLICATION_MODES:
        raise ValueError(
            f"mode must be one of {REPLICATION_MODES!r}, got {mode!r}")
    K = len(datasets)
    if len(labels) != K:
        raise ValueError(
            f"len(labels) ({len(labels)}) != len(datasets) ({K})")
    if K < 2:
        raise ValueError(
            f"Pairwise agreement needs at least 2 studies; got K={K}.")

    sig_signed = [_significant_signed(d, mode, threshold) for d in datasets]
    measured = [
        set(str(g) for g in d['names'].dropna().astype(str).tolist())
        if d is not None else set()
        for d in datasets
    ]

    out = np.full((K, K), np.nan, dtype=float)
    for i in range(K):
        sig_i = sig_signed[i]
        for j in range(K):
            if i == j:
                out[i, j] = 1.0 if sig_i else np.nan
                continue
            testable = [g for g in sig_i if g in measured[j]]
            if not testable:
                out[i, j] = np.nan
                continue
            sig_j = sig_signed[j]
            agreements = sum(1 for g in testable
                             if g in sig_j and sig_j[g] == sig_i[g])
            out[i, j] = agreements / len(testable)

    return pd.DataFrame(out, index=list(labels), columns=list(labels))


def _compute_all_modes(datasets: List[pd.DataFrame],
                       labels: List[str],
                       threshold: float = DEFAULT_FDR) -> dict:
    """Return '{mode: K x K DataFrame}' for all three modes."""
    return {
        mode: _compute_pairwise_agreement(
            datasets, labels, mode=mode, threshold=threshold)
        for mode in REPLICATION_MODES
    }


def _compute_gene_reproducibility(datasets: List[pd.DataFrame],
                                  labels: List[str],
                                  threshold: float = DEFAULT_FDR
                                  ) -> pd.DataFrame:
    """Per-gene reproducibility counts across studies.

    One row per gene observed in any study. 'n_studies_<mode>' reports
    the count of studies agreeing on the DOMINANT (majority) direction
    -- a gene called up in 4, down in 1 reports 4 (not 5). Reproducibility
    therefore includes both significance and direction concordance.
    """
    K = len(datasets)
    if len(labels) != K:
        raise ValueError(
            f"len(labels) ({len(labels)}) != len(datasets) ({K})")

    sig_maps = _all_mode_sig_maps(datasets, threshold)

    lfc_sum: dict[str, float] = {}
    lfc_count: dict[str, int] = {}
    all_genes: set = set()
    for d in datasets:
        if d is None or len(d) == 0:
            continue
        names = d['names'].astype(str).to_numpy()
        lfcs = d['logfoldchanges'].to_numpy(dtype=float)
        for n, lfc in zip(names, lfcs):
            if not n:
                continue
            all_genes.add(n)
            if np.isfinite(lfc):
                lfc_sum[n] = lfc_sum.get(n, 0.0) + float(lfc)
                lfc_count[n] = lfc_count.get(n, 0) + 1

    rows = []
    for g in sorted(all_genes):
        row = {
            'gene':        g,
            'mean_log2FC': (lfc_sum.get(g, 0.0) / lfc_count[g]
                            if lfc_count.get(g, 0) > 0 else float('nan')),
            'n_measured':  lfc_count.get(g, 0),
        }
        n_pos_strict = n_neg_strict = 0
        for mode in REPLICATION_MODES:
            n_pos = n_neg = 0
            for study_sigmap in sig_maps[mode]:
                s = study_sigmap.get(g)
                if s == 1:
                    n_pos += 1
                elif s == -1:
                    n_neg += 1
            if mode == 'strict':
                n_pos_strict, n_neg_strict = n_pos, n_neg
            if n_pos > n_neg:
                dominant, n_concordant = 1, n_pos
            elif n_neg > n_pos:
                dominant, n_concordant = -1, n_neg
            else:  # tie (including 0-vs-0)
                dominant, n_concordant = 0, n_pos
            sfx = _MODE_SUFFIX[mode]
            row[f'n_studies_{sfx}']          = n_concordant
            row[f'dominant_direction_{sfx}'] = dominant
        # direction_concordant is only shown for strict in the GUI, so
        # we only compute + store it for that mode.
        n_total_strict = n_pos_strict + n_neg_strict
        row['direction_concordant_strict'] = (
            n_total_strict > 0
            and (n_pos_strict == 0 or n_neg_strict == 0))
        rows.append(row)

    df = pd.DataFrame(rows)
    if 'n_studies_strict' in df.columns:
        df = df.sort_values(
            ['n_studies_strict', 'n_studies_loose', 'n_studies_direction'],
            ascending=False).reset_index(drop=True)
    return df


def _compute_orphan_stats(datasets: List[pd.DataFrame],
                          labels: List[str],
                          threshold: float = DEFAULT_FDR,
                          repro_df: Optional[pd.DataFrame] = None
                          ) -> pd.DataFrame:
    """Per-study orphan stats.

    For each study + mode: n_DEGs (sig in that study, direction matches
    cross-study majority), n_orphans (of those, called nowhere else with
    same direction), orphan_pct, n_in_geq_half (how many of the study's
    DEGs are also called in at least ceil(K/2) total studies).
    """
    K = len(datasets)
    if len(labels) != K:
        raise ValueError(
            f"len(labels) ({len(labels)}) != len(datasets) ({K})")

    if repro_df is None:
        repro_df = _compute_gene_reproducibility(
            datasets, labels, threshold=threshold)

    repro_lookup: dict[str, dict] = {
        row['gene']: row for _, row in repro_df.iterrows()
    }
    half_thresh = max(1, (K + 1) // 2)
    sig_maps = _all_mode_sig_maps(datasets, threshold)

    rows = []
    for i, lbl in enumerate(labels):
        row = {'study': lbl}
        for mode in REPLICATION_MODES:
            sfx = _MODE_SUFFIX[mode]
            study_sigmap = sig_maps[mode][i]
            n_DEGs = n_orphans = n_in_geq_half = 0
            for g, sign in study_sigmap.items():
                info = repro_lookup.get(g)
                if info is None:
                    continue
                dom = info[f'dominant_direction_{sfx}']
                n_studies_g = info[f'n_studies_{sfx}']
                # dom == 0 means tie -- still count the call
                if dom == 0 or sign == dom:
                    n_DEGs += 1
                    if n_studies_g <= 1:
                        n_orphans += 1
                    if n_studies_g >= half_thresh:
                        n_in_geq_half += 1
            row[f'n_DEGs_{sfx}']         = n_DEGs
            row[f'n_orphans_{sfx}']      = n_orphans
            row[f'orphan_pct_{sfx}']     = (
                100.0 * n_orphans / n_DEGs if n_DEGs > 0 else float('nan'))
            row[f'n_in_geq_half_{sfx}']  = n_in_geq_half
        rows.append(row)

    return pd.DataFrame(rows)


def _reproducibility_summary(repro_df: pd.DataFrame,
                             orphan_df: pd.DataFrame,
                             n_studies_total: int) -> dict:
    """Headline scalars per mode: mean orphan rate, max gene reproducibility,
    genes called in >= ceil(K/2) studies, and unique DEG count."""
    out = {
        'n_studies':      n_studies_total,
        'half_threshold': max(1, (n_studies_total + 1) // 2),
        'modes': {},
    }
    for mode in REPLICATION_MODES:
        sfx = _MODE_SUFFIX[mode]
        col = f'n_studies_{sfx}'
        if col not in repro_df.columns:
            continue
        n_studies_per_gene = repro_df[col].to_numpy(dtype=int)
        n_unique = int((n_studies_per_gene >= 1).sum())
        in_half = int((n_studies_per_gene >= out['half_threshold']).sum())
        max_n = int(n_studies_per_gene.max()) if len(n_studies_per_gene) else 0

        finite = orphan_df[f'orphan_pct_{sfx}'].dropna()
        mean_orphan = float(finite.mean()) if len(finite) > 0 else float('nan')

        out['modes'][sfx] = {
            'mean_orphan_pct':       mean_orphan,
            'max_n_studies':         max_n,
            'n_genes_in_geq_half':   in_half,
            'pct_genes_in_geq_half': (
                100.0 * in_half / n_unique if n_unique > 0 else 0.0),
            'n_unique_DEGs':         n_unique,
        }
    return out


def _pairwise_summary(matrix: pd.DataFrame) -> dict:
    """Off-diagonal mean/min/max agreement (diagonal is trivially 1.0)."""
    arr = matrix.to_numpy(dtype=float)
    K = arr.shape[0]
    off = [arr[i, j] for i in range(K) for j in range(K)
           if i != j and np.isfinite(arr[i, j])]
    if not off:
        return {'mean': float('nan'), 'min': float('nan'),
                'max': float('nan'), 'n_pairs': 0}
    return {
        'mean':    float(np.mean(off)),
        'min':     float(np.min(off)),
        'max':     float(np.max(off)),
        'n_pairs': len(off),
    }


def run_reproducibility(datasets: List[pd.DataFrame],
                        labels: List[str],
                        threshold: float = DEFAULT_FDR,
                        progress_cb=None) -> dict:
    """Full cross-dataset reproducibility analysis. Safe off the GUI thread.

    Returns ``{'matrices': {mode: K x K DataFrame}, 'repro_df': per-gene,
    'orphan_df': per-study, 'summary': scalar headline, 'pw_summaries':
    {strict, direction_only}, 'n_studies': K}``.
    """
    if progress_cb is not None:
        progress_cb('Computing pairwise matrices...')
    matrices = _compute_all_modes(datasets, labels, threshold=threshold)

    if progress_cb is not None:
        progress_cb('Computing per-gene reproducibility...')
    repro_df = _compute_gene_reproducibility(
        datasets, labels, threshold=threshold)

    if progress_cb is not None:
        progress_cb('Computing per-study orphan stats...')
    orphan_df = _compute_orphan_stats(
        datasets, labels, threshold=threshold, repro_df=repro_df)

    summary = _reproducibility_summary(
        repro_df, orphan_df, n_studies_total=len(datasets))
    pw_summaries = {
        m: _pairwise_summary(matrices[m])
        for m in ('strict', 'direction_only')
    }

    return {
        'matrices':     matrices,
        'repro_df':     repro_df,
        'orphan_df':    orphan_df,
        'summary':      summary,
        'pw_summaries': pw_summaries,
        'n_studies':    len(datasets),
    }
