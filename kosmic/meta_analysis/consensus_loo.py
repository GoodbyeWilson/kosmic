# Consensus Leave-One-Out Validation
# Held-out replication test for the consensus meta-analysis.
#
# For each fold j (j = 0..K-1):
#     1. Hold out study j.
#     2. Build the consensus from the remaining K-1 studies, using the
#        same pooling methods and the same calibration mode (analytical
#        or CC permutation) as the user's full-data run.
#     3. For each consensus gene g in C_(-j), check whether the held-out
#        study j independently calls g on its own using one of three
#        replication criteria:
#
#          strict           : pvals_adj < 0.05  AND same direction
#          loose            : pvals     < 0.05  AND same direction
#          direction_only   : same sign of log2FC, no p-value test
#
#     4. Report n_replicated / |C_(-j)| for that fold.
#
# Aggregate output: per-fold replication rates plus the headline number
# "average X% of consensus genes were independently replicated by the
# held-out study, across K folds (range Y-Z%)".
#
# This is **out-of-sample replication**, not a meta-analysis sensitivity
# LOO. The held-out study is never seen by the K-1 fold consensus, so
# the test is non-circular.
#
# Pure logic, no PyQt6 imports.

from __future__ import annotations

from typing import Callable, List, Optional

import numpy as np
import pandas as pd

from kosmic.meta_analysis.pooling_dispatch import (
    get_analytical_pool_fn,
    merge_consensus_methods,
)
from kosmic import DEFAULT_FDR


REPLICATION_MODES = ('strict', 'loose', 'direction_only')

REPLICATION_LABELS = {
    'strict':         'Strict (held-out FDR < 0.05 + same direction)',
    'loose':          'Loose (held-out p < 0.05 + same direction)',
    'direction_only': 'Direction-only (same sign of log2FC)',
}


# Held-out case/control AUC

def _module_score_sum_of_signed_z(held_out_pb: dict,
                                  consensus_genes: List[str],
                                  consensus_signs: np.ndarray) -> np.ndarray:
    """Signed z-score module score, one number per held-out sample.

    For each sample s and consensus gene g: z(expr[s, g]) * sign(pooled_logFC_g),
    with z computed per gene across samples in the held-out dataset.
    Up-genes contribute positively, down-genes negatively.

    Parameters
    ----------
    held_out_pb : dict
        Pseudobulk with keys 'expression' (n_samples, n_genes),
        'gene_names' (list), 'conditions' (n_samples,).
    consensus_genes : list of str
        Gene names from the K-1 consensus to use for scoring.
    consensus_signs : ndarray of +1/-1
        Pooled direction for each consensus gene (aligned with
        'consensus_genes'). Genes with sign 0 are treated as +1.

    Returns
    -------
    module_scores : ndarray (n_samples,)
    """
    expr = np.asarray(held_out_pb['expression'], dtype=float)
    gene_names = list(held_out_pb['gene_names'])
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    mask = np.array([g in gene_to_idx for g in consensus_genes])
    idxs = np.array([gene_to_idx[g] for g in consensus_genes if g in gene_to_idx])
    signs = consensus_signs[mask]
    signs = np.where(signs == 0, 1.0, signs)

    if idxs.size == 0:
        return np.zeros(expr.shape[0])

    sub = expr[:, idxs]
    mu = sub.mean(axis=0, keepdims=True)
    sd = sub.std(axis=0, ddof=0, keepdims=True)
    z = np.where(sd > 0, (sub - mu) / sd, 0.0)
    signed = z * signs[None, :]
    return signed.sum(axis=1)


def _module_score_ucell(held_out_pb: dict,
                        consensus_genes: List[str],
                        consensus_signs: np.ndarray) -> np.ndarray:
    """Rank-based module score matching SumRank paper's UCell approach.

    For each sample: rank all genes by expression (1 = lowest, N = highest),
    normalise to [0, 1]. Module score = mean normalised rank of up-genes
    minus mean normalised rank of down-genes. No absolute-expression
    information used -- robust to normalisation differences across
    datasets.

    Returns
    -------
    module_scores : ndarray (n_samples,)
    """
    from scipy.stats import rankdata

    expr = np.asarray(held_out_pb['expression'], dtype=float)
    gene_names = list(held_out_pb['gene_names'])
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    up_idx = [gene_to_idx[g] for g, s in zip(consensus_genes, consensus_signs)
              if g in gene_to_idx and s > 0]
    down_idx = [gene_to_idx[g] for g, s in zip(consensus_genes, consensus_signs)
                if g in gene_to_idx and s < 0]

    n_samples, n_genes = expr.shape
    if n_genes == 0:
        return np.zeros(n_samples)

    scores = np.zeros(n_samples)
    for i in range(n_samples):
        ranks = rankdata(expr[i, :]) / n_genes  # normalise to (0, 1]
        up_score = ranks[up_idx].mean() if up_idx else 0.5
        down_score = ranks[down_idx].mean() if down_idx else 0.5
        scores[i] = up_score - down_score
    return scores


def _resolve_pb_role_masks(held_out_pb: dict) -> tuple:
    """Return (is_disease, is_control) boolean masks for a pseudobulk.

    Roles come from the pseudobulk's 'role' column (populated by
    'load_pseudobulk_set'). Returns (None, None) if the pseudobulk
    has no role information.
    """
    role = held_out_pb.get('role')
    if role is None:
        return None, None
    role = np.asarray(role)
    return (role == 'disease'), (role == 'control')


def _held_out_auc(held_out_pb: Optional[dict],
                  consensus_df: pd.DataFrame) -> dict:
    """Compute held-out case/control AUC for a consensus gene list.

    Returns a dict with keys:
        'auc_sum_of_z':   float or nan -- primary (effect-size-aware)
        'auc_ucell':      float or nan -- SumRank-paper parity
        'n_disease':      int
        'n_control':      int
        'n_genes_used':   int -- consensus genes present in held-out
        'note':           str -- reason for nan if any

    Returns nans if pseudobulk is missing, if the consensus is empty,
    or if the held-out study has fewer than 2 samples per group.

    Per-sample roles come from the pseudobulk's 'role' column
    (populated by 'load_pseudobulk_set').
    """
    result = {
        'auc_sum_of_z': np.nan,
        'auc_ucell': np.nan,
        'n_disease': 0,
        'n_control': 0,
        'n_genes_used': 0,
        'note': '',
    }

    if held_out_pb is None:
        result['note'] = 'no pseudobulk (analytical-only LOO)'
        return result
    if consensus_df is None or len(consensus_df) == 0:
        result['note'] = 'empty consensus'
        return result

    is_disease, is_control = _resolve_pb_role_masks(held_out_pb)
    if is_disease is None:
        result['note'] = 'no role assignment in pseudobulk'
        return result

    n_dis = int(is_disease.sum())
    n_ctrl = int(is_control.sum())
    result['n_disease'] = n_dis
    result['n_control'] = n_ctrl
    if n_dis < 2 or n_ctrl < 2:
        result['note'] = f'insufficient samples (dis={n_dis}, ctrl={n_ctrl})'
        return result

    consensus_genes = consensus_df['names'].tolist()
    logfc = consensus_df['logfoldchanges'].to_numpy(dtype=float)
    signs = np.sign(logfc)

    gene_names_set = set(held_out_pb['gene_names'])
    n_used = sum(1 for g in consensus_genes if g in gene_names_set)
    result['n_genes_used'] = n_used
    if n_used == 0:
        result['note'] = 'no consensus genes measured in held-out'
        return result

    keep_mask = (is_disease | is_control)
    keep_idx = np.where(keep_mask)[0]
    y = is_disease[keep_idx].astype(int)

    pb_slim = {
        'expression': np.asarray(held_out_pb['expression'])[keep_idx, :],
        'gene_names': held_out_pb['gene_names'],
    }

    try:
        from sklearn.metrics import roc_auc_score
        scores_z = _module_score_sum_of_signed_z(pb_slim, consensus_genes, signs)
        if np.std(scores_z) > 0:
            result['auc_sum_of_z'] = float(roc_auc_score(y, scores_z))
        scores_u = _module_score_ucell(pb_slim, consensus_genes, signs)
        if np.std(scores_u) > 0:
            result['auc_ucell'] = float(roc_auc_score(y, scores_u))
    except Exception as e:
        result['note'] = f'AUC failed: {type(e).__name__}: {e}'

    return result


# Replication test

def _replication_test(consensus_df: pd.DataFrame,
                      held_out_df: pd.DataFrame,
                      mode: str,
                      fdr_threshold: float = DEFAULT_FDR) -> pd.DataFrame:
    """Per-gene replication table for one held-out study.

    Returns a DataFrame with one row per consensus gene:

        gene, consensus_log2FC, study_log2FC, study_pval,
        study_padj, replicated (bool), reason (str)

    'reason' is one of:
        'replicated'           -- passed the threshold
        'wrong_direction'      -- detected but opposite sign
        'not_significant'      -- detected, same direction, but p too high
        'missing'              -- gene not in held_out_df
    """
    if mode not in REPLICATION_MODES:
        raise ValueError(
            f"mode must be one of {REPLICATION_MODES!r}, got {mode!r}")

    _EMPTY_COLS = [
        'gene', 'consensus_log2FC', 'abs_consensus_log2FC',
        'consensus_fdr', 'heterogeneity_i2', 'n_studies', 'se_approx',
        'study_log2FC', 'study_pval', 'study_padj',
        'replicated', 'reason',
    ]
    if consensus_df is None or len(consensus_df) == 0:
        return pd.DataFrame({c: [] for c in _EMPTY_COLS})

    h = held_out_df.set_index('names', drop=False)

    rows = []
    for _, c in consensus_df.iterrows():
        gene = c['names']
        c_lfc = c.get('logfoldchanges', np.nan)

        if gene not in h.index:
            rows.append({
                'gene':              gene,
                'consensus_log2FC':  c_lfc,
                'study_log2FC':      np.nan,
                'study_pval':        np.nan,
                'study_padj':        np.nan,
                'replicated':        False,
                'reason':            'missing',
            })
            continue

        hr = h.loc[gene]
        if isinstance(hr, pd.DataFrame):  # gene appears more than once
            hr = hr.iloc[0]

        s_lfc  = hr.get('logfoldchanges', np.nan)
        s_p    = hr.get('pvals',         np.nan)
        s_padj = hr.get('pvals_adj',     np.nan)

        same_dir = (
            np.isfinite(c_lfc) and np.isfinite(s_lfc)
            and np.sign(c_lfc) == np.sign(s_lfc)
            and c_lfc != 0 and s_lfc != 0
        )

        if mode == 'direction_only':
            replicated = bool(same_dir)
            reason = 'replicated' if replicated else 'wrong_direction'
        elif mode == 'loose':
            sig = bool(np.isfinite(s_p) and s_p < fdr_threshold)
            replicated = bool(same_dir and sig)
            if replicated:
                reason = 'replicated'
            elif not same_dir:
                reason = 'wrong_direction'
            else:
                reason = 'not_significant'
        else:  # strict
            sig = bool(np.isfinite(s_padj) and s_padj < fdr_threshold)
            replicated = bool(same_dir and sig)
            if replicated:
                reason = 'replicated'
            elif not same_dir:
                reason = 'wrong_direction'
            else:
                reason = 'not_significant'

        # Features used downstream to answer "what predicts replication?".
        i2 = c.get('heterogeneity_i2', np.nan)
        n_studies = c.get('n_studies', np.nan)
        ci_lo = c.get('ci_lower', np.nan)
        ci_hi = c.get('ci_upper', np.nan)
        # SE from the normal-theory 95% CI (beta +/- 1.96*SE); NaN if absent.
        if (pd.notna(ci_lo) and pd.notna(ci_hi)
                and np.isfinite(ci_lo) and np.isfinite(ci_hi)):
            se_approx = float((ci_hi - ci_lo) / (2.0 * 1.96))
        else:
            se_approx = float('nan')

        rows.append({
            'gene':                  gene,
            'consensus_log2FC':      c_lfc,
            'abs_consensus_log2FC':  (abs(c_lfc) if np.isfinite(c_lfc)
                                      else float('nan')),
            'consensus_fdr':         c.get('fdr', float('nan')),
            'heterogeneity_i2':      float(i2) if pd.notna(i2)
                                     else float('nan'),
            'n_studies':             int(n_studies) if pd.notna(n_studies)
                                     else 0,
            'se_approx':             se_approx,
            'study_log2FC':          s_lfc,
            'study_pval':            s_p,
            'study_padj':            s_padj,
            'replicated':            replicated,
            'reason':                reason,
        })

    return pd.DataFrame(rows)


# Per-fold consensus rebuild

def _build_fold_consensus(fold_datasets: List[pd.DataFrame],
                          fold_pseudobulk: Optional[List[dict]],
                          method_keys: List[str],
                          calibration: str,
                          min_studies: int,
                          n_cc_perms: int,
                          de_method: str,
                          seed: int,
                          progress_cb: Optional[Callable[[str], None]] = None
                          ) -> pd.DataFrame:
    """Run the consensus pipeline on a single fold's K-1 datasets.

    Returns the merged consensus DataFrame (same shape as the full-data
    consensus). Mirrors GeneMAWorker.run() without the PyQt signals.
    """
    analytical_dfs = {}
    for key in method_keys:
        if progress_cb:
            progress_cb(f"  pool {key}...")
        pool_fn = get_analytical_pool_fn(key, min_studies)
        analytical_dfs[key] = pool_fn(fold_datasets, min_studies=min_studies)

    if calibration == 'cc_perm':
        if not fold_pseudobulk:
            raise ValueError(
                "CC perm calibration requires pseudobulk data; none "
                "provided for this fold.")
        from kosmic.meta_analysis.cc_permutation import consensus_cc_permutation
        if progress_cb:
            progress_cb(f"  CC perm ({n_cc_perms} perms)...")
        method_dfs = consensus_cc_permutation(
            fold_pseudobulk, analytical_dfs,
            n_perms=n_cc_perms,
            seed=seed,
            de_method=de_method,
            min_studies=min_studies,
        )
    else:
        method_dfs = analytical_dfs

    return merge_consensus_methods(method_keys, method_dfs)


# Main entry point

_BASE_METHOD_KEYS = ('dl', 'sumrank', 'gwop')

# Per base method, the corresponding Top 50% variant key.
_TOP50_VARIANT = {
    'dl':      'dl_top50',
    'sumrank': 'sumrank_top50',
}

BENCHMARK_LEVELS = ('quick', 'full', 'exhaustive')


def _label_for_methods(method_keys):
    """Human-readable label for a list of method keys."""
    if not method_keys:
        return '(empty)'
    if len(method_keys) == 1:
        return method_keys[0]
    return ' + '.join(method_keys)


def enumerate_method_configurations(level: str,
                                    base_methods=None,
                                    user_top50_per_method=None
                                    ) -> List[dict]:
    """Generate method configurations to benchmark.

    Parameters
    ----------
    level : 'quick' | 'full' | 'exhaustive'
        - 'quick': singles in both Top50 states + 31 subsets of the
          user's current per-method Top50 selection.
        - 'full': each method has 3 states (Off / On / On_Top50),
          giving 3^M - 1 non-empty configs.
        - 'exhaustive': each method+Top50 variant is independent;
          2^M - 1 configs.
    base_methods : iterable of str, optional
        Base method keys to consider. Default: all of '_BASE_METHOD_KEYS'.
    user_top50_per_method : dict, optional
        Required for 'quick': fixes the Top50 dimension before sweeping subsets.

    Returns
    -------
    list of dict
        Each entry: {'label', 'method_keys', 'n_methods',
        'states' (base -> '-' / 'all' / 'top50')}.
    """
    from itertools import combinations, product

    if level not in BENCHMARK_LEVELS:
        raise ValueError(
            f"level must be one of {BENCHMARK_LEVELS!r}, got {level!r}")

    base_methods = list(base_methods) if base_methods else list(_BASE_METHOD_KEYS)

    configs = []

    if level == 'quick':
        # Stage 1: each method alone in both Top50 states.
        for m in base_methods:
            configs.append({
                'label':       _label_for_methods([m]),
                'method_keys': [m],
                'n_methods':   1,
                'states':      {bm: ('all' if bm == m else '-')
                                for bm in base_methods},
            })
            top50 = _TOP50_VARIANT.get(m)
            if top50 is not None:
                configs.append({
                    'label':       _label_for_methods([top50]),
                    'method_keys': [top50],
                    'n_methods':   1,
                    'states':      {bm: ('top50' if bm == m else '-')
                                    for bm in base_methods},
                })

        # Stage 2: subsets of the user's chosen Top50 variants
        # (missing entries default to all-studies).
        user_top50_per_method = user_top50_per_method or {}
        chosen_keys = []
        chosen_state_map = {}
        for m in base_methods:
            if user_top50_per_method.get(m) and _TOP50_VARIANT.get(m):
                chosen_keys.append(_TOP50_VARIANT[m])
                chosen_state_map[m] = 'top50'
            else:
                chosen_keys.append(m)
                chosen_state_map[m] = 'all'

        # All non-empty subsets except singletons (already in stage 1).
        for r in range(2, len(chosen_keys) + 1):
            for combo_idx in combinations(range(len(chosen_keys)), r):
                method_keys = [chosen_keys[i] for i in combo_idx]
                states = {bm: ('-' if i not in combo_idx
                                else chosen_state_map[bm])
                          for i, bm in enumerate(base_methods)}
                configs.append({
                    'label':       _label_for_methods(method_keys),
                    'method_keys': method_keys,
                    'n_methods':   r,
                    'states':      states,
                })

    elif level == 'full':
        # 3 states per method (or 2 without a Top50 variant); skip all-Off.
        per_method_options = []
        for m in base_methods:
            opts = [(None, '-'), (m, 'all')]
            top50 = _TOP50_VARIANT.get(m)
            if top50 is not None:
                opts.append((top50, 'top50'))
            per_method_options.append(opts)

        for combo in product(*per_method_options):
            method_keys = [c[0] for c in combo if c[0] is not None]
            if not method_keys:
                continue
            states = {bm: combo[i][1]
                      for i, bm in enumerate(base_methods)}
            configs.append({
                'label':       _label_for_methods(method_keys),
                'method_keys': method_keys,
                'n_methods':   len(method_keys),
                'states':      states,
            })

    elif level == 'exhaustive':
        # 2^M - 1 non-empty subsets across all method+Top50 variants.
        all_variants = []
        variant_to_method = {}
        for m in base_methods:
            all_variants.append(m)
            variant_to_method[m] = (m, 'all')
            top50 = _TOP50_VARIANT.get(m)
            if top50 is not None:
                all_variants.append(top50)
                variant_to_method[top50] = (m, 'top50')

        for r in range(1, len(all_variants) + 1):
            for combo in combinations(all_variants, r):
                # Same base method twice (e.g. dl + dl_top50) -> state 'mixed'.
                states = {bm: '-' for bm in base_methods}
                for variant in combo:
                    base_m, st = variant_to_method[variant]
                    if states[base_m] != '-':
                        states[base_m] = 'mixed'
                    else:
                        states[base_m] = st
                configs.append({
                    'label':       _label_for_methods(list(combo)),
                    'method_keys': list(combo),
                    'n_methods':   r,
                    'states':      states,
                })

    return configs


def characterise_replicators(loo_result: dict) -> dict:
    """Pool per-fold per_gene_df tables, compare feature distributions
    across replication outcomes.

    Answers "what distinguishes replicated consensus genes from
    wrong-direction failures?" -- candidates for a post-hoc filter.

    Returns
    -------
    dict with 'pooled_df', 'feature_summary', 'n_per_reason',
    'n_folds', 'n_genes_total'.
    """
    folds = loo_result.get('folds', [])
    n_folds = len(folds)

    rows = []
    for f in folds:
        pg = f.get('per_gene_df')
        if pg is None or len(pg) == 0:
            continue
        for _, row in pg.iterrows():
            rows.append({
                'gene':                 row.get('gene', ''),
                'fold_index':           f['fold_index'],
                'study_held_out':       f['study_held_out'],
                'consensus_log2FC':     row.get('consensus_log2FC',
                                                 float('nan')),
                'abs_consensus_log2FC': row.get(
                    'abs_consensus_log2FC', float('nan')),
                'consensus_fdr':        row.get('consensus_fdr',
                                                 float('nan')),
                'heterogeneity_i2':     row.get('heterogeneity_i2',
                                                 float('nan')),
                'n_studies':            row.get('n_studies', 0),
                'se_approx':            row.get('se_approx',
                                                 float('nan')),
                'study_log2FC':         row.get('study_log2FC',
                                                 float('nan')),
                'study_pval':           row.get('study_pval',
                                                 float('nan')),
                'study_padj':           row.get('study_padj',
                                                 float('nan')),
                'replicated':           row.get('replicated', False),
                'reason':               row.get('reason', ''),
            })

    # Empty -> zero-row frame with expected cols, so the per-reason
    # summary still produces a 4-row table of NaNs.
    _POOLED_COLS = (
        'gene', 'fold_index', 'study_held_out',
        'consensus_log2FC', 'abs_consensus_log2FC', 'consensus_fdr',
        'heterogeneity_i2', 'n_studies', 'se_approx',
        'study_log2FC', 'study_pval', 'study_padj',
        'replicated', 'reason',
    )
    if rows:
        pooled = pd.DataFrame(rows)
    else:
        pooled = pd.DataFrame({c: [] for c in _POOLED_COLS})

    n_per_reason = (dict(pooled['reason'].value_counts())
                    if len(pooled) else {})

    feature_cols = ('abs_consensus_log2FC', 'heterogeneity_i2',
                    'n_studies', 'se_approx', 'consensus_fdr')
    summary_rows = []
    for reason in ('replicated', 'wrong_direction',
                   'not_significant', 'missing'):
        sub = pooled[pooled['reason'] == reason]
        n = len(sub)
        row = {'reason': reason, 'n': n}
        for col in feature_cols:
            if n == 0:
                row[f'{col}_mean']   = float('nan')
                row[f'{col}_median'] = float('nan')
                continue
            vals = sub[col].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            if len(vals) == 0:
                row[f'{col}_mean']   = float('nan')
                row[f'{col}_median'] = float('nan')
            else:
                row[f'{col}_mean']   = float(np.mean(vals))
                row[f'{col}_median'] = float(np.median(vals))
        summary_rows.append(row)

    feature_summary = pd.DataFrame(summary_rows)

    return {
        'pooled_df':       pooled,
        'feature_summary': feature_summary,
        'n_per_reason':    n_per_reason,
        'n_folds':         n_folds,
        'n_genes_total':   len(pooled),
    }


def characterise_per_study(loo_result: dict) -> pd.DataFrame:
    """Per-held-out-study breakdown of replication outcomes.

    For each of the K LOO folds, report how many consensus genes were
    replicated / wrong-direction / not-significant / missing in that
    fold's held-out study. Tells you whether specific studies are
    disproportionately driving wrong-direction failures -- e.g. a
    study with the DESeq2 transcriptome-size flip issue would show a
    much higher wrong_direction rate than the others.

    Columns: 'fold_index', 'study_held_out', 'n_consensus',
    'n_replicated', 'n_wrong_direction', 'n_not_significant',
    'n_missing', plus percentage equivalents.
    """
    folds = loo_result.get('folds', [])
    rows = []
    for f in folds:
        pg = f.get('per_gene_df')
        if pg is None:
            pg = pd.DataFrame()
        n = len(pg) if len(pg) > 0 else 0
        vc = (pg['reason'].value_counts().to_dict()
              if 'reason' in pg.columns else {})
        n_rep = int(vc.get('replicated', 0))
        n_wd  = int(vc.get('wrong_direction', 0))
        n_ns  = int(vc.get('not_significant', 0))
        n_ms  = int(vc.get('missing', 0))
        rows.append({
            'fold_index':          f.get('fold_index', -1),
            'study_held_out':      f.get('study_held_out', ''),
            'n_consensus':         n,
            'n_replicated':        n_rep,
            'n_wrong_direction':   n_wd,
            'n_not_significant':   n_ns,
            'n_missing':           n_ms,
            'pct_replicated':      (100.0 * n_rep / n) if n else 0.0,
            'pct_wrong_direction': (100.0 * n_wd / n) if n else 0.0,
            'pct_not_significant': (100.0 * n_ns / n) if n else 0.0,
            'pct_missing':         (100.0 * n_ms / n) if n else 0.0,
        })
    return pd.DataFrame(rows)


def characterise_per_gene(loo_result: dict,
                          sort_by: str = 'n_wrong_direction'
                          ) -> pd.DataFrame:
    """Per-gene breakdown: collapse across LOO folds.

    Identifies systematic-failure genes (wrong_direction in every fold) --
    concrete filter targets.

    Columns: 'gene', 'n_folds', 'n_replicated', 'n_wrong_direction',
    'n_not_significant', 'n_missing', 'pct_wrong_direction',
    'abs_log2FC_max/mean', 'i2_max/mean'. Sorted by 'sort_by' desc.
    """
    folds = loo_result.get('folds', [])
    gene_stats: dict = {}
    gene_lfc_accum: dict = {}
    gene_i2_accum: dict = {}
    for f in folds:
        pg = f.get('per_gene_df')
        if pg is None or len(pg) == 0:
            continue
        for _, row in pg.iterrows():
            gene = str(row.get('gene', ''))
            if not gene:
                continue
            st = gene_stats.setdefault(gene, {
                'gene':               gene,
                'n_folds':            0,
                'n_replicated':       0,
                'n_wrong_direction':  0,
                'n_not_significant':  0,
                'n_missing':          0,
                'abs_log2FC_max':     0.0,
                'i2_max':             0.0,
            })
            gene_lfc_accum.setdefault(gene, [])
            gene_i2_accum.setdefault(gene, [])

            st['n_folds'] += 1
            reason = row.get('reason', '')
            if reason == 'replicated':
                st['n_replicated'] += 1
            elif reason == 'wrong_direction':
                st['n_wrong_direction'] += 1
            elif reason == 'not_significant':
                st['n_not_significant'] += 1
            elif reason == 'missing':
                st['n_missing'] += 1

            lfc = row.get('abs_consensus_log2FC', float('nan'))
            if pd.notna(lfc) and np.isfinite(lfc):
                gene_lfc_accum[gene].append(float(lfc))
                if lfc > st['abs_log2FC_max']:
                    st['abs_log2FC_max'] = float(lfc)
            i2 = row.get('heterogeneity_i2', float('nan'))
            if pd.notna(i2) and np.isfinite(i2):
                gene_i2_accum[gene].append(float(i2))
                if i2 > st['i2_max']:
                    st['i2_max'] = float(i2)

    for g, st in gene_stats.items():
        lfcs = gene_lfc_accum.get(g, [])
        i2s = gene_i2_accum.get(g, [])
        st['abs_log2FC_mean'] = (float(np.mean(lfcs))
                                  if lfcs else float('nan'))
        st['i2_mean'] = (float(np.mean(i2s)) if i2s else float('nan'))

    df = pd.DataFrame(list(gene_stats.values()))
    if len(df):
        df['pct_wrong_direction'] = (
            100.0 * df['n_wrong_direction'] / df['n_folds'])
        df['pct_replicated'] = (
            100.0 * df['n_replicated'] / df['n_folds'])
        if sort_by in df.columns:
            df = df.sort_values(
                [sort_by, 'pct_wrong_direction'],
                ascending=[False, False]).reset_index(drop=True)
    return df


def run_consensus_loo(datasets: List[pd.DataFrame],
                      labels: List[str],
                      params: dict,
                      replication_mode: str = 'strict',
                      pseudobulk_data: Optional[List[dict]] = None,
                      consensus_fdr_threshold: float = DEFAULT_FDR,
                      replication_threshold: float = DEFAULT_FDR,
                      fold_progress_cb: Optional[Callable[[int, int, str], None]] = None,
                      stage_progress_cb: Optional[Callable[[str], None]] = None,
                      top_n: int = 1000,
                      stash_fold_details: bool = False,
                      ) -> dict:
    """Held-out replication LOO across all studies.

    For each fold j: build consensus on K-1 studies, test whether each
    consensus gene independently replicates in the held-out study j.

    Parameters
    ----------
    datasets : list of DataFrame
        Per-study DE results (must have names, logfoldchanges, pvals, pvals_adj).
    labels : list of str
        Study display names, parallel to 'datasets'.
    params : dict
        GeneMAWorker-shaped: method_keys, calibration, min_studies,
        cc_n_perms (if cc_perm), de_method, seed.
    replication_mode : 'strict' | 'loose' | 'direction_only'
    pseudobulk_data : list of dict, optional
        Required if params['calibration'] == 'cc_perm'.
    consensus_fdr_threshold : FDR for "gene in the K-1 fold consensus".
    replication_threshold   : FDR (strict) / p (loose) for held-out replication.

    Returns
    -------
    dict
        Top-level keys: replication_mode, n_studies, fdr_threshold,
        rep_threshold, folds, avg/min/max_replication_pct.
        Per-fold keys: fold_index, study_held_out, consensus_size,
        n_replicated, n_not_replicated, n_missing, replication_pct, per_gene_df.

    Raises
    ------
    ValueError
        If K < 3.
    """
    if replication_mode not in REPLICATION_MODES:
        raise ValueError(
            f"replication_mode must be one of {REPLICATION_MODES!r}, "
            f"got {replication_mode!r}")

    K = len(datasets)
    if K < 3:
        raise ValueError(
            f"LOO validation requires at least 3 studies; got K={K}.")
    if len(labels) != K:
        raise ValueError(
            f"len(labels) ({len(labels)}) != len(datasets) ({K})")

    method_keys  = params['method_keys']
    calibration  = params.get('calibration', 'analytical')
    min_studies  = params.get('min_studies', 3)
    n_cc_perms   = params.get('cc_n_perms', 1000)
    de_method    = params.get('de_method', 'deseq2_fast')
    seed         = params.get('seed', 42)

    if calibration == 'cc_perm' and pseudobulk_data is None:
        raise ValueError(
            "calibration='cc_perm' requires pseudobulk_data")
    if pseudobulk_data is not None and len(pseudobulk_data) != K:
        raise ValueError(
            f"len(pseudobulk_data) ({len(pseudobulk_data)}) != K ({K})")

    folds = []
    for j in range(K):
        held_out = datasets[j]
        held_out_label = labels[j]

        if fold_progress_cb:
            fold_progress_cb(j + 1, K, held_out_label)

        fold_datasets = [datasets[i] for i in range(K) if i != j]
        fold_pseudobulk = (
            [pseudobulk_data[i] for i in range(K) if i != j]
            if pseudobulk_data is not None else None)

        fold_consensus_df = _build_fold_consensus(
            fold_datasets=fold_datasets,
            fold_pseudobulk=fold_pseudobulk,
            method_keys=method_keys,
            calibration=calibration,
            min_studies=min_studies,
            n_cc_perms=n_cc_perms,
            de_method=de_method,
            seed=seed,
            progress_cb=stage_progress_cb,
        )

        consensus_genes = fold_consensus_df[
            fold_consensus_df['fdr'] < consensus_fdr_threshold
        ].copy()
        consensus_size = len(consensus_genes)

        # Partition: "obvious" = FDR-sig in >= 1 contributing study,
        # "meta-discovered" = sig in 0. The meta-discovered slice is
        # the genuine meta-analysis payoff.
        n_sig_per_gene = {}
        for ds in fold_datasets:
            if 'pvals_adj' in ds.columns:
                sig = ds.loc[ds['pvals_adj'] < replication_threshold, 'names']
            elif 'fdr' in ds.columns:
                sig = ds.loc[ds['fdr'] < replication_threshold, 'names']
            else:
                continue
            for g in sig:
                n_sig_per_gene[g] = n_sig_per_gene.get(g, 0) + 1

        # default-arg binding silences B023 (loop-variable closure capture).
        consensus_genes['n_contributing_sig'] = consensus_genes['names'].map(
            lambda g, _n=n_sig_per_gene: _n.get(g, 0)).astype(int)

        obvious_genes = consensus_genes[consensus_genes['n_contributing_sig'] >= 1]
        meta_disc_genes = consensus_genes[consensus_genes['n_contributing_sig'] == 0]

        # All three modes computed (cheap once C_(-j) and held-out DE
        # are in hand). Primary mode populates 'per_gene_df'; the others
        # populate the *_strict / *_loose / *_direction summary fields.
        per_mode_results = {
            mode: _replication_test(
                consensus_genes, held_out,
                mode=mode,
                fdr_threshold=replication_threshold,
            )
            for mode in REPLICATION_MODES
        }

        # default-arg captures held_out / replication_threshold so each
        # closure owns its own bindings (silences B023).
        def _partition_stats(genes_subset,
                             _held=held_out,
                             _thr=replication_threshold):
            out = {'n': len(genes_subset)}
            if len(genes_subset) == 0:
                for mode in REPLICATION_MODES:
                    out[f'replication_pct_{mode}'] = 0.0
                    out[f'n_replicated_{mode}']    = 0
                out['n_missing'] = 0
                out['n_testable'] = 0
                return out
            mode_rates = {}
            n_miss_p = 0
            for mode in REPLICATION_MODES:
                r = _replication_test(
                    genes_subset, _held,
                    mode=mode, fdr_threshold=_thr)
                n_rep = int(r['replicated'].sum())
                if mode == 'strict':
                    n_miss_p = int((r['reason'] == 'missing').sum())
                mode_rates[mode] = n_rep
            n_testable_p = max(len(genes_subset) - n_miss_p, 0)
            for mode in REPLICATION_MODES:
                out[f'n_replicated_{mode}'] = mode_rates[mode]
                out[f'replication_pct_{mode}'] = (
                    100.0 * mode_rates[mode] / n_testable_p
                    if n_testable_p > 0 else 0.0)
            out['n_missing'] = n_miss_p
            out['n_testable'] = n_testable_p
            return out

        partition_obvious = _partition_stats(obvious_genes)
        partition_meta    = _partition_stats(meta_disc_genes)

        # Breakdown by n_contributing_sig: 0 = meta-discovery,
        # 1 = single-signal amplification, 2+ = multi-study confirmation.
        partition_by_n = {}
        k_fold = K - 1
        for n_sig in range(k_fold + 1):
            subset = consensus_genes[consensus_genes['n_contributing_sig'] == n_sig]
            partition_by_n[n_sig] = _partition_stats(subset)
        per_gene_df = per_mode_results[replication_mode]

        # Missing count is mode-independent (presence in held-out study)
        n_missing = int(
            (per_mode_results['strict']['reason'] == 'missing').sum())

        # Testable denominator excludes missing genes -- a gene not
        # measured in the held-out study can't be tested at all, so
        # counting it as a failure misstates each mode's meaning.
        n_testable = max(consensus_size - n_missing, 0)

        per_mode_counts = {}
        per_mode_pcts = {}
        for mode in REPLICATION_MODES:
            n_rep = int(per_mode_results[mode]['replicated'].sum())
            per_mode_counts[mode] = n_rep
            per_mode_pcts[mode] = (
                100.0 * n_rep / n_testable
                if n_testable > 0 else 0.0)

        n_replicated     = per_mode_counts[replication_mode]
        n_not_replicated = consensus_size - n_replicated
        replication_pct  = per_mode_pcts[replication_mode]

        # Per-method standalone call sets (passing THAT method's own
        # FDR), partitioned into obvious vs meta-discovered. Surfaces
        # per-method meta-discoveries the intersection consensus hides.
        solo_per_method = {}
        for mkey in method_keys:
            fcol = f'fdr_{mkey}'
            if fcol not in fold_consensus_df.columns:
                continue
            solo = fold_consensus_df[
                fold_consensus_df[fcol] < consensus_fdr_threshold
            ].copy()
            n_total = len(solo)
            if n_total == 0:
                solo_per_method[mkey] = {
                    'n_total': 0,
                    'n_obvious': 0,
                    'n_meta_discovered': 0,
                    'partition_obvious': _partition_stats(solo),
                    'partition_meta_discovered': _partition_stats(solo),
                }
                continue
            solo['n_contributing_sig'] = solo['names'].map(
                lambda g, _n=n_sig_per_gene: _n.get(g, 0)).astype(int)
            solo_obvious = solo[solo['n_contributing_sig'] >= 1]
            solo_meta    = solo[solo['n_contributing_sig'] == 0]
            solo_per_method[mkey] = {
                'n_total':             n_total,
                'n_obvious':           len(solo_obvious),
                'n_meta_discovered':   len(solo_meta),
                'partition_obvious':   _partition_stats(solo_obvious),
                'partition_meta_discovered': _partition_stats(solo_meta),
            }

        # Top-N per method (fair cross-method comparison): each method's
        # top-N genes by its own p-value, replication tested on that
        # fixed-size slice -- removes the gene-set-size confound.
        top_n_per_method = {}
        for mkey in method_keys:
            pcol = f'pvals_{mkey}'
            if pcol not in fold_consensus_df.columns:
                continue
            ranked = fold_consensus_df.dropna(subset=[pcol]).sort_values(pcol)
            topn_genes = ranked.head(top_n).copy()
            topn_size = len(topn_genes)
            topn_res = _replication_test(
                topn_genes, held_out,
                mode=replication_mode,
                fdr_threshold=replication_threshold,
            )
            n_missing_m = int((topn_res['reason'] == 'missing').sum())
            n_testable_m = max(topn_size - n_missing_m, 0)
            n_rep_m = int(topn_res['replicated'].sum())
            n_wrong_m = int((topn_res['reason'] == 'wrong_direction').sum())
            n_notsig_m = int((topn_res['reason'] == 'not_significant').sum())
            top_n_per_method[mkey] = {
                'n':              topn_size,
                'n_testable':     n_testable_m,
                'n_replicated':   n_rep_m,
                'n_wrong_dir':    n_wrong_m,
                'n_not_sig':      n_notsig_m,
                'n_missing':      n_missing_m,
                'replication_pct': (100.0 * n_rep_m / n_testable_m
                                    if n_testable_m > 0 else 0.0),
                'wrong_dir_pct':   (100.0 * n_wrong_m / n_testable_m
                                    if n_testable_m > 0 else 0.0),
                'not_sig_pct':     (100.0 * n_notsig_m / n_testable_m
                                    if n_testable_m > 0 else 0.0),
            }

        # AUC via two module scores: sum-of-signed-z (primary) and
        # UCell rank (SumRank-paper parity). Requires pseudobulk;
        # silently nan for analytical-only runs.
        held_out_pb = (pseudobulk_data[j] if pseudobulk_data is not None
                       else None)
        auc_consensus = _held_out_auc(held_out_pb, consensus_genes)

        auc_per_method = {}
        for mkey in method_keys:
            pcol = f'pvals_{mkey}'
            if pcol not in fold_consensus_df.columns:
                continue
            ranked = fold_consensus_df.dropna(subset=[pcol]).sort_values(pcol)
            topn_genes = ranked.head(top_n)
            auc_per_method[mkey] = _held_out_auc(held_out_pb, topn_genes)

        fold_entry = {
            'fold_index':       j,
            'study_held_out':   held_out_label,
            'n_studies_in_fold': K - 1,
            'consensus_size':   consensus_size,
            'n_testable':       n_testable,
            'n_replicated':     n_replicated,
            'n_not_replicated': n_not_replicated,
            'n_missing':        n_missing,
            'replication_pct':  replication_pct,
            'n_replicated_strict':    per_mode_counts['strict'],
            'n_replicated_loose':     per_mode_counts['loose'],
            'n_replicated_direction': per_mode_counts['direction_only'],
            'replication_pct_strict':    per_mode_pcts['strict'],
            'replication_pct_loose':     per_mode_pcts['loose'],
            'replication_pct_direction': per_mode_pcts['direction_only'],
            'per_gene_df':      per_gene_df,
            'top_n_per_method': top_n_per_method,
            'top_n':            top_n,
            'partition_obvious':         partition_obvious,
            'partition_meta_discovered': partition_meta,
            'partition_by_n_contributing': partition_by_n,
            'auc_consensus':    auc_consensus,
            'auc_per_method':   auc_per_method,
            'solo_per_method':  solo_per_method,
        }

        # Stash the minimal held-out data + per-method gene sets so the
        # Consensus Evaluation page can compute direction replication
        # for arbitrary method combinations without re-running LOO.
        if stash_fold_details:
            called_by_method = {}
            for mkey in method_keys:
                fcol = f'fdr_{mkey}'
                if fcol not in fold_consensus_df.columns:
                    continue
                called_by_method[mkey] = [
                    str(g) for g in fold_consensus_df.loc[
                        fold_consensus_df[fcol] < consensus_fdr_threshold,
                        'names']
                ]
            keep_cons_cols = ['names', 'logfoldchanges']
            all_genes_with_lfc = fold_consensus_df[
                [c for c in keep_cons_cols
                 if c in fold_consensus_df.columns]].copy()
            keep_ho_cols = ['names', 'logfoldchanges', 'pvals', 'pvals_adj']
            held_out_min = held_out[
                [c for c in keep_ho_cols
                 if c in held_out.columns]].copy()
            fold_entry['_stash'] = {
                'called_genes_by_method': called_by_method,
                'all_genes_with_lfc':     all_genes_with_lfc,
                'held_out_min':           held_out_min,
                'n_sig_per_gene':         dict(n_sig_per_gene),
            }

        folds.append(fold_entry)

    valid = [f for f in folds if f['consensus_size'] > 0]

    def _agg(field):
        vals = [f[field] for f in valid]
        if not vals:
            return 0.0, 0.0, 0.0
        return (float(np.mean(vals)),
                float(np.min(vals)),
                float(np.max(vals)))

    avg_strict, min_strict, max_strict = _agg('replication_pct_strict')
    avg_loose,  min_loose,  max_loose  = _agg('replication_pct_loose')
    avg_dir,    min_dir,    max_dir    = _agg('replication_pct_direction')
    avg_pct,    min_pct,    max_pct    = _agg('replication_pct')

    def _agg_auc(key):
        vals = [f.get('auc_consensus', {}).get(key, np.nan) for f in valid]
        vals = [v for v in vals if np.isfinite(v)]
        if not vals:
            return np.nan, np.nan, np.nan
        return (float(np.mean(vals)),
                float(np.min(vals)),
                float(np.max(vals)))

    avg_auc_z, min_auc_z, max_auc_z = _agg_auc('auc_sum_of_z')
    avg_auc_u, min_auc_u, max_auc_u = _agg_auc('auc_ucell')

    # Per-method top-N: fixed gene-set size removes the consensus-size confound.
    top_n_summary = {}
    for mkey in method_keys:
        rep_vals, wrong_vals, notsig_vals, size_vals = [], [], [], []
        for f in folds:
            tn = f.get('top_n_per_method', {}).get(mkey)
            if tn is None or tn.get('n_testable', 0) == 0:
                continue
            rep_vals.append(tn['replication_pct'])
            wrong_vals.append(tn['wrong_dir_pct'])
            notsig_vals.append(tn['not_sig_pct'])
            size_vals.append(tn['n'])
        if not rep_vals:
            continue
        # AUC computable only when pseudobulk is present with valid group sizes.
        auc_z_vals, auc_u_vals = [], []
        for f in folds:
            am = f.get('auc_per_method', {}).get(mkey, {})
            z = am.get('auc_sum_of_z', np.nan)
            u = am.get('auc_ucell', np.nan)
            if np.isfinite(z):
                auc_z_vals.append(z)
            if np.isfinite(u):
                auc_u_vals.append(u)

        top_n_summary[mkey] = {
            'mean_replication_pct': float(np.mean(rep_vals)),
            'min_replication_pct':  float(np.min(rep_vals)),
            'max_replication_pct':  float(np.max(rep_vals)),
            'mean_wrong_dir_pct':   float(np.mean(wrong_vals)),
            'mean_not_sig_pct':     float(np.mean(notsig_vals)),
            'mean_n':               float(np.mean(size_vals)),
            'mean_auc_sum_of_z':    (float(np.mean(auc_z_vals))
                                     if auc_z_vals else np.nan),
            'mean_auc_ucell':       (float(np.mean(auc_u_vals))
                                     if auc_u_vals else np.nan),
        }

    # Per-method solo partition (Comparison-mode headline).
    solo_summary = {}
    for mkey in method_keys:
        total_vals, obv_vals, meta_vals = [], [], []
        meta_dir_vals, meta_strict_vals = [], []
        obv_dir_vals, obv_strict_vals = [], []
        for f in folds:
            sm = f.get('solo_per_method', {}).get(mkey)
            if sm is None:
                continue
            total_vals.append(sm.get('n_total', 0))
            obv_vals.append(sm.get('n_obvious', 0))
            meta_vals.append(sm.get('n_meta_discovered', 0))
            obv_p = sm.get('partition_obvious', {})
            meta_p = sm.get('partition_meta_discovered', {})
            if obv_p.get('n_testable', 0) > 0:
                obv_dir_vals.append(obv_p.get('replication_pct_direction_only', 0.0))
                obv_strict_vals.append(obv_p.get('replication_pct_strict', 0.0))
            if meta_p.get('n_testable', 0) > 0:
                meta_dir_vals.append(meta_p.get('replication_pct_direction_only', 0.0))
                meta_strict_vals.append(meta_p.get('replication_pct_strict', 0.0))
        if not total_vals:
            continue
        solo_summary[mkey] = {
            'mean_n_total':           float(np.mean(total_vals)),
            'mean_n_obvious':         float(np.mean(obv_vals)),
            'mean_n_meta_discovered': float(np.mean(meta_vals)),
            'mean_obvious_direction_pct': (
                float(np.mean(obv_dir_vals)) if obv_dir_vals else 0.0),
            'mean_obvious_strict_pct': (
                float(np.mean(obv_strict_vals)) if obv_strict_vals else 0.0),
            'mean_meta_direction_pct': (
                float(np.mean(meta_dir_vals)) if meta_dir_vals else 0.0),
            'mean_meta_strict_pct': (
                float(np.mean(meta_strict_vals)) if meta_strict_vals else 0.0),
        }

    def _agg_partition(key):
        rows = [f.get(key, {}) for f in valid]
        if not rows:
            return {}
        n_mean = float(np.mean([r.get('n', 0) for r in rows]))
        out = {'mean_n': n_mean}
        for mode in REPLICATION_MODES:
            pct_vals = [r.get(f'replication_pct_{mode}', 0.0) for r in rows
                        if r.get('n_testable', 0) > 0]
            if pct_vals:
                out[f'mean_replication_pct_{mode}'] = float(np.mean(pct_vals))
            else:
                out[f'mean_replication_pct_{mode}'] = 0.0
        return out

    partition_obvious_summary = _agg_partition('partition_obvious')
    partition_meta_summary    = _agg_partition('partition_meta_discovered')

    partition_by_n_summary = {}
    max_n = 0
    for f in valid:
        pbn = f.get('partition_by_n_contributing', {})
        if pbn:
            max_n = max(max_n, max(pbn.keys()))
    for n_sig in range(max_n + 1):
        rows = [f.get('partition_by_n_contributing', {}).get(n_sig, {})
                for f in valid]
        if not rows:
            continue
        n_mean = float(np.mean([r.get('n', 0) for r in rows]))
        out = {'mean_n': n_mean}
        for mode in REPLICATION_MODES:
            pct_vals = [r.get(f'replication_pct_{mode}', 0.0) for r in rows
                        if r.get('n_testable', 0) > 0]
            out[f'mean_replication_pct_{mode}'] = (
                float(np.mean(pct_vals)) if pct_vals else 0.0)
        partition_by_n_summary[n_sig] = out

    return {
        'replication_mode': replication_mode,
        'n_studies':        K,
        'fdr_threshold':    consensus_fdr_threshold,
        'rep_threshold':    replication_threshold,
        'folds':            folds,
        'avg_replication_pct': avg_pct,
        'min_replication_pct': min_pct,
        'max_replication_pct': max_pct,
        'avg_replication_pct_strict':    avg_strict,
        'min_replication_pct_strict':    min_strict,
        'max_replication_pct_strict':    max_strict,
        'avg_replication_pct_loose':     avg_loose,
        'min_replication_pct_loose':     min_loose,
        'max_replication_pct_loose':     max_loose,
        'avg_replication_pct_direction': avg_dir,
        'min_replication_pct_direction': min_dir,
        'max_replication_pct_direction': max_dir,
        # Held-out case/control AUC; nan if analytical-only LOO.
        'avg_auc_sum_of_z': avg_auc_z,
        'min_auc_sum_of_z': min_auc_z,
        'max_auc_sum_of_z': max_auc_z,
        'avg_auc_ucell':    avg_auc_u,
        'min_auc_ucell':    min_auc_u,
        'max_auc_ucell':    max_auc_u,
        'top_n':                  top_n,
        'top_n_per_method':       top_n_summary,
        # Honest per-method counts (each method's own FDR<0.05 call set,
        # not biased by intersection vetoes).
        'solo_per_method':        solo_summary,
        'partition_obvious':         partition_obvious_summary,
        'partition_meta_discovered': partition_meta_summary,
        'partition_by_n_contributing': partition_by_n_summary,
    }
