# Per-study direction of effect for pooled meta-analysis tables.
#
# Several pooling methods (Fisher, SumRank, gwOP, and to a lesser extent
# Stouffer) combine per-study evidence without regard to the sign of the
# effect, so a gene that is significantly up-regulated in one study and
# significantly down-regulated in another can still be called
# significant. 'annotate_direction' adds three columns to a pooled table
# so that such genes can be identified:
#
#   n_up                number of studies with study FDR < threshold and
#                       a positive log2 fold change
#   n_down              the same for a negative log2 fold change
#   direction_conflict  True when both n_up and n_down are non-zero
#
# The counts are computed from the per-study input tables, not from the
# 'study_effects' column of the pooled table, because the p-value stored
# there depends on the pooling method (SumRank records the raw p-value,
# the effect-size methods the adjusted one). A study is counted only if
# its table has a 'pvals_adj' column; 'studies_without_adjusted_p' lists
# the studies that do not, so that the caller can report them.
from __future__ import annotations

from typing import List, Optional, Sequence

import numpy as np
import pandas as pd

from kosmic import DEFAULT_FDR

DIRECTION_COLUMNS = ('n_up', 'n_down', 'direction_conflict')


def _study_name(df: pd.DataFrame, index: int) -> str:
    """Name a study table the way 'pooling_core' does."""
    if 'dataset' in df.columns and len(df) > 0:
        return str(df['dataset'].iloc[0])
    return f'Study_{index + 1}'


def studies_without_adjusted_p(datasets: Sequence[pd.DataFrame]) -> List[str]:
    """Return the names of the study tables that have no 'pvals_adj' column.

    These studies are left out of the 'n_up' / 'n_down' counts made by
    'annotate_direction'.

    Parameters
    ----------
    datasets : sequence of DataFrame
        The per-study tables passed to a pooling function.

    Returns
    -------
    list of str
        Study names, from the 'dataset' column where present and
        'Study_<i>' otherwise.
    """
    return [_study_name(df, i) for i, df in enumerate(datasets)
            if 'pvals_adj' not in df.columns]


def unadjusted_warning(datasets: Sequence[pd.DataFrame]) -> Optional[str]:
    """Describe the studies left out of the direction counts, or None.

    Parameters
    ----------
    datasets : sequence of DataFrame
        The per-study tables passed to a pooling function.

    Returns
    -------
    str or None
        A message for the output panel, or None when every study has
        'pvals_adj'.
    """
    missing = studies_without_adjusted_p(datasets)
    if not missing:
        return None
    if len(missing) == len(datasets):
        return ("No study has adjusted p-values (pvals_adj), so the "
                "Up / Down / Conflict columns are not computed.")
    return (f"No adjusted p-values (pvals_adj) in {', '.join(missing)}. "
            "These studies are not counted in the Up / Down / Conflict "
            "columns.")



def add_wald_pvalues(df: pd.DataFrame) -> pd.DataFrame:
    """Add Wald-test 'pvals' and within-study BH 'pvals_adj' to a summary.

    For tables that carry an effect and its standard error but no test
    result, such as the DESeq2 + VIF pathway summaries. The statistic is
    z = logfoldchanges / se against the standard normal, two-sided; the
    adjustment is Benjamini-Hochberg over the rows of this one study.

    Parameters
    ----------
    df : DataFrame
        One study's table with 'logfoldchanges' and 'se'.

    Returns
    -------
    DataFrame
        A copy with 'pvals' and 'pvals_adj'. Rows with a missing or
        non-positive SE get p = 1 and are left out of the BH denominator.
    """
    from scipy.stats import norm

    from kosmic.numerical import bh_fdr

    out = df.copy()
    lfc = pd.to_numeric(out['logfoldchanges'], errors='coerce').to_numpy(float)
    se = pd.to_numeric(out['se'], errors='coerce').to_numpy(float)
    valid = np.isfinite(lfc) & np.isfinite(se) & (se > 0)
    z = np.divide(lfc, se, out=np.zeros_like(lfc), where=valid)
    pvals = np.where(valid, 2 * norm.sf(np.abs(z)), np.nan)
    out['pvals'] = np.where(valid, np.clip(pvals, 1e-300, 1.0), 1.0)
    adj = bh_fdr(np.where(valid, out['pvals'], np.nan))
    out['pvals_adj'] = adj
    return out

def annotate_direction(meta_df: pd.DataFrame,
                       datasets: Sequence[pd.DataFrame],
                       fdr: float = DEFAULT_FDR) -> pd.DataFrame:
    """Add 'n_up', 'n_down' and 'direction_conflict' to a pooled table.

    Parameters
    ----------
    meta_df : DataFrame
        Pooled result with a 'names' column (genes or pathways).
    datasets : sequence of DataFrame
        The per-study tables that were pooled. Each needs 'names' and
        'logfoldchanges'; a table without 'pvals_adj' is not counted.
    fdr : float
        Study-level FDR threshold below which an effect is counted as
        significant.

    Returns
    -------
    DataFrame
        A copy of 'meta_df' with the three columns added. Rows whose name
        does not occur in any counted study get zero counts. When no study
        can be counted (none has 'pvals_adj'), the columns are not added,
        so that the absence of counts is not read as the absence of
        conflicts.
    """
    out = meta_df.copy()
    counted = [df for df in datasets
               if 'pvals_adj' in df.columns and 'logfoldchanges' in df.columns]
    if 'names' not in out.columns or not counted:
        return out

    names = pd.Index(out['names'])
    n_up = np.zeros(len(out), dtype=int)
    n_down = np.zeros(len(out), dtype=int)

    for df in counted:
        # Pooling assigns duplicate names last-wins, so do the same here.
        study = df.drop_duplicates('names', keep='last').set_index('names')
        padj = pd.to_numeric(study['pvals_adj'], errors='coerce')
        lfc = pd.to_numeric(study['logfoldchanges'], errors='coerce')
        sig = padj < fdr
        up = (sig & (lfc > 0)).reindex(names, fill_value=False).to_numpy()
        down = (sig & (lfc < 0)).reindex(names, fill_value=False).to_numpy()
        n_up += up.astype(int)
        n_down += down.astype(int)

    out['n_up'] = n_up
    out['n_down'] = n_down
    out['direction_conflict'] = (n_up > 0) & (n_down > 0)
    return out


def count_significant_conflicts(meta_df: pd.DataFrame,
                                fdr_col: str = 'fdr',
                                fdr: float = DEFAULT_FDR) -> tuple:
    """Count significant rows and, among them, direction conflicts.

    Parameters
    ----------
    meta_df : DataFrame
        Pooled table carrying 'fdr_col' and, optionally, the columns
        added by 'annotate_direction'.
    fdr_col : str
        Column holding the pooled FDR.
    fdr : float
        Pooled FDR threshold.

    Returns
    -------
    tuple of (int, int or None)
        The number of rows with 'fdr_col' below 'fdr', and how many of
        those have 'direction_conflict' set. The second value is None
        when the table has no 'direction_conflict' column: a result saved
        before the column existed, or one where no study had adjusted
        p-values.
    """
    if fdr_col not in meta_df.columns:
        return 0, None
    sig = pd.to_numeric(meta_df[fdr_col], errors='coerce') < fdr
    if 'direction_conflict' not in meta_df.columns:
        return int(sig.sum()), None
    conflict = meta_df['direction_conflict'].fillna(False).astype(bool)
    return int(sig.sum()), int((sig & conflict).sum())
