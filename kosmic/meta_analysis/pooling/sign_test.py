# Two-sided sign test on per-study effect directions.
#
# Under H0 of no effect, 'P(positive) = 0.5' in each study; the test
# statistic is the larger of 'n_pos' and 'n_neg', with a two-sided
# binomial tail.
#
# This is a non-parametric direction-only meta-analysis pool: it ignores
# effect magnitude and per-study standard errors, sensitive only to whether
# study log-fold-changes agree on sign. Useful as a robustness check beside
# DL / Stouffer / Fisher when SE estimates are noisy.
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import binom


def pool_sign_test(base_df: pd.DataFrame) -> pd.Series:
    """Compute per-row sign-test p-values from a 'study_effects' column.

    Parameters
    ----------
    base_df
        DataFrame whose rows are pathways/genes and which carries a
        'study_effects' column. Each entry is a list of dicts; each
        dict has at least 'logfc' (float). Studies with non-finite or
        zero 'logfc' are dropped.

    Returns
    -------
    pd.Series
        Index = 'base_df.index'; values = two-sided binomial p-values
        clipped to '[1e-300, 1.0]'. Rows with fewer than 2 valid studies
        return 'NaN'.
    """
    out = pd.Series(index=base_df.index, dtype=float)
    for name, row in base_df.iterrows():
        effects_data = row.get('study_effects', [])
        if not isinstance(effects_data, list) or len(effects_data) < 2:
            out[name] = np.nan
            continue
        signs = [
            np.sign(d.get('logfc', 0.0))
            for d in effects_data
            if np.isfinite(d.get('logfc', np.nan))
            and d.get('logfc', 0.0) != 0
        ]
        if len(signs) < 2:
            out[name] = np.nan
            continue
        n = len(signs)
        n_pos = sum(1 for s in signs if s > 0)
        k = max(n_pos, n - n_pos)
        # Two-sided binomial test: P(X >= k or X <= n-k) under p = 0.5.
        p = 2 * binom.sf(k - 1, n, 0.5)
        out[name] = float(np.clip(min(p, 1.0), 1e-300, 1.0))
    return out
