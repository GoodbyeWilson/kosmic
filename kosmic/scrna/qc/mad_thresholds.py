# MAD-based QC thresholds.
#
# Median Absolute Deviation (MAD) thresholds for the three standard scRNA QC
# metrics (genes/cell, total counts, MT %), computed in log1p space because
# all three are right-skewed count-derived quantities.
#
# The function returns thresholds in the original (un-logged) scale ready to
# feed into a filter step or display in QC spinboxes.
from __future__ import annotations

from typing import TYPE_CHECKING, Mapping

import numpy as np
from scipy.stats import median_abs_deviation

if TYPE_CHECKING:
    import pandas as pd


def compute_mad_thresholds(obs: Mapping[str, np.ndarray] | "pd.DataFrame",
                           nmads: float = 3.0) -> dict[str, dict[str, float]]:
    """Compute symmetric MAD thresholds for the three core QC metrics.

    Each metric is log1p-transformed before MAD is computed (right-skewed
    count data is poorly summarised by MAD on the raw scale). The returned
    bounds are in the original scale via 'expm1'. Lower bounds are
    floored at 0.

    Parameters
    ----------
    obs
        Anything that supports 'obs[col].values' -- typically
        'adata.obs'. Required columns: 'n_genes_by_counts',
        'total_counts', 'pct_counts_mt'. Missing columns yield
        'None' for that metric's entry.
    nmads
        Number of MADs from the median to set the bound. Default 3.

    Returns
    -------
    dict
        '{metric: {'min': float, 'max': float}}' for each available
        metric. 'pct_counts_mt' returns only 'max' (no biologically
        meaningful lower bound).
    """
    out: dict[str, dict[str, float]] = {}

    def _have(col: str) -> bool:
        try:
            return col in obs.columns  # type: ignore[attr-defined]
        except Exception:
            return col in obs

    def _values(col: str):
        try:
            return obs[col].values  # pandas-like
        except Exception:
            return np.asarray(obs[col])

    if _have('n_genes_by_counts'):
        log_v = np.log1p(_values('n_genes_by_counts'))
        med = float(np.median(log_v))
        mad = float(median_abs_deviation(log_v))
        out['n_genes_by_counts'] = {
            'min': float(max(0.0, np.expm1(med - nmads * mad))),
            'max': float(np.expm1(med + nmads * mad)),
        }

    if _have('total_counts'):
        log_v = np.log1p(_values('total_counts'))
        med = float(np.median(log_v))
        mad = float(median_abs_deviation(log_v))
        out['total_counts'] = {
            'min': float(max(0.0, np.expm1(med - nmads * mad))),
            'max': float(np.expm1(med + nmads * mad)),
        }

    if _have('pct_counts_mt'):
        log_v = np.log1p(_values('pct_counts_mt'))
        med = float(np.median(log_v))
        mad = float(median_abs_deviation(log_v))
        out['pct_counts_mt'] = {
            'max': float(np.expm1(med + nmads * mad)),
        }

    return out
