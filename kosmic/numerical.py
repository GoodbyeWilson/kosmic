"""Numerical guards over statsmodels / numpy."""
from __future__ import annotations

import numpy as np

__all__ = ["bh_fdr", "neg_log10"]


def bh_fdr(pvals, fill: float = 1.0) -> np.ndarray:
    """BH-adjusted p-values. Invalid entries (NaN/inf/<=0/>1) get 'fill'
    and are excluded from the BH denominator. Needed because
    statsmodels.multipletests does not mask invalid p-values and silently
    deflates every q-value when N is inflated."""
    from statsmodels.stats.multitest import multipletests
    arr = np.asarray(pvals, dtype=float)
    valid = np.isfinite(arr) & (arr > 0) & (arr <= 1)
    out = np.full(arr.shape, fill, dtype=float)
    if valid.sum() > 0:
        _, out[valid], _, _ = multipletests(arr[valid], method="fdr_bh")
    return out


def neg_log10(pvals, floor: float = 1e-300) -> np.ndarray:
    """-log10(p) with clipping so p=0 maps to a finite ceiling. Needed
    because -np.log10(0) returns inf, which crashes pyqtgraph and
    matplotlib axis autoscaling."""
    arr = np.asarray(pvals, dtype=float)
    return -np.log10(np.clip(arr, floor, 1.0))
