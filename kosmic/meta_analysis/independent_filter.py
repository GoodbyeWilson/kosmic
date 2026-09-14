# Cross-study Bourgon-style filter on the meta-analysis input universe.
#
# Drops genes whose aggregate expression across studies is below a user
# threshold *before* p-values are pooled. BH FDR then runs on the
# surviving universe; smaller universe -> threshold relaxes -> more
# power on the genes that matter, without inflating type I error
# (Bourgon's criterion: filter is on the marginal expression
# distribution, independent of the test statistic under the null).
from __future__ import annotations

import numpy as np
import pandas as pd

from kosmic import MIN_STUDIES

_REQUIRED_COLUMNS = ("names", "pct_disease", "pct_control")
_FILTER_STATS = ("mean_max_pct", "min_studies_above_pct")


def apply_filter(
    datasets: list[dict],
    threshold: float,
    stat: str = "mean_max_pct",
    *,
    min_studies: int = MIN_STUDIES,
) -> set[str]:
    """Return gene names passing a cross-study expression filter.

    Parameters
    ----------
    datasets : list of dict
        Each entry must have a ``df`` key holding a DataFrame with
        columns ``names``, ``pct_disease``, ``pct_control``. Other
        keys are ignored.
    threshold : float
        Per-gene aggregate must be >= ``threshold``. Pct columns are
        treated as fractions (0-1), matching the DE workspace's
        upstream ``annotate_pct_expressing``.
    stat : str
        How to aggregate across studies:

        * ``"mean_max_pct"`` -- mean (across studies in which the gene
          appears) of ``max(pct_disease, pct_control)``. Default.
        * ``"min_studies_above_pct"`` -- count of studies in which
          ``max(pct_disease, pct_control) >= threshold``. Gene passes
          if the count is at least ``min_studies``.

    min_studies : int
        Used by ``min_studies_above_pct``. Defaults to the package
        ``MIN_STUDIES`` constant.

    Returns
    -------
    set[str]
        Gene names passing the filter. Empty if no inputs supply the
        required columns.
    """
    if stat not in _FILTER_STATS:
        raise ValueError(
            f"Unknown stat {stat!r}. Expected one of {_FILTER_STATS}."
        )
    if threshold < 0:
        raise ValueError(f"threshold must be >= 0, got {threshold!r}.")

    # Collect per-gene per-study max(pct_disease, pct_control).
    per_gene_pcts: dict[str, list[float]] = {}
    for entry in datasets:
        df = entry.get("df")
        if df is None or df.empty:
            continue
        if any(col not in df.columns for col in _REQUIRED_COLUMNS):
            # Studies without pct columns can't contribute; skip rather
            # than raise, so a partial input still yields a usable
            # filter from whatever data IS annotated.
            continue
        names = df["names"].astype(str).values
        d = pd.to_numeric(df["pct_disease"], errors="coerce").fillna(0).values
        c = pd.to_numeric(df["pct_control"], errors="coerce").fillna(0).values
        max_pct = np.maximum(d, c)
        for gene, value in zip(names, max_pct):
            per_gene_pcts.setdefault(gene, []).append(float(value))

    if not per_gene_pcts:
        return set()

    if stat == "mean_max_pct":
        return {
            gene for gene, vals in per_gene_pcts.items()
            if (sum(vals) / len(vals)) >= threshold
        }

    # min_studies_above_pct
    return {
        gene for gene, vals in per_gene_pcts.items()
        if sum(1 for v in vals if v >= threshold) >= min_studies
    }


def filter_datasets_to_universe(
    datasets: list[dict],
    gene_universe: set[str] | None,
) -> list[dict]:
    """Return ``datasets`` with each df restricted to ``gene_universe``.

    Pass ``gene_universe=None`` to return the input unchanged. Each
    output dict mirrors the input with the same keys; only the ``df``
    is replaced by a filtered copy. Studies whose ``df`` lacks a
    ``names`` column are passed through (the filter is a no-op for
    them; pooling will fail later with a clearer error if it can't
    handle them).
    """
    if gene_universe is None:
        return datasets

    out: list[dict] = []
    for entry in datasets:
        df = entry.get("df")
        if df is None or "names" not in df.columns:
            out.append(entry)
            continue
        mask = df["names"].astype(str).isin(gene_universe)
        kept = df[mask].copy()
        out.append({**entry, "df": kept})
    return out
