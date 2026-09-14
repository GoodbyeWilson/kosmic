# Annotate consensus gene lists with Olink proteomics panel membership.
#
# Panel definitions are JSON in 'kosmic/reference/olink_panels/'
# (Olink_CVD_panels.json, Olink_Explore_panels.json); each value is a
# list of HGNC symbols.

from __future__ import annotations

import json
from functools import lru_cache
from importlib.resources import files
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd
from kosmic import DEFAULT_FDR

_DATA_DIR = files("kosmic.reference.olink_panels")

DEFAULT_PANELS = ("olink_cvd_iii", "olink_explore_3072")

PANEL_LABELS = {
    "olink_cvd_ii":              "Olink CVD II (92 proteins)",
    "olink_cvd_iii":             "Olink CVD III (91 proteins)",
    "olink_cvd_combined":        "Olink CVD II + III (102 proteins)",
    "olink_cardiometabolic":     "Olink Explore Cardiometabolic",
    "olink_cardiometabolic_ii":  "Olink Explore Cardiometabolic II",
    "olink_inflammation":        "Olink Explore Inflammation",
    "olink_inflammation_ii":     "Olink Explore Inflammation II",
    "olink_neurology":           "Olink Explore Neurology",
    "olink_neurology_ii":        "Olink Explore Neurology II",
    "olink_oncology":            "Olink Explore Oncology",
    "olink_oncology_ii":         "Olink Explore Oncology II",
    "olink_explore_3072":        "Olink Explore 3072 (all sub-panels, ~2925 unique)",
}


@lru_cache(maxsize=4)
def _load_panels(data_dir: Optional[Path] = None) -> dict[str, frozenset[str]]:
    base = Path(data_dir) if data_dir is not None else _DATA_DIR
    panels: dict[str, frozenset[str]] = {}
    for path, rename_union in (
        (base / "Olink_CVD_panels.json", None),
        (base / "Olink_Explore_panels.json", "explore_3072_combined"),
    ):
        if not path.exists():
            continue
        with open(path, encoding="utf-8") as f:
            raw = json.load(f)
        for k, genes in raw.items():
            key = "olink_explore_3072" if k == rename_union else f"olink_{k}"
            panels[key] = frozenset(str(g).upper() for g in genes if g)
    return panels


def _harmonise(genes: Iterable[str]) -> list[str]:
    """Uppercase + HGNC-alias-resolve, matching the consensus pipeline."""
    out = [str(g).upper() for g in genes]
    try:
        from kosmic.scrna.inspect.gene_names import load_hgnc_lookup
        lookup, _ = load_hgnc_lookup()
        out = [lookup.get(g, g) for g in out]
    except Exception:
        pass
    return out


def annotate_consensus(consensus_df: pd.DataFrame,
                       gene_col: str = "names",
                       panels: Optional[Iterable[str]] = None,
                       data_dir: Optional[Path] = None) -> pd.DataFrame:
    """Return a copy of 'consensus_df' with one 'on_<panel_key>' bool
    column per requested panel. If 'panels' is None, annotate every panel
    found on disk."""
    if gene_col not in consensus_df.columns:
        raise KeyError(f"gene_col '{gene_col}' not in consensus_df columns")

    all_panels = _load_panels(data_dir=data_dir)
    wanted = (list(all_panels) if panels is None
              else [p for p in panels if p in all_panels])

    out = consensus_df.copy()
    harmonised = _harmonise(out[gene_col].fillna("").astype(str).tolist())
    for key in wanted:
        members = all_panels[key]
        out[f"on_{key}"] = [g in members for g in harmonised]
    return out


def panel_summary(consensus_df: pd.DataFrame,
                  gene_col: str = "names",
                  fdr_col: Optional[str] = "fdr",
                  fdr_threshold: float = DEFAULT_FDR,
                  panels: Optional[Iterable[str]] = None,
                  data_dir: Optional[Path] = None) -> dict:
    """Per-panel coverage counts over the significant consensus subset.

    Returns ``{'n_consensus': int, 'panels': {key: {'label', 'n_on_panel',
    'pct'}}}`'. If 'fdr_col` is set, coverage is computed on rows where
    'fdr_col < fdr_threshold'; otherwise on all rows.
    """
    df = consensus_df
    if fdr_col is not None and fdr_col in df.columns:
        df = df[df[fdr_col] < fdr_threshold]

    all_panels = _load_panels(data_dir=data_dir)
    wanted = ([p for p in DEFAULT_PANELS if p in all_panels] if panels is None
              else [p for p in panels if p in all_panels])

    harmonised = set(_harmonise(df[gene_col].fillna("").astype(str).tolist()))
    harmonised.discard("")
    n_consensus = len(harmonised)

    out = {"n_consensus": n_consensus, "panels": {}}
    for key in wanted:
        n_on = len(harmonised & all_panels[key])
        out["panels"][key] = {
            "label":      PANEL_LABELS.get(key, key),
            "n_on_panel": n_on,
            "pct":        (100.0 * n_on / n_consensus) if n_consensus else 0.0,
        }
    return out
