# Preranked GSEA (fgsea-style) over the DESeq2 gene ranks.
#
# The independent gene-set test that corroborates the per-donor pathway
# scores, staying inside the DESeq2 framework: rank the genes by their
# DESeq2 statistic, then ask (threshold-free) whether each pathway is
# enriched toward the top or bottom of that ranking. Delegates to GSEApy's
# validated `prerank` (the fgsea/GSEA preranked algorithm).
#
# Pure logic (no Qt). Consumes the gene-level DE table; needs gene DE run.
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


FGSEA_RESULTS_FILE = 'fgsea_results.csv'
GSEA_LIBRARY_RESULTS_FILE = 'gsea_library_results.csv'


def select_top_nes(res, show_all_max=50, top_per_dir=20):
    """Adaptive gene-set selection for the NES bar chart.

    Shows every set when there are at most 'show_all_max' of them (so a small
    library like the 50 Hallmark sets is drawn whole); otherwise keeps the top
    'top_per_dir' by NES in each direction (for large libraries). Returns a
    DataFrame sorted by NES ascending (most negative first), so both the in-app
    plot and the matplotlib export select and order identically."""
    if res is None or len(res) == 0:
        return res
    df = res.copy()
    df['nes'] = pd.to_numeric(df['nes'], errors='coerce')
    df = df.dropna(subset=['nes'])
    if df.empty:
        return df
    if len(df) <= show_all_max:
        return df.sort_values('nes').reset_index(drop=True)
    down = df[df['nes'] < 0].sort_values('nes').head(top_per_dir)
    up = df[df['nes'] > 0].sort_values('nes', ascending=False).head(top_per_dir)
    return pd.concat([down, up]).sort_values('nes').reset_index(drop=True)


def save_fgsea_results(res, path):
    """Persist the fgsea summary table so figures reload the exact run instead of
    recomputing. Returns the written path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(path, index=False)
    return path


def load_fgsea_results(path):
    """Reload the fgsea summary written by 'save_fgsea_results', or None if the
    file is missing or lacks the expected columns."""
    path = Path(path)
    if not path.exists():
        return None
    df = pd.read_csv(path)
    if 'names' not in df.columns or 'nes' not in df.columns:
        return None
    return df


def run_fgsea(gene_de_results, pathway_gene_sets, rank_by='logfoldchanges',
              min_size=3, max_size=5000, permutation_num=1000, seed=0,
              threads=4, return_details=False):
    """Preranked GSEA of pathways against the ranked gene-level DE results.

    Parameters
    ----------
    gene_de_results : pandas.DataFrame
        Gene DE table; must have ``names`` and the ``rank_by`` column.
    pathway_gene_sets : dict
        ``{pathway_name: [gene_symbols]}``.
    rank_by : str
        Column used to rank genes (default ``'logfoldchanges'``).
    min_size, max_size : int
        Gene-set size bounds (genes present in the ranking). ``min_size=3``
        so small a-priori pathways are still tested (note: GSEA is noisy for
        very small sets).
    permutation_num : int
        Gene-set permutations for the null.
    seed : int
    threads : int
    return_details : bool
        If True, also return a per-pathway detail dict for plotting the
        running-enrichment (GSEA) curve: ``{term: {'RES', 'hits', 'es',
        'nes'}}`` (``RES`` = running enrichment score across the ranked
        genes; ``hits`` = rank positions of the set's genes).

    Returns
    -------
    pandas.DataFrame
        Columns: ``names``, ``nes``, ``es``, ``pvals`` (nominal),
        ``pvals_adj`` (FDR q), ``leading_edge``. One row per scored pathway.
        If ``return_details``, returns ``(summary_df, details_dict)``.
    """
    import gseapy as gp

    if gene_de_results is None or len(gene_de_results) == 0:
        raise ValueError("run_fgsea: gene DE results are empty. Run gene DE first.")
    df = gene_de_results
    if 'names' not in df.columns or rank_by not in df.columns:
        raise ValueError(
            f"run_fgsea: gene DE table needs 'names' and '{rank_by}' columns.")

    # Ranked list: gene -> score. Uppercase both sides (MSigDB convention)
    # so matching is robust; drop NaN scores and duplicate genes.
    rnk = pd.DataFrame({
        'gene': df['names'].astype(str).str.upper(),
        'score': pd.to_numeric(df[rank_by], errors='coerce'),
    }).dropna(subset=['score'])
    rnk = rnk[~rnk['gene'].duplicated()].sort_values('score', ascending=False)
    if rnk.empty:
        empty = pd.DataFrame(columns=['names', 'nes', 'es', 'pvals',
                                      'pvals_adj', 'leading_edge'])
        return (empty, {}) if return_details else empty

    gene_sets = {pw: [str(g).upper() for g in genes]
                 for pw, genes in pathway_gene_sets.items()}

    pre = gp.prerank(
        rnk=rnk, gene_sets=gene_sets,
        min_size=min_size, max_size=max_size,
        permutation_num=permutation_num,
        outdir=None, no_plot=True, seed=seed, threads=threads,
    )
    res = pre.res2d

    summary = pd.DataFrame({
        'names': res['Term'].astype(str),
        'nes': pd.to_numeric(res['NES'], errors='coerce'),
        'es': pd.to_numeric(res['ES'], errors='coerce'),
        'pvals': pd.to_numeric(res['NOM p-val'], errors='coerce'),
        'pvals_adj': pd.to_numeric(res['FDR q-val'], errors='coerce'),
        'leading_edge': res['Lead_genes'].astype(str),
    }).reset_index(drop=True)

    if not return_details:
        return summary

    details = {}
    for term, d in pre.results.items():
        details[term] = {
            'RES': np.asarray(d['RES'], dtype=float),
            'hits': [int(h) for h in d['hits']],
            'es': float(d['es']),
            'nes': float(d['nes']),
        }
    return summary, details


__all__ = ["run_fgsea", "save_fgsea_results", "load_fgsea_results",
           "select_top_nes", "FGSEA_RESULTS_FILE", "GSEA_LIBRARY_RESULTS_FILE"]
