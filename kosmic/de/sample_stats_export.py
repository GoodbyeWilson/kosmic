# Per-sample expression summary stats for the meta-analysis pipeline.
#
# Two entry points:
#
# * 'compute_pathway_stats' -- aggregates expression across each
#   pathway's gene set first (per sample), then summarises per condition.
# * 'compute_gene_level_stats' -- keeps per-gene expression per
#   sample, then summarises per condition.
#
# Both return '(norm_df, raw_df)' where 'norm' uses a log1p(CPM-like)
# transform when 'adata.raw' is available, and 'raw' uses the raw
# counts directly. Each row carries (n, mean, SD, SEM) for both the disease
# and control condition labels.
#
# The companion DE page writes these to ``<project>/results/de_analysis/
# de_stats/<GSE>_<level>Stats_<Norm|Raw>.csv``.
from __future__ import annotations

from typing import Iterable, List, Mapping, Tuple

import numpy as np
import pandas as pd
from scipy.stats import sem

from kosmic.scrna.counts import counts_adata, has_counts_layer
from kosmic.scrna.inspect.detection import detect_species, format_gene_for_species


def _format_genes_for_species(genes: Iterable[str], is_mouse: bool) -> List[str]:
    if not is_mouse:
        return list(genes)
    return [format_gene_for_species(g, 'mouse') for g in genes]


def _to_dense_float64(matrix) -> np.ndarray:
    if hasattr(matrix, 'toarray'):
        return matrix.toarray().astype(np.float64)
    return np.asarray(matrix, dtype=np.float64)


def _summary_stats(values: np.ndarray) -> Tuple[int, float, float, float]:
    n = len(values)
    mean = float(np.mean(values)) if n > 0 else 0.0
    if n > 1:
        sd = float(np.std(values, ddof=1))
        se = float(sem(values))
    else:
        sd = 0.0
        se = 0.0
    return n, mean, sd, se


def _row_summary(base: dict, name_col: str, dis_stats: Tuple[int, float, float, float],
                 ctrl_stats: Tuple[int, float, float, float],
                 disease_label: str, control_label: str) -> dict:
    dn, dm, dsd, dse = dis_stats
    cn, cm, csd, cse = ctrl_stats
    return {
        **base,
        f'N_{disease_label}':       dn,
        f'{disease_label}_Mean':     dm,
        f'{disease_label}_SD':       dsd,
        f'{disease_label}_SEM':      dse,
        f'N_{control_label}':        cn,
        f'{control_label}_Mean':     cm,
        f'{control_label}_SD':       csd,
        f'{control_label}_SEM':      cse,
    }


def _per_sample_pathway_means(adata, pathway_gene_sets: Mapping[str, List[str]],
                              sample_df: pd.DataFrame, sample_col: str,
                              ) -> dict[str, list[dict]]:
    """For each pathway, compute (raw_mean, norm_mean) per sample.

    Returns '{pathway_name: [{condition, raw, norm}, ...]}'. Pathways
    with zero gene coverage are dropped.
    """
    counts = counts_adata(adata, copy=False)
    from_counts = has_counts_layer(adata)
    ref_var_names = list(counts.var_names)
    is_mouse = detect_species(ref_var_names) == 'mouse'
    coverage_var_names = set(ref_var_names)

    out: dict[str, list[dict]] = {}
    for pathway_name, genes in pathway_gene_sets.items():
        genes_fmt = _format_genes_for_species(genes, is_mouse)
        available_genes = [g for g in genes_fmt if g in coverage_var_names]
        if not available_genes:
            continue

        sample_expressions: list[dict] = []
        for _, row in sample_df.iterrows():
            sample_id = row['sample']
            condition = row['condition']
            sample_mask = adata.obs[sample_col] == sample_id

            sample_cells = counts[sample_mask, available_genes]

            if sample_cells.n_obs == 0:
                continue

            raw_expr = _to_dense_float64(sample_cells.X)
            raw_sample_mean = float(np.mean(np.mean(raw_expr, axis=1)))

            if from_counts:
                cell_totals = raw_expr.sum(axis=1, keepdims=True)
                cell_totals[cell_totals == 0] = 1
                norm_expr = np.log1p(raw_expr / cell_totals * 1e4)
            else:
                norm_expr = raw_expr
            norm_sample_mean = float(np.mean(np.mean(norm_expr, axis=1)))

            sample_expressions.append({
                'condition': condition,
                'norm': norm_sample_mean,
                'raw': raw_sample_mean,
            })

        if sample_expressions:
            out[pathway_name] = sample_expressions
    return out


def compute_pathway_stats(adata, sample_df: pd.DataFrame,
                          sample_col: str, condition_col: str,
                          control_label: str, disease_label: str,
                          pathway_gene_sets: Mapping[str, List[str]],
                          gse_id: str,
                          ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Pathway-level per-sample expression stats."""
    pathway_samples = _per_sample_pathway_means(
        adata, pathway_gene_sets, sample_df, sample_col)

    norm_rows: list[dict] = []
    raw_rows: list[dict] = []
    for pathway_name, sample_exprs in pathway_samples.items():
        expr_df = pd.DataFrame(sample_exprs)

        cond_norm: dict[str, np.ndarray] = {}
        cond_raw: dict[str, np.ndarray] = {}
        for cond in (control_label, disease_label):
            cond_df = expr_df[expr_df['condition'] == cond]
            if cond_df.empty:
                continue
            cond_norm[cond] = cond_df['norm'].values
            cond_raw[cond] = cond_df['raw'].values

        if (control_label not in cond_norm
                or disease_label not in cond_norm):
            continue

        base = {'GSE': gse_id, 'Pathway': pathway_name}
        norm_rows.append(_row_summary(
            base, 'Pathway',
            _summary_stats(cond_norm[disease_label]),
            _summary_stats(cond_norm[control_label]),
            disease_label, control_label,
        ))
        raw_rows.append(_row_summary(
            base, 'Pathway',
            _summary_stats(cond_raw[disease_label]),
            _summary_stats(cond_raw[control_label]),
            disease_label, control_label,
        ))

    return pd.DataFrame(norm_rows), pd.DataFrame(raw_rows)


def _per_sample_per_gene_means(adata, gene_list: List[str],
                               sample_df: pd.DataFrame, sample_col: str,
                               ) -> pd.DataFrame:
    """Return long-form '(sample, condition, gene, norm_mean, raw_mean)'
    dataframe for the requested gene list."""
    counts = counts_adata(adata, copy=False)
    from_counts = has_counts_layer(adata)
    rows: list[dict] = []
    for _, row in sample_df.iterrows():
        sample_id = row['sample']
        condition = row['condition']
        sample_mask = adata.obs[sample_col] == sample_id

        genes_to_use = [g for g in gene_list if g in counts.var_names]
        if not genes_to_use:
            continue
        sample_cells = counts[sample_mask, genes_to_use]

        if sample_cells.n_obs == 0:
            continue

        raw_expr = _to_dense_float64(sample_cells.X)
        if from_counts:
            cell_totals = raw_expr.sum(axis=1, keepdims=True)
            cell_totals[cell_totals == 0] = 1
            norm_expr = np.log1p(raw_expr / cell_totals * 1e4)
        else:
            norm_expr = raw_expr

        for i, gene in enumerate(genes_to_use):
            rows.append({
                'sample': sample_id,
                'condition': condition,
                'gene': gene,
                'norm_mean': float(np.mean(norm_expr[:, i])),
                'raw_mean': float(np.mean(raw_expr[:, i])),
            })
    return pd.DataFrame(rows)


def compute_gene_level_stats(adata, sample_df: pd.DataFrame,
                             sample_col: str, condition_col: str,
                             control_label: str, disease_label: str,
                             pathway_gene_sets: Mapping[str, List[str]],
                             gse_id: str,
                             ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Gene-level per-sample expression stats restricted to pathway genes."""
    ref_var_names = list(counts_adata(adata, copy=False).var_names)
    is_mouse = detect_species(ref_var_names) == 'mouse'
    ref_var_set = set(ref_var_names)

    all_pathway_genes: dict[str, list[str]] = {}
    for pathway_name, genes in pathway_gene_sets.items():
        genes_fmt = _format_genes_for_species(genes, is_mouse)
        for gene in genes_fmt:
            if gene in ref_var_set:
                all_pathway_genes.setdefault(gene, []).append(pathway_name)

    if not all_pathway_genes:
        return pd.DataFrame(), pd.DataFrame()

    gene_list = list(all_pathway_genes.keys())
    long_df = _per_sample_per_gene_means(adata, gene_list, sample_df, sample_col)
    if long_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    norm_rows: list[dict] = []
    raw_rows: list[dict] = []
    for gene in gene_list:
        gene_data = long_df[long_df['gene'] == gene]
        ctrl_data = gene_data[gene_data['condition'] == control_label]
        dis_data = gene_data[gene_data['condition'] == disease_label]
        if ctrl_data.empty or dis_data.empty:
            continue
        pathways_for_gene = '|'.join(all_pathway_genes.get(gene, ['Unknown']))
        base = {'GSE': gse_id, 'Gene': gene, 'Pathway': pathways_for_gene}

        norm_rows.append(_row_summary(
            base, 'Gene',
            _summary_stats(dis_data['norm_mean'].values),
            _summary_stats(ctrl_data['norm_mean'].values),
            disease_label, control_label,
        ))
        raw_rows.append(_row_summary(
            base, 'Gene',
            _summary_stats(dis_data['raw_mean'].values),
            _summary_stats(ctrl_data['raw_mean'].values),
            disease_label, control_label,
        ))

    return pd.DataFrame(norm_rows), pd.DataFrame(raw_rows)
