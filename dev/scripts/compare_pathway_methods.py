"""
Simulation: Compare pathway-level meta-analysis methods
=========================================================
Realistic synthetic scRNA-seq with REAL gene names and REAL pathway databases
(KEGG 2021 + MSigDB Hallmark). 5 KEGG pathways get true effects, rest are null.

Methods:
  1. Mean expression + DerSimonian-Laird
  2. AUCell (decoupler) + DerSimonian-Laird
  3. DESeq2 (pydeseq2) gene-level → DL within pathway → DL across datasets
  4. GSEA (gseapy prerank) + DerSimonian-Laird on NES
  5. GSEA SumRank (scRNA_Reproducibility paper method)

Usage:
  conda activate kirk
  python scripts/compare_pathway_methods.py
"""

import numpy as np
import pandas as pd
import anndata as ad
import gseapy as gp
from scipy import stats
from scipy.stats import rankdata
from multiprocessing import Pool, cpu_count
import time
import os
import warnings
warnings.filterwarnings('ignore')
os.environ['NUMBA_NUM_THREADS'] = '1'

# ---------------------------------------------------------------------------
# 1. Load real pathway databases and set up simulation
# ---------------------------------------------------------------------------

print("Loading pathway databases...", flush=True)
KEGG = gp.get_library('KEGG_2021_Human')
HALLMARK = gp.get_library('MSigDB_Hallmark_2020')

# Combine into one database for GSEA
ALL_GENE_SETS = {}
ALL_GENE_SETS.update({f'KEGG_{k}': v for k, v in KEGG.items()})
ALL_GENE_SETS.update({f'HM_{k}': v for k, v in HALLMARK.items()})
print(f"  {len(KEGG)} KEGG + {len(HALLMARK)} Hallmark = {len(ALL_GENE_SETS)} total pathways", flush=True)

# Collect all genes that appear in any pathway
ALL_PATHWAY_GENES = set()
for genes in ALL_GENE_SETS.values():
    ALL_PATHWAY_GENES.update(genes)

# Use ~15,000 genes: all pathway genes + random "background" genes to fill out
N_PATHWAY_GENES = len(ALL_PATHWAY_GENES)
N_BACKGROUND_GENES = 10000
GENE_NAMES = sorted(ALL_PATHWAY_GENES) + [f'BG_{i}' for i in range(N_BACKGROUND_GENES)]
N_GENES = len(GENE_NAMES)
GENE_TO_IDX = {g: i for i, g in enumerate(GENE_NAMES)}
print(f"  {N_PATHWAY_GENES} pathway genes + {N_BACKGROUND_GENES} background = {N_GENES} total genes", flush=True)

# Convert pathway gene sets to index arrays
PATHWAY_IDX = {}
for pw_name, gene_list in ALL_GENE_SETS.items():
    idxs = [GENE_TO_IDX[g] for g in gene_list if g in GENE_TO_IDX]
    if len(idxs) >= 5:
        PATHWAY_IDX[pw_name] = idxs

# Select 5 KEGG pathways to have TRUE effects
TRUE_EFFECT_PATHWAYS = {
    'KEGG_Glycolysis / Gluconeogenesis': 0.5,    # moderate up
    'KEGG_Citrate cycle (TCA cycle)':    -0.4,   # moderate down
    'KEGG_Oxidative phosphorylation':    0.7,    # strong up
    'KEGG_Fatty acid degradation':       0.3,    # mild up
    'KEGG_Pyruvate metabolism':          0.15,    # weak signal
}

# Also track some null pathways for reporting
NULL_REPORT_PATHWAYS = [
    'KEGG_ABC transporters',
    'KEGG_Calcium signaling pathway',
    'KEGG_Cell cycle',
]

# Pathways to report in results table
REPORT_PATHWAYS = list(TRUE_EFFECT_PATHWAYS.keys()) + NULL_REPORT_PATHWAYS

print(f"  True effect pathways: {len(TRUE_EFFECT_PATHWAYS)}", flush=True)
for pw, eff in TRUE_EFFECT_PATHWAYS.items():
    n_genes = len(PATHWAY_IDX.get(pw, []))
    print(f"    {pw}: log2FC={eff:+.2f}, {n_genes} genes", flush=True)
print()

# Simulation parameters
N_DATASETS = 8
SAMPLES_PER_GROUP_RANGE = (3, 8)
CELLS_PER_SAMPLE_RANGE = (50, 500)
NB_DISPERSION_RANGE = (0.05, 0.8)
DROPOUT_RATE = 0.05
BETWEEN_STUDY_TAU = 0.15
WITHIN_PATHWAY_GENE_SD = 0.3
WITHIN_PATHWAY_CORRELATION = 0.3


# ---------------------------------------------------------------------------
# 2. Data generation
# ---------------------------------------------------------------------------

def _build_correlation_matrix(pw_size, rho):
    C = np.full((pw_size, pw_size), rho)
    np.fill_diagonal(C, 1.0)
    return C


def simulate_dataset(dataset_idx, rng=None):
    if rng is None:
        rng = np.random.default_rng()

    # Per-gene true log2FC (only pathway genes with true effects get non-zero)
    gene_true_log2fc = np.zeros(N_GENES)

    for pw_name, base_effect in TRUE_EFFECT_PATHWAYS.items():
        if pw_name not in PATHWAY_IDX:
            continue
        pw_genes = PATHWAY_IDX[pw_name]
        pw_size = len(pw_genes)

        # This dataset's effect for this pathway
        dataset_effect = base_effect + rng.normal(0, BETWEEN_STUDY_TAU)

        # Correlated gene effects within pathway
        C = _build_correlation_matrix(pw_size, WITHIN_PATHWAY_CORRELATION)
        L = np.linalg.cholesky(C)
        z = rng.normal(0, 1, pw_size)
        correlated_noise = L @ z * WITHIN_PATHWAY_GENE_SD

        for i, g in enumerate(pw_genes):
            # Add effect (genes can be in multiple pathways — effects accumulate)
            gene_true_log2fc[g] += dataset_effect + correlated_noise[i]

    # Baseline expression calibrated to real scRNA-seq
    log_baseline = rng.normal(-4.0, 2.0, N_GENES)
    baseline_mu = np.exp(log_baseline)
    baseline_mu = np.clip(baseline_mu, 0.001, 200)

    gene_dispersion = np.exp(rng.uniform(
        np.log(NB_DISPERSION_RANGE[0]),
        np.log(NB_DISPERSION_RANGE[1]),
        N_GENES
    ))

    n_per_group = rng.integers(*SAMPLES_PER_GROUP_RANGE, endpoint=True)

    pseudobulk_sum_data = []
    pseudobulk_mean_data = []
    cell_level_data = []
    conditions = []
    sample_ids = []
    cells_per_sample_list = []

    for group in ['control', 'disease']:
        for s in range(n_per_group):
            sample_id = f'D{dataset_idx}_{group}_{s}'
            n_cells = rng.integers(*CELLS_PER_SAMPLE_RANGE)
            cells_per_sample_list.append(n_cells)

            mu = baseline_mu.copy()
            if group == 'disease':
                mu = mu * (2.0 ** gene_true_log2fc)

            r = gene_dispersion
            p = r / (r + mu)
            p = np.clip(p, 1e-10, 1 - 1e-10)

            cells = rng.negative_binomial(r, p, size=(n_cells, N_GENES)).astype(np.float64)
            dropout_mask = rng.random((n_cells, N_GENES)) < DROPOUT_RATE
            cells[dropout_mask] = 0

            cell_level_data.append(cells)
            pseudobulk_sum_data.append(cells.sum(axis=0))
            pseudobulk_mean_data.append(cells.mean(axis=0))
            conditions.append(group)
            sample_ids.append(sample_id)

    pseudobulk_sum = np.array(pseudobulk_sum_data, dtype=np.float64)
    pseudobulk_mean = np.array(pseudobulk_mean_data, dtype=np.float64)
    cell_level = np.vstack(cell_level_data)

    cell_conditions = np.empty(cell_level.shape[0], dtype='U10')
    cell_sample_ids = np.empty(cell_level.shape[0], dtype='U30')
    offset = 0
    for i, (cond, sid) in enumerate(zip(conditions, sample_ids)):
        nc = cells_per_sample_list[i]
        cell_conditions[offset:offset + nc] = cond
        cell_sample_ids[offset:offset + nc] = sid
        offset += nc

    return {
        'pseudobulk_sum': pseudobulk_sum,
        'pseudobulk_mean': pseudobulk_mean,
        'cell_level': cell_level,
        'cell_conditions': cell_conditions,
        'cell_sample_ids': cell_sample_ids,
        'conditions': conditions,
        'sample_ids': sample_ids,
        'cells_per_sample': cells_per_sample_list,
    }


# ---------------------------------------------------------------------------
# 3. DESeq2 (pydeseq2)
# ---------------------------------------------------------------------------

def run_deseq2(dataset):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts_df = pd.DataFrame(
        dataset['pseudobulk_sum'].astype(int),
        columns=GENE_NAMES,
        index=dataset['sample_ids'],
    )
    meta_df = pd.DataFrame(
        {'condition': dataset['conditions']},
        index=dataset['sample_ids'],
    )

    nonzero_mask = counts_df.sum(axis=0) > 0
    counts_filt = counts_df.loc[:, nonzero_mask]

    try:
        dds = DeseqDataSet(counts=counts_filt, metadata=meta_df, design='~condition', quiet=True)
        dds.deseq2()
        stat = DeseqStats(dds, contrast=['condition', 'disease', 'control'], quiet=True)
        stat.summary()
        res = stat.results_df

        out = pd.DataFrame({
            'gene': GENE_NAMES,
            'gene_idx': np.arange(N_GENES),
            'log2fc': 0.0,
            'se': 1.0,
            'pval': 1.0,
        })
        for gene_name in res.index:
            if gene_name in GENE_TO_IDX:
                idx = GENE_TO_IDX[gene_name]
                row = res.loc[gene_name]
                lfc = row['log2FoldChange']
                se = row['lfcSE']
                pv = row['pvalue']
                if pd.notna(lfc) and pd.notna(se) and pd.notna(pv):
                    out.loc[idx, 'log2fc'] = lfc
                    out.loc[idx, 'se'] = max(se, 0.001)
                    out.loc[idx, 'pval'] = pv
        return out

    except Exception:
        return _run_de_ttest(dataset)


def _run_de_ttest(dataset):
    pb = dataset['pseudobulk_sum']
    conds = np.array(dataset['conditions'])
    d_mask = conds == 'disease'
    c_mask = conds == 'control'

    lib_sizes = pb.sum(axis=1, keepdims=True)
    lib_sizes[lib_sizes == 0] = 1
    median_lib = np.median(lib_sizes)
    pb_log = np.log2(pb / lib_sizes * median_lib + 1)

    d_vals = pb_log[d_mask]
    c_vals = pb_log[c_mask]
    nd, nc = d_vals.shape[0], c_vals.shape[0]

    log2fc = d_vals.mean(axis=0) - c_vals.mean(axis=0)
    d_var = d_vals.var(axis=0, ddof=1)
    c_var = c_vals.var(axis=0, ddof=1)
    se = np.sqrt(d_var / nd + c_var / nc)
    se = np.maximum(se, 0.001)

    t_stat = log2fc / se
    num = (d_var / nd + c_var / nc) ** 2
    denom = (d_var / nd) ** 2 / max(nd - 1, 1) + (c_var / nc) ** 2 / max(nc - 1, 1)
    denom = np.maximum(denom, 1e-30)
    df = np.maximum(num / denom, 1.0)

    pvals = 2 * stats.t.sf(np.abs(t_stat), df)
    pvals = np.where(np.isnan(pvals), 1.0, pvals)

    return pd.DataFrame({
        'gene': GENE_NAMES,
        'gene_idx': np.arange(N_GENES),
        'log2fc': log2fc,
        'se': se,
        'pval': pvals,
    })


# ---------------------------------------------------------------------------
# 4. Method 1: Mean expression + DL
# ---------------------------------------------------------------------------

def method1_mean_expression(datasets):
    results = {}
    for pw_name in REPORT_PATHWAYS:
        if pw_name not in PATHWAY_IDX:
            results[pw_name] = {'effect': 0, 'se': 1, 'pval': 1, 'I2': 0}
            continue
        pw_idx = np.array(PATHWAY_IDX[pw_name])
        dataset_effects = []
        dataset_ses = []

        for ds in datasets:
            pb_log = np.log2(ds['pseudobulk_mean'][:, pw_idx] + 1)
            pw_expr = pb_log.mean(axis=1)

            conds = np.array(ds['conditions'])
            d_vals = pw_expr[conds == 'disease']
            c_vals = pw_expr[conds == 'control']

            effect = d_vals.mean() - c_vals.mean()
            se = np.sqrt(d_vals.var(ddof=1) / len(d_vals) + c_vals.var(ddof=1) / len(c_vals))
            dataset_effects.append(effect)
            dataset_ses.append(max(se, 0.001))

        pooled_effect, pooled_se, pval, i2 = dersimonian_laird(
            np.array(dataset_effects), np.array(dataset_ses)
        )
        results[pw_name] = {'effect': pooled_effect, 'se': pooled_se, 'pval': pval, 'I2': i2}
    return results


# ---------------------------------------------------------------------------
# 5. Method 2: AUCell via decoupler
# ---------------------------------------------------------------------------

def method2_aucell(datasets):
    import decoupler as dc

    # Build net for our report pathways only (speed)
    net_rows = []
    for pw_name in REPORT_PATHWAYS:
        if pw_name not in ALL_GENE_SETS:
            # Try without prefix
            raw_name = pw_name.replace('KEGG_', '').replace('HM_', '')
            gene_list = KEGG.get(raw_name, HALLMARK.get(raw_name, []))
        else:
            gene_list = ALL_GENE_SETS[pw_name]
        for g in gene_list:
            if g in GENE_TO_IDX:
                net_rows.append({'source': pw_name, 'target': g, 'weight': 1.0})
    net = pd.DataFrame(net_rows)

    results = {}
    for pw_name in REPORT_PATHWAYS:
        dataset_effects = []
        dataset_ses = []

        for ds in datasets:
            # Build AnnData and run AUCell (cache per dataset)
            if '_aucell_scores' not in ds:
                adata = ad.AnnData(
                    X=ds['cell_level'].astype(np.float32),
                    obs=pd.DataFrame({
                        'condition': ds['cell_conditions'],
                        'sample': ds['cell_sample_ids'],
                    }, index=[f'cell_{i}' for i in range(ds['cell_level'].shape[0])]),
                    var=pd.DataFrame(index=GENE_NAMES),
                )
                dc.mt.aucell(adata, net, verbose=False)
                ds['_aucell_scores'] = adata.obsm['score_aucell']

            scores_df = ds['_aucell_scores']
            if pw_name not in scores_df.columns:
                dataset_effects.append(0.0)
                dataset_ses.append(1.0)
                continue

            scores = scores_df[pw_name].values
            sample_ids = ds['cell_sample_ids']
            unique_samples = ds['sample_ids']
            sample_scores = np.empty(len(unique_samples))
            sample_conds = np.array(ds['conditions'])

            for s_i, sid in enumerate(unique_samples):
                mask = sample_ids == sid
                sample_scores[s_i] = scores[mask].mean()

            d_scores = sample_scores[sample_conds == 'disease']
            c_scores = sample_scores[sample_conds == 'control']

            effect = d_scores.mean() - c_scores.mean()
            se = np.sqrt(d_scores.var(ddof=1) / len(d_scores) + c_scores.var(ddof=1) / len(c_scores))
            dataset_effects.append(effect)
            dataset_ses.append(max(se, 0.001))

        pooled_effect, pooled_se, pval, i2 = dersimonian_laird(
            np.array(dataset_effects), np.array(dataset_ses)
        )
        results[pw_name] = {'effect': pooled_effect, 'se': pooled_se, 'pval': pval, 'I2': i2}

    return results


# ---------------------------------------------------------------------------
# 6. Method 3: DESeq2 gene-level → DL within pathway → DL across datasets
# ---------------------------------------------------------------------------

def method3_deseq2_dl_pooling(de_results_list):
    results = {}
    for pw_name in REPORT_PATHWAYS:
        if pw_name not in PATHWAY_IDX:
            results[pw_name] = {'effect': 0, 'se': 1, 'pval': 1, 'I2': 0}
            continue
        pw_set = set(PATHWAY_IDX[pw_name])
        dataset_pw_effects = []
        dataset_pw_ses = []

        for de_df in de_results_list:
            pw_mask = de_df['gene_idx'].isin(pw_set)
            pw_de = de_df.loc[pw_mask]
            pw_de = pw_de[(pw_de['pval'] < 1.0) & (pw_de['se'] < 10)]
            if len(pw_de) < 2:
                continue

            pw_effect, pw_se, _, _ = dersimonian_laird(
                pw_de['log2fc'].values, pw_de['se'].values
            )
            dataset_pw_effects.append(pw_effect)
            dataset_pw_ses.append(pw_se)

        if len(dataset_pw_effects) < 2:
            results[pw_name] = {'effect': 0, 'se': 1, 'pval': 1, 'I2': 0}
            continue

        pooled_effect, pooled_se, pval, i2 = dersimonian_laird(
            np.array(dataset_pw_effects), np.array(dataset_pw_ses)
        )
        results[pw_name] = {'effect': pooled_effect, 'se': pooled_se, 'pval': pval, 'I2': i2}
    return results


# ---------------------------------------------------------------------------
# 7. Run GSEA once per dataset (shared by Methods 4 and 5)
# ---------------------------------------------------------------------------

def _run_gsea_all_datasets(de_results_list):
    """Run gseapy prerank once per dataset against all pathways.
    Returns:
      nes_matrix: (n_datasets x n_pathways) NES values
      pval_matrix: (n_datasets x n_pathways) nominal p-values
      es_matrix: (n_datasets x n_pathways) enrichment scores (for sign)
      all_pw_names: list of pathway names
    """
    all_pw_names = list(ALL_GENE_SETS.keys())
    n_all = len(all_pw_names)
    n_datasets = len(de_results_list)
    nes_matrix = np.zeros((n_datasets, n_all))
    pval_matrix = np.ones((n_datasets, n_all))
    es_matrix = np.zeros((n_datasets, n_all))
    pw_name_to_idx = {pw: i for i, pw in enumerate(all_pw_names)}

    from kosmic.numerical import neg_log10
    for d, de_df in enumerate(de_results_list):
        ranking = de_df.copy()
        ranking['stat'] = neg_log10(ranking['pval']) * np.sign(ranking['log2fc'])
        rank_series = pd.Series(
            ranking['stat'].values,
            index=[GENE_NAMES[int(i)] for i in ranking['gene_idx'].values]
        ).sort_values(ascending=False)

        try:
            res = gp.prerank(
                rnk=rank_series,
                gene_sets=ALL_GENE_SETS,
                min_size=5,
                max_size=500,
                permutation_num=1000,
                seed=42,
                no_plot=True,
                verbose=False,
            )
            for _, row in res.res2d.iterrows():
                term = row['Term']
                if term in pw_name_to_idx:
                    idx = pw_name_to_idx[term]
                    nes_matrix[d, idx] = float(row['NES'])
                    es_matrix[d, idx] = float(row['ES'])
                    pv = float(row['NOM p-val'])
                    pval_matrix[d, idx] = max(pv, 1e-300)
        except Exception:
            pass

    return nes_matrix, pval_matrix, es_matrix, all_pw_names


# ---------------------------------------------------------------------------
# 8. Method 4: GSEA + DerSimonian-Laird on NES
# ---------------------------------------------------------------------------

def method4_gsea_dl(nes_matrix, all_pw_names):
    """Pool NES across datasets using DerSimonian-Laird (SE=1 by convention)."""
    pw_name_to_idx = {pw: i for i, pw in enumerate(all_pw_names)}
    results = {}
    for pw_name in REPORT_PATHWAYS:
        if pw_name not in pw_name_to_idx:
            results[pw_name] = {'effect': 0, 'se': 1, 'pval': 1, 'I2': 0}
            continue
        p_i = pw_name_to_idx[pw_name]
        nes_values = nes_matrix[:, p_i]
        ses = np.ones(len(nes_values))
        pooled_effect, pooled_se, pval, i2 = dersimonian_laird(nes_values, ses)
        results[pw_name] = {'effect': pooled_effect, 'se': pooled_se, 'pval': pval, 'I2': i2}
    return results


# ---------------------------------------------------------------------------
# 9. Method 5: GSEA SumRank — faithful to scRNA_Reproducibility code
# ---------------------------------------------------------------------------

def method5_gsea_sumrank(nes_matrix, pval_matrix, es_matrix, all_pw_names):
    """Faithful implementation of SumRank_GSEAPathways.R:
    1. Per dataset: signed -log10(GSEA p-value) per pathway (negative if ES < 0)
    2. Rank pathways descending by signed stat (rank 1 = most significant upregulated)
    3. Restrict to pathways present in all datasets, re-rank
    4. Normalize ranks to [0,1]: (rank-1)/(N-1)
    5. Sum across datasets (using all datasets = 100%)
    6. Cap sum at n_datasets/2
    7. Irwin-Hall CDF for p-value
    8. Test both directions (up and down) for two-sided result
    """
    n_datasets, n_all = nes_matrix.shape
    pw_name_to_idx = {pw: i for i, pw in enumerate(all_pw_names)}

    # Step 1: Compute signed -log10(GSEA p-value) per pathway per dataset
    # Matches their code: SignedNegLogPVal = -log10(pvalue), negative if ES < 0
    from kosmic.numerical import neg_log10
    signed_stat = np.zeros_like(pval_matrix)
    for d in range(n_datasets):
        neglogp = neg_log10(pval_matrix[d, :])
        signs = np.sign(es_matrix[d, :])
        signs[signs == 0] = 1  # treat zero ES as positive
        signed_stat[d, :] = neglogp * signs

    # Step 2: Find pathways present in all datasets (non-zero ES in all)
    # A pathway is "present" if GSEA returned a result for it
    present_mask = np.all(np.abs(es_matrix) > 0, axis=0)
    # Also need at least some pathways
    present_indices = np.where(present_mask)[0]
    n_present = len(present_indices)
    if n_present < 2:
        # Fallback: use all pathways
        present_indices = np.arange(n_all)
        n_present = n_all

    # Step 3: Rank pathways within each dataset by signed stat (descending)
    # Rank 1 = highest signed stat = most significant upregulated
    # Matches their: PathwaysDataframe$SignedNegLogPVal_rank = 1:nrow(...)
    rank_matrix = np.zeros((n_datasets, n_present))
    for d in range(n_datasets):
        vals = signed_stat[d, present_indices]
        # rankdata gives rank 1 = smallest, so negate for descending
        rank_matrix[d, :] = rankdata(-vals)

    # Step 4: Normalize to [0,1]: (rank-1)/(N-1)
    # Matches their: (currentTest$SignedNegLogPVal_rank-1)/(nrow(FinalPathwayList_Table)-1)
    norm_ranks = (rank_matrix - 1) / max(n_present - 1, 1)

    # Step 5: Sum across all datasets (ProportionofTopDatasets = 1.0)
    # With 8 datasets, using 100% — no top-K% filtering needed
    sum_ranks = norm_ranks.sum(axis=0)

    # Build lookup from present index position to original pathway index
    present_pw_lookup = {present_indices[k]: k for k in range(n_present)}

    # Step 6 & 7: Cap and compute p-values
    results = {}
    n_used = n_datasets  # all datasets used (100%)

    for pw_name in REPORT_PATHWAYS:
        if pw_name not in pw_name_to_idx:
            results[pw_name] = {'effect': 0, 'se': 0, 'pval': 1, 'I2': 0}
            continue
        orig_idx = pw_name_to_idx[pw_name]
        if orig_idx not in present_pw_lookup:
            results[pw_name] = {'effect': 0, 'se': 0, 'pval': 1, 'I2': 0}
            continue

        k = present_pw_lookup[orig_idx]
        sr = sum_ranks[k]

        # Test UPREGULATION: small sum = consistently top-ranked (upregulated)
        # Cap at n_used/2 (their code: if sum > n/2, set to n/2)
        sr_up = min(sr, n_used / 2)
        p_up = _irwin_hall_cdf(sr_up, n_used)

        # Test DOWNREGULATION: large sum = consistently bottom-ranked
        # Equivalent to flipping ranks and testing the other tail
        sr_down = min(n_used - sr, n_used / 2)
        p_down = _irwin_hall_cdf(sr_down, n_used)

        # Two-sided: take better direction, Bonferroni for 2 tests
        pval = min(p_up, p_down) * 2
        pval = min(pval, 1.0)

        results[pw_name] = {
            'effect': float(nes_matrix[:, orig_idx].mean()),
            'se': 0, 'pval': pval, 'I2': 0,
        }
    return results


def _irwin_hall_cdf(x, n):
    from math import factorial, floor
    if x <= 0:
        return 0.0
    if x >= n:
        return 1.0
    cdf = 0.0
    for k in range(int(floor(x)) + 1):
        sign = (-1) ** k
        binom_coeff = factorial(n) / (factorial(k) * factorial(n - k))
        cdf += sign * binom_coeff * (x - k) ** n
    cdf /= factorial(n)
    return max(0.0, min(1.0, cdf))


# ---------------------------------------------------------------------------
# DerSimonian-Laird
# ---------------------------------------------------------------------------

def dersimonian_laird(effects, ses):
    k = len(effects)
    if k < 2:
        return effects[0], ses[0], 1.0, 0.0

    weights = 1.0 / (ses ** 2)
    fe_effect = np.sum(weights * effects) / np.sum(weights)
    Q = np.sum(weights * (effects - fe_effect) ** 2)
    c = np.sum(weights) - np.sum(weights ** 2) / np.sum(weights)
    tau2 = max(0, (Q - (k - 1)) / c)

    re_weights = 1.0 / (ses ** 2 + tau2)
    pooled_effect = np.sum(re_weights * effects) / np.sum(re_weights)
    pooled_se = np.sqrt(1.0 / np.sum(re_weights))

    z = pooled_effect / pooled_se
    pval = 2 * (1 - stats.norm.cdf(abs(z)))

    i2 = max(0, (Q - (k - 1)) / Q * 100) if Q > 0 else 0
    return pooled_effect, pooled_se, pval, i2


# ---------------------------------------------------------------------------
# Single run
# ---------------------------------------------------------------------------

def _run_single(run_idx):
    rng = np.random.default_rng(run_idx * 1000 + 42)

    datasets = [simulate_dataset(i, rng) for i in range(N_DATASETS)]
    de_results_list = [run_deseq2(ds) for ds in datasets]

    r1 = method1_mean_expression(datasets)
    r2 = method2_aucell(datasets)
    r3 = method3_deseq2_dl_pooling(de_results_list)

    # Run GSEA once, shared by Methods 4 and 5
    nes_matrix, pval_matrix, es_matrix, all_pw_names = _run_gsea_all_datasets(de_results_list)
    r4 = method4_gsea_dl(nes_matrix, all_pw_names)
    r5 = method5_gsea_sumrank(nes_matrix, pval_matrix, es_matrix, all_pw_names)

    method_names = ['1_MeanExpr', '2_AUCell', '3_DESeq2_DL', '4_GSEA_DL', '5_GSEA_SumRank']
    all_results = {'1_MeanExpr': r1, '2_AUCell': r2, '3_DESeq2_DL': r3, '4_GSEA_DL': r4, '5_GSEA_SumRank': r5}

    pvals = {}
    for m in method_names:
        for pw in REPORT_PATHWAYS:
            pvals[(m, pw)] = all_results[m][pw]['pval']

    return pvals


# ---------------------------------------------------------------------------
# Run simulation
# ---------------------------------------------------------------------------

def run_simulation(n_runs=50, n_workers=None):
    if n_workers is None:
        n_workers = max(1, cpu_count() - 1)

    method_names = ['1_MeanExpr', '2_AUCell', '3_DESeq2_DL', '4_GSEA_DL', '5_GSEA_SumRank']

    detections = {m: {pw: 0 for pw in REPORT_PATHWAYS} for m in method_names}
    pval_sums = {m: {pw: 0.0 for pw in REPORT_PATHWAYS} for m in method_names}

    print(f"  Using {n_workers} worker processes", flush=True)

    completed = 0
    with Pool(n_workers) as pool:
        for result in pool.imap_unordered(_run_single, range(n_runs)):
            completed += 1
            if completed % 5 == 0 or completed == n_runs:
                print(f"  Completed {completed}/{n_runs}", flush=True)

            for m in method_names:
                for pw in REPORT_PATHWAYS:
                    pval = result[(m, pw)]
                    pval_sums[m][pw] += pval
                    if pval < 0.05:
                        detections[m][pw] += 1

    return detections, pval_sums, n_runs


def print_results(detections, pval_sums, n_runs):
    method_names = ['1_MeanExpr', '2_AUCell', '3_DESeq2_DL', '4_GSEA_DL', '5_GSEA_SumRank']

    print("\n" + "=" * 120)
    print("SIMULATION RESULTS: Power (% of runs with p < 0.05)")
    print(f"  {n_runs} runs, {N_DATASETS} datasets/run, {N_GENES} genes")
    print(f"  {SAMPLES_PER_GROUP_RANGE[0]}-{SAMPLES_PER_GROUP_RANGE[1]} samples/group, "
          f"{CELLS_PER_SAMPLE_RANGE[0]}-{CELLS_PER_SAMPLE_RANGE[1]} cells/sample")
    print(f"  DE: pydeseq2 | AUCell: decoupler | GSEA: gseapy (370 real pathways)")
    print("=" * 120)

    header = f"{'Pathway':<45} {'Effect':>6} {'Size':>4}"
    for m in method_names:
        header += f" {m:>14}"
    print(header)
    print("-" * 120)

    for pw in REPORT_PATHWAYS:
        true_eff = TRUE_EFFECT_PATHWAYS.get(pw, 0.0)
        pw_size = len(PATHWAY_IDX.get(pw, []))
        short_name = pw.replace('KEGG_', '').replace('HM_', '')[:42]
        row = f"{short_name:<45} {true_eff:>+6.2f} {pw_size:>4}"
        for m in method_names:
            power = detections[m][pw] / n_runs * 100
            row += f" {power:>13.1f}%"
        print(row)

    print("-" * 120)

    print("\n" + "=" * 120)
    print("SUMMARY")
    print("=" * 120)
    null_pws = [pw for pw in REPORT_PATHWAYS if pw not in TRUE_EFFECT_PATHWAYS]
    real_pws = [pw for pw in REPORT_PATHWAYS if pw in TRUE_EFFECT_PATHWAYS]

    for m in method_names:
        true_pos = sum(detections[m][pw] for pw in real_pws) / (len(real_pws) * n_runs) * 100
        false_pos = sum(detections[m][pw] for pw in null_pws) / (len(null_pws) * n_runs) * 100
        print(f"  {m:<18}  Sensitivity: {true_pos:5.1f}%   FPR: {false_pos:5.1f}%")


if __name__ == '__main__':
    t0 = time.time()
    print(f"\nStarting simulation...\n")

    detections, pval_sums, n_runs = run_simulation(n_runs=50)
    print_results(detections, pval_sums, n_runs)

    elapsed = time.time() - t0
    print(f"\nTotal time: {elapsed/60:.1f} minutes ({elapsed:.0f}s)")
