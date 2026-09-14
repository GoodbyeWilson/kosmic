"""
Proteomics Validation
======================
Validates a gene-level consensus DE result against an external proteomics
dataset (PRIDE PXD043768 -- Behounek et al. 2025, Circ Heart Fail).

Pipeline:
    1. Download MaxQuant proteinGroups.txt + sample description from PRIDE
       (cached locally on first use).
    2. Parse proteinGroups: drop reverse/contaminant/only-by-site rows,
       extract gene names from fasta headers (GN= field), keep LFQ
       intensity columns.
    3. Compute differential abundance (Welch t-test on log2-LFQ) for
       HFrEF (combined RVD + non-RVD) vs Control donor LV samples.
    4. Match to a transcriptomic consensus gene list and report
       concordance metrics.
"""

from __future__ import annotations

import os
import re
import urllib.request
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from kosmic import DEFAULT_FDR

# PRIDE PXD043768 file URLs (LV proteinGroups + sample description)
_PXD = "PXD043768"
_BASE_URL = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD043768"
_LV_PROTEINGROUPS_URL = f"{_BASE_URL}/Left_Ventricle_proteinGroups.txt"
_RV_PROTEINGROUPS_URL = f"{_BASE_URL}/Right_Ventricle_proteinGroups.txt"
_DESCRIPTION_URL = f"{_BASE_URL}/RAWs_Description.xlsx"

_CACHE_SUBDIR = "proteomics_validation"


def _cache_dir() -> Path:
    """Local cache directory for downloaded proteomics files."""
    base = Path(os.environ.get("LOCALAPPDATA",
                               Path.home() / ".cache")) / "KOSMIC" / _CACHE_SUBDIR
    base.mkdir(parents=True, exist_ok=True)
    return base


def _download_if_missing(url: str, filename: str,
                         progress_cb=None) -> Path:
    """Download a file to the cache if not already present."""
    target = _cache_dir() / filename
    if target.exists() and target.stat().st_size > 0:
        return target

    if progress_cb is not None:
        progress_cb(f"Downloading {filename}...")

    tmp = target.with_suffix(target.suffix + ".tmp")
    try:
        urllib.request.urlretrieve(url, tmp)
        tmp.replace(target)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise
    return target


def _parse_protein_groups(path: Path,
                          apply_hgnc: bool = True) -> pd.DataFrame:
    """Parse a MaxQuant proteinGroups.txt file.

    Returns
    -------
    DataFrame with columns:
        gene : str             -- HGNC-approved gene symbol from
                                  fasta header GN=, harmonised via the
                                  bundled HGNC table so it matches the
                                  symbols used in KOSMIC DE results.
        sample_<NN> : float    -- LFQ intensity per sample (one column
                                  per sample, NN matches the original
                                  '01'..'30' suffix)
    """
    keep = ['Protein IDs', 'Fasta headers',
            'Reverse', 'Potential contaminant',
            'Only identified by site']

    head = pd.read_csv(path, sep='\t', nrows=0)
    lfq_cols = [c for c in head.columns if c.startswith('LFQ intensity ')]
    keep_full = keep + lfq_cols

    df = pd.read_csv(path, sep='\t', usecols=keep_full, low_memory=False)

    df = df[(df['Reverse'] != '+') &
            (df['Potential contaminant'] != '+') &
            (df['Only identified by site'] != '+')]

    gn_pattern = re.compile(r'GN=(\S+)')

    def _gene(header):
        if not isinstance(header, str):
            return None
        m = gn_pattern.search(header)
        return m.group(1) if m else None

    df = df.assign(gene=df['Fasta headers'].apply(_gene))
    df = df[df['gene'].notna()]

    # Harmonise gene names against the bundled HGNC table so they match
    # the symbols used by KOSMIC's DE pipeline (which itself harmonises
    # AnnData gene names on import).
    if apply_hgnc:
        try:
            from kosmic.scrna.inspect.gene_names import load_hgnc_lookup
            lookup, _ = load_hgnc_lookup()
            df['gene'] = df['gene'].map(lambda g: lookup.get(g, g))
        except Exception:
            pass

    rename_map = {c: c.replace('LFQ intensity ', 'sample_')
                  for c in lfq_cols}
    df = df.rename(columns=rename_map)

    out = df[['gene'] + list(rename_map.values())]

    # If multiple proteinGroups rows now collapse to the same gene
    # (e.g. one isoform mapped to a previous symbol that resolves to
    # the same approved symbol as another row), aggregate by SUM of
    # LFQ intensity per gene.
    if out['gene'].duplicated().any():
        out = out.groupby('gene', as_index=False).sum(numeric_only=True)
    return out


def _parse_description(path: Path) -> pd.DataFrame:
    """Parse RAWs_Description.xlsx to map sample IDs to conditions.

    Returns
    -------
    DataFrame with columns: raw_number (int), ventricle (str), group (str)
    """
    df = pd.read_excel(path)

    # Left ventricle entries (cols Raw_number, Ventricle, Group)
    lv = df[['Raw_number', 'Ventricle', 'Group']].dropna()
    lv = lv.rename(columns={
        'Raw_number': 'raw_number',
        'Ventricle': 'ventricle',
        'Group': 'group',
    })

    # Right ventricle entries (cols Raw_number.1, Ventricle.1, Group.1)
    rv = df[['Raw_number.1', 'Ventricle.1', 'Group.1']].dropna()
    rv = rv.rename(columns={
        'Raw_number.1': 'raw_number',
        'Ventricle.1': 'ventricle',
        'Group.1': 'group',
    })

    out = pd.concat([lv, rv], ignore_index=True)
    out['raw_number'] = out['raw_number'].astype(int)
    return out


def fetch_and_parse_lv(progress_cb=None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Download (if needed) and parse the LV proteinGroups + sample map.

    Returns (protein_df, sample_df).
    """
    pg_path = _download_if_missing(
        _LV_PROTEINGROUPS_URL, "Left_Ventricle_proteinGroups.txt", progress_cb)
    desc_path = _download_if_missing(
        _DESCRIPTION_URL, "RAWs_Description.xlsx", progress_cb)

    if progress_cb is not None:
        progress_cb("Parsing proteinGroups...")
    proteins = _parse_protein_groups(pg_path)

    if progress_cb is not None:
        progress_cb("Parsing sample description...")
    samples = _parse_description(desc_path)

    return proteins, samples


def compute_protein_de(proteins: pd.DataFrame,
                       samples: pd.DataFrame,
                       ventricle: str = 'Left',
                       disease_groups=('RVD', 'noRVD'),
                       control_group: str = 'Control',
                       min_detection_frac: float = 1.0,
                       impute: bool = False,
                       median_normalize: bool = True) -> pd.DataFrame:
    """Run a Welch t-test per protein on log2-LFQ values.

    Parameters
    ----------
    proteins : DataFrame
        Output of `fetch_and_parse_lv()`.
    samples : DataFrame
        Sample-to-condition mapping (raw_number, ventricle, group).
    ventricle : str
    disease_groups : tuple[str]
    control_group : str
    min_detection_frac : float
        Minimum fraction of samples per arm in which the protein must
        be detected (LFQ > 0) for it to be included in the test.
        Default 0.7 means a protein needs detection in at least 70%
        of disease samples AND 70% of control samples.
    impute : bool
        If True, impute remaining missing log2-LFQ values with a small
        constant below the lowest detected value in that sample
        (left-censored MNAR imputation, MaxQuant/Perseus convention).

    Returns
    -------
    DataFrame with one row per protein passing the detection filter.
    """
    from scipy.stats import t as t_dist
    from kosmic.numerical import bh_fdr

    samples = samples[samples['ventricle'] == ventricle].copy()
    samples['sample_col'] = samples['raw_number'].apply(
        lambda n: f"sample_{int(n):02d}")

    disease_cols = samples[
        samples['group'].isin(disease_groups)]['sample_col'].tolist()
    control_cols = samples[
        samples['group'] == control_group]['sample_col'].tolist()
    disease_cols = [c for c in disease_cols if c in proteins.columns]
    control_cols = [c for c in control_cols if c in proteins.columns]
    if not disease_cols or not control_cols:
        raise ValueError(
            f"No matching sample columns found for ventricle={ventricle}")

    n_d_total = len(disease_cols)
    n_c_total = len(control_cols)

    d_mat = proteins[disease_cols].values.astype(float)
    c_mat = proteins[control_cols].values.astype(float)

    # Detection filter: require min_detection_frac in EACH arm
    n_d_detected = (d_mat > 0).sum(axis=1)
    n_c_detected = (c_mat > 0).sum(axis=1)
    keep_mask = ((n_d_detected >= min_detection_frac * n_d_total)
                 & (n_c_detected >= min_detection_frac * n_c_total))

    if not keep_mask.any():
        return pd.DataFrame(
            columns=['gene', 'log2FC', 'pval', 'fdr',
                     'n_disease', 'n_control'])

    proteins_kept = proteins[keep_mask].reset_index(drop=True)
    d_mat = d_mat[keep_mask]
    c_mat = c_mat[keep_mask]

    # log2 transform with NaN for zeros
    d_log = np.log2(np.where(d_mat > 0, d_mat, np.nan))
    c_log = np.log2(np.where(c_mat > 0, c_mat, np.nan))

    if median_normalize:
        # Per-sample median centring (subtract sample median from every
        # protein in that sample) so any small loading differences are
        # removed before computing log2FC.
        with np.errstate(invalid='ignore'):
            d_med = np.nanmedian(d_log, axis=0)
            c_med = np.nanmedian(c_log, axis=0)
        # Use the global grand median as reference so we don't introduce
        # an arbitrary shift between groups
        all_meds = np.concatenate([d_med, c_med])
        grand = float(np.nanmedian(all_meds[np.isfinite(all_meds)]))
        d_log = d_log - (d_med - grand)
        c_log = c_log - (c_med - grand)

    if impute:
        # Sample-min imputation: for each sample, find the minimum
        # detected log2-LFQ value, replace NaNs with (min - 1) to model
        # left-censored MNAR (Perseus default uses width=0.3, downshift=1.8;
        # we use a simpler 1-bit shift below the per-sample minimum).
        full = np.concatenate([d_log, c_log], axis=1)
        # Per-column (sample) min of detected values
        with np.errstate(invalid='ignore'):
            sample_min = np.nanmin(full, axis=0)
        sample_min = np.where(np.isfinite(sample_min), sample_min, 0.0)
        impute_val = sample_min - 1.0  # one log2 unit below the min

        # Impute disease columns
        for j in range(d_log.shape[1]):
            col = d_log[:, j]
            mask = ~np.isfinite(col)
            d_log[mask, j] = impute_val[j]
        # Impute control columns (sample_min indices are concatenated)
        for j in range(c_log.shape[1]):
            col = c_log[:, j]
            mask = ~np.isfinite(col)
            c_log[mask, j] = impute_val[d_log.shape[1] + j]

    n_d = (~np.isnan(d_log)).sum(axis=1)
    n_c = (~np.isnan(c_log)).sum(axis=1)

    mean_d = np.nanmean(d_log, axis=1)
    mean_c = np.nanmean(c_log, axis=1)
    var_d = np.nanvar(d_log, axis=1, ddof=1)
    var_c = np.nanvar(c_log, axis=1, ddof=1)

    log2fc = mean_d - mean_c

    se_sq = (np.where(n_d > 0, var_d / np.maximum(n_d, 1), 0) +
             np.where(n_c > 0, var_c / np.maximum(n_c, 1), 0))
    se = np.sqrt(np.clip(se_sq, 1e-20, None))

    valid = (n_d >= 2) & (n_c >= 2) & np.isfinite(log2fc) & (se > 0)

    t_stat = np.where(valid, log2fc / se, 0.0)

    with np.errstate(divide='ignore', invalid='ignore'):
        num = se_sq ** 2
        denom = ((var_d / np.maximum(n_d, 1)) ** 2
                 / np.maximum(n_d - 1, 1)
                 + (var_c / np.maximum(n_c, 1)) ** 2
                 / np.maximum(n_c - 1, 1))
        denom = np.where(denom > 0, denom, 1e-20)
        df = num / denom
        df = np.where(np.isfinite(df) & (df > 0), df, 1.0)

    pvals = np.ones(len(log2fc))
    pvals[valid] = 2 * t_dist.sf(np.abs(t_stat[valid]), df[valid])
    pvals = np.clip(pvals, 1e-300, 1.0)

    fdr = np.ones(len(pvals))
    if valid.any():
        fdr[valid] = bh_fdr(pvals[valid])

    return pd.DataFrame({
        'gene': proteins_kept['gene'].values,
        'log2FC': log2fc,
        'pval': pvals,
        'fdr': fdr,
        'n_disease': n_d,
        'n_control': n_c,
    })


def validate_consensus(consensus_df: pd.DataFrame,
                       protein_de: pd.DataFrame,
                       gene_col: str = 'names',
                       logfc_col: str = 'logfoldchanges',
                       fdr_col: str = 'fdr',
                       fdr_threshold: float = DEFAULT_FDR) -> dict:
    """Match a transcript consensus DataFrame to protein DE results.

    Parameters
    ----------
    consensus_df : DataFrame
        Gene-level consensus output (must have gene name, log2FC, FDR cols).
    protein_de : DataFrame
        Output of `compute_protein_de()`.
    gene_col, logfc_col, fdr_col : str
    fdr_threshold : float

    Returns
    -------
    dict with keys:
        merged : DataFrame   -- per-gene table joined on gene symbol
        n_consensus : int    -- # genes significant in transcript consensus
        n_overlap : int      -- # consensus genes with proteomic measurement
        n_concordant : int   -- # overlapping genes with same direction
        concordance_pct : float
        n_validated : int    -- # consensus genes also significant in
                                proteomics (FDR < threshold, same direction)
        validation_pct : float
        spearman_r : float   -- correlation across overlap (all directions)
        spearman_p : float
    """
    from scipy.stats import spearmanr

    transcript = consensus_df[[gene_col, logfc_col, fdr_col]].copy()
    transcript = transcript.rename(columns={
        gene_col: 'gene',
        logfc_col: 'transcript_log2FC',
        fdr_col: 'transcript_fdr',
    })
    transcript['gene'] = transcript['gene'].astype(str).str.upper()
    # Harmonise the consensus gene names too, in case any were aliases
    try:
        from kosmic.scrna.inspect.gene_names import load_hgnc_lookup
        lookup, _ = load_hgnc_lookup()
        # Lookup is case-sensitive on the HGNC table but symbols are
        # uppercase already after .str.upper()
        transcript['gene'] = transcript['gene'].map(
            lambda g: lookup.get(g, g))
    except Exception:
        pass
    transcript = transcript.drop_duplicates(subset='gene')

    protein = protein_de.rename(columns={
        'log2FC': 'protein_log2FC',
        'pval': 'protein_pval',
        'fdr': 'protein_fdr',
    })[['gene', 'protein_log2FC', 'protein_pval', 'protein_fdr']].copy()
    protein['gene'] = protein['gene'].astype(str).str.upper()
    protein = protein.drop_duplicates(subset='gene')

    merged = transcript.merge(protein, on='gene', how='left')

    # Same direction = product of log2FCs > 0
    merged['same_direction'] = (
        np.sign(merged['transcript_log2FC'])
        == np.sign(merged['protein_log2FC']))
    merged.loc[merged['protein_log2FC'].isna(), 'same_direction'] = False

    # Subsets
    consensus_mask = merged['transcript_fdr'] < fdr_threshold
    overlap_mask = consensus_mask & merged['protein_log2FC'].notna()
    concordant_mask = overlap_mask & merged['same_direction']
    validated_mask = (
        concordant_mask
        & (merged['protein_fdr'] < fdr_threshold)
    )

    n_consensus = int(consensus_mask.sum())
    n_overlap = int(overlap_mask.sum())
    n_concordant = int(concordant_mask.sum())
    n_validated = int(validated_mask.sum())

    # Spearman across all genes that have both measurements
    full_overlap = merged.dropna(
        subset=['transcript_log2FC', 'protein_log2FC'])
    if len(full_overlap) >= 3:
        rho, sp = spearmanr(
            full_overlap['transcript_log2FC'],
            full_overlap['protein_log2FC'])
        rho = float(rho) if np.isfinite(rho) else float('nan')
        sp = float(sp) if np.isfinite(sp) else float('nan')
    else:
        rho, sp = float('nan'), float('nan')

    # Stratified concordance: among consensus genes that ALSO have a
    # meaningful protein effect, what fraction agree in direction?
    strat = {}
    for prot_thr in (0.1, 0.2, 0.3, 0.5):
        meaningful_mask = (
            overlap_mask
            & (merged['protein_log2FC'].abs() > prot_thr))
        n_meaningful = int(meaningful_mask.sum())
        n_meaningful_concordant = int(
            (meaningful_mask & merged['same_direction']).sum())
        strat[f'|protLFC|>{prot_thr}'] = {
            'n': n_meaningful,
            'n_concordant': n_meaningful_concordant,
            'pct': (100.0 * n_meaningful_concordant / n_meaningful
                    if n_meaningful else 0.0),
        }

    # Concordance among consensus genes that pass protein FDR<0.05 too
    sig_proteins_mask = (
        overlap_mask
        & (merged['protein_fdr'] < fdr_threshold))
    n_sig_overlap = int(sig_proteins_mask.sum())
    n_sig_concordant = int(
        (sig_proteins_mask & merged['same_direction']).sum())
    strat['protein_FDR<0.05'] = {
        'n': n_sig_overlap,
        'n_concordant': n_sig_concordant,
        'pct': (100.0 * n_sig_concordant / n_sig_overlap
                if n_sig_overlap else 0.0),
    }

    return {
        'merged': merged,
        'n_consensus': n_consensus,
        'n_overlap': n_overlap,
        'n_concordant': n_concordant,
        'concordance_pct': (
            100.0 * n_concordant / n_overlap if n_overlap else 0.0),
        'n_validated': n_validated,
        'validation_pct': (
            100.0 * n_validated / n_overlap if n_overlap else 0.0),
        'spearman_r': rho,
        'spearman_p': sp,
        'fdr_threshold': fdr_threshold,
        'stratified': strat,
    }


def run_validation(consensus_df: pd.DataFrame,
                   ventricle: str = 'Left',
                   progress_cb=None,
                   **validate_kwargs) -> dict:
    """End-to-end: download (if needed) -> protein DE -> match to consensus."""
    proteins, samples = fetch_and_parse_lv(progress_cb=progress_cb)
    if progress_cb is not None:
        progress_cb(f"Computing protein DE ({ventricle} ventricle)...")
    protein_de = compute_protein_de(proteins, samples, ventricle=ventricle)
    if progress_cb is not None:
        progress_cb("Matching to consensus...")
    return validate_consensus(consensus_df, protein_de, **validate_kwargs)
