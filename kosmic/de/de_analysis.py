# Differential Expression Analysis
# Pseudobulk DE analysis for scRNA-seq data.
# Supports Welch's t-test (with optional EB moderation) and DESeq2 methods.

import numpy as np
import pandas as pd

from kosmic.scrna.inspect.detection import detect_species, format_gene_for_species
from kosmic import (
    DEFAULT_FDR, DEFAULT_LFC_THRESHOLD,
    PATHWAY_GENE_DETECTION_PCT, PATHWAY_MIN_GENES, PATHWAY_MIN_COVERAGE,
    DE_MIN_CELLS, DE_MIN_EXPRESSING_SAMPLES, DE_MAX_SE, DE_SE_FLOOR,
    DE_FILTER_MIN_COUNT, DE_FILTER_MIN_TOTAL_COUNT, DE_FILTER_MIN_SAMPLES,
    DE_DETECTION_MIN_PCT, DE_DETECTION_ON, DE_DETECTION_MIN_DONOR_FRAC,
    DE_DESEQ2_INDEPENDENT_FILTER, DE_DESEQ2_COOKS_FILTER,
)

# Default detection pre-filter fraction; 0 disables (config 'de.detection_on').
DE_DETECTION_DEFAULT = DE_DETECTION_MIN_PCT if DE_DETECTION_ON else 0.0


def prepare_gene_coverage(pathway_gene_sets, var_names, species='human'):
    """Check which pathway genes are available in the dataset.

    Parameters
    ----------
    pathway_gene_sets : dict
        {pathway_name: [gene_list]}.
    var_names : set or list
        Available gene names in the dataset.
    species : str
        'human' or 'mouse' for gene name conversion.

    Returns
    -------
    tuple of (metabolic_genes, pathway_coverage)
        metabolic_genes: list of unique available genes (species-converted).
        pathway_coverage: dict with coverage stats per pathway.
    """
    var_names = set(var_names)
    is_mouse = species == 'mouse'

    all_genes = []
    pathway_coverage = {}

    for pathway_name, genes in pathway_gene_sets.items():
        if is_mouse:
            genes_fmt = [format_gene_for_species(g, 'mouse') for g in genes]
        else:
            genes_fmt = genes

        available = [g for g in genes_fmt if g in var_names]
        coverage_pct = len(available) / len(genes) * 100 if genes else 0

        pathway_coverage[pathway_name] = {
            'total_genes': len(genes),
            'available_genes': len(available),
            'coverage_pct': coverage_pct,
            'genes': available,
        }
        all_genes.extend(available)

    # Deduplicate preserving order
    metabolic_genes = list(dict.fromkeys(all_genes))
    return metabolic_genes, pathway_coverage


def create_pseudobulk(adata, genes, sample_col, condition_col, min_cells=DE_MIN_CELLS,
                      aggregate='mean', counts_layer=None, covariates=None):
    """Create pseudobulk expression profiles by sample.

    Uses raw counts (adata.raw if available) without per-cell normalisation.

    Parameters
    ----------
    adata : anndata.AnnData
        Input data. Uses adata.raw if available for raw counts.
    genes : list of str
        Genes to include.
    sample_col, condition_col : str
        Columns in obs for sample grouping and condition labels.
    min_cells : int
        Minimum cells per sample to include.
    aggregate : 'mean' | 'sum'
        'mean' for t-test, 'sum' for DESeq2.
    counts_layer : str, optional
        Name of a counts layer to use as the count source (e.g.
        'decontX_counts' for DecontX-corrected counts). When given, that
        layer -- keyed on 'adata.var_names' -- overrides the default
        counts / raw / X resolution. Genes absent from the layer are
        dropped.
    covariates : sequence of str, optional
        Extra 'obs' columns to carry through onto 'sample_df', one value
        per sample, for use in the DE design. 'study' is the one this
        exists for: pooling several cohorts into a single model needs the
        cohort offset in the design, or it lands in the residual and
        inflates dispersion. A column that varies *within* a sample is
        not a sample-level property and is skipped, since taking the
        first cell's value would silently invent one.

    Returns
    -------
    pseudobulk_matrix : np.ndarray
        Shape '(n_samples, n_genes)'.
    sample_df : pd.DataFrame
        Columns: 'sample', 'condition', 'n_cells'. Adds a 'role'
        column ('disease' / 'control' / 'exclude') when
        'adata.obs._role' is present, so DE / CC perm consumers
        read roles directly without re-resolving from condition labels.
    genes_used : list of str
        Gene names in the matrix (subset of 'genes' present in adata).
    """
    # Prefer layers['counts'], then adata.raw, then adata.X (raw retains
    # the full gene set if HVG filtering trimmed adata.var_names).
    all_var = set(adata.raw.var_names) if adata.raw is not None else set(adata.var_names)
    genes_used = [g for g in genes if g in all_var]
    if not genes_used:
        return np.array([]), pd.DataFrame(), []

    if counts_layer is not None:
        # Explicit count source (e.g. DecontX-corrected counts). The layer
        # is keyed on the current var_names, so restrict to genes present
        # there rather than the raw superset.
        layer_genes = [g for g in genes_used if g in adata.var_names]
        if not layer_genes:
            return np.array([]), pd.DataFrame(), []
        if counts_layer not in adata.layers:
            raise ValueError(
                f"Counts layer '{counts_layer}' not found in adata.layers.")
        adata_for_de = adata[:, layer_genes].copy()
        adata_for_de.X = adata_for_de.layers[counts_layer].copy()
        genes_used = layer_genes
    elif 'counts' in adata.layers and all(g in adata.var_names for g in genes_used):
        # True raw counts stored before normalisation -- best source.
        adata_for_de = adata[:, genes_used].copy()
        adata_for_de.X = adata_for_de.layers['counts'].copy()
    elif adata.raw is not None:
        raw_genes = [g for g in genes_used if g in adata.raw.var_names]
        if raw_genes:
            adata_for_de = adata.raw[:, raw_genes].to_adata()
            genes_used = raw_genes
        else:
            adata_for_de = adata[:, genes_used].copy()
    else:
        adata_for_de = adata[:, genes_used].copy()

    pseudobulk_data = []
    sample_metadata = []

    # All cells of one sample share the same role -- take the first cell's _role.
    has_role = '_role' in adata_for_de.obs.columns

    samples = adata_for_de.obs[sample_col].unique()

    for sample_id in samples:
        mask = adata_for_de.obs[sample_col] == sample_id
        sample_cells = adata_for_de[mask]

        if sample_cells.n_obs < min_cells:
            continue

        if hasattr(sample_cells.X, 'toarray'):
            X = sample_cells.X.toarray()
        else:
            X = np.asarray(sample_cells.X)

        if aggregate == 'sum':
            expression = np.sum(X, axis=0)
        else:
            expression = np.mean(X, axis=0)

        expression = np.asarray(expression).flatten()
        pseudobulk_data.append(expression)

        meta = {
            'sample': sample_id,
            'condition': sample_cells.obs[condition_col].iloc[0],
            'n_cells': sample_cells.n_obs,
        }
        if has_role:
            meta['role'] = str(sample_cells.obs['_role'].iloc[0])
        for col in (covariates or ()):
            if col in sample_cells.obs.columns:
                values = sample_cells.obs[col].astype(str)
                # One value per sample or nothing: a covariate that
                # differs between a sample's own cells is not a
                # sample-level property.
                if values.nunique() == 1:
                    meta[col] = str(values.iloc[0])
        sample_metadata.append(meta)

    if not pseudobulk_data:
        return np.array([]), pd.DataFrame(), genes_used

    pseudobulk_matrix = np.array(pseudobulk_data)
    sample_df = pd.DataFrame(sample_metadata)

    return pseudobulk_matrix, sample_df, genes_used


def _pathway_retention(adata, pathway_gene_sets, gene_detection_pct,
                       min_coverage, min_genes, counts_layer=None):
    """Gate pathways by gene detection + coverage.

    Two passes: species-format each pathway's genes and keep those present in
    the count matrix, then drop genes not detected in >= 'gene_detection_pct'
    of cells (disease OR control). A pathway is retained only if at least
    'min_genes' of its annotated genes survive and they cover >= 'min_coverage'
    of the annotated set.

    Returns
    -------
    (pathway_gene_map, report_rows, base_var_names)
        pathway_gene_map : {pathway: [retained_genes]} for pathways that pass.
        report_rows : one dict per pathway (kept and dropped) with keys
            'pathway', 'annotated', 'present', 'detected', 'coverage', 'status'.
        base_var_names : the gene universe used, for downstream normalisation.
    """
    from kosmic.scrna.inspect.detection import detect_species, format_gene_for_species

    # With an explicit counts layer the gene universe is the current var_names
    # (where the layer lives), not the raw superset.
    if counts_layer is not None:
        base_var_names = list(adata.var_names)
    else:
        base_var_names = list(adata.raw.var_names if adata.raw is not None else adata.var_names)
    var_names = set(base_var_names)
    species = detect_species(list(var_names))

    # Pass 1: species-format + presence in the count matrix.
    pathway_present = {}  # pathway -> (n_annotated, [present_genes])
    candidate_genes = set()
    for pw_name, genes in pathway_gene_sets.items():
        genes_fmt = [format_gene_for_species(g, species) for g in genes]
        present = [g for g in genes_fmt if g in var_names]
        pathway_present[pw_name] = (len(genes), present)
        candidate_genes.update(present)

    # Pass 2: detection filter -- keep genes detected in >= gene_detection_pct
    # of cells in disease OR control (the same "% of cells" metric as gene DE),
    # so a pathway score reflects the biology measurable in this cell type.
    if gene_detection_pct and gene_detection_pct > 0 and candidate_genes:
        retained_set = set(genes_passing_detection(
            adata, list(candidate_genes), gene_detection_pct,
            counts_layer=counts_layer))
    else:
        retained_set = set(candidate_genes)

    # Coverage / min-genes gate on RETAINED genes; denominator is the full
    # annotated set.
    pathway_gene_map = {}
    report_rows = []
    for pw_name, (n_annotated, present) in pathway_present.items():
        retained = [g for g in present if g in retained_set]
        coverage = len(retained) / n_annotated if n_annotated else 0.0
        if len(retained) < min_genes:
            status = f"dropped: {len(retained)} detected < {min_genes} required"
        elif coverage < min_coverage:
            status = f"dropped: coverage {coverage:.0%} < {min_coverage:.0%}"
        else:
            status = "scored"
            pathway_gene_map[pw_name] = retained
        report_rows.append({
            'pathway': pw_name,
            'annotated': n_annotated,
            'present': len(present),
            'detected': len(retained),
            'coverage': coverage,
            'status': status,
        })
    return pathway_gene_map, report_rows, base_var_names


def pathway_coverage_report(adata, pathway_gene_sets, *,
                            gene_detection_pct=PATHWAY_GENE_DETECTION_PCT,
                            min_coverage=PATHWAY_MIN_COVERAGE,
                            min_genes=PATHWAY_MIN_GENES, counts_layer=None):
    """Per-pathway detection/coverage table -- which pathways score, and why.

    One row per pathway: annotated / present / detected gene counts, coverage
    fraction (detected / annotated), and a status ('scored' or the reason it
    was dropped). Mirrors exactly the gate 'compute_pathway_scores' applies.
    """
    _, report_rows, _ = _pathway_retention(
        adata, pathway_gene_sets, gene_detection_pct, min_coverage, min_genes,
        counts_layer=counts_layer)
    return pd.DataFrame(report_rows, columns=[
        'pathway', 'annotated', 'present', 'detected', 'coverage', 'status'])


def compute_pathway_scores(adata, pathway_gene_sets, sample_col, condition_col,
                           min_cells=DE_MIN_CELLS,
                           gene_detection_pct=PATHWAY_GENE_DETECTION_PCT,
                           min_coverage=PATHWAY_MIN_COVERAGE,
                           min_genes=PATHWAY_MIN_GENES,
                           normalization='cpm', counts_layer=None,
                           standardize_genes=False):
    """Compute per-donor pathway scores from raw counts.

    For each pathway, builds pseudobulk profiles then averages (equal
    weight) across the pathway's genes per donor. Returns a
    (samples x pathways) matrix passable to run_pseudobulk_de().

    'normalization' sets the per-donor transform:
    'cpm' (full-genome library size, log2), 'deseq2' (median-of-ratios
    size factors then log2 -- the same estimator the DESeq2 DE uses), or
    'vst' (DESeq2 variance-stabilising transform via pydeseq2, fit on the
    full-genome pseudobulk). All three yield a donor-level score whose
    disease-vs-control difference is a pathway log2FC.

    Parameters
    ----------
    adata : anndata.AnnData
        Input data.
    pathway_gene_sets : dict
        {pathway_name: [gene_list]}.
    sample_col, condition_col : str
        Metadata columns.
    min_cells : int
        Minimum cells per sample.
    gene_detection_pct : float
        A pathway gene is retained only if detected (non-zero count) in at
        least this fraction of cells in disease OR control. 0 disables the
        detection filter. Defaults to 'PATHWAY_GENE_DETECTION_PCT'. This is
        the same "% of cells" metric as the gene-level DE detection filter.
    min_coverage : float
        Minimum fraction of a pathway's annotated genes that must be retained
        (present AND detected) for the pathway to be scored (0-1). Defaults
        to 'PATHWAY_MIN_COVERAGE'.
    min_genes : int
        Minimum absolute number of retained pathway genes required.
        Defaults to 'PATHWAY_MIN_GENES'.
    counts_layer : str, optional
        Counts layer to score (e.g. 'decontX_counts'); when set, that
        layer -- keyed on 'adata.var_names' -- is the count source and the
        gene universe, in place of the counts / raw / X default.
    standardize_genes : bool
        Z-score each gene across donors before averaging the module, so every
        subunit contributes equally and high-abundance genes cannot dominate.
        The score is then in standardised units (a robustness view), not log2
        expression -- its group difference is not a log2FC.

    Returns
    -------
    score_matrix : np.ndarray
        (n_samples, n_pathways).
    sample_df : pd.DataFrame
        Sample metadata with 'sample', 'condition', 'n_cells'.
    pathway_names : list of str
        Pathway names matching columns of score_matrix.
    """
    pathway_gene_map, _, base_var_names = _pathway_retention(
        adata, pathway_gene_sets, gene_detection_pct, min_coverage, min_genes,
        counts_layer=counts_layer)

    if not pathway_gene_map:
        return np.array([]), pd.DataFrame(), []

    # Deduplicate the retained-gene union while preserving order.
    all_genes = list(dict.fromkeys(
        g for genes in pathway_gene_map.values() for g in genes))

    # Create pseudobulk (sum counts) for pathway genes
    pb_matrix, sample_df, genes_used = create_pseudobulk(
        adata, all_genes, sample_col, condition_col,
        min_cells=min_cells, aggregate='sum', counts_layer=counts_layer,
    )

    if pb_matrix.size == 0:
        return np.array([]), pd.DataFrame(), []

    # Full-genome pseudobulk drives normalisation (library sizes / size
    # factors / VST dispersion trend).
    pb_full, _, full_genes_used = create_pseudobulk(
        adata, base_var_names, sample_col, condition_col,
        min_cells=min_cells, aggregate='sum', counts_layer=counts_layer,
    )

    # Each branch yields a per-donor gene matrix 'values' plus the gene
    # index into its columns; the pathway loop is shared.
    if normalization == 'vst':
        values = _vst_transform(pb_full, full_genes_used, sample_df, condition_col)
        gene_to_idx = {g: i for i, g in enumerate(full_genes_used)}
    elif normalization == 'deseq2':
        # Median-of-ratios size factors, matching the gene-level DESeq2 DE.
        from kosmic.meta_analysis.deseq2_fast import _compute_size_factors
        size_factors = _compute_size_factors(pb_full).reshape(-1, 1)
        size_factors = np.clip(size_factors, 1e-8, None)
        values = np.log2(pb_matrix / size_factors + 1)
        gene_to_idx = {g: i for i, g in enumerate(genes_used)}
    else:
        # CPM-normalise pathway genes using full-genome library sizes.
        lib_sizes = pb_full.sum(axis=1, keepdims=True)
        lib_sizes = np.clip(lib_sizes, 1, None)
        values = np.log2(pb_matrix / lib_sizes * 1e6 + 1)
        gene_to_idx = {g: i for i, g in enumerate(genes_used)}

    # Optionally z-score each gene across donors so every subunit contributes
    # equally and high-abundance genes cannot dominate the module average. The
    # score is then standardised, not log2 -- its group difference is not a
    # log2FC (see the 'deseq2_vst_zscore' scoring method).
    if standardize_genes:
        mu = values.mean(axis=0, keepdims=True)
        sd = values.std(axis=0, ddof=0, keepdims=True)
        sd = np.where(sd == 0, 1.0, sd)
        values = (values - mu) / sd

    # Equal-weight mean of the pathway's genes per donor.
    pathway_names = []
    score_columns = []
    for pw_name, pw_genes in pathway_gene_map.items():
        col_indices = [gene_to_idx[g] for g in pw_genes if g in gene_to_idx]
        if not col_indices:
            continue
        pw_score = values[:, col_indices].mean(axis=1)
        score_columns.append(pw_score)
        pathway_names.append(pw_name)

    if not score_columns:
        return np.array([]), pd.DataFrame(), []

    score_matrix = np.column_stack(score_columns)
    return score_matrix, sample_df, pathway_names


def _vst_transform(pb_full, gene_names, sample_df, condition_col):
    """Blind DESeq2 variance-stabilising transform of a pseudobulk matrix.

    Returns a (donors x genes) array on a log2-like scale. The fit is blind
    to the design ('use_design=False'), so the transform does not depend on
    the disease/control contrast it is later tested on.
    """
    from pydeseq2.dds import DeseqDataSet
    counts_df = pd.DataFrame(
        np.rint(pb_full).astype(int),
        index=[str(s) for s in sample_df['sample'].values],
        columns=list(gene_names))
    cond = [str(c) for c in sample_df['condition'].values]
    if len(set(cond)) < 2:
        # Blind VST ignores the design; supply a valid two-level factor so
        # the dataset builds even for a single-condition edge case.
        cond = ['a' if i % 2 == 0 else 'b' for i in range(len(cond))]
    metadata = pd.DataFrame({'condition': cond}, index=counts_df.index)
    dds = DeseqDataSet(counts=counts_df, metadata=metadata,
                       design='~condition', quiet=True)
    dds.vst_fit(use_design=False)
    return np.asarray(dds.vst_transform())


def pseudobulk_expression_matrix(adata, sample_col, condition_col,
                                 normalization='vst', counts_layer=None,
                                 min_cells=DE_MIN_CELLS):
    """Per-donor normalised expression matrix (donors x genes) for plots.

    Pseudobulk-sums each donor, then normalises to match the DE method the
    values are shown alongside, so a heatmap or expression bar reflects the
    same quantity the test used rather than a mismatched transform:

      - 'vst'    : DESeq2 blind variance-stabilising transform (DESeq2 engine).
      - 'deseq2' : median-of-ratios size-factor log2 (the DESeq2 DE basis).
      - 'cpm'    : log2 CPM on full-genome library sizes (the Welch engine).

    Donors whose '_role' is 'exclude' are left out. These values are
    shown beside DE results computed from the disease/control contrast,
    so carrying an excluded arm into the plot puts donors on the chart
    that contributed nothing to the numbers next to it -- a third group
    in the expression bars, extra columns in the heatmap.

    Returns a DataFrame of donor-level values with leading 'sample' and
    'condition' columns; empty if no donor clears 'min_cells'.
    """
    if counts_layer is not None:
        base_var_names = list(adata.var_names)
    else:
        base_var_names = list(
            adata.raw.var_names if adata.raw is not None else adata.var_names)
    pb_full, sample_df, genes_used = create_pseudobulk(
        adata, base_var_names, sample_col, condition_col,
        min_cells=min_cells, aggregate='sum', counts_layer=counts_layer)
    if pb_full.size == 0 or not genes_used:
        return pd.DataFrame()

    if normalization == 'vst':
        values = _vst_transform(pb_full, genes_used, sample_df, condition_col)
    elif normalization == 'deseq2':
        from kosmic.meta_analysis.deseq2_fast import _compute_size_factors
        sf = np.clip(_compute_size_factors(pb_full).reshape(-1, 1), 1e-8, None)
        values = np.log2(pb_full / sf + 1)
    else:  # 'cpm'
        lib = np.clip(pb_full.sum(axis=1, keepdims=True), 1, None)
        values = np.log2(pb_full / lib * 1e6 + 1)

    df = pd.DataFrame(np.asarray(values), columns=list(genes_used))
    df.insert(0, 'condition', sample_df['condition'].values)
    df.insert(0, 'sample', sample_df['sample'].values)

    # Drop the excluded arm. Roles are per donor, so this is exact.
    if 'role' in sample_df.columns:
        keep = sample_df['role'].astype(str).str.lower().isin(
            ('disease', 'control')).values
        if keep.any():
            df = df[keep].reset_index(drop=True)
    return df


def run_pseudobulk_de(pseudobulk_matrix, sample_df, gene_names,
                      min_expressing_samples=DE_MIN_EXPRESSING_SAMPLES, method='ttest',
                      max_se=DE_MAX_SE, moderate=True,
                      pre_transformed=False, normalize='cpm', se_floor=DE_SE_FLOOR,
                      filter_min_count=DE_FILTER_MIN_COUNT,
                      filter_min_samples=DE_FILTER_MIN_SAMPLES,
                      study_col=None):
    """Run pseudobulk differential expression analysis (Welch's t-test on log2-CPM).

    Sample filtering is by sample_df['role'] ('disease' / 'control' /
    'exclude'), populated by create_pseudobulk from adata.obs['_role'].

    Parameters
    ----------
    pseudobulk_matrix : numpy.ndarray
        '(n_samples, n_genes)'. Sum counts per sample, unless
        'pre_transformed=True' in which case already log2-CPM.
    sample_df : pandas.DataFrame
        Must have a 'role' column.
    gene_names : list of str
        Gene names matching matrix columns.
    min_expressing_samples : int
        Min samples with expression > 0 to test a gene.
    method : str
        'ttest' (Welch's t-test).
    max_se : float
        Standard-error cap.
    moderate : bool
        Apply limma-style empirical Bayes variance moderation.
    pre_transformed : bool
        Set True when the matrix is already on the log2 scale.
    normalize : 'cpm' | 'log2'
        'cpm': CPM-normalise then log2(CPM+1) (for sum-count input).
        'log2': log2(x+1) only (for mean-count input).

    Returns
    -------
    pandas.DataFrame
        Columns: names, logfoldchanges, pvals, pvals_adj, disease_mean,
        control_mean, disease_var, control_var, disease_samples,
        control_samples, cohens_d, se, ci_lower, ci_upper, se_diff,
        statistic.
    """
    from scipy.stats import ttest_ind
    from kosmic.numerical import bh_fdr

    # Normalise + log2-transform per sample.
    if pre_transformed:
        log2_cpm = pseudobulk_matrix
    elif normalize == 'log2':
        log2_cpm = np.log2(pseudobulk_matrix + 1)
    else:
        lib_sizes = pseudobulk_matrix.sum(axis=1, keepdims=True)
        lib_sizes = np.clip(lib_sizes, 1, None)
        cpm_matrix = pseudobulk_matrix / lib_sizes * 1e6
        log2_cpm = np.log2(cpm_matrix + 1)

    if 'role' not in sample_df.columns:
        raise ValueError(
            "sample_df is missing the 'role' column. Set roles in the "
            "Inspect tab (Sample Setup) so adata.obs['_role'] is "
            "populated, then re-create the pseudobulk.")
    disease_mask = (sample_df['role'] == 'disease').values
    control_mask = (sample_df['role'] == 'control').values

    # Sample-level gene filter. Raw counts use the expression-prevalence rule
    # (donor-aware, library-size-adjusted); pre-transformed inputs (e.g. pathway
    # scores) are not counts, so they keep a simple presence threshold.
    if pre_transformed:
        genes_to_test = []
        for idx, gene in enumerate(gene_names):
            n_expr_d = (pseudobulk_matrix[disease_mask, idx] > 0).sum()
            n_expr_c = (pseudobulk_matrix[control_mask, idx] > 0).sum()
            if (n_expr_d >= min_expressing_samples
                    or n_expr_c >= min_expressing_samples):
                genes_to_test.append((idx, gene))
    else:
        keep = filter_by_expression(
            pseudobulk_matrix, sample_df['role'].astype(str).values,
            min_count=filter_min_count, min_samples=filter_min_samples,
            studies=_study_values(sample_df, study_col))
        genes_to_test = [(idx, gene) for idx, gene in enumerate(gene_names)
                         if keep[idx]]

    if not genes_to_test:
        return pd.DataFrame()

    # Run DE tests on the log2-CPM matrix.
    results = []
    for gene_idx, gene_name in genes_to_test:
        disease_vals = log2_cpm[disease_mask, gene_idx]
        control_vals = log2_cpm[control_mask, gene_idx]

        if len(disease_vals) < 2 or len(control_vals) < 2:
            continue

        disease_mean = np.mean(disease_vals)
        control_mean = np.mean(control_vals)

        fold_change = disease_mean - control_mean  # log2-CPM mean difference

        # Welch's t-test on log2-CPM values.
        try:
            statistic, pvalue = ttest_ind(disease_vals, control_vals, equal_var=False)
        except ValueError:
            pvalue, statistic = 1.0, 0.0  # n<2 or both groups zero-variance

        disease_var = np.var(disease_vals, ddof=1)
        control_var = np.var(control_vals, ddof=1)

        nd, nc = len(disease_vals), len(control_vals)
        pooled_std = np.sqrt(
            ((nd - 1) * disease_var + (nc - 1) * control_var) /
            max(nd + nc - 2, 1)
        )

        cohens_d = fold_change / pooled_std if pooled_std > 0 else 0.0

        # SE directly from t-test on log2-CPM (no delta method needed).
        se_logfc = np.sqrt(disease_var / nd + control_var / nc)
        if se_logfc > max_se or not np.isfinite(se_logfc):
            se_logfc = np.sqrt(1 / max(nd, 1) + 1 / max(nc, 1))
        se_logfc = np.clip(se_logfc, se_floor, max_se)

        ci_lower = fold_change - 1.96 * se_logfc
        ci_upper = fold_change + 1.96 * se_logfc
        se_diff = pooled_std * np.sqrt(1 / nd + 1 / nc) if pooled_std > 0 else 0

        results.append({
            'names': gene_name,
            'logfoldchanges': fold_change,
            'pvals': pvalue,
            'disease_mean': disease_mean,
            'control_mean': control_mean,
            'disease_var': disease_var,
            'control_var': control_var,
            'disease_samples': nd,
            'control_samples': nc,
            'cohens_d': cohens_d,
            'se': se_logfc,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'se_diff': se_diff,
            'statistic': statistic,
        })

    if not results:
        return pd.DataFrame()

    de_results = pd.DataFrame(results)

    # Unmoderated SE preserved for meta-analysis: moderated SEs are
    # deflated and cause overconfident pooling (same issue as DESeq2 shrinkage).
    unshrunk_se = de_results['se'].copy()

    # --- Empirical Bayes variance moderation (limma-style) ---
    if moderate:
        from scipy.stats import t as t_dist

        dv = de_results['disease_var'].values
        cv = de_results['control_var'].values
        nd = de_results['disease_samples'].values
        nc = de_results['control_samples'].values
        df_gene = nd + nc - 2

        # Pooled residual variance (sum-of-squares would inflate the prior).
        s2 = ((nd - 1) * dv + (nc - 1) * cv) / np.maximum(df_gene, 1)

        s2_pos = s2[s2 > 0]
        if len(s2_pos) >= 3:
            log_s2 = np.log(s2_pos)
            s0_sq = np.exp(np.median(log_s2))
            var_log_s2 = np.var(log_s2, ddof=1)
            from scipy.special import polygamma
            mean_trigamma = np.mean([polygamma(1, d / 2) for d in df_gene[s2 > 0]])
            trigamma_d0_half = max(var_log_s2 - mean_trigamma, 0.01)
            d0 = max(2.0 / trigamma_d0_half, 1.0)
            d0 = min(d0, 50.0)  # cap to prevent over-shrinkage

            s2_mod = (df_gene * s2 + d0 * s0_sq) / (df_gene + d0)
            df_mod = df_gene + d0

            for i in range(len(de_results)):
                if s2[i] <= 0 or nd[i] < 2 or nc[i] < 2:
                    continue
                se_pooled = np.sqrt(s2_mod[i] * (1.0 / nd[i] + 1.0 / nc[i]))
                if se_pooled > 0:
                    d_mean = de_results.iloc[i]['disease_mean']
                    c_mean = de_results.iloc[i]['control_mean']
                    t_mod = (d_mean - c_mean) / se_pooled
                    p_mod = 2.0 * t_dist.sf(abs(t_mod), df_mod[i])
                    de_results.at[de_results.index[i], 'pvals'] = p_mod
                    de_results.at[de_results.index[i], 'statistic'] = t_mod

                    # Update SE of logFC using moderation ratio
                    var_ratio = s2_mod[i] / s2[i]
                    se_orig = de_results.iloc[i]['se']
                    se_mod = se_orig * np.sqrt(var_ratio)
                    se_mod = np.clip(se_mod, 0.01, max_se)
                    lfc = de_results.iloc[i]['logfoldchanges']
                    de_results.at[de_results.index[i], 'se'] = se_mod
                    de_results.at[de_results.index[i], 'ci_lower'] = lfc - 1.96 * se_mod
                    de_results.at[de_results.index[i], 'ci_upper'] = lfc + 1.96 * se_mod

    # 'se' for meta-analysis stays unmoderated; moderated kept as a separate column.
    de_results['se_moderated'] = de_results['se'].copy()
    de_results['se'] = unshrunk_se

    # FDR correction
    de_results['pvals_adj'] = bh_fdr(de_results['pvals'])

    return de_results


def annotate_de_results(de_results, pathway_coverage):
    """Add derived columns and pathway annotations to DE results.

    Parameters
    ----------
    de_results : pandas.DataFrame
        Raw DE results from run_pseudobulk_de.
    pathway_coverage : dict
        From prepare_gene_coverage.

    Returns
    -------
    tuple of (de_results, significant_genes)
        de_results: annotated DataFrame with abs_logfoldchange, -log10_pval, pathways.
        significant_genes: rows with FDR<0.05 whose estimated
            |logFC| also exceeds 0.25. The FDR applies to the first
            condition only -- the fold-change cut is applied to the
            point estimate and carries no error rate of its own.
    """
    if de_results.empty:
        return de_results, pd.DataFrame()

    from kosmic.numerical import neg_log10
    de_results = de_results.copy()
    de_results['abs_logfoldchange'] = de_results['logfoldchanges'].abs()
    de_results['-log10_pval'] = neg_log10(de_results['pvals_adj'])

    # Map genes to pathways.
    pathway_annotation = {}
    for pathway_name, coverage in pathway_coverage.items():
        for gene in coverage['genes']:
            if gene in pathway_annotation:
                pathway_annotation[gene] += f";{pathway_name}"
            else:
                pathway_annotation[gene] = pathway_name

    de_results['pathways'] = de_results['names'].map(pathway_annotation).fillna('')

    significant_genes = de_results[
        (de_results['pvals_adj'] < DEFAULT_FDR) &
        (de_results['abs_logfoldchange'] > DEFAULT_LFC_THRESHOLD)
    ].copy()

    return de_results, significant_genes


def significant_subset(de_results, fdr=DEFAULT_FDR,
                       lfc=DEFAULT_LFC_THRESHOLD):
    """Rows of a DE-results frame passing both thresholds.

    Used for the significant-genes CSV every run writes, so the in-app
    count and the exported file are the same selection rather than two
    reimplementations of it. Prefers the precomputed
    'abs_logfoldchange' column, falling back to the absolute value of
    'logfoldchanges' for frames that predate it. A missing p-value
    (independent-filtered or a Cook's outlier) is not significant.
    """
    if de_results is None or de_results.empty:
        return de_results
    if 'pvals_adj' not in de_results.columns:
        return de_results.iloc[0:0]
    if 'abs_logfoldchange' in de_results.columns:
        magnitude = de_results['abs_logfoldchange']
    else:
        magnitude = de_results['logfoldchanges'].abs()
    keep = (de_results['pvals_adj'] < fdr) & (magnitude > lfc)
    return de_results[keep.fillna(False)].copy()


def tested_gene_set(de_results):
    """Gene names that survived DE filtering, or None when unavailable.

    Returns the 'filter_status == tested' names from a DE-results frame (all
    names if the frame has no 'filter_status' column), or None if 'de_results'
    is None or has no 'names' column. Per-gene figures restrict to this set so a
    gene the DE table filtered out (e.g. one expressed in a single donor) never
    appears -- applied in shared code so the in-app and exported plots agree.
    """
    if de_results is None:
        return None
    cols = getattr(de_results, 'columns', [])
    if 'names' not in cols:
        return None
    dr = de_results
    if 'filter_status' in cols:
        dr = dr[dr['filter_status'] == 'tested']
    return set(dr['names'].astype(str))


def required_donor_count(n_disease, n_control, min_samples=None,
                         large_n=10, min_prop=0.7):
    """How many donors must carry a gene for it to be tested.

    The smaller arm's size, so a gene present in only one arm still
    qualifies, eased for large cohorts by 'large_n + (n - large_n) *
    min_prop' -- which is why a 12 v 12 design asks for 11.4 rather than
    12. 'min_samples' overrides it outright.

    Public because the DE page states this number in its filter summary.
    Computing it twice, once for the model and once for the label, is
    how a page ends up describing a rule it is not running.
    """
    if min_samples is not None and min_samples > 0:
        required = float(min_samples)
    else:
        sizes = [int(n) for n in (n_disease, n_control) if n and n > 0]
        if not sizes:
            return None
        required = float(min(sizes))
    if required > large_n:
        required = large_n + (required - large_n) * min_prop
    return required


def filter_by_expression(counts, roles, *, min_count=DE_FILTER_MIN_COUNT,
                         min_total_count=DE_FILTER_MIN_TOTAL_COUNT,
                         min_samples=None, large_n=10, min_prop=0.7,
                         studies=None):
    """Sample-level expression-prevalence gene filter for pseudobulk counts.

    This is edgeR's 'filterByExpr', and it is the gene filter for
    pseudobulk DE: once cells are aggregated into one profile per donor
    you are doing bulk RNA-seq, and this is what bulk RNA-seq uses.

    'counts' is a samples x genes matrix of raw pseudobulk counts; 'roles' is a
    per-sample array of 'disease' / 'control' / other. A gene is kept when its
    CPM clears a library-size-adjusted cutoff (derived from 'min_count' against
    the median library size) in at least the smallest group's number of samples,
    and its summed count clears 'min_total_count'. This is donor-aware: the
    sample count is what enforces cross-donor reproducibility.

    'min_samples' overrides the auto smallest-group size (None = auto). The
    'large_n' / 'min_prop' rule caps the required sample count so very large
    cohorts don't demand detection in every sample. Returns a boolean mask over
    genes.

    'studies' names each sample's cohort. When given, the whole rule is
    applied within each study and a gene must pass in every one -- which
    is what running each dataset separately would do, and the only
    honest rule for an atlas. Without it, a pooled cohort hides the
    problem: on 12 disease + 12 control donors the requirement is 11.4
    samples, so a gene present in exactly one study's 12 donors and
    absent from the other clears it, and gets tested as though both
    cohorts could see it. Library size, group sizes and the resulting
    cutoff are all computed within each study.
    """
    counts = np.asarray(counts, dtype=float)
    if studies is not None:
        studies_arr = np.asarray(studies).astype(str)
        if counts.ndim != 2 or counts.shape[0] == 0:
            return np.zeros(
                counts.shape[1] if counts.ndim == 2 else 0, dtype=bool)
        levels = pd.unique(studies_arr)
        if len(levels) > 1:
            keep = np.ones(counts.shape[1], dtype=bool)
            roles_arr = np.asarray(roles).astype(str)
            for level in levels:
                rows = studies_arr == level
                if not rows.any():
                    continue
                keep &= filter_by_expression(
                    counts[rows], roles_arr[rows], min_count=min_count,
                    min_total_count=min_total_count, min_samples=min_samples,
                    large_n=large_n, min_prop=min_prop)
            return keep

    if counts.ndim != 2 or counts.shape[0] == 0:
        return np.zeros(counts.shape[1] if counts.ndim == 2 else 0, dtype=bool)

    lib = np.clip(counts.sum(axis=1), 1.0, None)
    median_lib = float(np.median(lib))
    cpm = counts / lib[:, None] * 1e6
    cpm_cutoff = min_count / median_lib * 1e6

    roles_arr = np.asarray(roles).astype(str)
    required = required_donor_count(
        int((roles_arr == 'disease').sum()),
        int((roles_arr == 'control').sum()),
        min_samples=min_samples, large_n=large_n, min_prop=min_prop)
    if required is None:
        required = float(counts.shape[0])

    tol = 1e-14
    keep_cpm = (cpm >= cpm_cutoff).sum(axis=0) >= (required - tol)
    keep_total = counts.sum(axis=0) >= (min_total_count - tol)
    return keep_cpm & keep_total


def _detection_counts_by_donor(adata, counts_layer, donor_values, cell_mask):
    """Per-donor non-zero cell counts, and cells per donor.

    Returns '(donors, nonzero_counts, n_cells, var_names)' where
    'nonzero_counts' is 'n_donors x n_genes'. Everything the detection
    rules need is an aggregation of this, so the matrix is walked once
    per donor rather than once per (study, arm, donor) combination.
    """
    from scipy import sparse as sp

    if counts_layer is not None and counts_layer in adata.layers:
        X, var_names = adata.layers[counts_layer], list(adata.var_names)
    elif adata.raw is not None:
        X, var_names = adata.raw.X, list(adata.raw.var_names)
    else:
        X, var_names = adata.X, list(adata.var_names)

    binary = (X != 0)
    donors = [d for d in pd.unique(donor_values) if d is not None]
    counts = np.zeros((len(donors), len(var_names)), dtype=np.float64)
    n_cells = np.zeros(len(donors), dtype=np.float64)
    for k, donor in enumerate(donors):
        m = (donor_values == donor)
        if cell_mask is not None:
            m = m & cell_mask
        n_cells[k] = m.sum()
        if not n_cells[k]:
            continue
        if sp.issparse(binary):
            counts[k] = np.asarray(binary[m].sum(axis=0)).ravel()
        else:
            counts[k] = np.asarray(binary[m]).sum(axis=0)
    return donors, counts, n_cells, var_names


def genes_passing_detection(adata, genes, min_pct, counts_layer=None, *,
                            donor_col=None, min_donor_frac=0.0,
                            study_col=None, cell_mask=None):
    """Genes detected (count > 0) widely enough to be worth testing.

    The base rule keeps a gene when it is non-zero in at least 'min_pct'
    of the cells of one arm (disease or control from 'obs['_role']'),
    every donor of that arm pooled together. That pooled figure is a
    prevalence over *cells*, which is the wrong unit in two ways, both
    fixable here:

    'donor_col' + 'min_donor_frac'
        Pooling lets a gene pass on the strength of a few donors: one
        donor expressing it in 30% of their cells outvotes eleven donors
        at 1%. Supplying a donor column computes the rate *within* each
        donor and requires at least 'min_donor_frac' of an arm's donors
        to clear 'min_pct' individually. 0 disables it (pooled rule).

    'study_col'
        On a combined object the pooled rate blends cohorts, so a gene
        visible in one study and absent from another passes on the
        average of the two. Supplying a study column applies the whole
        rule within each study and requires every study to pass, which
        is the per-dataset filter people assume is running. Studies
        contributing no cells to either arm are skipped rather than
        failing the gene.

    'cell_mask'
        Restrict the calculation to these cells. The pipeline passes the
        donors that survive 'min_cells', so cells that never reach the
        model do not shape the gene list either.

    'min_pct' is a fraction; 0 disables the filter. Detection is read
    from the same counts the DE uses (layer / raw / X). Returns the
    input list unchanged when disabled or roles are missing; genes
    absent from the detection matrix are kept (they cannot be assessed).
    """
    if not min_pct or min_pct <= 0 or '_role' not in adata.obs.columns:
        return list(genes)

    roles = adata.obs['_role'].astype(str).str.lower().to_numpy()
    is_arm = {'disease': roles == 'disease', 'control': roles == 'control'}
    if cell_mask is not None:
        cell_mask = np.asarray(cell_mask, dtype=bool)
        is_arm = {k: v & cell_mask for k, v in is_arm.items()}
    if not is_arm['disease'].any() or not is_arm['control'].any():
        return list(genes)

    stratify_donors = bool(donor_col) and min_donor_frac and min_donor_frac > 0
    if stratify_donors and donor_col not in adata.obs.columns:
        stratify_donors = False
    if study_col and study_col not in adata.obs.columns:
        study_col = None

    studies = (adata.obs[study_col].astype(str).to_numpy()
               if study_col else np.full(adata.n_obs, '_all', dtype=object))

    # Grouping key per cell. A donor sits inside exactly one arm and one
    # study, so grouping by donor is safe for every rule: the pooled rate
    # is just the cell-weighted sum over that slice's donors. Without a
    # donor column, group by arm x study instead -- anything coarser
    # would let a bucket span two slices and average them together.
    if donor_col and donor_col in adata.obs.columns:
        donor_values = adata.obs[donor_col].astype(str).to_numpy()
    else:
        donor_values = np.array(
            [f"{r}|{t}" for r, t in zip(roles, studies)], dtype=object)
        stratify_donors = False

    donors, counts, n_cells, var_names = _detection_counts_by_donor(
        adata, counts_layer, donor_values, cell_mask)
    donor_index = {d: k for k, d in enumerate(donors)}

    def arm_passes(cells):
        """Boolean per gene: does this arm, within these cells, pass?"""
        if not cells.any():
            return None
        rows, weights = [], []
        for donor in pd.unique(donor_values[cells]):
            k = donor_index.get(donor)
            if k is None or not n_cells[k]:
                continue
            rows.append(k)
            weights.append(n_cells[k])
        if not rows:
            return None
        rows = np.asarray(rows)
        weights = np.asarray(weights, dtype=np.float64)
        if stratify_donors:
            rates = counts[rows] / weights[:, None]
            n_clearing = (rates >= min_pct).sum(axis=0)
            required = max(1, int(np.ceil(min_donor_frac * len(rows))))
            return n_clearing >= required
        pooled = counts[rows].sum(axis=0) / weights.sum()
        return pooled >= min_pct

    passing = None
    for study in pd.unique(studies):
        in_study = studies == study
        per_arm = [arm_passes(is_arm[a] & in_study) for a in ('disease', 'control')]
        per_arm = [p for p in per_arm if p is not None]
        if not per_arm:
            # No cells from this study in either arm -- it has no opinion.
            continue
        # A gene needs one arm to express it; both arms failing is a fail.
        study_pass = np.logical_or.reduce(per_arm)
        passing = study_pass if passing is None else (passing & study_pass)

    if passing is None:
        return list(genes)

    max_det = dict(zip(var_names, passing))
    return [g for g in genes if max_det.get(g, True)]


def donors_meeting_min_cells(adata, sample_col, min_cells):
    """Boolean cell mask for donors contributing at least *min_cells* cells.

    'create_pseudobulk' drops thin donors, so their cells never reach the
    model. Detection is computed over this mask so the gene list is
    decided by the same donors the contrast is.
    """
    if not min_cells or sample_col not in adata.obs.columns:
        return None
    sizes = adata.obs.groupby(sample_col, observed=True).size()
    keep = set(sizes[sizes >= min_cells].index)
    if len(keep) == len(sizes):
        return None
    return adata.obs[sample_col].isin(keep).to_numpy()


def run_de_pipeline(adata, sample_col, condition_col,
                    pathway_gene_sets, min_cells=DE_MIN_CELLS, min_expressing_samples=DE_MIN_EXPRESSING_SAMPLES,
                    de_method='ttest', max_se=DE_MAX_SE, full_genome=False,
                    moderate=False, unit='sample', counts_layer=None,
                    detection_min_pct=DE_DETECTION_DEFAULT,
                    detection_min_donor_frac=DE_DETECTION_MIN_DONOR_FRAC,
                    detection_study_col=None,
                    fdr_genes=None,
                    deseq2_independent_filter=DE_DESEQ2_INDEPENDENT_FILTER,
                    deseq2_cooks_filter=DE_DESEQ2_COOKS_FILTER,
                    filter_min_count=DE_FILTER_MIN_COUNT,
                    filter_min_samples=DE_FILTER_MIN_SAMPLES,
                    covariates=None,
                    progress_callback=None):
    """End-to-end DE pipeline.

    Disease/control assignment comes from 'adata.obs['_role']'
    (set in the Inspect tab); the role column propagates through
    'create_pseudobulk' so engines read it directly.

    'covariates' names extra sample-level 'obs' columns to adjust for --
    'study' when running one model over pooled cohorts. They enter the
    DESeq2 design before 'condition', so the contrast stays the condition
    coefficient. Ignored by the t-test engine, which has no design.
    Covariates that cannot be estimated (constant, unique per sample, or
    collinear with condition) are dropped, and the design actually used
    is recorded on the results' 'attrs'.

    Parameters
    ----------
    adata : anndata.AnnData
        Input data. Must have 'obs['_role']' populated.
    sample_col : str
        Sample grouping column.
    condition_col : str
        Condition column (display/stratification only -- actual split is by role).
    pathway_gene_sets : dict
        '{pathway: [genes]}'.
    min_cells : int
        Min cells per sample.
    min_expressing_samples : int
        Min expressing samples per gene.
    de_method : 'ttest' | 'ttest_raw' | 'deseq2'
        Engine for pseudobulk mode (unit='sample'). Ignored when unit='cell'.
    max_se : float
        Standard-error cap.
    full_genome : bool
        Test all genes. Pathway membership is still annotated on the results.
    detection_min_donor_frac : float
        Fraction of an arm's donors that must individually clear
        'detection_min_pct'. Without it the detection rate is pooled over
        an arm's cells, so a single donor expressing a gene widely can
        carry it past the filter. 0 restores the pooled rule.
    detection_study_col : str, optional
        Apply the detection rule within each level of this 'obs' column
        and require every level to pass. Set it to the study column on a
        combined object: pooled detection otherwise blends cohorts, and
        a gene visible in one study and absent from another passes on
        the average of the two.
    fdr_genes : set of str, optional
        Restrict the BH correction to this committed gene list (hypothesis
        mode). The fit + normalisation stay genome-wide; only the FDR
        denominator changes, and results are subset to these genes. When None,
        FDR is computed over every tested gene (discovery mode).
    moderate : bool
        Apply EB variance moderation (t-test only).
    unit : 'sample' | 'cell'
        Unit of replication. 'sample' (default) aggregates cells into one
        pseudobulk profile per 'sample_col' and tests across samples --
        the statistically correct unit. 'cell' skips aggregation and runs
        scanpy 'rank_genes_groups' (Wilcoxon) treating each cell as an
        observation; this is the scanpy marker-gene workflow applied to the
        disease/control contrast. Cell mode ignores donor grouping, so it is
        pseudoreplicated and inflates significance -- exposed as an explicit
        opt-in for studies with too few samples for pseudobulk.
    progress_callback : callable(msg), optional

    Returns
    -------
    de_results, significant_genes, pathway_coverage, sample_df
    """
    if '_role' not in adata.obs.columns:
        # Reconstruct _role from a saved role_map if one exists (h5ad written
        # mid-flow may have role_map in uns but no materialised _role column).
        role_map = adata.uns.get('role_map') if hasattr(adata, 'uns') else None
        if isinstance(role_map, dict) and role_map and condition_col in adata.obs.columns:
            from kosmic.scrna.inspect.roles import resolve_roles
            cond_values = adata.obs[condition_col].values
            is_disease, two_group = resolve_roles(cond_values, role_map)
            roles = np.where(
                is_disease, 'disease',
                np.where(two_group, 'control', 'exclude'))
            adata.obs['_role'] = pd.Categorical(
                roles, categories=['control', 'disease', 'exclude'])
        else:
            raise ValueError(
                "adata.obs is missing the '_role' column. Set Disease / "
                "Control / Exclude in the Inspect tab (Sample Setup), hit "
                "Save, then re-run DE.")
    def _progress(msg):
        if progress_callback:
            progress_callback(msg)

    species = detect_species(list(adata.var_names))

    # Always compute coverage so pathway annotation is available downstream.
    coverage_var_names = set(adata.raw.var_names) if adata.raw is not None else set(adata.var_names)
    metabolic_genes, pathway_coverage = prepare_gene_coverage(
        pathway_gene_sets, coverage_var_names, species
    )

    if unit == 'cell':
        return _run_cell_level_de_pipeline(
            adata, metabolic_genes, pathway_coverage,
            full_genome=full_genome, progress_callback=progress_callback)

    # A specified counts layer (e.g. DecontX) lives on the current
    # var_names, so the genome-wide gene universe comes from there, not raw.
    def _genome_var_names():
        if counts_layer is not None:
            return list(adata.var_names)
        return list(adata.raw.var_names) if adata.raw is not None else list(adata.var_names)

    if de_method == 'deseq2' and not full_genome:
        # DESeq2 must see the full transcriptome -- median-of-ratios
        # normalisation assumes most genes don't change, which is violated
        # by curated subsets where most pathway members shift one way.
        # Run genome-wide, then filter results to pathway genes.
        genes_to_test = _genome_var_names()
        _progress(f"DESeq2: running on full transcriptome ({len(genes_to_test):,} genes) "
                  f"for proper normalization, then filtering to {len(metabolic_genes)} pathway genes")
    elif full_genome:
        genes_to_test = _genome_var_names()
        _progress(f"Full-genome DE: testing {len(genes_to_test):,} genes")
    else:
        genes_to_test = metabolic_genes
        _progress(f"Pathway DE: testing {len(genes_to_test)} pathway genes")

    if not genes_to_test:
        return pd.DataFrame(), pd.DataFrame(), pathway_coverage, pd.DataFrame()

    # Genes considered before any filtering (the dataset / set size), for the
    # results-count readout.
    n_genes_input = len(genes_to_test)

    # Biological pre-filter: drop the near-zero detection floor before DESeq2
    # estimates dispersion / size factors, so noise genes never enter the model.
    # One resolution of the study column for the whole run: None unless
    # the object really holds more than one cohort.
    study_col = detection_study_col
    if study_col and study_col in adata.obs.columns:
        if adata.obs[study_col].astype(str).nunique() < 2:
            study_col = None
    else:
        study_col = None
    detection_study_col = study_col
    if study_col:
        _progress(
            f"Gene filter stratified by '{study_col}': a gene must pass in "
            f"every study, not on the pooled cohort")

    if detection_min_pct and detection_min_pct > 0:
        n_before = len(genes_to_test)
        # Donors below min_cells are dropped by create_pseudobulk, so their
        # cells must not shape the gene list either -- otherwise the filter
        # is decided partly by cells the model never sees.
        cell_mask = donors_meeting_min_cells(adata, sample_col, min_cells)
        genes_to_test = genes_passing_detection(
            adata, genes_to_test, detection_min_pct, counts_layer=counts_layer,
            donor_col=sample_col, min_donor_frac=detection_min_donor_frac,
            study_col=study_col, cell_mask=cell_mask)
        scope = f"detected in >={detection_min_pct:.0%} of cells"
        if detection_min_donor_frac and detection_min_donor_frac > 0:
            scope += f" in >={detection_min_donor_frac:.0%} of a group's donors"
        else:
            scope += " in >=1 condition"
        if study_col:
            scope += f", in every '{study_col}'"
        if cell_mask is not None:
            scope += f" (donors with >={min_cells} cells only)"
        _progress(
            f"Detection pre-filter: kept {len(genes_to_test):,} of {n_before:,} "
            f"genes ({scope}); dropped {n_before - len(genes_to_test):,}")
        if not genes_to_test:
            return pd.DataFrame(), pd.DataFrame(), pathway_coverage, pd.DataFrame()

    # Pseudobulk: sum counts for CPM/DESeq2, mean for raw t-test.
    use_mean = (de_method == 'ttest_raw')
    aggregate = 'mean' if use_mean else 'sum'
    _progress("Creating pseudobulk profiles...")
    # The study column has to reach the gene filter, which lives past
    # pseudobulk aggregation -- so it travels as a sample-level covariate
    # whether or not it is in the design.
    pb_covariates = list(covariates or ())
    if detection_study_col and detection_study_col not in pb_covariates:
        pb_covariates.append(detection_study_col)
    pseudobulk_matrix, sample_df, genes_used = create_pseudobulk(
        adata, genes_to_test, sample_col, condition_col,
        min_cells=min_cells, aggregate=aggregate, counts_layer=counts_layer,
        covariates=pb_covariates,
    )

    if pseudobulk_matrix.size == 0:
        return pd.DataFrame(), pd.DataFrame(), pathway_coverage, sample_df

    # Log the DE engine + version so provenance is captured in the run log.
    deseq2_filter_counts = None
    if de_method == 'deseq2':
        import pydeseq2
        _progress(
            f"DE engine: pydeseq2 v{pydeseq2.__version__} (DESeq2) on "
            f"{len(genes_used):,} genes...")
        de_results = _run_deseq2(
            pseudobulk_matrix, sample_df, genes_used,
            min_expressing_samples=min_expressing_samples,
            independent_filter=deseq2_independent_filter,
            cooks_filter=deseq2_cooks_filter,
            filter_min_count=filter_min_count,
            filter_min_samples=filter_min_samples,
            covariates=covariates,
            study_col=detection_study_col,
        )
        if not de_results.empty and de_results.attrs.get('design'):
            _progress(f"DESeq2 design: {de_results.attrs['design']}")
            dropped = [c for c in (covariates or ())
                       if c not in de_results.attrs.get('covariates_used', [])]
            if dropped:
                _progress("  dropped as not estimable: " + ", ".join(dropped))
        fcounts = de_results.attrs.get('deseq2_filter_counts', {}) if not de_results.empty else {}
        deseq2_filter_counts = fcounts
        n_indep = fcounts.get('independent_filter', 0)
        n_cooks = fcounts.get('cooks_outlier', 0)
        if n_indep or n_cooks:
            _progress(
                f"DESeq2 filtering: {n_indep:,} independent-filtered (low mean "
                f"count), {n_cooks:,} Cook's outliers -- FDR left blank for "
                f"these, not collapsed to 1.0")
        # If we ran genome-wide for normalisation, filter results to pathway genes.
        if not full_genome and not de_results.empty:
            metabolic_set = set(metabolic_genes)
            de_results = de_results[de_results['names'].isin(metabolic_set)].copy()
            de_results.reset_index(drop=True, inplace=True)
            _progress(f"Filtered to {len(de_results)} pathway genes")
    else:
        normalize = 'log2' if de_method == 'ttest_raw' else 'cpm'
        _progress(
            f"DE engine: Welch t-test (normalize={normalize}, "
            f"moderate={moderate}) on {len(genes_used):,} genes...")
        de_results = run_pseudobulk_de(
            pseudobulk_matrix, sample_df, genes_used,
            min_expressing_samples=min_expressing_samples,
            method='ttest', max_se=max_se,
            moderate=moderate, normalize=normalize,
            filter_min_count=filter_min_count,
            filter_min_samples=filter_min_samples,
            study_col=detection_study_col,
        )

    # Hypothesis mode: correct over the committed gene list only. The fit and
    # normalisation above stayed genome-wide; here we just swap the BH
    # denominator so a small a-priori gene set is not penalised by the
    # transcriptome-wide multiple-testing burden. The full genome-wide ranking
    # is preserved (see below) for consumers that need it -- notably preranked
    # GSEA, whose null is the whole transcriptome.
    genome_wide_ranking = None
    if fdr_genes is not None and not de_results.empty:
        from kosmic.numerical import bh_fdr
        # The committed list may be in a different species namespace than the
        # data (e.g. human pathway symbols against a mouse dataset). Convert
        # to the dataset's namespace before matching, mirroring
        # prepare_gene_coverage / _pathway_retention.
        fdr_genes = {format_gene_for_species(g, species) for g in fdr_genes}
        genome_wide_ranking = de_results.copy()
        de_results = de_results[de_results['names'].isin(fdr_genes)].copy()
        de_results.reset_index(drop=True, inplace=True)
        if 'pvals' in de_results.columns and not de_results.empty:
            de_results['pvals_adj'] = bh_fdr(de_results['pvals'].astype(float).values)
            # Committed-set BH re-tests genes DESeq2 independent-filtered
            # genome-wide, so retag from the fresh result: a finite p-value is
            # now 'tested'; a Cook's outlier (NaN p-value) stays flagged.
            if 'filter_status' in de_results.columns:
                finite_p = np.isfinite(de_results['pvals'].astype(float).values)
                de_results['filter_status'] = np.where(
                    finite_p, 'tested', 'cooks_outlier')
        _progress(
            f"Hypothesis-mode FDR: BH over {len(de_results):,} committed genes "
            f"(genome-wide normalization retained)")

    # Annotate (adds pathway column + significance flags).
    _progress("Annotating results...")
    de_results, significant_genes = annotate_de_results(de_results, pathway_coverage)

    # "% cells expressing" separates real biological hits from sparse-noise
    # artefacts in snRNA-seq where dropout is brutal.
    _progress("Annotating cell-level expression rates...")
    de_results = annotate_pct_expressing(de_results, adata)
    significant_genes = annotate_pct_expressing(significant_genes, adata)

    # Stash the pre-restriction genome-wide ranking so genome-wide consumers
    # (e.g. fgsea) are not limited to the hypothesis gene set.
    if genome_wide_ranking is not None:
        genome_wide_ranking, _ = annotate_de_results(genome_wide_ranking, pathway_coverage)
        de_results.attrs['genome_wide_ranking'] = genome_wide_ranking

    # Re-attach DESeq2 filter counts (attrs are lost through the subset/annotate
    # steps above) so the GUI can surface how many genes DESeq2 excluded.
    if deseq2_filter_counts is not None:
        de_results.attrs['deseq2_filter_counts'] = deseq2_filter_counts
    de_results.attrs['n_genes_input'] = n_genes_input

    return de_results, significant_genes, pathway_coverage, sample_df


def run_cell_level_de(adata):
    """Cell-level DE via scanpy 'rank_genes_groups' (Wilcoxon).

    Treats each cell as an observation -- no per-sample pseudobulk
    aggregation. This is the same scanpy machinery KOSMIC uses for cluster
    marker dotplots ('compute_cluster_marker_genes'), applied to the
    disease vs control contrast read from 'adata.obs['_role']'.

    Input normalisation
    -------------------
    scanpy 'rank_genes_groups' expects log-normalised expression. KOSMIC's
    'adata.raw' holds *raw counts* (normalize.py snapshots raw counts before
    normalising), so we rebuild a log1p(CPM) matrix from counts here rather
    than passing 'use_raw=True'. Running Wilcoxon on un-normalised counts
    would leave the test confounded by per-cell sequencing depth -- biasing
    every gene toward whichever arm is sequenced deeper (typically showing up
    as spurious global up-regulation) -- and would misapply scanpy's logFC
    formula (mean-of-logs + expm1), which is only meaningful on log data.

    Statistically this is still pseudoreplication: cells within a donor are
    correlated, so the effective n is inflated and p-values are
    anticonservative relative to pseudobulk. Use only as an explicit
    opt-in for datasets with too few samples for pseudobulk.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have 'obs['_role']' with 'disease' / 'control' labels.

    Returns
    -------
    pandas.DataFrame
        Genome-wide results with columns: names, scores, logfoldchanges,
        pvals, pvals_adj, disease_samples, control_samples (the last two are
        cell counts per arm). Empty frame if either arm has no cells.
    """
    import scanpy as sc

    if '_role' not in adata.obs.columns:
        raise ValueError(
            "adata.obs is missing the '_role' column. Set Disease / "
            "Control / Exclude in the Inspect tab (Sample Setup), hit "
            "Save, then re-run DE.")

    roles = adata.obs['_role'].astype(str)
    keep = roles.isin(['disease', 'control']).values
    n_d = int((roles == 'disease').sum())
    n_c = int((roles == 'control').sum())
    if n_d < 1 or n_c < 1:
        return pd.DataFrame()

    # Source raw counts on the full gene set, then CPM + log1p so the
    # Wilcoxon test sees depth-normalised log expression (see docstring).
    if adata.raw is not None:
        work = adata.raw.to_adata()
        from_counts = True
    elif 'counts' in adata.layers:
        work = adata.copy()
        work.X = work.layers['counts'].copy()
        from_counts = True
    else:
        # No stored counts -- assume adata.X is already log-normalised.
        work = adata.copy()
        from_counts = False

    work = work[keep].copy()
    work.raw = None
    if from_counts:
        sc.pp.normalize_total(work, target_sum=1e4)
        sc.pp.log1p(work)

    work.obs['_role'] = pd.Categorical(
        work.obs['_role'].astype(str), categories=['control', 'disease'])

    sc.tl.rank_genes_groups(
        work, groupby='_role', groups=['disease'], reference='control',
        method='wilcoxon', use_raw=False, pts=False,
    )

    res = work.uns['rank_genes_groups']
    grp = 'disease'
    de_results = pd.DataFrame({
        'names': list(res['names'][grp]),
        'scores': np.asarray(res['scores'][grp], dtype=float),
        'logfoldchanges': np.asarray(res['logfoldchanges'][grp], dtype=float),
        'pvals': np.asarray(res['pvals'][grp], dtype=float),
        'pvals_adj': np.asarray(res['pvals_adj'][grp], dtype=float),
    })
    de_results['disease_samples'] = n_d
    de_results['control_samples'] = n_c
    return de_results


def _run_cell_level_de_pipeline(adata, metabolic_genes, pathway_coverage,
                                full_genome=False, progress_callback=None):
    """Cell-level branch of run_de_pipeline. Returns the same 4-tuple."""
    def _progress(msg):
        if progress_callback:
            progress_callback(msg)

    _progress("Cell-level DE: scanpy rank_genes_groups (Wilcoxon), "
              "disease vs control cells (pseudoreplicated -- ignores donor)...")
    de_results = run_cell_level_de(adata)

    roles = adata.obs['_role'].astype(str)
    sample_df = pd.DataFrame({
        'sample': ['disease_cells', 'control_cells'],
        'condition': ['disease', 'control'],
        'role': ['disease', 'control'],
        'n_cells': [int((roles == 'disease').sum()),
                    int((roles == 'control').sum())],
    })

    if de_results.empty:
        return de_results, pd.DataFrame(), pathway_coverage, sample_df

    # rank_genes_groups corrects FDR genome-wide; filter to pathway genes
    # afterwards (mirrors the DESeq2 path) unless the user asked for all genes.
    if not full_genome:
        de_results = de_results[
            de_results['names'].isin(set(metabolic_genes))
        ].reset_index(drop=True)
        _progress(f"Filtered to {len(de_results)} pathway genes")

    _progress("Annotating results...")
    de_results, significant_genes = annotate_de_results(de_results, pathway_coverage)

    _progress("Annotating cell-level expression rates...")
    de_results = annotate_pct_expressing(de_results, adata)
    significant_genes = annotate_pct_expressing(significant_genes, adata)

    return de_results, significant_genes, pathway_coverage, sample_df


def annotate_pct_expressing(de_results, adata):
    """Add 'pct_disease' / 'pct_control' columns to a DE results frame.

    Each value is the fraction of cells (cell-level, not pseudobulk) in
    that role group with a non-zero count for the gene. Vectorised over
    the entire matrix in a single pass; per-gene lookup follows.

    Reads role assignment from 'adata.obs['_role']'. When the column
    is missing, the input frame is returned unchanged (no-op).
    """
    if de_results is None or de_results.empty:
        return de_results
    if '_role' not in adata.obs.columns:
        return de_results
    if 'names' not in de_results.columns:
        return de_results

    from scipy import sparse as sp

    is_disease = (adata.obs['_role'].astype(str) == 'disease').values
    is_control = (adata.obs['_role'].astype(str) == 'control').values
    n_d = int(is_disease.sum())
    n_c = int(is_control.sum())

    X = adata.X
    if sp.issparse(X):
        # (X != 0) on a sparse matrix returns sparse; sum over axis 0 gives per-gene counts.
        binary = (X != 0).astype(np.float32)
        d_counts = (np.asarray(binary[is_disease].sum(axis=0)).flatten()
                    if n_d > 0 else np.zeros(X.shape[1]))
        c_counts = (np.asarray(binary[is_control].sum(axis=0)).flatten()
                    if n_c > 0 else np.zeros(X.shape[1]))
    else:
        X_arr = np.asarray(X)
        d_counts = (X_arr[is_disease] > 0).sum(axis=0) if n_d > 0 else np.zeros(X_arr.shape[1])
        c_counts = (X_arr[is_control] > 0).sum(axis=0) if n_c > 0 else np.zeros(X_arr.shape[1])

    pct_d_per_var = d_counts / max(n_d, 1)
    pct_c_per_var = c_counts / max(n_c, 1)

    var_to_idx = {str(g): i for i, g in enumerate(adata.var_names)}
    names_arr = de_results['names'].astype(str).values
    pct_disease = np.array([
        pct_d_per_var[var_to_idx[g]] if g in var_to_idx else np.nan
        for g in names_arr
    ])
    pct_control = np.array([
        pct_c_per_var[var_to_idx[g]] if g in var_to_idx else np.nan
        for g in names_arr
    ])

    out = de_results.copy()
    out['pct_disease'] = pct_disease
    out['pct_control'] = pct_control
    return out


def usable_covariates(sample_df, covariates, conditions):
    """Covariates that can actually go in the design, in order.

    A covariate is dropped when it is absent, constant (contributes no
    information and makes the model matrix singular), unique per sample
    (perfectly fits the data, leaving nothing for condition), or
    collinear with condition -- which is the one that matters: if every
    sample of a cohort shares one condition, 'study' and 'condition'
    explain the same split and the contrast is not estimable. Better to
    drop it and say so than to hand DESeq2 a rank-deficient design.
    """
    import numpy as np

    usable = []
    n = len(sample_df)
    for col in (covariates or ()):
        if col not in sample_df.columns:
            continue
        values = sample_df[col].astype(str)
        n_levels = values.nunique()
        if n_levels < 2 or n_levels >= n:
            continue
        # Collinear with condition? Every level maps to one condition.
        by_level = {lv: set(np.asarray(conditions)[(values == lv).values])
                    for lv in values.unique()}
        if all(len(conds) < 2 for conds in by_level.values()):
            continue
        usable.append(col)
    return usable


#: Historical private alias -- the rule is now public because the DE page
#: needs it too: what the GUI offers as a covariate and what the model
#: will accept have to be the same list, or a user ticks something that
#: is silently dropped at run time.
_usable_covariates = usable_covariates


def _study_values(sample_df, study_col):
    """Per-sample study labels for the gene filter, or None.

    None whenever there is nothing to stratify -- no column, or one
    cohort -- so a single study takes the plain pooled path.
    """
    if not study_col or study_col not in sample_df.columns:
        return None
    values = sample_df[study_col].astype(str).values
    return values if len(pd.unique(values)) > 1 else None


def _run_deseq2(pseudobulk_matrix, sample_df, gene_names,
                min_expressing_samples=DE_MIN_EXPRESSING_SAMPLES,
                independent_filter=DE_DESEQ2_INDEPENDENT_FILTER,
                cooks_filter=DE_DESEQ2_COOKS_FILTER,
                filter_min_count=DE_FILTER_MIN_COUNT,
                filter_min_samples=DE_FILTER_MIN_SAMPLES,
                covariates=None, study_col=None):
    """Run DESeq2 on pseudobulk sum counts.

    'independent_filter' / 'cooks_filter' pass through to pydeseq2's
    'DeseqStats'. When a gene is dropped by either, its 'padj' (and, for a
    Cook's outlier, 'pvalue') comes back NaN -- kept as NaN here, not collapsed
    to 1.0, and tagged in a 'filter_status' column so the exclusion is visible.

    Parameters
    ----------
    pseudobulk_matrix : numpy.ndarray
        '(n_samples, n_genes)' sum counts.
    sample_df : pandas.DataFrame
        Sample metadata; must include a 'role' column.
    gene_names : list of str
    min_expressing_samples : int

    Returns
    -------
    pandas.DataFrame
        DE results in the standard format.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    if 'role' not in sample_df.columns:
        raise ValueError(
            "sample_df is missing the 'role' column. Set roles in the "
            "Inspect tab and regenerate the pseudobulk.")
    # Drop 'exclude' / unassigned samples.
    keep_mask = sample_df['role'].isin(['disease', 'control']).values
    sample_df = sample_df[keep_mask].reset_index(drop=True)
    pseudobulk_matrix = pseudobulk_matrix[keep_mask]

    if len(sample_df) < 2:
        return pd.DataFrame()

    # PyDESeq2 needs unique string sample IDs.
    sample_ids = sample_df['sample'].astype(str).values
    if len(set(sample_ids)) != len(sample_ids):
        sample_ids = np.array([f"{s}_{i}" for i, s in enumerate(sample_ids)])

    counts_df = pd.DataFrame(
        pseudobulk_matrix.round().astype(int),
        index=sample_ids,
        columns=gene_names,
    )

    # Sample-level gene filter (expression-prevalence, donor-aware).
    keep_mask = filter_by_expression(
        counts_df.values, sample_df['role'].astype(str).values,
        min_count=filter_min_count, min_samples=filter_min_samples,
        studies=_study_values(sample_df, study_col))
    genes_passing = [g for g, k in zip(counts_df.columns, keep_mask) if k]
    if not genes_passing:
        return pd.DataFrame()
    counts_df = counts_df[genes_passing]

    # Also drop zero-total genes (PyDESeq2 errors on them).
    nonzero_genes = counts_df.columns[counts_df.sum(axis=0) > 0]
    counts_df = counts_df[nonzero_genes]
    if counts_df.empty:
        return pd.DataFrame()

    # Emit 'disease' / 'control' strings for the PyDESeq2 design matrix.
    synth_conds = sample_df['role'].astype(str).values
    metadata_df = pd.DataFrame(
        {'condition': pd.Categorical(
            synth_conds, categories=['control', 'disease'],
        )},
        index=sample_ids,
    )

    # Covariates precede 'condition' in the formula, so the contrast is
    # still the condition coefficient -- adjusted for them rather than
    # confounded with them. Pooling cohorts into one model without
    # 'study' here leaves the cohort offset in the residual, inflating
    # dispersion and shrinking every effect toward zero.
    used_covariates = usable_covariates(sample_df, covariates, synth_conds)
    for col in used_covariates:
        metadata_df[col] = pd.Categorical(sample_df[col].astype(str).values)
    design = "~" + " + ".join([*used_covariates, "condition"])

    dds = DeseqDataSet(
        counts=counts_df, metadata=metadata_df,
        design=design,
        quiet=True,
    )
    dds.deseq2()

    contrast = ["condition", "disease", "control"]
    ds = DeseqStats(dds, contrast=contrast, alpha=DEFAULT_FDR,
                    independent_filter=independent_filter,
                    cooks_filter=cooks_filter)
    ds.summary()

    # Unshrunk estimates kept for meta-analysis (shrinkage distorts pooled SE).
    unshrunk_lfc = ds.results_df['log2FoldChange'].copy()
    unshrunk_se = ds.results_df['lfcSE'].copy()

    # Apply LFC shrinkage on the non-intercept coefficient.
    _shrinkage_status = "not attempted"
    try:
        # The *condition* coefficient specifically. Taking the first
        # non-intercept column worked only while condition was the sole
        # term; with a covariate in the design that is 'study', and
        # shrinking it instead would leave the reported fold changes
        # unshrunk while announcing that they had been.
        coeff_name = None
        if hasattr(ds, 'LFC') and hasattr(ds.LFC, 'columns'):
            condition_coeffs = [c for c in ds.LFC.columns
                                if c.lower().startswith('condition')]
            if condition_coeffs:
                coeff_name = condition_coeffs[0]

        if coeff_name is None:
            coeff_name = "condition[T.disease]"

        ds.lfc_shrink(coeff=coeff_name)

        lfc_after = ds.results_df['log2FoldChange']
        max_change = float((unshrunk_lfc - lfc_after).abs().max())
        if max_change < 1e-10:
            _shrinkage_status = f"no effect (coeff='{coeff_name}')"
        else:
            _shrinkage_status = f"applied (coeff='{coeff_name}', max_change={max_change:.4f})"
    except (ValueError, RuntimeError, ImportError, TypeError, KeyError) as exc:
        _shrinkage_status = f"failed (coeff='{coeff_name}'): {exc}"

    results = ds.results_df.copy()

    disease_mask = metadata_df['condition'] == 'disease'
    control_mask = metadata_df['condition'] == 'control'

    de_list = []
    for gene in results.index:
        row = results.loc[gene]
        d_vals = counts_df.loc[disease_mask.values, gene].values.astype(float)
        c_vals = counts_df.loc[control_mask.values, gene].values.astype(float)

        d_mean, c_mean = np.mean(d_vals), np.mean(c_vals)
        d_var = np.var(d_vals, ddof=1) if len(d_vals) > 1 else 0.0
        c_var = np.var(c_vals, ddof=1) if len(c_vals) > 1 else 0.0

        nd, nc = len(d_vals), len(c_vals)
        pooled_std = np.sqrt(((nd - 1) * d_var + (nc - 1) * c_var) / max(nd + nc - 2, 1))
        cohens_d = (d_mean - c_mean) / pooled_std if pooled_std > 0 else 0.0

        # Unshrunk: meta-analysis pooling. Shrunk: single-study visualisation.
        lfc_raw = unshrunk_lfc[gene] if np.isfinite(unshrunk_lfc[gene]) else 0.0
        se_raw = unshrunk_se[gene] if np.isfinite(unshrunk_se[gene]) else np.nan
        lfc_shrunk = row['log2FoldChange'] if np.isfinite(row['log2FoldChange']) else 0.0
        se_shrunk = row['lfcSE'] if np.isfinite(row['lfcSE']) else np.nan
        # Keep DESeq2's NaN verdicts instead of collapsing to 1.0: an
        # independent-filtered gene (padj NaN, pvalue finite) or a Cook's
        # outlier (pvalue NaN) is not the same as a tested-but-non-significant
        # gene. 'filter_status' records which so the exclusion is visible.
        pval = row['pvalue']
        padj = row['padj']
        if np.isfinite(padj):
            filter_status = 'tested'
        elif np.isfinite(pval):
            filter_status = 'independent_filter'
        else:
            filter_status = 'cooks_outlier'
        stat = row['stat'] if np.isfinite(row['stat']) else 0.0
        se_diff = pooled_std * np.sqrt(1 / max(nd, 1) + 1 / max(nc, 1)) if pooled_std > 0 else 0.0

        de_list.append({
            'names': gene,
            'logfoldchanges': lfc_raw,
            'pvals': pval,
            'pvals_adj': padj,
            'disease_mean': d_mean,
            'control_mean': c_mean,
            'disease_var': d_var,
            'control_var': c_var,
            'disease_samples': nd,
            'control_samples': nc,
            'cohens_d': cohens_d,
            'se': se_raw,
            'ci_lower': lfc_raw - 1.96 * se_raw,
            'ci_upper': lfc_raw + 1.96 * se_raw,
            'se_diff': se_diff,
            'statistic': stat,
            'logfoldchanges_shrunk': lfc_shrunk,
            'se_shrunk': se_shrunk,
            'filter_status': filter_status,
        })

    result_df = pd.DataFrame(de_list)
    result_df.attrs['shrinkage_status'] = _shrinkage_status
    # What was actually fitted, not what was asked for -- a covariate can
    # be dropped as inestimable, and the methods should say which.
    result_df.attrs['design'] = design
    result_df.attrs['covariates_used'] = list(used_covariates)
    result_df.attrs['deseq2_filter_counts'] = (
        result_df['filter_status'].value_counts().to_dict())
    result_df.attrs['deseq2_filter_settings'] = {
        'independent_filter': bool(independent_filter),
        'cooks_filter': bool(cooks_filter),
    }
    return result_df
