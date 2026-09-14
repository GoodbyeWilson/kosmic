# Batch Column Detection
# Auto-detect batch/sample and condition columns in AnnData obs, plus
# helpers for the Inspect tab: priority-sort obs columns by biological
# relevance, designate roles by name heuristics, summarise per-sample
# metadata.
from typing import Optional

import pandas as pd


def detect_batch_column(adata) -> Optional[str]:
    """Auto-detect the most likely batch/sample column in obs."""
    # 'study' leads: on a combined master (concat_studies writes that
    # column) the study effect is the one that has to come out. Falling
    # through to 'sample' there would correct 159 donors while leaving the
    # five cohorts intact -- exactly backwards.
    priority = ['study', 'Study',
                'sample', 'Sample', 'donor', 'Donor', 'patient', 'Patient',
                'batch', 'Batch', 'orig.ident', 'library_id']

    for col in priority:
        if col in adata.obs.columns:
            n = adata.obs[col].nunique()
            if 2 <= n <= 100:
                return col

    for col in adata.obs.columns:
        if col in ('leiden', 'clusters', 'seurat_clusters', 'cell_type',
                    'louvain', 'n_genes', 'n_counts', 'total_counts',
                    'n_genes_by_counts', 'pct_counts_mt', 'doublet_score'):
            continue
        try:
            n = adata.obs[col].nunique()
        except TypeError:
            continue  # unhashable dtype (e.g. list-valued cells)
        if 2 <= n <= 20:
            return col

    return None


def detect_condition_column(adata, batch_key: Optional[str] = None) -> Optional[str]:
    """Auto-detect the most likely condition/disease column in obs.

    Looks for columns with 2-5 unique values that aren't the batch column.
    """
    priority = ['condition', 'Condition', 'disease', 'Disease',
                'group', 'Group', 'treatment', 'Treatment',
                'status', 'Status', 'phenotype', 'Phenotype',
                'diagnosis', 'Diagnosis']

    for col in priority:
        if col == batch_key:
            continue
        if col in adata.obs.columns:
            n = adata.obs[col].nunique()
            if 2 <= n <= 10:
                return col

    # Skip known non-condition columns
    skip = {'leiden', 'clusters', 'seurat_clusters', 'cell_type', 'louvain',
            'n_genes', 'n_counts', 'total_counts', 'n_genes_by_counts',
            'pct_counts_mt', 'doublet_score', batch_key}

    for col in adata.obs.columns:
        if col in skip:
            continue
        try:
            n = adata.obs[col].nunique()
        except TypeError:
            continue  # unhashable dtype (e.g. list-valued cells)
        if 2 <= n <= 5:
            return col

    return None


# Inspect-tab config: pattern dictionaries used by the Inspect tab to
# sort and auto-detect obs columns. The Qt-side consumers live in the
# GUI; the patterns themselves live here so they're centrally defined.

# Standardised obs column names always sort to the top.
STANDARDISED_OBS_NAMES = ('condition', 'sample', 'cell_type', '_role')

# Biological / clinical column-name patterns -- middle tier of the
# metadata-table sort, also drives Designate-as auto-detection.
USER_OBS_PATTERNS = (
    'condition', 'group', 'disease', 'status', 'treatment',
    'diagnosis', 'phenotype', 'etiology', 'cohort',
    'donor', 'patient', 'subject', 'sample', 'biosample',
    'individual', 'specimen', 'orig.ident', 'orig_ident', 'library',
    'cell_type', 'celltype', 'cell.type', 'annotation', 'cluster',
    'sex', 'gender', 'age', 'race', 'ethnicity', 'bmi',
    'batch', 'tissue', 'source', 'subtype', 'region',
)

# Scanpy-style auto-computed column-name patterns -- bottom tier.
COMPUTED_OBS_PATTERNS = (
    'n_genes', 'n_counts', 'total_counts', 'pct_counts',
    'doublet', '_mt', 'mt_', '_mito', 'leiden', 'louvain',
    'seurat', 'percent', 'log1p', 'highly_variable',
    'umap_', 'pca_', '_outlier',
)


_SAMPLE_COL_CANDIDATES = ('sample', 'donor_id', 'donor', 'patient',
                          'patient_id', 'subject', 'subject_id',
                          'orig.ident', 'library_id')


#: obs columns that mean "which study did this cell come from", best
#: first. Deliberately excludes 'batch', which in a single study usually
#: means a technical run rather than a separate cohort.
_STUDY_COL_CANDIDATES = ('study', 'dataset', 'cohort', 'accession', 'gse')


def detect_study_column(adata) -> Optional[str]:
    """The obs column identifying the source study, or None.

    A column only counts when it has more than one value: a combined
    object that happens to carry a constant 'study' label is a single
    study, and treating it as an atlas would put an unfittable term in
    the design.
    """
    if adata is None:
        return None
    for cand in _STUDY_COL_CANDIDATES:
        if cand in adata.obs.columns:
            if adata.obs[cand].astype(str).nunique() > 1:
                return cand
    return None


def cohort_shape(adata, sample_col: Optional[str] = None,
                 study_col: Optional[str] = None) -> dict:
    """Is this one study with several donors, or several studies pooled?

    That single fact decides most of a DE setup -- whether 'study'
    belongs in the design, whether the detection filter has to be
    applied per study, and whether an unexplained batch effect is
    expected or alarming -- so it is worth stating rather than leaving
    the user to infer it from a column list.

    Returns a dict with 'is_atlas', 'study_col', 'n_studies',
    'n_donors', 'n_disease', 'n_control', 'n_excluded', 'n_cells' and
    'studies' (one row per study: 'study', 'n_donors', 'n_disease',
    'n_control', 'n_cells'). Every field degrades to 0 / empty rather
    than raising when the columns are missing.
    """
    result = {
        'is_atlas': False, 'study_col': None, 'n_studies': 0,
        'n_donors': 0, 'n_disease': 0, 'n_control': 0, 'n_excluded': 0,
        'n_cells': 0, 'sample_col': None, 'studies': [],
    }
    if adata is None:
        return result

    obs = adata.obs
    result['n_cells'] = int(adata.n_obs)

    if sample_col is None:
        for cand in _SAMPLE_COL_CANDIDATES:
            if cand in obs.columns:
                sample_col = cand
                break
    if sample_col is None or sample_col not in obs.columns:
        return result
    result['sample_col'] = sample_col

    if study_col is None:
        study_col = detect_study_column(adata)
    if study_col is not None and study_col not in obs.columns:
        study_col = None

    frame = pd.DataFrame({'sample': obs[sample_col].astype(str).values})
    if '_role' in obs.columns:
        frame['role'] = obs['_role'].astype(str).str.lower().values
    else:
        frame['role'] = 'unassigned'
    if study_col:
        frame['study'] = obs[study_col].astype(str).values

    # One row per donor. A donor belongs to one arm and one study, so
    # taking the first cell's value is exact, not an approximation.
    per_donor = frame.groupby('sample', observed=True).first()
    counts = frame.groupby('sample', observed=True).size()

    result['n_donors'] = int(len(per_donor))
    result['n_disease'] = int((per_donor['role'] == 'disease').sum())
    result['n_control'] = int((per_donor['role'] == 'control').sum())
    result['n_excluded'] = int((per_donor['role'] == 'exclude').sum())

    if study_col:
        result['study_col'] = study_col
        rows = []
        for study, group in per_donor.groupby('study', observed=True):
            rows.append({
                'study': str(study),
                'n_donors': int(len(group)),
                'n_disease': int((group['role'] == 'disease').sum()),
                'n_control': int((group['role'] == 'control').sum()),
                'n_cells': int(counts[group.index].sum()),
            })
        rows.sort(key=lambda r: r['study'])
        result['studies'] = rows
        result['n_studies'] = len(rows)
        result['is_atlas'] = len(rows) > 1
    else:
        result['n_studies'] = 1

    return result


def cohort_headline(shape: dict) -> tuple:
    """Two short lines describing a cohort, for a summary card.

    ('Atlas -- 2 studies', '24 donors: 12 disease / 12 control')
    ('Single study', '14 donors: 6 disease / 8 control')
    """
    if not shape or not shape.get('n_donors'):
        return ('No donors resolved', 'Set the sample column in Setup')

    if shape['is_atlas']:
        names = [row['study'] for row in shape['studies']]
        # Name them while the list is short enough to read; the point is
        # to recognise the cohort at a glance, not to enumerate it.
        if 0 < len(names) <= 3:
            first = "Atlas \u2014 " + " + ".join(names)
        else:
            first = f"Atlas \u2014 {shape['n_studies']} studies"
    else:
        first = "Single study"

    second = f"{shape['n_donors']} donors: "
    if shape['n_disease'] or shape['n_control']:
        second += f"{shape['n_disease']} disease / {shape['n_control']} control"
    else:
        second += "no roles set"
    if shape['n_excluded']:
        second += f" ({shape['n_excluded']} excluded)"
    return (first, second)


def sample_overview(adata, sample_col: Optional[str] = None) -> dict:
    """Headline counts plus a per-sample QC table for the Samples overview.

    Returns {'n_cells', 'n_genes', 'n_samples', 'n_conditions', 'table'}.
    'table' has one row per sample: 'sample', 'n_cells', and, when
    derivable, 'condition' (from '_role' else 'condition'),
    'median_genes' (per-cell genes detected), 'pct_mito' (mean per-cell
    mitochondrial %), and 'pct_doublets' (when Scrublet has run).
    Metrics use the scanpy QC obs columns when present and fall back to
    cheap sparse ops on X; every field degrades gracefully.
    """
    import numpy as np

    result = {'n_cells': 0, 'n_genes': 0, 'n_samples': 0,
              'n_conditions': 0, 'table': pd.DataFrame()}
    if adata is None:
        return result

    obs = adata.obs
    result['n_cells'] = int(adata.n_obs)
    result['n_genes'] = int(adata.n_vars)

    condition_col = ('_role' if '_role' in obs.columns
                     else 'condition' if 'condition' in obs.columns
                     else None)
    if condition_col:
        result['n_conditions'] = int(obs[condition_col].nunique())

    if sample_col is None:
        for cand in _SAMPLE_COL_CANDIDATES:
            if cand in obs.columns:
                sample_col = cand
                break
    if sample_col is None or sample_col not in obs.columns:
        return result

    per_cell = pd.DataFrame({'sample': obs[sample_col].astype(str).values})

    # Genes detected per cell
    if 'n_genes_by_counts' in obs.columns:
        per_cell['genes'] = obs['n_genes_by_counts'].to_numpy(dtype=float)
    else:
        try:
            import scipy.sparse as sp
            X = adata.X
            per_cell['genes'] = (
                np.asarray(X.getnnz(axis=1), dtype=float) if sp.issparse(X)
                else (np.asarray(X) > 0).sum(axis=1).astype(float))
        except (AttributeError, TypeError, MemoryError):
            pass

    # Mitochondrial % per cell
    if 'pct_counts_mt' in obs.columns:
        per_cell['mito'] = obs['pct_counts_mt'].to_numpy(dtype=float)
    else:
        try:
            mask = np.asarray(
                adata.var_names.str.upper().str.startswith('MT-'))
            if mask.any():
                X = adata.X
                total = np.asarray(X.sum(axis=1)).ravel().astype(float)
                mt = np.asarray(
                    X[:, np.where(mask)[0]].sum(axis=1)).ravel().astype(float)
                with np.errstate(divide='ignore', invalid='ignore'):
                    per_cell['mito'] = np.where(
                        total > 0, mt / total * 100, 0.0)
        except (AttributeError, TypeError, MemoryError):
            pass

    if condition_col:
        per_cell['condition'] = obs[condition_col].astype(str).values
    if 'predicted_doublet' in obs.columns:
        per_cell['doublet'] = obs['predicted_doublet'].astype(bool).values

    grouped = per_cell.groupby('sample', observed=True)
    table = pd.DataFrame({'n_cells': grouped.size()})
    if 'condition' in per_cell.columns:
        table['condition'] = grouped['condition'].agg(
            lambda s: s.mode().iat[0] if len(s.mode()) else s.iat[0])
    if 'genes' in per_cell.columns:
        table['median_genes'] = grouped['genes'].median()
    if 'mito' in per_cell.columns:
        table['pct_mito'] = grouped['mito'].mean()
    if 'doublet' in per_cell.columns:
        table['pct_doublets'] = grouped['doublet'].mean() * 100

    table = table.reset_index().rename(columns={'index': 'sample'})
    if table.columns[0] != 'sample':
        table = table.rename(columns={table.columns[0]: 'sample'})
    sort_cols = (['condition', 'sample'] if 'condition' in table.columns
                 else ['sample'])
    table = table.sort_values(sort_cols).reset_index(drop=True)

    result['n_samples'] = len(table)
    result['table'] = table
    return result


def summarise_samples(adata, sample_col: Optional[str] = None) -> pd.DataFrame:
    """One row per unique sample with consensus metadata + cell count.

    Returns an empty DataFrame if no usable sample column is available.
    """
    if adata is None:
        return pd.DataFrame()

    obs = adata.obs
    if sample_col is None:
        for cand in ('sample', 'donor_id', 'donor', 'patient',
                     'patient_id', 'subject', 'subject_id', 'orig.ident',
                     'library_id'):
            if cand in obs.columns:
                sample_col = cand
                break
    if sample_col is None or sample_col not in obs.columns:
        return pd.DataFrame()

    first = obs.groupby(sample_col, observed=True).first()
    n_cells = obs.groupby(sample_col, observed=True).size()
    first.insert(0, 'n_cells', n_cells)
    first = first.reset_index()
    cols = [sample_col, 'n_cells'] + [
        c for c in first.columns if c not in (sample_col, 'n_cells')]
    return first[cols]
