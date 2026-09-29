# Batch differential expression across cell types.
#
# One cell type at a time is the statistically correct unit for this kind
# of study: a disease effect in endothelium and a disease effect in
# cardiomyocytes are different contrasts, and pooling them into a single
# pseudobulk per donor averages them away. The consequence is a lot of
# repetitive subsetting -- for a meta-vs-mega design, every cell type has
# to be run once on the pooled master (mega) and once per study (meta),
# which is dozens of manual runs.
#
# This module does that loop. It implements no statistics of its own:
# each cell type is handed to 'run_de_pipeline' exactly as the DE page
# would hand it a manually-subset dataset, and the results are written to
# the standard per-study layout so the Meta workspace discovers them
# without extra plumbing.
#
# Naming: results land as '{accession}_{slug}_DE_{method}.csv' beside the
# whole-dataset result. 'discover_de_results' reads the part before
# '_DE_' as the accession, so each cell type appears in the Meta
# workspace as its own entry ('GSE292067_Endothelial') -- which is what
# you want, since pooling happens across studies within a cell type,
# never across cell types.
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional, Sequence

import numpy as np
import pandas as pd

from kosmic import (
    DEFAULT_FDR, DEFAULT_LFC_THRESHOLD, DE_MIN_CELLS, DE_MIN_COUNTS,
)
from kosmic.paths import (
    de_result_path, de_stats_dir, pseudobulk_path, significant_path,
)


# DESeq2 needs at least this many donors per arm to fit a dispersion; below
# it the run is skipped with a reason rather than producing a fit nothing
# should be read off.
MIN_SAMPLES_PER_ARM = 2


# Method label used in the DE CSV filename. The meta-analysis discovery
# regex only recognises this fixed vocabulary, so it lives in one place
# rather than being spelled out at each call site.
_METHOD_LABELS = {
    ('ttest_raw', False): 'welch_raw',
    ('ttest', False): 'welch_cpm',
    ('ttest', True): 'welch_cpm_eb',
    ('deseq2', False): 'deseq2',
    ('deseq2', True): 'deseq2',
}


def method_label(de_method: str, moderate: bool = False) -> str:
    """Filename token for a DE method: ('ttest', True) -> 'welch_cpm_eb'."""
    return _METHOD_LABELS.get((de_method, bool(moderate)), de_method)


def slugify_cell_type(name: str) -> str:
    """Filename-safe token for a cell type.

    Spaces and punctuation collapse to underscores. A literal '_DE_' is
    neutralised because 'discover_de_results' splits the filename on the
    first occurrence of it, so a cell type containing 'DE' would
    otherwise truncate the accession.
    """
    slug = re.sub(r'[^0-9A-Za-z]+', '_', str(name)).strip('_')
    slug = re.sub(r'_{2,}', '_', slug)
    slug = slug.replace('_DE_', '_De_')
    if slug.startswith('DE_'):
        slug = 'De_' + slug[3:]
    if slug.endswith('_DE'):
        slug = slug[:-3] + '_De'
    return slug or 'unnamed'


@dataclass
class CellTypePlan:
    """What a batch run would do for one cell type, before running it."""

    cell_type: str
    slug: str
    n_cells: int
    n_disease_samples: int
    n_control_samples: int
    eligible: bool
    reason: str = ''

    @property
    def n_samples(self) -> int:
        return self.n_disease_samples + self.n_control_samples


@dataclass
class CellTypeRun:
    """What a batch run actually did for one cell type."""

    cell_type: str
    slug: str
    accession: str
    status: str  # 'ok' | 'skipped' | 'failed'
    message: str = ''
    n_cells: int = 0
    n_samples: int = 0
    n_genes_tested: int = 0
    n_significant: int = 0
    de_path: Optional[Path] = None
    significant_path: Optional[Path] = None
    pseudobulk_path: Optional[Path] = None
    design: str = ''
    de_results: Optional[pd.DataFrame] = field(default=None, repr=False)


def _role_series(adata) -> pd.Series:
    if '_role' not in adata.obs.columns:
        raise ValueError(
            "adata.obs is missing the '_role' column. Set Disease / Control / "
            "Exclude in the Inspect tab (Sample Setup), hit Save, then re-run.")
    return adata.obs['_role'].astype(str).str.lower()


def included_roles(roles: pd.Series) -> pd.Series:
    """Boolean mask for the cells that can enter a contrast."""
    return roles.isin(('disease', 'control'))


def plan_cell_types(adata, cell_type_col: str, sample_col: str,
                    cell_types: Optional[Sequence[str]] = None,
                    min_cells: int = DE_MIN_CELLS,
                    min_samples_per_arm: int = MIN_SAMPLES_PER_ARM,
                    min_counts: int = DE_MIN_COUNTS,
                    counts_layer: Optional[str] = None,
                    ) -> list[CellTypePlan]:
    """Report, per cell type, whether a per-type DE run is worth attempting.

    Donor counts mirror what 'create_pseudobulk' will actually keep: a
    donor contributing fewer than *min_cells* cells or fewer than
    *min_counts* transcripts of this type is dropped there, so it is not
    counted here either. A cell type is eligible when both arms retain
    at least *min_samples_per_arm* donors.

    Cells whose role is 'exclude' are ignored throughout -- they never
    enter a contrast, so they must not make a cell type look testable.
    """
    if cell_type_col not in adata.obs.columns:
        raise ValueError(f"Cell-type column '{cell_type_col}' not in obs.")
    if sample_col not in adata.obs.columns:
        raise ValueError(f"Sample column '{sample_col}' not in obs.")

    from kosmic.de.de_analysis import cell_depths

    roles = _role_series(adata)
    labels = adata.obs[cell_type_col].astype(str)
    samples = adata.obs[sample_col].astype(str)
    depths = cell_depths(adata, counts_layer) if min_counts else None

    included = included_roles(roles)
    # Unlabelled cells: 'nan' under pandas 2's astype(str), NaN under pandas 3.
    observed = [t for t in pd.unique(labels[included])
                if pd.notna(t) and t not in ('nan', 'None')]
    wanted = list(cell_types) if cell_types is not None else sorted(observed)

    plans: list[CellTypePlan] = []
    for cell_type in wanted:
        in_type = included & (labels == str(cell_type))
        n_cells = int(in_type.sum())
        if not n_cells:
            plans.append(CellTypePlan(
                str(cell_type), slugify_cell_type(cell_type), 0, 0, 0,
                eligible=False, reason='no cells with this label'))
            continue

        per_cell = pd.DataFrame({
            'sample': samples[in_type].values,
            'role': roles[in_type].values,
            'depth': depths[in_type.values] if depths is not None else 0.0,
        })
        by_donor = per_cell.groupby(['sample', 'role'], observed=True)['depth']
        counts = by_donor.size()
        passing = counts >= min_cells
        if min_counts:
            passing &= by_donor.sum() >= min_counts
        kept = counts[passing]
        n_disease = int(sum(1 for (_s, r) in kept.index if r == 'disease'))
        n_control = int(sum(1 for (_s, r) in kept.index if r == 'control'))

        eligible = (n_disease >= min_samples_per_arm
                    and n_control >= min_samples_per_arm)
        reason = ''
        if not eligible:
            floor = f">={min_cells} cells"
            if min_counts:
                floor += f" and >={min_counts:,} transcripts"
            reason = (f"needs >={min_samples_per_arm} donors per arm with "
                      f"{floor}; has {n_disease} disease / "
                      f"{n_control} control")
        plans.append(CellTypePlan(
            str(cell_type), slugify_cell_type(cell_type), n_cells,
            n_disease, n_control, eligible, reason))

    return plans


def pseudobulk_frame(matrix: np.ndarray, sample_df: pd.DataFrame,
                     genes: Sequence[str]) -> pd.DataFrame:
    """Shape a pseudobulk matrix into the CSV the meta-analysis reads back.

    Index is the donor id; the leading columns are 'condition', 'n_cells',
    'total_counts' (the donor's summed transcripts over all genes, when
    'create_pseudobulk' supplied it) and -- when the source h5ad had
    roles set -- 'role'. 'role' is the
    canonical disease/control assignment: the consensus case-control
    permutation refuses to run without it rather than guess which raw
    condition label is the disease arm.
    """
    df = pd.DataFrame(
        matrix, columns=[str(g) for g in genes],
        index=[str(s) for s in sample_df['sample'].values])
    df.insert(0, 'condition', [str(c) for c in sample_df['condition'].values])
    df.insert(1, 'n_cells', sample_df['n_cells'].values.astype(int))
    pos = 2
    if 'total_counts' in sample_df.columns:
        df.insert(pos, 'total_counts',
                  sample_df['total_counts'].values.astype(int))
        pos += 1
    if 'role' in sample_df.columns:
        df.insert(pos, 'role', [str(r) for r in sample_df['role'].values])
    return df


def run_de_by_cell_type(adata, cell_type_col: str, sample_col: str,
                        condition_col: str, output_dir,
                        accession: str,
                        cell_types: Optional[Sequence[str]] = None,
                        de_method: str = 'deseq2',
                        moderate: bool = False,
                        pathway_gene_sets: Optional[dict] = None,
                        fdr_genes=None,
                        covariates: Optional[Sequence[str]] = None,
                        min_cells: int = DE_MIN_CELLS,
                        min_samples_per_arm: int = MIN_SAMPLES_PER_ARM,
                        min_counts: int = DE_MIN_COUNTS,
                        fdr: float = DEFAULT_FDR,
                        lfc: float = DEFAULT_LFC_THRESHOLD,
                        write: bool = True,
                        progress_callback: Optional[Callable[[str], None]] = None,
                        cancelled: Optional[Callable[[], bool]] = None,
                        **pipeline_kwargs) -> list[CellTypeRun]:
    """Run the DE pipeline once per cell type and write per-type results.

    Each cell type is subset out of *adata* and handed to
    'run_de_pipeline' unchanged, so a batch run and the equivalent manual
    run on a hand-subset dataset produce the same numbers. Ineligible
    cell types (see 'plan_cell_types') are skipped, and a failure in one
    type does not abort the rest -- both are reported in the returned
    rows.

    *covariates* is what makes the mega arm a mega arm: pass ['study']
    when running on a combined master so the design becomes
    '~study + condition'. Leave it empty for the per-study meta arm.

    *pathway_gene_sets* annotates the results' 'pathways' column.
    *fdr_genes* is hypothesis mode: the fit and normalisation stay
    genome-wide -- DESeq2 needs the whole transcriptome for its size
    factors and dispersion trend -- and only the BH denominator narrows
    to the committed list, with results subset to it. Passing neither
    runs discovery, which is what this did for every mode before they
    were threaded through.

    Extra keyword arguments are forwarded to 'run_de_pipeline', so
    detection thresholds and DESeq2 filter settings can be kept identical
    across arms.
    """
    from kosmic.de.de_analysis import (
        create_pseudobulk, run_de_pipeline, significant_subset,
    )

    # write=False is a dry run -- computing results without a place to
    # put them is legitimate, so do not demand a path for one.
    output_dir = Path(output_dir) if output_dir is not None else None
    if write and output_dir is None:
        raise ValueError("output_dir is required when write=True.")
    label = method_label(de_method, moderate)

    def _progress(msg):
        if progress_callback:
            progress_callback(msg)

    plans = plan_cell_types(
        adata, cell_type_col, sample_col, cell_types=cell_types,
        min_cells=min_cells, min_samples_per_arm=min_samples_per_arm,
        min_counts=min_counts,
        counts_layer=pipeline_kwargs.get('counts_layer'))

    if write:
        de_stats_dir(output_dir).mkdir(parents=True, exist_ok=True)

    roles = _role_series(adata)
    labels = adata.obs[cell_type_col].astype(str)
    runs: list[CellTypeRun] = []

    for i, plan in enumerate(plans, start=1):
        run_accession = f"{accession}_{plan.slug}"
        head = f"[{i}/{len(plans)}] {plan.cell_type}"

        if cancelled is not None and cancelled():
            runs.append(CellTypeRun(
                plan.cell_type, plan.slug, run_accession, 'skipped',
                'cancelled before this cell type ran'))
            continue

        if not plan.eligible:
            _progress(f"{head}: skipped -- {plan.reason}")
            runs.append(CellTypeRun(
                plan.cell_type, plan.slug, run_accession, 'skipped',
                plan.reason, n_cells=plan.n_cells))
            continue

        _progress(f"{head}: {plan.n_cells:,} cells, "
                  f"{plan.n_disease_samples} disease / "
                  f"{plan.n_control_samples} control donors")

        # 'exclude' cells are dropped here rather than left to the
        # pipeline, so the detection pre-filter sees only cells that can
        # influence the contrast.
        mask = (labels == plan.cell_type) & included_roles(roles)
        subset = adata[mask.values].copy()

        try:
            de_results, _significant, _coverage, sample_df = run_de_pipeline(
                subset, sample_col, condition_col, pathway_gene_sets or {},
                min_cells=min_cells, min_counts=min_counts,
                de_method=de_method, full_genome=True,
                moderate=moderate, covariates=covariates,
                fdr_genes=fdr_genes,
                progress_callback=lambda m, h=head: _progress(f"{h}: {m}"),
                **pipeline_kwargs)
        except Exception as exc:  # noqa: BLE001 - one bad type must not stop the batch
            _progress(f"{head}: FAILED -- {exc}")
            runs.append(CellTypeRun(
                plan.cell_type, plan.slug, run_accession, 'failed',
                str(exc), n_cells=plan.n_cells))
            continue

        if de_results is None or de_results.empty:
            _progress(f"{head}: no genes tested")
            runs.append(CellTypeRun(
                plan.cell_type, plan.slug, run_accession, 'skipped',
                'no genes survived filtering', n_cells=plan.n_cells))
            continue

        n_sig = int(((de_results['pvals_adj'] < fdr)
                     & (de_results['abs_logfoldchange'] > lfc)).sum())
        run = CellTypeRun(
            plan.cell_type, plan.slug, run_accession, 'ok',
            n_cells=plan.n_cells, n_samples=len(sample_df),
            n_genes_tested=len(de_results), n_significant=n_sig,
            design=str(de_results.attrs.get('design', '')),
            de_results=de_results)

        if write:
            run.de_path = de_result_path(output_dir, run_accession, label)
            de_results.to_csv(run.de_path, index=False)

            # The significant subset alongside the full table -- the same
            # selection the counts above report, so nobody has to re-derive
            # it with slightly different thresholds.
            run.significant_path = significant_path(output_dir, run_accession)
            significant_subset(de_results, fdr, lfc).to_csv(
                run.significant_path, index=False)

            gene_names = [str(g) for g in de_results['names'].tolist()]
            pb_matrix, pb_sample_df, pb_genes = create_pseudobulk(
                subset, gene_names, sample_col, condition_col,
                min_cells=min_cells, min_counts=min_counts, aggregate='sum')
            if pb_matrix.size:
                run.pseudobulk_path = pseudobulk_path(output_dir, run_accession)
                pseudobulk_frame(pb_matrix, pb_sample_df, pb_genes).to_csv(
                    run.pseudobulk_path)

        _progress(f"{head}: {n_sig:,} significant of {len(de_results):,} tested"
                  + (f" ({run.design})" if run.design else ''))
        runs.append(run)
        del subset

    return runs



def batch_summary(runs: Sequence[CellTypeRun]) -> pd.DataFrame:
    """One row per cell type: what ran, what it found, where it landed."""
    return pd.DataFrame([{
        'cell_type': r.cell_type,
        'accession': r.accession,
        'status': r.status,
        'n_cells': r.n_cells,
        'n_samples': r.n_samples,
        'n_genes_tested': r.n_genes_tested,
        'n_significant': r.n_significant,
        'design': r.design,
        'message': r.message,
        'de_path': str(r.de_path) if r.de_path else '',
    } for r in runs])


def common_tested_genes(runs_or_frames) -> list[str]:
    """Genes tested by every run -- the only fair set to compare arms on.

    Each arm applies its own detection pre-filter, so their gene
    universes differ. Comparing without intersecting first conflates
    "not significant here" with "never tested here".
    """
    sets = []
    for item in runs_or_frames:
        df = item.de_results if isinstance(item, CellTypeRun) else item
        if df is None or getattr(df, 'empty', True) or 'names' not in df.columns:
            continue
        sets.append({str(g) for g in df['names']})
    if not sets:
        return []
    return sorted(set.intersection(*sets))
