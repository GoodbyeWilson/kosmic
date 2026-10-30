# What the DE steps save to a study's folder, read back for figures.
#
# Figure Export draws a study's DE figures from these files when the DE
# workspace does not hold the study (ADR-014). Every path is rebuilt from
# the study folder and a file name; a path stored inside a saved file
# (the batch summary's 'de_path' column, say) names the computer the
# analysis ran on and is never used, so a copied project folder reads the
# same.
#
# Most outputs were already saved by the steps that make them: the gene
# DE table per cell type and method, the pathway DE table, the per-donor
# pathway scores and the fgsea and library GSEA summaries. Three were not,
# and are saved from now on by the functions here:
#
# - each gene DE run's settings, as a JSON record beside its table
#   ('{accession}_DE_{method}.json'): the count source, the columns and
#   labels of the contrast, the gene sets selected and the thresholds;
# - an over-representation (GO / Enrichr) run's table;
# - the running-enrichment curves of a library GSEA run, which the
#   mountain figure draws.
from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from kosmic.paths import de_pathway_scoring_dir, de_stats_dir

#: An over-representation run's table, in the DE workspace's normalised form.
ENRICHMENT_RESULTS_FILE = 'enrichment_results.csv'

#: Pathway scoring levels and their file-name prefixes.
PATHWAY_LEVELS = {'': 'Pathways', 'parent_': 'Pathway groups'}


# -- Gene DE runs ------------------------------------------------------------------

@dataclass(frozen=True)
class GeneDERun:
    """One saved gene DE table.

    ``name`` is the file's accession part ('Chaffin_Endothelial_Cell');
    ``cell_type`` the cell type it was run on, or None for a whole-study
    run (or one whose cell type cannot be told from its name); ``method``
    the method token ('deseq2', 'welch_cpm_eb', ...).
    """

    name: str
    cell_type: Optional[str]
    method: str
    path: Path

    @property
    def label(self) -> str:
        return f"{self.cell_type or self.name} ({self.method})"

    @property
    def record_path(self) -> Path:
        return run_record_path(self.path)


def list_gene_de_runs(study_dir) -> list[GeneDERun]:
    """The study's saved gene DE tables, whole-study runs first, then by
    cell type and method.

    The cell type of a per-cell-type run is taken from its run record,
    else from the batch summary written beside it
    ('{accession}_by_cell_type_summary.csv'), else from the file name: the
    part after the study folder's name, compared without case (the Stroke
    project has both 'GSE174574_...' and 'gse174574_...').
    """
    return gene_de_runs_in(de_stats_dir(study_dir), Path(study_dir).name)


def gene_de_runs_in(stats, study: str) -> list[GeneDERun]:
    """:func:`list_gene_de_runs` for a statistics folder, with 'study' the
    name a per-cell-type file starts with."""
    from kosmic.meta_analysis.io import _DE_RESULTS_PATTERN

    stats = Path(stats)
    if not stats.is_dir():
        return []
    labels = _summary_cell_types(stats)
    runs = []
    for path in sorted(stats.glob('*_DE_*.csv')):
        m = _DE_RESULTS_PATTERN.match(path.name)
        if not m:
            continue
        name, method = m.group(1), m.group(2)
        record = read_run_record(path) or {}
        if 'cell_type' in record:
            cell_type = record['cell_type'] or None
        elif name in labels:
            cell_type = labels[name]
        else:
            cell_type = _cell_type_from_name(name, study)
        runs.append(GeneDERun(name=name, cell_type=cell_type, method=method, path=path))
    runs.sort(key=lambda r: (r.cell_type is not None, (r.cell_type or '').lower(),
                             r.method, r.name))
    return runs


def _summary_cell_types(stats: Path) -> dict:
    """{run accession: cell type} from the batch summaries in 'stats'."""
    out = {}
    for path in stats.glob('*_by_cell_type_summary.csv'):
        try:
            df = pd.read_csv(path, usecols=['cell_type', 'accession'])
        except (OSError, ValueError):
            continue
        out.update(zip(df['accession'].astype(str), df['cell_type'].astype(str)))
    return out


def _cell_type_from_name(name: str, study: str) -> Optional[str]:
    if name.lower() == study.lower():
        return None
    prefix = study.lower() + '_'
    if name.lower().startswith(prefix):
        return name[len(prefix):].replace('_', ' ')
    return None


def read_gene_de(run: GeneDERun) -> Optional[pd.DataFrame]:
    """The run's DE table, or None when it is missing, unreadable or has no
    'names' column."""
    try:
        df = pd.read_csv(run.path)
    except (OSError, ValueError):
        return None
    return df if 'names' in df.columns and not df.empty else None


# -- Run records ---------------------------------------------------------------

def run_record_path(de_path) -> Path:
    """'{accession}_DE_{method}.json' beside '{accession}_DE_{method}.csv'."""
    return Path(de_path).with_suffix('.json')


def build_run_record(adata, *, method: str, sample_col: str, condition_col: str,
                     disease_label: str, control_label: str,
                     cell_type: Optional[str] = None, cell_type_col: Optional[str] = None,
                     min_cells: Optional[int] = None, min_counts: Optional[int] = None,
                     covariates=(), pathway_gene_sets: Optional[dict] = None,
                     geneset_label: Optional[str] = None,
                     h5ad_path=None) -> dict:
    """The settings of a gene DE run, as saved beside its table.

    ``counts_source`` is the count matrix the run used: 'raw' for the
    original counts or a layer name ('decontX_counts'), as
    :func:`kosmic.scrna.counts.loaded_count_source` reports it. Only the
    file name of ``h5ad_path`` is kept.
    """
    from kosmic.scrna.counts import loaded_count_source

    return {
        'saved': datetime.now().isoformat(timespec='seconds'),
        'method': method,
        'cell_type': cell_type,
        'cell_type_column': cell_type_col,
        'sample_column': sample_col,
        'condition_column': condition_col,
        'disease_label': disease_label,
        'control_label': control_label,
        'counts_source': loaded_count_source(adata) if adata is not None else 'raw',
        'min_cells': min_cells,
        'min_counts': min_counts,
        'covariates': list(covariates or ()),
        'geneset_label': geneset_label or '',
        'pathway_gene_sets': {str(k): [str(g) for g in v]
                              for k, v in (pathway_gene_sets or {}).items()},
        'dataset_file': Path(str(h5ad_path)).name if h5ad_path else None,
    }


def save_run_record(de_path, record: dict) -> Path:
    path = run_record_path(de_path)
    path.write_text(json.dumps(record, indent=2), encoding='utf-8')
    return path


def read_run_record(de_path) -> Optional[dict]:
    """The run's record, or None for a run saved before records were kept
    (or one whose record cannot be read)."""
    path = run_record_path(de_path)
    try:
        data = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


# -- Pathway scoring -------------------------------------------------------------

def pathway_levels(study_dir) -> list[str]:
    """The pathway levels with a saved pathway DE table: '' (pathways)
    and 'parent_' (pathway groups), in that order."""
    d = de_pathway_scoring_dir(study_dir)
    return [lv for lv in PATHWAY_LEVELS if (d / f'{lv}pathway_de_results.csv').is_file()]


def read_pathway_de(study_dir, level: str = '') -> Optional[pd.DataFrame]:
    path = de_pathway_scoring_dir(study_dir) / f'{level}pathway_de_results.csv'
    try:
        return pd.read_csv(path)
    except (OSError, ValueError):
        return None


def read_per_donor_scores(study_dir, level: str = '') -> Optional[dict]:
    from kosmic.de.pathway_scoring import PER_DONOR_SCORES_FILE, load_per_donor_scores
    return load_per_donor_scores(
        de_pathway_scoring_dir(study_dir) / f'{level}{PER_DONOR_SCORES_FILE}')


def read_fgsea(study_dir) -> Optional[pd.DataFrame]:
    from kosmic.de.fgsea import FGSEA_RESULTS_FILE, load_fgsea_results
    return load_fgsea_results(de_pathway_scoring_dir(study_dir) / FGSEA_RESULTS_FILE)


def gsea_library_path(study_dir) -> Path:
    from kosmic.de.fgsea import GSEA_LIBRARY_RESULTS_FILE
    return de_pathway_scoring_dir(study_dir) / GSEA_LIBRARY_RESULTS_FILE


def read_gsea_library(study_dir) -> Optional[pd.DataFrame]:
    from kosmic.de.fgsea import load_fgsea_results
    return load_fgsea_results(gsea_library_path(study_dir))


# -- Library GSEA running-enrichment curves --------------------------------------

def gsea_details_path(results_path) -> Path:
    """'<results stem>_details.npz' beside a library GSEA summary."""
    p = Path(results_path)
    return p.with_name(f'{p.stem}_details.npz')


def save_gsea_details(results_path, details: dict, library_name: str = '') -> Path:
    """Save the running-enrichment curves of a library GSEA run (the
    ``details`` of :func:`kosmic.de.fgsea.run_fgsea`) beside its summary."""
    terms = list(details)
    arrays = {'terms': np.array(terms, dtype=str),
              'library': np.array(library_name or ''),
              'es': np.array([float(details[t]['es']) for t in terms]),
              'nes': np.array([float(details[t]['nes']) for t in terms])}
    for i, t in enumerate(terms):
        arrays[f'res_{i}'] = np.asarray(details[t]['RES'], dtype=float)
        arrays[f'hits_{i}'] = np.asarray(details[t]['hits'], dtype=np.int64)
    path = gsea_details_path(results_path)
    with open(path, 'wb') as fh:
        np.savez_compressed(fh, **arrays)
    return path


def read_gsea_details(results_path) -> tuple[Optional[dict], str]:
    """(details, library name) saved beside a library GSEA summary;
    (None, '') for a run saved before the curves were kept."""
    path = gsea_details_path(results_path)
    try:
        with np.load(path, allow_pickle=False) as z:
            terms = [str(t) for t in z['terms']]
            details = {t: {'RES': z[f'res_{i}'], 'hits': [int(h) for h in z[f'hits_{i}']],
                           'es': float(z['es'][i]), 'nes': float(z['nes'][i])}
                       for i, t in enumerate(terms)}
            return details, str(z['library'])
    except (OSError, KeyError, ValueError):
        return None, ''


# -- Over-representation ----------------------------------------------------------

def enrichment_path(folder) -> Path:
    return Path(folder) / ENRICHMENT_RESULTS_FILE


def save_enrichment(folder, df: pd.DataFrame) -> Path:
    path = enrichment_path(folder)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def read_enrichment(folder) -> Optional[pd.DataFrame]:
    """The over-representation table saved in 'folder', or None."""
    try:
        df = pd.read_csv(enrichment_path(folder))
    except (OSError, ValueError):
        return None
    return df if 'Term' in df.columns and not df.empty else None
