# I/O helpers for the meta-analysis subsystem.
#
# Three entry points, all pure-logic (no PyQt6 imports):
#
# * 'discover_de_results' — walk a project folder, return a structured
#   list of every DE-related CSV found. The bridge between the scRNA side
#   (which writes CSVs in a standard layout) and the meta-analysis side
#   (which needs to locate inputs).
# * 'load_de_results' — load one DE CSV into a canonical-schema
#   DataFrame (standardised column names, 'dataset' column added).
# * 'load_pseudobulk_set' — load a list of pseudobulk CSVs into the
#   dict structure consumed by 'cc_permutation' and consensus LOO.
#
# 'split_cell_types' names the cell type of each per-cell-type result, so
# the Select Studies page can offer one cell type across studies.
#
# Plus 'require_se' — small validator for "is this DataFrame usable
# for meta-analysis pooling?" Fails loudly with a helpful remediation
# message rather than silently fabricating weights.
from __future__ import annotations

import re
from pathlib import Path
from typing import Callable, List, Optional, Sequence

import pandas as pd

from kosmic.paths import de_stats_dir, de_pathway_scoring_dir


# Filesystem discovery: find DE CSVs written by the scRNA pipeline

_GENE_STATS_PATTERN = re.compile(r"^(.+?)_GeneStats_(Norm|Raw)\.csv$")
_DE_RESULTS_PATTERN = re.compile(
    r"^(.+?)_DE_(welch_raw|welch_cpm|welch_cpm_eb|welch_ttest|welch_ttest_eb|deseq2)\.csv$"
)
_PATHWAY_DE_PATTERN = re.compile(r"^(?:parent_)?pathway_de_results\.csv$")


def discover_de_results(project_folder: Path) -> list[dict]:
    """Walk immediate subdirectories of *project_folder* for DE CSVs.

    Expected layout per analysis folder::

        {project_folder}/{analysis_name}/results/de_analysis/statistics/
            {accession}_GeneStats_{Norm|Raw}.csv
            {accession}_DE_{method}.csv
        {project_folder}/{analysis_name}/results/de_analysis/pathway_scoring/
            pathway_de_results.csv

    Returns
    -------
    list[dict]
        One dict per file with keys:
        'accession' (str), 'file_type' ('"gene"', '"de_results"',
        or '"pathway_de"'), 'variant' ('"Norm"'/'"Raw"' for gene
        stats, DE method name for de_results, '"pathway_de"'), 'path'
        (Path), 'analysis_folder' (str), and 'de_method' (str) for
        de_results entries.
    """
    project_folder = Path(project_folder)
    if not project_folder.is_dir():
        return []

    results: list[dict] = []

    for child in sorted(project_folder.iterdir()):
        if not child.is_dir():
            continue

        # Gene-stats + per-method DE tables in statistics/
        stats_dir = de_stats_dir(child)
        if stats_dir.is_dir():
            for csv_file in sorted(stats_dir.iterdir()):
                if not csv_file.is_file() or csv_file.suffix.lower() != ".csv":
                    continue

                name = csv_file.name
                gene_match = _GENE_STATS_PATTERN.match(name)
                de_match = _DE_RESULTS_PATTERN.match(name)

                if gene_match:
                    results.append({
                        "accession": gene_match.group(1),
                        "file_type": "gene",
                        "variant": gene_match.group(2),
                        "path": csv_file.resolve(),
                        "analysis_folder": child.name,
                    })
                elif de_match:
                    results.append({
                        "accession": de_match.group(1),
                        "file_type": "de_results",
                        "variant": de_match.group(2),
                        "de_method": de_match.group(2),
                        "path": csv_file.resolve(),
                        "analysis_folder": child.name,
                    })

        # Pathway-level DE results in pathway_scoring/
        scoring_dir = de_pathway_scoring_dir(child)
        if scoring_dir.is_dir():
            # Try to infer accession from sibling stats files; fall back
            # to the folder name.
            accession = child.name
            for entry in results:
                if entry["analysis_folder"] == child.name and entry.get("accession"):
                    accession = entry["accession"]
                    break

            for csv_file in sorted(scoring_dir.iterdir()):
                if not csv_file.is_file() or csv_file.suffix.lower() != ".csv":
                    continue
                if _PATHWAY_DE_PATTERN.match(csv_file.name):
                    results.append({
                        "accession": accession,
                        "file_type": "pathway_de",
                        "variant": "pathway_de",
                        "path": csv_file.resolve(),
                        "analysis_folder": child.name,
                    })

    return results


# Choosing which discovered files a run actually uses

def select_de_entries(entries: Sequence[dict],
                      preferred_methods: Sequence[str],
                      allowed: Optional[Sequence[str]] = None,
                      ) -> tuple[list[dict], list[str], list[str]]:
    """Pick one DE-results entry per accession, honouring a study selection.

    'discover_de_results' walks the whole project folder, so it returns
    every study it can see -- including ones the user unticked, and
    per-cell-type results sitting beside a whole-dataset result. Pooling
    is only meaningful over the studies actually chosen, so *allowed*
    gates the accessions that survive. 'None' (or an empty selection)
    means no gate has been set yet and everything discovered is used.

    Within an accession, the first method in *preferred_methods* that
    exists wins. If none does, the accession still runs on whatever it
    has and is reported in *fallback*, since a study analysed with a
    different engine is usually better than a silently missing study.

    Returns
    -------
    (chosen, fallback, skipped)
        'chosen' is one entry per surviving accession, ordered by
        accession. 'fallback' holds '"accession:method"' strings for
        accessions that had no preferred method. 'skipped' holds
        accessions dropped because they were not in *allowed*.
    """
    allowed_set = {str(a) for a in allowed} if allowed else None

    by_accession: dict[str, list[dict]] = {}
    for entry in entries:
        if entry.get("file_type") != "de_results":
            continue
        by_accession.setdefault(entry["accession"], []).append(entry)

    chosen: list[dict] = []
    fallback: list[str] = []
    skipped: list[str] = []

    for accession in sorted(by_accession):
        if allowed_set is not None and accession not in allowed_set:
            skipped.append(accession)
            continue
        candidates = by_accession[accession]
        match = None
        for method in preferred_methods:
            match = next((e for e in candidates
                          if e.get("de_method") == method), None)
            if match is not None:
                break
        if match is None:
            match = candidates[0]
            fallback.append(f"{accession}:{match.get('de_method', '?')}")
        chosen.append(match)

    return chosen, fallback, skipped


def split_cell_types(accessions_by_folder: dict[str, Sequence[str]]
                     ) -> dict[str, Optional[str]]:
    """Name the cell type of each accession, or None for a whole-study result.

    Per-cell-type DE writes '{accession}_{cell type}' beside the study's
    own result ('kosmic.de.batch'), and cell-type names contain
    underscores themselves ('CD8_T_Cell'), so the name cannot be split on
    '_'. Within one study folder, an accession is a cell type of that
    study when it is the folder name, or another accession of the same
    folder, followed by '_'; the longest such prefix wins.

    Parameters
    ----------
    accessions_by_folder : dict
        '{study folder name: [accession, ...]}', as 'discover_de_results'
        reports them ('analysis_folder' and 'accession').

    Returns
    -------
    dict
        '{accession: cell type or None}'.
    """
    out: dict[str, Optional[str]] = {}
    for folder, accessions in accessions_by_folder.items():
        prefixes = {folder, *accessions}
        for acc in accessions:
            owners = [p for p in prefixes
                      if p != acc and acc.startswith(p + '_')]
            if owners:
                owner = max(owners, key=len)
                out[acc] = acc[len(owner) + 1:]
            else:
                out[acc] = None
    return out


# Per-study DE CSV loader

def load_de_results(filepath: str,
                    dataset_name: Optional[str] = None) -> pd.DataFrame:
    """Load a DE-results CSV with canonicalised column names.

    Renames common column variants ('gene', 'logFC', 'pvalue',
    'padj', 'lfcSE', …) to the project-canonical names ('names',
    'logfoldchanges', 'pvals', 'pvals_adj', 'se', 'pathways').
    A 'dataset' column is added, taking the file stem by default.
    """
    df = pd.read_csv(filepath)

    if dataset_name is None:
        dataset_name = Path(filepath).stem
    df['dataset'] = dataset_name

    column_mapping = {
        'gene': 'names', 'gene_name': 'names', 'Gene': 'names',
        'logFC': 'logfoldchanges', 'log2FoldChange': 'logfoldchanges',
        'logfoldchange': 'logfoldchanges',
        'pvalue': 'pvals', 'p_value': 'pvals', 'pval': 'pvals',
        'padj': 'pvals_adj', 'p_adj': 'pvals_adj', 'FDR': 'pvals_adj',
        'qvalue': 'pvals_adj',
        'pathway': 'pathways', 'Pathway': 'pathways',
        # Standard error — common variants from DESeq2 (lfcSE), edgeR/limma
        # (stderr/std_error), upstream KOSMIC outputs (se), manual CSVs.
        'SE': 'se', 'stderr': 'se', 'std_error': 'se', 'lfcSE': 'se',
        'stdError': 'se',
    }
    for old, new in column_mapping.items():
        if old in df.columns and new not in df.columns:
            df.rename(columns={old: new}, inplace=True)

    required = ['names', 'logfoldchanges']
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns: {missing}\n"
            f"Available: {list(df.columns)}")

    return df


def require_se(df: pd.DataFrame, filepath: str = '<unknown>') -> None:
    """Validate that 'df' has a usable 'se' column for meta-analysis.

    Raises 'ValueError' with a helpful message if SE is missing or entirely
    NaN. Meta-analysis pooling needs per-gene SE to compute inverse-variance
    weights; there is no reliable way to impute it from p-values or effect
    sizes alone, so we fail loudly rather than fabricate weights.
    """
    if 'se' not in df.columns:
        raise ValueError(
            f"{filepath}: DE table is missing an 'se' (standard error) column. "
            f"Meta-analysis pooling requires per-gene SE. If your CSV has SE "
            f"under a different header, rename it to one of: "
            f"se, SE, stderr, std_error, lfcSE, stdError."
        )
    if df['se'].isna().all():
        raise ValueError(
            f"{filepath}: 'se' column is present but entirely NaN. "
            f"Meta-analysis pooling requires per-gene SE.")


# Pseudobulk-set loader (for CC permutation + consensus LOO)

def load_pseudobulk_set(
    pb_paths: Sequence,
    labels: Optional[Sequence[str]] = None,
    progress_cb: Optional[Callable[[str], None]] = None,
    skip_errors: bool = False,
) -> List[dict]:
    """Load a list of pseudobulk CSVs into the dict structure consumed by
    'kosmic.meta_analysis.cc_permutation' routines.

    Parameters
    ----------
    pb_paths : sequence of str or Path
        Paths to pseudobulk CSVs (one per study).
    labels : sequence of str, optional
        Per-study labels. Defaults to 'Study_{i}' for unlabelled slots.
    progress_cb : callable, optional
        Invoked with a short message per study — both on successful load
        ("Loaded pseudobulk: NAME (N samples, settings: found)") and (if
        'skip_errors=True') on a skipped failure.
    skip_errors : bool
        If 'True', per-path load failures are reported via 'progress_cb'
        (if set) and the path is skipped. If 'False', the first failure
        re-raises.

    Returns
    -------
    list of dict
        One dict per successfully-loaded path, with keys: 'expression'
        (ndarray), 'conditions' (ndarray of str), 'gene_names' (list
        of str), 'n_cells' (ndarray), 'dataset_name' (str), 'role'
        (ndarray of str: 'disease' / 'control' / 'exclude' per sample).

    Notes
    -----
    The 'role' array is the canonical disease/control assignment
    downstream code consumes. It comes from the pseudobulk CSV's own
    'role' column, written by 'create_pseudobulk' from the source
    h5ad's 'adata.obs['_role']'. If the column is absent, 'role' is
    'None' and downstream code raises a clear error pointing the user
    back to the Inspect tab to set roles.
    """
    from kosmic.meta_analysis.cc_permutation import load_pseudobulk

    pb_data: List[dict] = []
    labels_seq = list(labels) if labels else []

    for i, pb_path in enumerate(pb_paths):
        try:
            expr, conds, genes, samples, n_cells, roles = load_pseudobulk(pb_path)
        except Exception as e:  # noqa: BLE001
            if skip_errors:
                if progress_cb is not None:
                    progress_cb(f"Failed to load {pb_path}: {e}")
                continue
            raise

        ds_name = (labels_seq[i] if i < len(labels_seq)
                   else f'Study_{i}')
        pb_data.append({
            'expression': expr,
            'conditions': conds,
            'gene_names': genes,
            'n_cells': n_cells,
            'dataset_name': ds_name,
            'role': roles,
        })

        if progress_cb is not None:
            progress_cb(
                f"  Loaded pseudobulk: {ds_name} ({len(samples)} samples)")

    return pb_data


# Per-result settings sidecar (companion JSON next to an exported meta CSV)

def write_meta_settings_sidecar(path, *, tool, analysis_type, parameters,
                                results_summary) -> None:
    """Write a settings sidecar JSON alongside an exported meta-analysis CSV.

    Captures the run parameters, a short result summary and the pooling-stack
    software versions, so each saved CSV carries its exact provenance.
    """
    import json
    import sys
    from datetime import datetime

    def _ver(pkg):
        try:
            import importlib
            return getattr(importlib.import_module(pkg), '__version__', None)
        except ImportError:
            return None

    payload = {
        'timestamp': datetime.now().isoformat(timespec='seconds'),
        'tool': tool,
        'analysis_type': analysis_type,
        'parameters': parameters,
        'results_summary': results_summary,
        'software_versions': {
            'python': sys.version.split()[0],
            'numpy': _ver('numpy'),
            'pandas': _ver('pandas'),
            'scipy': _ver('scipy'),
            'statsmodels': _ver('statsmodels'),
        },
    }
    with open(path, 'w') as f:
        json.dump(payload, f, indent=2)


# Provenance aggregation: combined cross-study methods

ALL_CELLS = 'all_cells'
MIXED = 'mixed'


def meta_selection_folder(cell_types: Sequence[Optional[str]]) -> str:
    """Name the output folder for a meta-analysis selection (ADR-007).

    Parameters
    ----------
    cell_types : sequence of str or None
        The cell type of each selected project result, None for a
        whole-study result (see 'split_cell_types'). External imports are
        left out by the caller: they have no cell type.

    Returns
    -------
    str
        The cell type when every result is of one cell type; 'all_cells'
        when all are whole-study results or there are none; 'mixed' when
        the selection spans more than one.
    """
    kinds = set(cell_types)
    if not kinds or kinds == {None}:
        return ALL_CELLS
    if len(kinds) == 1:
        return next(iter(kinds))
    return MIXED


def mixed_selection_warning(selection: Optional[str]) -> Optional[str]:
    """The output-panel warning for a selection spanning cell types, or None."""
    if selection != MIXED:
        return None
    return ("Warning: the selected results span more than one cell type, so "
            "they were pooled together and saved under meta_analysis/mixed/. "
            "Choose one cell type on Select Studies to pool it on its own.")


def meta_provenance_dir(project_folder, selection: Optional[str] = None) -> Path:
    """Directory holding the meta-analysis provenance record of a selection.

    Without 'selection', the top-level folder earlier versions used.
    """
    from kosmic.paths import meta_output_dir
    return meta_output_dir(project_folder, selection)


def gather_study_provenance(project_folder, study_names=None, *,
                            only_named: bool = False) -> list:
    """(label, provenance_record) for each study folder with a sidecar.

    Provenance lives next to the study's h5ad in 'processed_data/'. Walks the
    immediate subfolders of the project (skipping '_'-prefixed and the
    meta-analysis output). When 'study_names' is given, only folders whose name
    is in that set are returned, in that order; otherwise all found, sorted.

    'only_named' drops the studies that were not asked for. Off by
    default, which keeps the project-wide view; the methods document
    turns it on, because a methods section should describe the studies
    that went into the analysis and not every folder that happens to
    carry a sidecar.
    """
    from kosmic import provenance
    from kosmic.paths import processed_data_dir

    project_folder = Path(project_folder)
    if not project_folder.is_dir():
        return []

    found = {}
    for child in sorted(project_folder.iterdir()):
        if (not child.is_dir() or child.name.startswith('_')
                or child.name == 'meta_analysis'):
            continue
        rec = provenance.load(processed_data_dir(child))
        if rec is not None:
            found[child.name] = rec

    if study_names:
        ordered = [(n, found[n]) for n in study_names if n in found]
        if only_named:
            return ordered
        # Include any provenance-bearing studies not named, after the named set.
        ordered += [(n, r) for n, r in found.items() if n not in set(study_names)]
        return ordered
    return list(found.items())


def _meta_staleness_note(studies, meta_rec) -> str:
    """Warn when a study's data changed since the last meta run.

    Compares each study's current data token to the token the most recent
    meta stage recorded consuming. '' when nothing is stale or unknowable.
    """
    if not meta_rec:
        return ''
    stages = meta_rec.get('stages') or []
    if not stages:
        return ''
    recorded = (stages[-1].get('params') or {}).get('study_tokens') or {}
    if not recorded:
        return ''
    current = {label: (rec.get('fingerprint') or {}).get('token')
               for label, rec in studies}
    changed = [label for label, tok in recorded.items()
               if label in current and current[label] and tok
               and current[label] != tok]
    if not changed:
        return ''
    return ("NOTE: these studies changed since the last meta-analysis run -- "
            "re-run to refresh: " + ", ".join(sorted(changed)) + "\n\n")


def record_meta_stage(project_folder, stage: str, params: dict,
                      selection: Optional[str] = None) -> bool:
    """Record one meta-analysis stage into the selection's sidecar.

    'selection' is the output folder of what was pooled
    ('meta_selection_folder'); each cell type keeps its own record
    (ADR-007).

    Every step of the discovery flow writes here -- pooling, enrichment,
    validation -- so the Methods page can render the whole run from one
    record instead of each page keeping its own account of what it did.

    Only the latest run of each stage is kept. A meta-analysis run
    overwrites its output file, so appending one entry per click built a
    history describing results that no longer existed -- 31 recorded
    pooling runs against 8 files on disk. Each result also carries its
    own sidecar next to the CSV, which records the outcome and library
    versions; this record is the current state, not an audit log.

    Best-effort: returns False rather than raising, because failing to
    write a methods note must never lose an analysis result.
    """
    if not project_folder:
        return False
    try:
        from kosmic import provenance
        name = Path(project_folder).name
        if selection:
            name = f"{name} ({selection})"
        provenance.record_stage(
            meta_provenance_dir(project_folder, selection),
            name, stage, params, replace=True)
        return True
    except Exception:
        return False


def build_combined_methods(project_folder, study_names=None, *,
                           selection: Optional[str] = None,
                           full: bool = False) -> str:
    """Render a combined methods document: each study's provenance followed by
    the meta-analysis pooling step. '' when no provenance exists anywhere.

    'study_names' are study folder names. The meta-analysis record is the
    one of 'selection' (ADR-007), falling back to the top-level record
    that earlier versions wrote when the selection has none yet.

    Reports the final configuration of each stage and summarises values
    too long to read. 'full=True' gives the unabridged audit trail.
    """
    from kosmic import provenance

    studies = gather_study_provenance(project_folder, study_names,
                                      only_named=bool(study_names))
    meta_rec = None
    if selection:
        meta_rec = provenance.load(meta_provenance_dir(project_folder, selection))
    if meta_rec is None:
        meta_rec = provenance.load(meta_provenance_dir(project_folder))
    doc = provenance.render_combined(studies, meta_rec, full=full)
    if not doc:
        return ''
    return _meta_staleness_note(studies, meta_rec) + doc
