# The genes a meta-analysis tests, and the rules that chose them.
#
# Every rule here is optional and off by default, so a meta-analysis that
# sets none of them pools what it always pooled. They exist so that a gene
# set chosen for a paper (current symbols only, the expression filter, no
# gene that harmonisation merged differently in different studies) is chosen
# by KOSMIC itself and recorded, not by a script beside it.
#
# The set is saved beside the run's results, so that the merged
# (pooled-pseudobulk) arm can correct its FDR over exactly the genes the
# meta-analysis tested (kosmic.de.batch.run_de_by_cell_type, fdr_genes).
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

GENE_SET_SUFFIX = "_genes.csv"


def current_symbols(species: str = "human") -> set[str]:
    """Current approved symbols: HGNC for human, MGI for mouse."""
    from kosmic.scrna.inspect.gene_names import load_locus_groups
    return set(load_locus_groups(species=species))


def inconsistent_merge_genes(project_dir, accessions) -> tuple[set[str], list[str]]:
    """Genes that gene-name harmonisation built from a different number of
    features in different studies (``kosmic.combine.gene_merge``).

    The same check the Combine step offers when building an atlas, run on the
    studies of this meta-analysis. Returns ``(genes, unknown)``, ``unknown``
    being the studies with no harmonisation record, which cannot be checked.
    """
    from kosmic.combine.gene_merge import project_merge_asymmetry
    result = project_merge_asymmetry(project_dir, list(accessions))
    return set(result["genes"]), list(result["unknown"])


def study_folders(datasets: list[dict], project_dir) -> dict[str, Optional[str]]:
    """Dataset name -> the project's study folder its DE table came from.

    None for a dataset from outside the project (an imported table), whose
    harmonisation record KOSMIC cannot see.
    """
    root = Path(project_dir).resolve() if project_dir else None
    out = {}
    for entry in datasets:
        folder = None
        path = entry.get("path")
        if root is not None and path and not entry.get("external"):
            try:
                folder = Path(path).resolve().relative_to(root).parts[0]
            except (ValueError, IndexError):
                folder = None
        out[entry.get("name", "?")] = folder
    return out


def build_gene_set(datasets: list[dict], *,
                   min_studies: int,
                   filter_threshold: Optional[float] = None,
                   filter_stat: Optional[str] = None,
                   symbols: Optional[set[str]] = None,
                   exclude: Optional[Iterable[str]] = None,
                   exclude_reason: Optional[str] = None,
                   ) -> tuple[Optional[set[str]], dict]:
    """Apply the selected gene-set rules to the studies' DE tables.

    Parameters
    ----------
    datasets : list of dict
        Each with a ``df`` holding a ``names`` column (as for pooling).
    min_studies : int
        Passed to the expression filter (its per-study count rule).
    filter_threshold : float, optional
        The meta-analysis expression filter
        (``kosmic.meta_analysis.independent_filter.apply_filter``); None
        leaves it off.
    filter_stat : str, optional
        The filter's statistic; defaults to ``kosmic.META_IF_FILTER_STAT``.
    symbols : set of str, optional
        Keep only these names, e.g. :func:`current_symbols`.
    exclude : iterable of str, optional
        Leave these genes out, e.g. :func:`inconsistent_merge_genes`.
    exclude_reason : str, optional
        What ``exclude`` is, for the record.

    Returns
    -------
    (set or None, dict)
        The gene set, or None when no rule is set (the datasets are then
        pooled as they are), and a record of each rule with how many genes
        it removed, in the order applied.
    """
    reported = set()
    for entry in datasets:
        df = entry.get("df")
        if df is not None and "names" in df.columns:
            reported.update(df["names"].astype(str))

    steps = []
    genes = set(reported)
    if filter_threshold is not None:
        from kosmic import META_IF_FILTER_STAT
        from kosmic.meta_analysis.independent_filter import apply_filter
        stat = filter_stat or META_IF_FILTER_STAT
        passing = apply_filter(datasets, threshold=filter_threshold, stat=stat,
                               min_studies=min_studies)
        steps.append(_step("expression filter", genes, genes & passing,
                           threshold=filter_threshold, stat=stat,
                           min_studies=min_studies))
        genes &= passing
    if symbols is not None:
        steps.append(_step("current symbols only", genes, genes & symbols))
        genes &= symbols
    if exclude is not None:
        exclude = set(exclude)
        steps.append(_step("excluded genes", genes, genes - exclude,
                           reason=exclude_reason, n_listed=len(exclude)))
        genes -= exclude

    record = {"n_genes_reported": len(reported), "rules": steps,
              "n_genes": len(genes) if steps else len(reported)}
    return (genes if steps else None), record


def _step(rule, before, after, **extra):
    return {"rule": rule, "n_before": len(before),
            "n_removed": len(before) - len(after), **extra}


def genes_tested(method_dfs: dict, method_keys: Iterable[str]) -> list[str]:
    """Genes every pooling method returned a result for: what the run tested.

    Rule methods (vote counting) report every gene they saw and are left out.
    """
    from kosmic.meta_analysis.pooling_dispatch import RULE_METHODS
    sets = [set(method_dfs[k]["names"].astype(str)) for k in method_keys
            if k in method_dfs and k not in RULE_METHODS]
    if not sets:
        sets = [set(method_dfs[k]["names"].astype(str)) for k in method_keys
                if k in method_dfs]
    return sorted(set.intersection(*sets)) if sets else []


def gene_set_path(results_csv) -> Path:
    """``consensus_reml_cc.csv`` -> ``consensus_reml_cc_genes.csv``."""
    p = Path(results_csv)
    return p.with_name(p.stem + GENE_SET_SUFFIX)


def save_gene_set(path, genes: Iterable[str]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.Series(sorted(genes), name="names").to_csv(path, index=False)
    return path


def read_gene_set(path) -> set[str]:
    return set(pd.read_csv(path, dtype=str, keep_default_na=False)["names"])


def newest_gene_set(folder) -> Optional[Path]:
    """The newest saved gene set in a meta-analysis selection folder."""
    folder = Path(folder)
    if not folder.is_dir():
        return None
    found = [p for p in folder.glob(f"consensus_*{GENE_SET_SUFFIX}") if p.is_file()]
    return max(found, key=lambda p: p.stat().st_mtime) if found else None


def meta_gene_sets(project_dir, counts_source, cell_types) -> dict[str, Optional[Path]]:
    """Cell type -> the newest gene set its meta-analysis saved, or None.

    Looks in ``meta_analysis/<counts>/<cell type>/``, the folder a
    meta-analysis of that cell type's per-study results writes to: the
    counts are those a DE run on ``counts_source`` writes under
    (``kosmic.paths.counts_folder_name``), and the cell type folder is the
    batch DE's name for it.
    """
    from kosmic.de.batch import slugify_cell_type
    from kosmic.meta_analysis.io import meta_selection_path
    from kosmic.paths import counts_folder_name, meta_output_dir
    counts = counts_folder_name(counts_source)
    return {ct: newest_gene_set(meta_output_dir(
                project_dir, meta_selection_path(counts, slugify_cell_type(ct))))
            for ct in cell_types}


__all__ = ["GENE_SET_SUFFIX", "build_gene_set", "meta_gene_sets", "current_symbols", "gene_set_path",
           "inconsistent_merge_genes", "newest_gene_set", "read_gene_set",
           "save_gene_set", "study_folders", "genes_tested"]
