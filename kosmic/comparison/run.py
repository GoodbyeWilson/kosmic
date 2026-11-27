"""Finding what can be compared, and running a comparison (ADR-017).

:func:`available_selections` lists the project's counts and cell types
with a per-study DESeq2 run in at least two studies, the shared atlas
excluded. :func:`check` resolves the merged design for one of them and
finds the atlas's saved merged run. :func:`run_comparison` collects every
method's result, builds the comparison and saves it.

What happens to merged pseudobulk is chosen by ``merged_action``:

- ``'auto'``: the saved run when it matches the design, a new fit when
  there is no saved run; a saved run that does not match is an error,
  because the user has to choose (the page asks first);
- ``'refit'``: fit again with the matched design (Resolve), replacing the
  saved run;
- ``'saved'``: use the saved run as it is (Run anyway); the differences
  are recorded in the comparison's parameters.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import pandas as pd

from kosmic.comparison.agreement import Rule, gene_characteristics
from kosmic.comparison.bundle import ComparisonBundle
from kosmic.comparison.design import Design, load_design
from kosmic.comparison.ma import MA_METHODS, load_ma
from kosmic.comparison.merged import (
    DEFAULT_FILTERS_NOTE, SavedMerged, find_saved_merged, run_merged,
)
from kosmic.comparison.standard import METHOD_LABELS, to_standard

ALL_METHODS = ('pb_merged', *MA_METHODS)
MERGED_ACTIONS = ('auto', 'refit', 'saved')


@dataclass(frozen=True)
class Selection:
    """One counts source and cell type with per-study runs to compare.

    ``runs`` maps each study to its gene DE table; ``names`` are the
    tables' dataset names ('Chaffin_Ventricular_Cardiomyocyte'), which the
    Meta-Analysis workspace records as ``study_names``.
    """

    counts: str
    cell_type: str
    path: str                 # '<counts>/<cell type folder>'
    runs: tuple               # ((study, Path), ...)
    names: tuple

    @property
    def de_paths(self) -> dict:
        return dict(self.runs)

    @property
    def studies(self) -> list:
        return [s for s, _ in self.runs]

    @property
    def label(self) -> str:
        return (f"{self.cell_type}, {self.counts.replace('_', ' ')} "
                f"({len(self.runs)} studies)")

    def restricted(self, studies) -> 'Selection':
        """The same selection with only ``studies``."""
        keep = [i for i, (s, _) in enumerate(self.runs) if s in set(studies)]
        return Selection(self.counts, self.cell_type, self.path,
                         tuple(self.runs[i] for i in keep), tuple(self.names[i] for i in keep))


def available_selections(project, method: str = 'deseq2') -> list:
    """The project's comparable selections, by cell type then counts."""
    from kosmic.de.batch import slugify_cell_type
    from kosmic.de.saved import list_gene_de_runs
    from kosmic.meta_analysis.io import meta_selection_folder, meta_selection_path
    from kosmic.paths import MASTER_ACCESSION, list_studies, study_dir

    groups: dict = {}
    for study in list_studies(project):
        if study == MASTER_ACCESSION:
            continue
        for run in list_gene_de_runs(study_dir(project, study)):
            if run.method != method or not run.cell_type or not run.counts:
                continue
            key = (run.counts, run.cell_type)
            groups.setdefault(key, {}).setdefault(study, run)
    out = []
    for (counts, cell_type), runs in groups.items():
        if len(runs) < 2:
            continue
        studies = sorted(runs)
        out.append(Selection(
            counts=counts, cell_type=cell_type,
            # The Meta-Analysis workspace names the folder by the tables'
            # cell type part, the slug ('Ventricular_Cardiomyocyte').
            path=meta_selection_path(counts, meta_selection_folder(
                [slugify_cell_type(cell_type)])),
            runs=tuple((s, runs[s].path) for s in studies),
            names=tuple(runs[s].name for s in studies)))
    out.sort(key=lambda s: (s.cell_type.lower(), s.counts))
    return out


@dataclass
class Check:
    """The merged design for a selection and the atlas's saved run."""

    design: Design
    saved: Optional[SavedMerged]

    @property
    def needs_choice(self) -> bool:
        """Whether the user has to choose between Resolve and Run anyway."""
        return self.saved is not None and not self.saved.matches


def check(project, selection: Selection) -> Check:
    """Resolve the merged design and find the saved merged run."""
    design = load_design(project, selection.de_paths)
    return Check(design, find_saved_merged(project, design))


def run_comparison(project, selection: Selection, methods, rule: Rule, *,
                   merged_action: str = 'auto', checked: Optional[Check] = None,
                   de_method: str = 'deseq2',
                   progress: Optional[Callable[[str], None]] = None) -> ComparisonBundle:
    """Collect the methods' results, compare them and save the comparison."""
    if merged_action not in MERGED_ACTIONS:
        raise ValueError(f"merged_action must be one of {MERGED_ACTIONS}")
    methods = [m for m in ALL_METHODS if m in set(methods)]
    if len(methods) < 2:
        raise ValueError("Choose at least two methods to compare.")
    say = progress or (lambda _m: None)
    checked = checked or check(project, selection)
    design = checked.design
    tables, sources = [], {}
    merged_info: dict = {}

    if 'pb_merged' in methods:
        saved = checked.saved
        if merged_action == 'refit' or saved is None:
            if saved is None:
                say("No saved merged pseudobulk run for this selection: fitting it.")
            merged = run_merged(project, design, progress=say)
            sources['pb_merged'] = {'source': 'fitted', 'path': str(merged.path),
                                    'note': ('Fitted with the matched design, replacing '
                                             'the saved run.' if saved is not None else
                                             'Fitted: there was no saved run.')}
        elif saved.matches:
            merged = saved
            sources['pb_merged'] = {'source': 'saved', 'path': str(saved.path),
                                    'note': 'Read from the saved run; its design matches.'}
        elif merged_action == 'saved':
            merged = saved
            sources['pb_merged'] = {'source': 'saved', 'path': str(saved.path),
                                    'note': 'Read from the saved run although its design '
                                            'differs (Run anyway).'}
            merged_info['design_differences'] = list(saved.differences)
        else:
            raise ValueError("The saved merged run's design differs from the studies'; "
                             "choose to refit it or to use it anyway.")
        merged_info['saved_design'] = (merged.record or {}).get('design')
        if design.filters is None:
            merged_info['filters_note'] = DEFAULT_FILTERS_NOTE
        tables.append(to_standard('pb_merged', merged.read(), rule.alpha))

    ma_ids = [m for m in methods if m in MA_METHODS]
    ma_settings = {}
    reml = None
    if ma_ids:
        results = load_ma(project, selection.path, list(selection.names), ma_ids,
                          de_method=de_method, counts=selection.counts,
                          studies=selection.studies, progress=say)
        for mid, res in results.items():
            sources[mid] = {'source': res.source, 'note': res.note,
                            'path': str(res.path) if res.path else None}
            ma_settings = res.settings
            tables.append(to_standard(mid, res.table, rule.alpha))
        if 'ma_reml' in results:
            reml = results['ma_reml'].table
    genes = _gene_characteristics(project, selection, de_method, reml, rule, say)

    params = {'counts': selection.counts, 'cell_type': selection.cell_type,
              'studies': selection.studies, 'dataset_names': list(selection.names),
              'design_formula': design.formula, 'design_mapping': design.mapping,
              'design_warnings': list(design.warnings), 'merged_action': merged_action,
              'sources': sources, 'meta_analysis_settings': ma_settings,
              'method_labels': {m: METHOD_LABELS[m] for m in methods}, **merged_info}
    bundle = ComparisonBundle.build(pd.concat(tables, ignore_index=True), rule, params,
                                    methods, genes=genes)
    folder = bundle.save(project, selection.path)
    say(f"Comparison saved in {folder}")
    return bundle


def _gene_characteristics(project, selection: Selection, de_method: str, reml, rule: Rule,
                          say) -> Optional[pd.DataFrame]:
    """Per gene, from the studies' own DE tables and REML
    (``agreement.gene_characteristics``). REML is pooled with the comparison's
    meta-analysis settings when it was not among the methods compared; None
    when the studies' tables cannot be read."""
    from kosmic.meta_analysis.saved import study_tables

    tables, missing = study_tables(project, list(selection.names), de_method=de_method,
                                   counts=selection.counts)
    if missing or not tables:
        say("Per-gene study summaries left out: no readable DE table for "
            + ', '.join(missing or ['any study']) + '.')
        return None
    if reml is None:
        reml = load_ma(project, selection.path, list(selection.names), ['ma_reml'],
                       de_method=de_method, counts=selection.counts,
                       studies=selection.studies)['ma_reml'].table
    return gene_characteristics(tables, reml, rule.alpha)
