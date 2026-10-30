"""A project's saved meta-analysis runs, read back for Figure Export (ADR-014).

Every pooling run is saved in its selection's folder,
``meta_analysis/{selection}/``:

- a gene-level run as ``consensus_{methods}_{calibration}.csv`` with a
  ``.json`` sidecar whose parameters name the DE method and the studies
  pooled (``study_names``, each a study's DE accession such as
  ``Chaffin_Endothelial_Cell``);
- a pathway-level run as ``pathway_ma_{method}.csv`` and ``.json``;
- the Enrichment step's runs on the pooled gene list:
  ``enrichment_results.csv`` (over-representation) and one
  ``gsea_library_results_{library}.csv`` per library, with its
  running-enrichment curves beside it (:mod:`kosmic.de.saved`).

:func:`list_meta_runs` lists them and :func:`figure_inputs` builds, from
one run, the mapping the Meta-Analysis workspace's ``figure_inputs()``
gives Figure Export, so the Meta figures draw a saved run the same way
as one held in memory. Each study's DE table is loaded as the workspace
loads it (:func:`kosmic.meta_analysis.io.load_study_de`). The mapping
reads each part when it is first asked for, and :meth:`MetaFigureInputs.available`
says whether a part exists without reading a table.
"""
from __future__ import annotations

import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

import pandas as pd

from kosmic.paths import meta_output_dir

_GENE_RUN = re.compile(r'^consensus_(.+)_(analytical|cc)\.csv$')
_PATHWAY_RUN = re.compile(r'^pathway_ma_(.+)\.csv$')
_GSEA_LIBRARY = re.compile(r'^gsea_library_results_(.+)\.csv$')
_CALIBRATION_TEXT = {'analytical': 'analytical', 'cc': 'calibrated'}


@dataclass(frozen=True)
class MetaRun:
    """One saved meta-analysis run: a selection's gene-level run and the
    newest pathway-level run beside it (either may be absent)."""

    project: Path
    selection: str
    gene_path: Optional[Path]
    pathway_path: Optional[Path]

    @property
    def folder(self) -> Path:
        return meta_output_dir(self.project, self.selection)

    @property
    def label(self) -> str:
        """'Endothelial Cell: reml + vote count (analytical)'."""
        sel = self.selection.replace('_', ' ')
        if self.gene_path is None:
            m = _PATHWAY_RUN.match(self.pathway_path.name)
            return f"{sel}: pathways, {m.group(1).replace('_', ' ')}"
        m = _GENE_RUN.match(self.gene_path.name)
        methods = m.group(1).replace('vote_count', 'vote count').replace('_', ' + ')
        return f"{sel}: {methods} ({_CALIBRATION_TEXT[m.group(2)]})"


def _newest(paths) -> Optional[Path]:
    paths = list(paths)
    return max(paths, key=lambda p: p.stat().st_mtime) if paths else None


def list_meta_runs(project_folder) -> list[MetaRun]:
    """The project's saved meta-analysis runs, by selection and then newest
    first. A selection with pathway runs only gives one run."""
    root = meta_output_dir(project_folder)
    if not root.is_dir():
        return []
    runs = []
    for folder in sorted(p for p in root.iterdir() if p.is_dir()):
        files = [p for p in folder.iterdir() if p.is_file()]
        genes = sorted((p for p in files if _GENE_RUN.match(p.name)),
                       key=lambda p: p.stat().st_mtime, reverse=True)
        pathway = _newest(p for p in files if _PATHWAY_RUN.match(p.name))
        for gene in genes or [None]:
            if gene is None and pathway is None:
                continue
            runs.append(MetaRun(Path(project_folder), folder.name, gene, pathway))
    return runs


def read_parameters(csv_path) -> dict:
    """The parameters in a run's ``.json`` sidecar ({} without one)."""
    try:
        data = json.loads(Path(csv_path).with_suffix('.json').read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return {}
    return data.get('parameters') or {}


def study_tables(project_folder, names, *, level: str = 'gene',
                 de_method: Optional[str] = None) -> tuple[list, list]:
    """The DE tables of the studies a run pooled.

    Parameters
    ----------
    project_folder : path
    names : list of str
        The run's ``study_names``.
    level : {'gene', 'pathway'}
    de_method : str, optional
        The gene-level run's DE method; a study without a table of that
        method uses another, as the Meta-Analysis workspace does.

    Returns
    -------
    (tables, missing)
        ``tables`` is ``[(name, DataFrame)]`` in the run's order;
        ``missing`` the names with no readable table in the project (an
        external CSV imported for the run, or a study since removed).
    """
    loaders, missing = _study_loaders(project_folder, names, level, de_method)
    tables = []
    for name in names:
        if name in missing:
            continue
        try:
            tables.append((name, loaders[name]()))
        except (OSError, ValueError):
            missing.append(name)
    return tables, missing


def _study_loaders(project_folder, names, level, de_method):
    """({name: function loading its table}, names without a table), from
    the project's file listing only."""
    from kosmic.meta_analysis.io import (
        discover_de_results, load_de_results, load_study_de, require_se,
        select_de_entries,
    )
    discovered = discover_de_results(Path(project_folder))
    found: dict = {}
    if level == 'gene':
        chosen, _fallback, _skipped = select_de_entries(
            discovered, [de_method] if de_method else [], allowed=names)
        for entry in chosen:
            found[entry['accession']] = (
                lambda e=entry: load_study_de(e['path'], e.get('de_method')))
    else:
        by_name: dict = {}
        for entry in discovered:
            if entry['file_type'] == 'pathway_de':
                by_name.setdefault(entry['accession'], []).append(entry)
        for name, entries in by_name.items():
            # The pathway level, not pathway groups, as the workspace pools.
            base = [e for e in entries if 'parent_' not in e['path'].name]
            entry = (base or entries)[0]

            def _load(e=entry, n=name):
                df = load_de_results(str(e['path']), dataset_name=n)
                require_se(df, str(e['path']))
                return df
            found[name] = _load
    return found, [n for n in names if n not in found]


class MetaFigureInputs(Mapping):
    """The Meta-Analysis workspace's ``figure_inputs()`` keys for a saved
    run. Each value is read the first time it is asked for."""

    def __init__(self, run: MetaRun):
        self.run = run
        self._values: dict = {
            'selection': run.selection, 'project_folder': str(run.project)}
        self._loaders: dict[str, Callable] = {
            'gene_results': lambda: self._csv(run.gene_path),
            'pathway_results': lambda: self._csv(run.pathway_path),
            'gene_studies': lambda: self._studies('gene')[0],
            'pathway_studies': lambda: self._studies('pathway')[0],
            'missing_gene_studies': lambda: self._run_studies('gene')[1],
            'missing_pathway_studies': lambda: self._run_studies('pathway')[1],
            'enrichment_results': self._enrichment,
            'gsea_library_results': lambda: self._csv(self._gsea_path()),
            'gsea_library_details': lambda: self._gsea_details()[0],
            'gsea_library_name': lambda: self._gsea_details()[1],
            'enrichment_ranking': lambda: self['gene_results'],
        }
        self._study_cache: dict = {}
        self._loader_cache: dict = {}

    # -- Mapping ------------------------------------------------------------------

    def __getitem__(self, key):
        if key not in self._values:
            if key not in self._loaders:
                raise KeyError(key)
            self._values[key] = self._loaders[key]()
        return self._values[key]

    def __iter__(self):
        return iter([*self._values.keys() - self._loaders.keys(), *self._loaders])

    def __len__(self):
        return len(self._values.keys() | self._loaders.keys())

    # -- Cheap availability ----------------------------------------------------

    def available(self, key: str) -> bool:
        """Whether 'key' has a value, from which files exist."""
        run = self.run
        if key in ('gene_results', 'gene_studies', 'enrichment_ranking'):
            return run.gene_path is not None and run.gene_path.is_file()
        if key in ('pathway_results', 'pathway_studies'):
            return run.pathway_path is not None and run.pathway_path.is_file()
        if key == 'enrichment_results':
            from kosmic.de.saved import enrichment_path
            return enrichment_path(run.folder).is_file()
        if key == 'gsea_library_results':
            return self._gsea_path() is not None
        if key == 'gsea_library_details':
            from kosmic.de.saved import gsea_details_path
            path = self._gsea_path()
            return path is not None and gsea_details_path(path).is_file()
        return self.get(key) is not None

    # -- Readers -----------------------------------------------------------------

    @staticmethod
    def _csv(path) -> Optional[pd.DataFrame]:
        if path is None:
            return None
        try:
            return pd.read_csv(path)
        except (OSError, ValueError):
            return None

    def _run_studies(self, level: str):
        """(names, {name: loader}, missing names) of the run's studies at
        'level', from the run's sidecar and the project's file listing."""
        if level not in self._loader_cache:
            path = self.run.gene_path if level == 'gene' else self.run.pathway_path
            params = read_parameters(path) if path is not None else {}
            names = [str(n) for n in params.get('study_names') or []]
            loaders, missing = (_study_loaders(self.run.project, names, level,
                                               params.get('de_method'))
                                if names else ({}, []))
            self._loader_cache[level] = (names, loaders, missing)
        names, loaders, missing = self._loader_cache[level]
        return loaders, missing, names

    def _studies(self, level: str):
        """([(name, table)], missing names) of the run's studies."""
        if level not in self._study_cache:
            loaders, missing, names = self._run_studies(level)
            missing = list(missing)
            tables = []
            for name in names:
                if name in missing:
                    continue
                try:
                    tables.append((name, loaders[name]()))
                except (OSError, ValueError):
                    missing.append(name)
            self._study_cache[level] = (tables, missing)
        return self._study_cache[level]

    def _enrichment(self):
        from kosmic.de.saved import read_enrichment
        return read_enrichment(self.run.folder)

    def _gsea_path(self) -> Optional[Path]:
        folder = self.run.folder
        if not folder.is_dir():
            return None
        return _newest(p for p in folder.iterdir() if _GSEA_LIBRARY.match(p.name))

    def _gsea_details(self):
        path = self._gsea_path()
        if path is None:
            return None, None
        from kosmic.de.saved import read_gsea_details
        details, library = read_gsea_details(path)
        return details, library or _GSEA_LIBRARY.match(path.name).group(1)


def figure_inputs(run: MetaRun) -> MetaFigureInputs:
    """The Meta figures' inputs for a saved run (see the module header)."""
    return MetaFigureInputs(run)
