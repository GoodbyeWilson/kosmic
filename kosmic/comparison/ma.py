"""The four meta-analysis methods of the comparison (ADR-017 §3).

REML, SumRank, Stouffer and vote counting are taken from a saved
single-method run of the Meta-Analysis workspace when one exists for the
selection's studies with the same settings, and are otherwise pooled
again from the studies' saved DE tables with the workspace's own
functions: the expression filter (``independent_filter.apply_filter``),
``pooling_dispatch.get_analytical_pool_fn`` and Benjamini-Hochberg
(``numerical.bh_fdr``), as ``GeneMAWorker`` runs them.

All four use one set of settings: those of the newest saved run of the
selection's studies (minimum studies, the expression filter, the gene-set
settings and, for vote counting, the minimum votes), else the
Meta-Analysis defaults. A saved run made with other settings is not read,
so the four methods are always pooled the same way. The gene set is built
as the Meta-Analysis workspace builds it
(``kosmic.meta_analysis.gene_set.build_gene_set``): the expression filter,
then current gene symbols only, then the genes harmonisation merged
inconsistently across the studies left out, each when the run used it.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import pandas as pd

#: Comparison method id -> pooling key (``pooling_dispatch``).
MA_METHODS = {'ma_reml': 'reml', 'ma_sumrank': 'sumrank', 'ma_stouffer': 'stouffer',
              'ma_vote': 'vote_count'}

_SINGLE_RUN = re.compile(r'^consensus_([a-z_]+?)_analytical\.csv$')


@dataclass(frozen=True)
class MASettings:
    """How the four methods are pooled."""

    min_studies: int
    filter_active: bool
    filter_threshold: Optional[float]
    filter_stat: Optional[str]
    min_votes: Optional[int] = None
    source: str = 'defaults'   # 'defaults' or the saved run's file name
    current_symbols_only: bool = False
    drop_inconsistent_merges: bool = False

    def as_dict(self) -> dict:
        return {'min_studies': self.min_studies,
                'independent_filter': {'active': self.filter_active,
                                       'threshold': self.filter_threshold,
                                       'stat': self.filter_stat},
                'current_symbols_only': self.current_symbols_only,
                'drop_inconsistent_merges': self.drop_inconsistent_merges,
                'min_votes': self.min_votes, 'settings_from': self.source}


@dataclass
class MAResult:
    """One method's pooled table and where it came from."""

    method: str
    key: str
    table: pd.DataFrame
    source: str          # 'saved' or 'pooled'
    note: str = ''
    path: Optional[Path] = None
    settings: dict = field(default_factory=dict)


def default_settings() -> MASettings:
    from kosmic import (
        META_IF_DEFAULT_PCT, META_IF_FILTER_STAT, META_IF_ON_BY_DEFAULT, MIN_STUDIES,
    )
    active = bool(META_IF_ON_BY_DEFAULT)
    return MASettings(min_studies=int(MIN_STUDIES), filter_active=active,
                      filter_threshold=(META_IF_DEFAULT_PCT / 100.0) if active else None,
                      filter_stat=META_IF_FILTER_STAT if active else None)


def _settings_from(params: dict, source: str) -> Optional[MASettings]:
    if 'min_studies' not in params:
        return None
    filt = params.get('independent_filter') or {}
    return MASettings(min_studies=int(params['min_studies']),
                      filter_active=bool(filt.get('active')),
                      filter_threshold=filt.get('threshold'),
                      filter_stat=filt.get('stat'),
                      min_votes=params.get('min_votes'), source=source,
                      current_symbols_only=bool(params.get('current_symbols_only')),
                      drop_inconsistent_merges=bool(params.get('drop_inconsistent_merges')))


def _same_pooling(a: MASettings, b: MASettings, key: str, n_studies: int) -> bool:
    """Whether two runs pooled ``key`` the same way. No minimum votes means
    every study (``vote_count_fast``), so it equals ``n_studies``."""
    same = (a.min_studies == b.min_studies and a.filter_active == b.filter_active
            and a.filter_threshold == b.filter_threshold and a.filter_stat == b.filter_stat
            and a.current_symbols_only == b.current_symbols_only
            and a.drop_inconsistent_merges == b.drop_inconsistent_merges)
    votes = (a.min_votes or n_studies) == (b.min_votes or n_studies)
    return same and (key != 'vote_count' or votes)


def saved_runs(project, selection: str, names) -> list:
    """``[(path, parameters)]`` of the selection's saved analytical runs that
    pooled exactly ``names``, newest first."""
    from kosmic.meta_analysis.saved import read_parameters
    from kosmic.paths import meta_output_dir

    folder = meta_output_dir(project, selection)
    if not folder.is_dir():
        return []
    want = sorted(map(str, names))
    out = []
    for path in sorted(folder.glob('consensus_*_analytical.csv'),
                       key=lambda p: p.stat().st_mtime, reverse=True):
        params = read_parameters(path)
        if sorted(map(str, params.get('study_names') or [])) == want:
            out.append((path, params))
    return out


def choose_settings(runs: list) -> MASettings:
    """The settings of the newest saved run, else the defaults; the minimum
    votes from the newest run that included vote counting."""
    from dataclasses import replace
    settings = next((s for s in (_settings_from(p, path.name) for path, p in runs)
                     if s is not None), None) or default_settings()
    votes = next((p.get('min_votes') for _path, p in runs
                  if 'vote_count' in (p.get('method_keys') or [])), None)
    return replace(settings, min_votes=votes)


def _read_saved(path: Path, key: str) -> pd.DataFrame:
    """A saved single-method run's table, with 'fdr' that method's BH FDR."""
    from kosmic.numerical import bh_fdr
    df = pd.read_csv(path)
    if f'fdr_{key}' in df.columns:
        df['fdr'] = df[f'fdr_{key}']
    elif 'pvals_pooled' in df.columns:
        df['fdr'] = bh_fdr(df['pvals_pooled'].values)
    return df


def pool(tables: list, key: str, settings: MASettings,
         exclude: Optional[set] = None) -> pd.DataFrame:
    """Pool the studies' DE tables with one method, as the Meta-Analysis
    workspace does: the gene set, pooling function, BH FDR.

    ``tables`` is ``[(name, DataFrame)]`` (``saved.study_tables``);
    ``exclude`` the genes merged inconsistently across the studies, used
    when ``settings.drop_inconsistent_merges``.
    """
    from kosmic.meta_analysis.gene_set import build_gene_set, current_symbols
    from kosmic.meta_analysis.independent_filter import filter_datasets_to_universe
    from kosmic.meta_analysis.pooling_dispatch import RULE_METHODS, get_analytical_pool_fn
    from kosmic.numerical import bh_fdr
    from kosmic.scrna.inspect.detection import detect_species

    datasets = [{'name': n, 'df': df} for n, df in tables]
    symbols = None
    if settings.current_symbols_only:
        names = [g for _n, df in tables for g in df['names'].astype(str)[:200]]
        symbols = current_symbols(detect_species(names))
    filtered = settings.filter_active and settings.filter_threshold is not None
    genes, _record = build_gene_set(
        datasets, min_studies=settings.min_studies,
        filter_threshold=settings.filter_threshold if filtered else None,
        filter_stat=settings.filter_stat or 'mean_max_pct',
        symbols=symbols,
        exclude=exclude if settings.drop_inconsistent_merges else None,
        exclude_reason='merged inconsistently across studies')
    if genes is not None:
        datasets = filter_datasets_to_universe(datasets, genes)
    min_votes = settings.min_votes if key == 'vote_count' else None
    fn = get_analytical_pool_fn(key, settings.min_studies, min_votes=min_votes)
    df = fn([d['df'] for d in datasets], min_studies=settings.min_studies)
    if key not in RULE_METHODS and 'pvals_pooled' in df.columns:
        df['fdr'] = bh_fdr(df['pvals_pooled'].values)
    return df


def load_ma(project, selection: str, names, method_ids, *, de_method: Optional[str] = None,
            counts: Optional[str] = None, studies=None,
            progress: Optional[Callable[[str], None]] = None) -> dict:
    """``{method id: MAResult}`` for the selected meta-analysis methods.

    ``names`` are the pooled datasets' names ('Chaffin_Ventricular_Cardiomyocyte'),
    as the Meta-Analysis workspace records them in ``study_names``;
    ``studies`` their study folders, needed when the settings leave out
    the genes merged inconsistently across studies.
    Raises ValueError when a study's DE table cannot be read for a method
    that has to be pooled, or when the inconsistent merges cannot be
    checked because ``studies`` is not given.
    """
    from kosmic.meta_analysis.saved import study_tables

    say = progress or (lambda _m: None)
    runs = saved_runs(project, selection, names)
    settings = choose_settings(runs)
    single = {}
    for path, params in runs:
        m = _SINGLE_RUN.match(path.name)
        if m and m.group(1) in MA_METHODS.values() and m.group(1) not in single:
            single[m.group(1)] = (path, params)

    tables = None
    exclude = None
    merge_note = ''
    out = {}
    for mid in method_ids:
        key = MA_METHODS[mid]
        saved = single.get(key)
        if saved is not None:
            saved_settings = _settings_from(saved[1], saved[0].name)
            if saved_settings is not None and _same_pooling(saved_settings, settings, key,
                                                            len(names)):
                say(f"{key}: read from {saved[0].name}")
                out[mid] = MAResult(mid, key, _read_saved(saved[0], key), 'saved',
                                    f"Read from the saved run {saved[0].name}.", saved[0],
                                    settings.as_dict())
                continue
            note = (f"Pooled again: the saved run {saved[0].name} used other settings "
                    f"than {settings.source}.")
        else:
            note = "Pooled again: no saved single-method run for these studies."
        if tables is None:
            tables, missing = study_tables(project, list(names), de_method=de_method,
                                           counts=counts)
            if missing:
                raise ValueError("No readable DE table for: " + ', '.join(missing))
            if settings.drop_inconsistent_merges:
                if not studies:
                    raise ValueError("The study folders are needed to leave out genes "
                                     "merged inconsistently across studies.")
                from kosmic.meta_analysis.gene_set import inconsistent_merge_genes
                exclude, unknown = inconsistent_merge_genes(project, list(studies))
                if unknown:
                    merge_note = (" Inconsistent merges could not be checked for (no "
                                  "harmonisation record): " + ', '.join(unknown) + '.')
        say(f"{key}: pooling {len(tables)} studies...")
        out[mid] = MAResult(mid, key, pool(tables, key, settings, exclude), 'pooled',
                            note + merge_note, None, settings.as_dict())
    return out
