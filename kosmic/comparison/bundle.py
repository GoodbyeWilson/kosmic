"""A comparison's results, saved in the project (ADR-017 §6).

A :class:`ComparisonBundle` holds the standard table of every method
(:mod:`kosmic.comparison.standard`) with each gene's ``in_universe`` and
``significant`` under the rule used, the call sets, the pairwise measures
and the parameters. It is saved in
``methods_comparison/<counts>/<selection>/`` (``paths.comparison_output_dir``):

- ``results.csv``: the long table, every gene of every method, so the
  comparison can be cut under another rule without refitting;
- ``pairwise.csv``: :func:`kosmic.comparison.agreement.pairwise`;
- ``params.json``: the methods, the rule, the merged design and its
  warnings, where each result came from, and what the comparison was made
  of (selection, studies, KOSMIC version).

Saving also records a ``methods_comparison`` stage in the folder's
``provenance.json``, replacing the previous one, because a new run
overwrites the saved comparison.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd

from kosmic.comparison.agreement import Rule, apply_rule, call_sets, methods_in, pairwise

RESULTS_FILE = 'results.csv'
PAIRWISE_FILE = 'pairwise.csv'
PARAMS_FILE = 'params.json'
PROVENANCE_STAGE = 'methods_comparison'


def _kosmic_version() -> Optional[str]:
    from importlib.metadata import PackageNotFoundError, version
    try:
        return version('kosmic')
    except PackageNotFoundError:
        return None


@dataclass
class ComparisonBundle:
    """What the Methods Comparison workspace hands to Figures."""

    results: pd.DataFrame
    rule: Rule
    params: dict = field(default_factory=dict)
    folder: Optional[Path] = None

    @classmethod
    def build(cls, standard: pd.DataFrame, rule: Rule, params: Optional[dict] = None,
              methods=None) -> 'ComparisonBundle':
        """A bundle from the concatenated standard tables of the methods."""
        return cls(apply_rule(standard, rule, methods), rule, dict(params or {}))

    def recut(self, rule: Rule, methods=None) -> 'ComparisonBundle':
        """The same results under another rule or set of methods."""
        base = self.results.drop(columns=['in_universe', 'significant'], errors='ignore')
        return ComparisonBundle(apply_rule(base, rule, methods), rule, dict(self.params),
                                self.folder)

    @property
    def methods(self) -> list:
        return methods_in(self.results)

    @property
    def sets(self) -> dict:
        return call_sets(self.results)

    @property
    def pairwise(self) -> pd.DataFrame:
        return pairwise(self.results)

    @property
    def n_universe(self) -> int:
        return int(self.results.loc[self.results['in_universe'], 'gene'].nunique())

    def full_params(self) -> dict:
        return {**self.params, 'methods': self.methods, 'rule': self.rule.as_dict(),
                'rule_text': self.rule.text, 'n_genes_compared': self.n_universe,
                'n_called': {m: len(s) for m, s in self.sets.items()},
                'saved': datetime.now().isoformat(timespec='seconds'),
                'kosmic_version': _kosmic_version()}

    def save(self, project, selection: str) -> Path:
        """Write the bundle to the selection's folder and record the stage."""
        from kosmic import provenance
        from kosmic.paths import comparison_output_dir

        folder = comparison_output_dir(project, selection)
        folder.mkdir(parents=True, exist_ok=True)
        self.results.to_csv(folder / RESULTS_FILE, index=False)
        self.pairwise.to_csv(folder / PAIRWISE_FILE, index=False)
        params = {**self.full_params(), 'selection': selection}
        (folder / PARAMS_FILE).write_text(json.dumps(params, indent=2, default=str),
                                          encoding='utf-8')
        try:
            provenance.record_stage(folder, f"{Path(project).name} ({selection})",
                                    PROVENANCE_STAGE, params, replace=True)
        except Exception:   # noqa: BLE001 - a methods note must never lose the result
            pass
        self.folder = folder
        self.params = params
        return folder


def load(folder) -> Optional[ComparisonBundle]:
    """A saved bundle, or None when the folder does not hold a readable one."""
    folder = Path(folder)
    try:
        results = pd.read_csv(folder / RESULTS_FILE)
        params = json.loads((folder / PARAMS_FILE).read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return None
    for col in ('tested', 'called', 'in_universe', 'significant'):
        if col in results.columns:
            results[col] = results[col].astype(str).str.lower().isin(['true', '1'])
    results['rank'] = pd.to_numeric(results['rank'], errors='coerce').astype('Int64')
    rule = Rule(**{k: v for k, v in (params.get('rule') or {}).items()
                   if k in ('rule', 'alpha', 'top_n')})
    return ComparisonBundle(results, rule, params, folder)


def list_saved(project) -> list:
    """``[(selection, folder)]`` of the project's saved comparisons, newest first."""
    from kosmic.paths import comparison_output_dir

    root = comparison_output_dir(project)
    if not root.is_dir():
        return []
    found = [(p.parent, p.stat().st_mtime) for p in root.glob(f'*/*/{PARAMS_FILE}')]
    found.sort(key=lambda x: x[1], reverse=True)
    return [(f"{f.parent.name}/{f.name}", f) for f, _t in found]
