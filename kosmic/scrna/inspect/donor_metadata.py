# Per-donor metadata from a table (a paper's supplementary table, a
# clinical sheet) joined onto obs by sample id.
#
# Deposited objects often carry only the donor id and the arm; age, sex
# and aetiology live in the paper. This module reads such a table, works
# out which of its columns identifies the study's samples, reports how
# many samples it matches, and writes the chosen columns into obs, one
# value per cell. Only the chosen columns: a supplementary table has a
# dozen clinical fields no design will use, and every column added is one
# more entry in every metadata dropdown.
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd

#: Column names, lower-cased, that carry a role the analysis uses. These
#: are pre-selected in the import dialog; everything else is opt-in.
SEX_LIKE = ('sex', 'gender', 'donor_sex')
AGE_LIKE = ('age', 'age_years', 'donor_age', 'age_at_death', 'age at death')


def read_table(path) -> pd.DataFrame:
    """Read a CSV, TSV or Excel table into a DataFrame of strings.

    Values are kept as text: '9' and 'NICM (chemo)' are both just labels
    until the user decides what a column is for, and a numeric age
    survives the round trip through obs either way.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in ('.xlsx', '.xls'):
        df = pd.read_excel(path, dtype=str)
    elif suffix == '.tsv' or suffix == '.txt':
        df = pd.read_csv(path, sep='\t', dtype=str)
    else:
        df = pd.read_csv(path, sep=None, engine='python', dtype=str)
    df.columns = [str(c).strip() for c in df.columns]
    return df.apply(lambda col: col.str.strip() if col.dtype == object else col)


def match_column(table: pd.DataFrame, sample_ids: Sequence[str]) -> Optional[str]:
    """The table column whose values best cover the study's sample ids.

    Returns None when no column matches even one id.
    """
    ids = set(map(str, sample_ids))
    best, best_n = None, 0
    for col in table.columns:
        n = len(ids & set(table[col].dropna().astype(str)))
        if n > best_n:
            best, best_n = col, n
    return best


@dataclass
class MatchReport:
    key_column: str
    n_samples: int
    n_matched: int
    unmatched_samples: list = field(default_factory=list)   # in obs, not in table
    unmatched_rows: list = field(default_factory=list)      # in table, not in obs
    duplicate_keys: list = field(default_factory=list)      # table rows sharing a key


def match_report(table: pd.DataFrame, key_column: str,
                 sample_ids: Sequence[str]) -> MatchReport:
    """How the table's key column lines up with the study's sample ids."""
    ids = [str(s) for s in pd.unique(pd.Series(list(sample_ids)))]
    keys = table[key_column].dropna().astype(str)
    key_set = set(keys)
    dup = sorted(keys[keys.duplicated()].unique())
    return MatchReport(
        key_column=key_column,
        n_samples=len(ids),
        n_matched=sum(1 for s in ids if s in key_set),
        unmatched_samples=[s for s in ids if s not in key_set],
        unmatched_rows=sorted(k for k in key_set if k not in set(ids)),
        duplicate_keys=dup,
    )


def suggested_columns(table: pd.DataFrame, key_column: str,
                      existing: Sequence[str] = ()) -> list[str]:
    """Columns to pre-select: sex-like and age-like, not already in obs."""
    existing_lower = {str(c).lower() for c in existing}
    out = []
    for col in table.columns:
        if col == key_column or col.lower() in existing_lower:
            continue
        low = col.lower()
        if low in SEX_LIKE or low in AGE_LIKE:
            out.append(col)
    return out


def column_summary(table: pd.DataFrame, col: str, n_examples: int = 3) -> tuple[int, str]:
    """(number of distinct values, a few example values) for the dialog."""
    values = table[col].dropna().astype(str)
    values = values[values != '']
    distinct = values.unique()
    return len(distinct), ', '.join(distinct[:n_examples])


def safe_obs_name(name: str, existing: Sequence[str] = ()) -> str:
    """An obs column name from a table header: lower-case, spaces and
    punctuation to underscores, unique against ``existing``."""
    import re
    base = re.sub(r'[^0-9a-zA-Z]+', '_', str(name)).strip('_').lower() or 'column'
    out, i = base, 2
    existing = set(map(str, existing))
    while out in existing:
        out = f"{base}_{i}"
        i += 1
    return out


def apply_metadata(adata, sample_col: str, table: pd.DataFrame, key_column: str,
                   columns: Sequence[str], rename: Optional[dict] = None) -> dict:
    """Write ``columns`` from ``table`` into ``adata.obs`` by sample id.

    Cells whose sample has no row in the table get an empty value, never
    an error: papers and deposits disagree on the odd id. A table with
    duplicate keys uses the first row for each. Returns
    ``{obs_name: table_column}`` for what was written.

    ``rename`` maps table column -> obs column name; unspecified columns
    get :func:`safe_obs_name`.
    """
    rename = dict(rename or {})
    keyed = table.dropna(subset=[key_column]).copy()
    keyed[key_column] = keyed[key_column].astype(str)
    keyed = keyed.drop_duplicates(subset=[key_column], keep='first').set_index(key_column)
    samples = adata.obs[sample_col].astype(str)
    written = {}
    for col in columns:
        obs_name = rename.get(col) or safe_obs_name(col, list(adata.obs.columns) + list(written))
        values = samples.map(keyed[col]).fillna('')
        adata.obs[obs_name] = pd.Categorical(values.astype(str))
        written[obs_name] = col
    return written
