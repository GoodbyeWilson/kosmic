# Sex inference per sample from XIST and Y-chromosome gene expression.
#
# XIST is expressed in every female cell and in no male cell; the Y-linked
# genes (DDX3Y, UTY, KDM5D, RPS4Y1, EIF1AY, USP9Y) the reverse. Summed over
# a donor's cells the two signals separate completely, so the call is a
# genotype read from expression rather than a statistic. It stands on its
# own where the depositor recorded no sex, and where they did it checks
# the metadata: a sample whose call disagrees with its record is mislabelled
# or swapped, and a sample carrying both signals is a mixture (ambient
# contamination, or a multiplexed library that was not demultiplexed).
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

#: obs columns a depositor might have used for sex, best first.
SEX_COLUMN_CANDIDATES = (
    'sex', 'Sex', 'SEX', 'gender', 'Gender', 'donor_sex',
)

XIST = {'human': 'XIST', 'mouse': 'Xist'}
Y_GENES = {
    'human': ('DDX3Y', 'UTY', 'KDM5D', 'RPS4Y1', 'EIF1AY', 'USP9Y'),
    'mouse': ('Ddx3y', 'Uty', 'Kdm5d', 'Eif2s3y'),
}

#: Below this many counts per million on both signals a sample gets no
#: call: too few cells, or a matrix the genes are missing from.
MIN_SIGNAL_CPM = 20.0
#: A minor signal above this fraction of the major one marks a mixture.
#: Clean donors sit at 0.001-0.01; TWCM-11-104 (Koenig) is at 0.29.
MIXED_FRACTION = 0.10

FEMALE, MALE, UNCLEAR = 'female', 'male', 'unclear'


def normalise_sex_value(value) -> Optional[str]:
    """Map a recorded sex label to 'female' / 'male', or None if neither."""
    if value is None:
        return None
    v = str(value).strip().lower()
    if v in ('f', 'female', 'woman', 'women', 'xx'):
        return FEMALE
    if v in ('m', 'male', 'man', 'men', 'xy'):
        return MALE
    return None


def recorded_sex_column(obs) -> Optional[str]:
    """The first sex-like column present in ``obs``, or None."""
    for col in SEX_COLUMN_CANDIDATES:
        if col in obs.columns:
            return col
    return None


def sex_genes(var_names, species: Optional[str] = None):
    """The XIST and Y-gene names present in ``var_names`` for the species
    (detected from the names when not given). Returns (xist, y_genes);
    ``xist`` is None when absent."""
    names = set(map(str, var_names))
    if species is None:
        from kosmic.scrna.inspect.detection import detect_species
        species = detect_species(list(var_names))
    species = 'mouse' if species == 'mouse' else 'human'
    xist = XIST[species] if XIST[species] in names else None
    y = [g for g in Y_GENES[species] if g in names]
    return xist, y


def infer_sex(adata, sample_col: str, species: Optional[str] = None,
              counts_layer: Optional[str] = None) -> pd.DataFrame:
    """Call the sex of every sample from its summed XIST and Y-gene counts.

    Parameters
    ----------
    adata : anndata.AnnData
        Counts in ``X`` (or ``layers[counts_layer]``). Normalised values
        also work, since only the ratio of the two signals is used.
    sample_col : str
        obs column identifying the donor or sample.
    species : 'human' | 'mouse', optional
        Detected from the gene names when not given.
    counts_layer : str, optional

    Returns
    -------
    pd.DataFrame
        Indexed by sample, columns ``n_cells``, ``xist_cpm``, ``y_cpm``
        (counts per million of the sample's total), ``sex_inferred``
        (``'female'`` / ``'male'`` / ``'unclear'``) and ``sex_flag``
        (``''``, ``'mixed signal'`` when the minor signal exceeds
        :data:`MIXED_FRACTION` of the major, or ``'no signal'``).
        Empty when neither XIST nor any Y gene is in the matrix.
    """
    import scipy.sparse as sp

    xist, y = sex_genes(adata.var_names, species)
    if xist is None and not y:
        return pd.DataFrame(columns=['n_cells', 'xist_cpm', 'y_cpm',
                                     'sex_inferred', 'sex_flag'])
    X = adata.layers[counts_layer] if counts_layer else adata.X
    names = list(map(str, adata.var_names))
    idx = [names.index(g) for g in ([xist] if xist else []) + y]
    sub = X[:, idx]
    sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
    total = np.asarray(X.sum(axis=1)).ravel().astype(float)

    per_cell = pd.DataFrame({
        'sample': adata.obs[sample_col].astype(str).values,
        'xist': sub[:, 0] if xist else 0.0,
        'y': sub[:, 1 if xist else 0:].sum(axis=1) if y else 0.0,
        'total': total,
    })
    g = per_cell.groupby('sample', observed=True).sum(numeric_only=True)
    n_cells = per_cell.groupby('sample', observed=True).size()
    with np.errstate(divide='ignore', invalid='ignore'):
        xist_cpm = np.where(g['total'] > 0, g['xist'] / g['total'] * 1e6, 0.0)
        y_cpm = np.where(g['total'] > 0, g['y'] / g['total'] * 1e6, 0.0)

    major = np.maximum(xist_cpm, y_cpm)
    minor = np.minimum(xist_cpm, y_cpm)
    call = np.where(xist_cpm > y_cpm, FEMALE, MALE).astype(object)
    flag = np.full(len(g), '', dtype=object)
    no_signal = major < MIN_SIGNAL_CPM
    call[no_signal] = UNCLEAR
    flag[no_signal] = 'no signal'
    with np.errstate(divide='ignore', invalid='ignore'):
        mixed = (~no_signal) & (minor / np.where(major > 0, major, np.nan) > MIXED_FRACTION)
    flag[mixed] = 'mixed signal'

    return pd.DataFrame({
        'n_cells': n_cells.reindex(g.index).astype(int).values,
        'xist_cpm': np.round(xist_cpm, 1),
        'y_cpm': np.round(y_cpm, 1),
        'sex_inferred': call,
        'sex_flag': flag,
    }, index=g.index)


def compare_with_recorded(table: pd.DataFrame, obs, sample_col: str,
                          sex_col: Optional[str] = None) -> pd.DataFrame:
    """Add ``sex_recorded`` and ``sex_check`` to an :func:`infer_sex` table.

    ``sex_col`` defaults to the first of :data:`SEX_COLUMN_CANDIDATES` in
    ``obs``. ``sex_check`` is ``'mismatch'`` where the recorded value is
    female/male and differs from a confident call, ``'not recorded'``
    where there is no recorded value, and ``''`` otherwise. A sample
    whose recorded sex varies between its cells is reported as
    ``'inconsistent record'``.
    """
    out = table.copy()
    sex_col = sex_col or recorded_sex_column(obs)
    if sex_col is None or sex_col not in obs.columns:
        out['sex_recorded'] = None
        out['sex_check'] = 'not recorded'
        return out
    per = (pd.DataFrame({'sample': obs[sample_col].astype(str).values,
                         'sex': obs[sex_col].map(normalise_sex_value).values})
           .groupby('sample', observed=True)['sex']
           .agg(lambda s: sorted(set(x for x in s if x is not None))))
    recorded, check = [], []
    for sample, call in zip(out.index, out['sex_inferred']):
        values = per.get(sample, [])
        if len(values) > 1:
            recorded.append('/'.join(values))
            check.append('inconsistent record')
        elif not values:
            recorded.append(None)
            check.append('not recorded')
        else:
            recorded.append(values[0])
            check.append('mismatch' if call in (FEMALE, MALE) and call != values[0] else '')
    out['sex_recorded'] = recorded
    out['sex_check'] = check
    return out


def write_inferred_sex(adata, sample_col: str, table: pd.DataFrame,
                       column: str = 'sex_inferred') -> None:
    """Write each sample's call into ``adata.obs[column]`` (one value per cell)."""
    calls = table['sex_inferred'].to_dict()
    adata.obs[column] = pd.Categorical(
        adata.obs[sample_col].astype(str).map(calls).fillna(UNCLEAR))
