# Comparing clusters with the labels the authors deposited.
#
# Many public objects carry the authors' own cell-type calls. The Annotate
# tab keeps them in obs['cell_type_author'] before writing its own labels,
# and this module summarises how the clusters line up with them: per
# cluster, which author label dominates and what share of the cluster it
# covers; overall, what share of author-labelled cells carry their
# cluster's dominant label (cluster purity against the authors). No mapping
# between the two vocabularies is assumed, so the figure says how cleanly
# the clustering reproduces the authors' partition, not whether the names
# agree.
#
# Cells the authors left unlabelled (NaN, empty, 'unknown', 'NA' and the
# like) are counted separately and never as disagreement.
#
# A combined dataset holds the labels of several studies, each in its own
# authors' vocabulary ('CM' in one, 'Cardiomyocyte' in another), so no
# single name can dominate a cluster even when every study agrees with it.
# Given each cell's study, agreement is counted within each study: a cell
# agrees when it carries the most common label of its own study in its
# cluster.
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

UNLABELLED = frozenset({'', 'nan', 'none', 'na', 'n/a', 'unknown', 'unassigned',
                        'unlabelled', 'unlabeled', 'undefined'})


def author_label_series(values) -> np.ndarray:
    """Author labels as a string array with unlabelled cells set to ''.

    Accepts any obs column (categorical, object, numeric); labels are
    compared case-insensitively against UNLABELLED.
    """
    arr = pd.Series(values).astype(object)
    out = np.where(arr.isna(), '', arr.astype(str).values).astype(object)
    low = np.char.lower(np.char.strip(out.astype(str)))
    out[np.isin(low, list(UNLABELLED))] = ''
    return out


@dataclass
class AuthorBreakdown:
    """Author labels within one cluster.

    'n_agree' is the number of labelled cells carrying the most common
    label, counted within each study when 'by_study' is set.
    """
    counts: pd.Series          # label -> cells, descending, unlabelled removed
    n_unlabelled: int
    n_agree: int = 0
    by_study: Optional[dict] = None   # study -> (top label, its cells, labelled cells)

    @property
    def n_labelled(self) -> int:
        return int(self.counts.sum())

    @property
    def top(self) -> tuple[str, int]:
        if self.counts.empty:
            return '', 0
        return str(self.counts.index[0]), int(self.counts.iloc[0])

    @property
    def agreement(self) -> float:
        """Share of the labelled cells that carry the most common label."""
        return self.n_agree / self.n_labelled if self.n_labelled else float('nan')

    def describe(self, max_rows: int = 6) -> str:
        """Multi-line text: the labels and their shares, per study when grouped."""
        if self.by_study:
            lines = [f"{study}: {lb} {100 * k / n:.0f}% of {n:,} labelled"
                     for study, (lb, k, n) in self.by_study.items()]
        else:
            n = self.n_labelled
            lines = [f"{lb}: {c:,} ({100 * c / n:.1f}%)"
                     for lb, c in self.counts.head(max_rows).items()]
            rest = self.counts.iloc[max_rows:]
            if len(rest):
                lines.append(f"{len(rest)} other labels: {int(rest.sum()):,}")
        if self.n_unlabelled:
            lines.append(f"no author label: {self.n_unlabelled:,}")
        return "\n".join(lines)


def author_breakdown(labels: np.ndarray, studies=None) -> AuthorBreakdown:
    """Count author labels in one cluster (labels from author_label_series).

    'studies', each cell's study, makes agreement count within each study;
    without it, or with a single study, the cluster is one group.
    """
    arr = np.asarray(labels, dtype=object)
    labelled = arr != ''
    counts = pd.Series(arr[labelled]).value_counts()
    b = AuthorBreakdown(counts=counts, n_unlabelled=int((~labelled).sum()),
                        n_agree=int(counts.iloc[0]) if len(counts) else 0)
    if studies is not None:
        st = np.asarray(studies, dtype=object)[labelled].astype(str)
        if len(set(st)) > 1:
            per = pd.Series(arr[labelled]).groupby(st).value_counts()
            b.by_study = {}
            for study, sub in per.groupby(level=0):
                sub = sub.droplevel(0).sort_values(ascending=False)
                b.by_study[study] = (str(sub.index[0]), int(sub.iloc[0]), int(sub.sum()))
            b.n_agree = sum(k for _, k, _ in b.by_study.values())
    return b


def study_groups(obs: pd.DataFrame, col: str = 'study'):
    """Each cell's study when 'obs' combines more than one, else None."""
    if col in obs.columns and obs[col].nunique() > 1:
        return obs[col].astype(str).to_numpy()
    return None


def cluster_purity(clusters, author_values, studies=None) -> dict:
    """Overall agreement between a clustering and the authors' labels.

    Returns n_labelled, n_agree (cells carrying their cluster's dominant
    author label, within their own study when 'studies' is given),
    n_unlabelled and purity (n_agree / n_labelled, NaN when nothing is
    labelled). Cells with no cluster (NaN) are ignored.
    """
    cl = pd.Series(clusters).astype(object)
    labels = author_label_series(author_values)
    st = None if studies is None else np.asarray(studies, dtype=object)
    has_cluster = ~cl.isna().values
    n_labelled = n_agree = n_unlabelled = 0
    for _, idx in pd.Series(np.flatnonzero(has_cluster)).groupby(
            cl.values[has_cluster].astype(str)).groups.items():
        idx = np.asarray(idx)
        b = author_breakdown(labels[idx], None if st is None else st[idx])
        n_labelled += b.n_labelled
        n_agree += b.n_agree
        n_unlabelled += b.n_unlabelled
    purity = n_agree / n_labelled if n_labelled else float('nan')
    return {'n_labelled': n_labelled, 'n_agree': n_agree,
            'n_unlabelled': n_unlabelled, 'purity': purity}
