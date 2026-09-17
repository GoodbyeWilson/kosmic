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
from __future__ import annotations

from dataclasses import dataclass

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
    """Author labels within one cluster."""
    counts: pd.Series          # label -> cells, descending, unlabelled removed
    n_unlabelled: int

    @property
    def n_labelled(self) -> int:
        return int(self.counts.sum())

    @property
    def top(self) -> tuple[str, int]:
        if self.counts.empty:
            return '', 0
        return str(self.counts.index[0]), int(self.counts.iloc[0])

    def describe(self, max_rows: int = 6) -> str:
        """Multi-line text: each author label with its cell count and share."""
        n = self.n_labelled
        lines = [f"{lb}: {c:,} ({100 * c / n:.1f}%)"
                 for lb, c in self.counts.head(max_rows).items()]
        rest = self.counts.iloc[max_rows:]
        if len(rest):
            lines.append(f"{len(rest)} other labels: {int(rest.sum()):,}")
        if self.n_unlabelled:
            lines.append(f"no author label: {self.n_unlabelled:,}")
        return "\n".join(lines)


def author_breakdown(labels: np.ndarray) -> AuthorBreakdown:
    """Count author labels in one cluster (labels from author_label_series)."""
    arr = np.asarray(labels, dtype=object)
    unlabelled = int((arr == '').sum())
    counts = pd.Series(arr[arr != '']).value_counts()
    return AuthorBreakdown(counts=counts, n_unlabelled=unlabelled)


def cluster_purity(clusters, author_values) -> dict:
    """Overall agreement between a clustering and the authors' labels.

    Returns n_labelled, n_agree (cells carrying their cluster's dominant
    author label), n_unlabelled and purity (n_agree / n_labelled, NaN when
    nothing is labelled). Cells with no cluster (NaN) are ignored.
    """
    cl = pd.Series(clusters).astype(object)
    labels = author_label_series(author_values)
    has_cluster = ~cl.isna().values
    n_labelled = n_agree = n_unlabelled = 0
    for _, idx in pd.Series(np.flatnonzero(has_cluster)).groupby(
            cl.values[has_cluster].astype(str)).groups.items():
        b = author_breakdown(labels[np.asarray(idx)])
        n_labelled += b.n_labelled
        n_agree += b.top[1]
        n_unlabelled += b.n_unlabelled
    purity = n_agree / n_labelled if n_labelled else float('nan')
    return {'n_labelled': n_labelled, 'n_agree': n_agree,
            'n_unlabelled': n_unlabelled, 'purity': purity}
