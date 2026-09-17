"""Cluster-vs-author-label summaries (kosmic.scrna.annotate.author_labels)."""
import numpy as np
import pandas as pd

from kosmic.scrna.annotate.author_labels import (
    author_breakdown, author_label_series, cluster_purity,
)


def test_unlabelled_values_are_blanked():
    s = pd.Categorical(['CM', 'unknown', None, 'NA', '', 'Fib', 'Unassigned'])
    out = author_label_series(s)
    assert list(out) == ['CM', '', '', '', '', 'Fib', '']


def test_breakdown_top_and_describe():
    b = author_breakdown(np.array(['CM'] * 90 + ['Fib'] * 10 + [''] * 5, dtype=object))
    assert b.n_labelled == 100 and b.n_unlabelled == 5
    assert b.top == ('CM', 90)
    text = b.describe()
    assert 'CM: 90 (90.0%)' in text and 'no author label: 5' in text


def test_breakdown_all_unlabelled():
    b = author_breakdown(np.array(['', ''], dtype=object))
    assert b.n_labelled == 0 and b.top == ('', 0)


def test_cluster_purity_ignores_unclustered_and_unlabelled():
    clusters = ['0'] * 4 + ['1'] * 4 + [np.nan] * 2
    authors = ['CM', 'CM', 'CM', 'Fib',      # cluster 0: 3 of 4 agree
               'Fib', 'Fib', 'unknown', None,  # cluster 1: 2 of 2 labelled agree
               'CM', 'CM']                    # no cluster: ignored
    r = cluster_purity(clusters, authors)
    assert r['n_labelled'] == 6 and r['n_agree'] == 5 and r['n_unlabelled'] == 2
    assert abs(r['purity'] - 5 / 6) < 1e-12


def test_cluster_purity_nothing_labelled():
    r = cluster_purity(['0', '0'], [None, 'unknown'])
    assert r['n_labelled'] == 0 and np.isnan(r['purity'])
