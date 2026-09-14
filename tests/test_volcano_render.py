"""Sentinel: the volcano plot must render even when the DE frame has a
non-unique index or a duplicated column -- either used to raise
'truth value of a Series/DataFrame is ambiguous' and blanked the figure export."""
from __future__ import annotations

import matplotlib
matplotlib.use('Agg')
import pandas as pd

from kosmic.visualisation.de.volcano import create_volcano_plot


def _frame():
    return pd.DataFrame({
        'names': ['ENO1', 'PKM', 'LDHA', 'SDHB', 'CS', 'IDH1'],
        'logfoldchanges': [-0.42, 0.31, 0.8, -0.19, 0.05, -0.6],
        'pvals': [0.007, 0.4, 0.001, 0.09, 0.7, 0.002],
        'pvals_adj': [0.014, 0.47, 0.004, 0.13, 0.72, 0.006],
    })


def test_volcano_renders_with_clean_frame():
    fig = create_volcano_plot(_frame(), 'DCM', 'Donor')
    assert fig is not None


def test_volcano_survives_non_unique_index():
    df = _frame()
    df.index = ['ENO1'] * len(df)  # degenerate index from a names-indexed frame
    fig = create_volcano_plot(df, 'DCM', 'Donor')
    assert fig is not None


def test_volcano_survives_duplicate_column():
    df = _frame()
    df = pd.concat([df, df[['pvals_adj']]], axis=1)  # duplicated 'pvals_adj'
    fig = create_volcano_plot(df, 'DCM', 'Donor')
    assert fig is not None


def test_volcano_survives_dataframe_valued_attrs():
    # Hypothesis-mode DE stashes a frame in attrs['genome_wide_ranking']; pandas
    # propagates and compares attrs during nsmallest, which used to blow up.
    df = _frame()
    df.attrs['genome_wide_ranking'] = df.copy()
    fig = create_volcano_plot(df, 'DCM', 'Donor', genes_of_interest={'LDHA', 'PKM'})
    assert fig is not None
