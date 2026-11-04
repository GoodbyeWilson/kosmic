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


def test_points_at_p_of_one_sit_inside_the_axes():
    """A gene with P = 1 is drawn at 0. With the axis starting at 0, half of
    each such point was clipped, and showed below the axis in a vector
    editor that released the clipping mask."""
    df = _frame()
    df.loc[4, 'pvals_adj'] = 1.0
    for y_limit in (None, 5):
        ax = create_volcano_plot(df, 'DCM', 'Donor', y_limit=y_limit).axes[0]
        bottom, top = ax.get_ylim()
        assert -0.05 * top < bottom < 0


def test_unnamed_genes_are_recognised():
    from kosmic.visualisation.de.volcano import is_unnamed_gene
    for name in ('AC004936.1', 'AL008633.1', 'ENSG00000268744', 'LINC00670',
                 'CDC42EP3-AS1', 'CRIM1-DT', 'RP11-34P13.3', 'RP11-442N24--B.1'):
        assert is_unnamed_gene(name), name
    for name in ('ABCG2', 'NKX2-5', 'HLA-A', 'MT-ND1', 'C1orf232', 'LINC-PINT', 'AQP3'):
        assert not is_unnamed_gene(name), name


def test_named_labels_only_skips_unnamed_genes():
    """Without them the next named genes are labelled; the genes still plot."""
    from kosmic.visualisation.de.volcano import _select_genes_to_label
    df = pd.DataFrame({
        'names': ['AC004936.1', 'LINC00670', 'ABCG2', 'EXT1', 'CDC42EP3-AS1', 'CHDH'],
        'logfoldchanges': [3.0, 2.5, 2.0, 1.5, 1.0, -2.0],
        'pvals_adj': [1e-30, 1e-25, 1e-20, 1e-15, 1e-10, 1e-12],
    })
    every = _select_genes_to_label(df, 'pvals_adj', 0.05, 0, None, 4, True)
    named = _select_genes_to_label(df, 'pvals_adj', 0.05, 0, None, 4, False)
    assert list(df.loc[every, 'names']) == ['AC004936.1', 'LINC00670', 'CHDH']
    assert list(df.loc[named, 'names']) == ['ABCG2', 'EXT1', 'CHDH']
    # Genes asked for by name are labelled either way.
    asked = _select_genes_to_label(df, 'pvals_adj', 0.05, 0, ['LINC00670'], 4, False)
    assert df.loc[asked[0], 'names'] == 'LINC00670'
