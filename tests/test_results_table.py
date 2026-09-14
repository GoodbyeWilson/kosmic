"""Sentinel for the shared results-table cell formatter.

'format_value' is the single source of truth for both the item-based
ResultsTable and the virtualized ResultsTableView, so their cell text and
sort order stay identical. This guards the display text + numeric sort key
(magnitude, not lexical; N/A sorts last).
"""
from __future__ import annotations

import numpy as np

from kosmic.gui.shared.widgets.results_table import Column, format_value


def test_format_value_text_and_sort_keys():
    s = Column('Gene', 'names', 's')
    assert format_value('TNNT2', s) == ('TNNT2', None)      # strings sort by text
    assert format_value(None, s) == ('', None)
    assert format_value(float('nan'), s) == ('', None)

    e = Column('FDR', 'pvals_adj', '.2e')
    text, key = format_value(1.43e-23, e)
    assert text == '1.43e-23' and key == 1.43e-23           # numeric key = magnitude
    # N/A gets +inf so it sorts last, not as text.
    assert format_value(np.nan, e) == ('N/A', float('inf'))

    f = Column('log2FC', 'logfoldchanges', '.3f')
    assert format_value(-0.6582, f) == ('-0.658', -0.6582)

    b = Column('Sig', '_significant', 'bool')
    assert format_value(True, b) == ('Yes', 1.0)
    assert format_value(False, b) == ('', 0.0)

    d = Column('n', 'n', 'd')
    assert format_value(12.0, d) == ('12', 12.0)


def test_numeric_sort_key_beats_lexical():
    """1e-10 must sort before 1e-9 by magnitude (the bug the sort key fixes)."""
    e = Column('p', 'p', '.2e')
    _, k_small = format_value(1e-10, e)
    _, k_big = format_value(1e-9, e)
    assert k_small < k_big                       # magnitude: e-10 before e-9
    # Rendered text sorts them the opposite (wrong) way lexically.
    assert format_value(1e-10, e)[0] > format_value(1e-9, e)[0]
