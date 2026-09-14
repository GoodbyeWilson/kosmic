# DataFrame-backed sortable results table.
#
# Replaces ~14 hand-rolled 'setHorizontalHeaderLabels' + 'for row in df.iterrows()'
# loops scattered across the GUI.  Usage:
#
#     cols = [Column('Gene', 'names', 's'),
#             Column('log2FC', 'logfoldchanges', '.3f'),
#             Column('FDR', 'fdr', '.2e')]
#     table = ResultsTable()
#     table.set_schema(cols)
#     table.set_data(df)
#
# The per-cell 'fmt' string drives both the display text and the hidden
# 'UserRole' numeric payload that makes the column sort numerically
# instead of lexicographically ("1.00e-10" < "9.00e-01" when sorted as
# text).  'N/A' values sort to the bottom via '+inf'.
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np
import pandas as pd
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QTableWidgetItem

from kosmic.gui.shared.widgets.copyable_table import CopyableTableWidget


@dataclass(frozen=True)
class Column:
    """One column definition for 'ResultsTable'.

    Attributes
    ----------
    header     : Display label for the column header.
    key        : DataFrame column name.  Missing keys render as empty/N/A.
    fmt        : One of:
                 - 's'    string (as-is; NaN renders empty)
                 - 'd'    integer (NaN renders 'na_text')
                 - 'bool' Yes/blank by default; customise with
                   'true_text' / 'false_text'.  Sortable via UserRole 1/0.
                 - any Python format spec for floats, e.g. '.3f', '.2e',
                   '.1%', '.1f'.  NaN/non-finite render 'na_text'
                   and sort last.
    tooltip    : Optional tooltip text for the column header.
    true_text  : Text for True cells when fmt='bool' (default 'Yes').
    false_text : Text for False cells when fmt='bool' (default '').
    na_text    : Text shown for NaN / non-finite numeric cells (default
                 'N/A').  Some tables use '-'.
    """
    header: str
    key: str
    fmt: str
    tooltip: str = ''
    true_text: str = 'Yes'
    false_text: str = ''
    na_text: str = 'N/A'


class _SortItem(QTableWidgetItem):
    """QTableWidgetItem that sorts by its UserRole numeric payload.

    'QTableWidget.setSortingEnabled(True)' calls '__lt__', which by
    default compares display text.  That breaks numeric columns
    ("1.00e-10" sorts after "9.00e-01" lexicographically).  This
    subclass compares UserRole floats instead, so migrated tables
    finally sort numerically.  Items without a UserRole payload
    fall back to text comparison (string columns behave unchanged).
    """

    def __lt__(self, other):  # noqa: D401
        try:
            a = self.data(Qt.ItemDataRole.UserRole)
            b = other.data(Qt.ItemDataRole.UserRole)
        except AttributeError:
            return super().__lt__(other)
        if a is None and b is None:
            return super().__lt__(other)
        if a is None:
            return False  # items without payload sort last
        if b is None:
            return True
        try:
            return float(a) < float(b)
        except (TypeError, ValueError):
            return super().__lt__(other)


class ResultsTable(CopyableTableWidget):
    """Schema + DataFrame driven sortable results table.

    - 'set_schema' installs column definitions and headers.
    - 'set_data' populates rows from a DataFrame, respecting the
      fmt of each column.  Cells are non-editable.  Sorting is disabled
      during population then re-enabled, and columns auto-size on finish.
    - Numeric + bool columns sort by magnitude, not by rendered text
      ('"1.00e-10"' sorts above '"9.00e-01"').  'N/A' sorts last.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._cols: tuple[Column, ...] = ()

    # Public API
    def set_schema(self, cols: Sequence[Column]) -> None:
        self._cols = tuple(cols)
        self.setColumnCount(len(self._cols))
        self.setHorizontalHeaderLabels([c.header for c in self._cols])
        for i, c in enumerate(self._cols):
            if c.tooltip:
                hdr = self.horizontalHeaderItem(i)
                if hdr is not None:
                    hdr.setToolTip(c.tooltip)

    def set_data(
        self,
        df: pd.DataFrame,
        decorate: Optional[
            Callable[[QTableWidgetItem, pd.Series, int], None]
        ] = None,
    ) -> None:
        """Populate rows from 'df'.

        'decorate': optional per-cell hook, called after each cell is
        built with '(item, row_series, col_index)'.  Callers use this
        to apply bold font, foreground colour, cell-specific tooltip,
        extra UserRole payload for drill-down, etc. without subclassing.
        """
        if not self._cols:
            raise RuntimeError(
                "ResultsTable.set_data() called before set_schema()")
        self.setSortingEnabled(False)
        self.setRowCount(0)
        self.setRowCount(len(df))
        for row_i, (_, row) in enumerate(df.iterrows()):
            for col_i, col in enumerate(self._cols):
                item = self._make_item(row.get(col.key), col)
                if decorate is not None:
                    decorate(item, row, col_i)
                self.setItem(row_i, col_i, item)
        self.setSortingEnabled(True)
        self.resizeColumnsToContents()

    # Cell construction
    @staticmethod
    def _make_item(val, col: Column) -> QTableWidgetItem:
        """Build one cell per the column's fmt spec.

        All cells are non-editable.  Numeric and bool cells store a sort
        key in UserRole (via '_SortItem') so sorting respects magnitude.
        """
        text, sort_key = format_value(val, col)
        if sort_key is None:
            item = QTableWidgetItem()
        else:
            item = _SortItem()
            item.setData(Qt.ItemDataRole.UserRole, sort_key)
        item.setText(text)
        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
        return item


# Helpers
def format_value(val, col: Column) -> tuple[str, Optional[float]]:
    """Return '(display_text, sort_key)' for a cell per the column fmt.

    'sort_key' is a float for numeric / bool columns (float('inf') for
    N/A, so they sort last) or None for plain strings (sort by text).
    Shared by the item table and the model-backed 'ResultsTableView' so
    both format cells identically.
    """
    fmt = col.fmt
    if fmt == 's':
        return ('' if _isna(val) else str(val), None)
    if fmt == 'bool':
        b = bool(val) and not _isna(val)
        return (col.true_text if b else col.false_text, 1.0 if b else 0.0)
    # Numeric: 'd' | ',d' | '.Nf' | '.Ne' | '.N%' | etc. '*d' specs
    # require int input, so cast; other specs take float.
    if _is_finite_number(val):
        f = float(val)
        text = f'{int(f):{fmt}}' if fmt.endswith('d') else f'{f:{fmt}}'
        return (text, f)
    return (col.na_text, float('inf'))
def _isna(val) -> bool:
    """Scalar NaN/None test. Fast paths for the common cell types so
    per-cell formatting doesn't pay 'pd.isna' overhead on big tables."""
    if val is None:
        return True
    if isinstance(val, (float, np.floating)):
        return val != val  # NaN
    if isinstance(val, (int, bool, str, np.integer)):
        return False
    try:
        return bool(pd.isna(val))
    except (TypeError, ValueError):
        return False


def _is_finite_number(val) -> bool:
    if val is None:
        return False
    if isinstance(val, bool):  # bool is a subclass of int; treat as non-numeric
        return False
    if isinstance(val, (int, np.integer)):
        return True
    if isinstance(val, (float, np.floating)):
        return np.isfinite(val)
    return False


__all__ = ["Column", "ResultsTable", "format_value"]
