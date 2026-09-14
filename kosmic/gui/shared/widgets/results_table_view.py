# Model-backed results table for large DataFrames.
#
# A 'QTableView' + 'QAbstractTableModel' twin of 'ResultsTable' that only
# renders visible rows, so a 26k-gene DE table shows instantly instead of
# building one QTableWidgetItem per cell. Shares 'Column' and
# 'format_value' with the item table, so cell text and numeric ordering are
# identical. Use it for the few genuinely large tables; the item-based
# 'ResultsTable' stays fine for the many small ones.
#
# Cells are formatted once in 'set_data' into per-column lists, so 'data()'
# is a plain lookup. Sorting reorders those lists via a single numpy
# argsort on the precomputed keys -- no per-comparison Python callback,
# which a QSortFilterProxyModel would otherwise incur ~n*log(n) times.

from __future__ import annotations

from typing import List, Optional, Sequence

import numpy as np
import pandas as pd
from PyQt6.QtCore import QAbstractTableModel, QModelIndex, Qt
from PyQt6.QtGui import QKeySequence
from PyQt6.QtWidgets import QApplication, QTableView

from kosmic.gui.shared.widgets.results_table import Column, format_value


class _DataFrameModel(QAbstractTableModel):
    """Serves precomputed cell text; sorts by precomputed keys via argsort."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._cols: tuple[Column, ...] = ()
        self._nrows = 0
        self._text: List[List[str]] = []            # [col][row] display text
        self._keys: List[Optional[np.ndarray]] = []  # [col] numeric sort key, or None

    def set_schema(self, cols: Sequence[Column]) -> None:
        self.beginResetModel()
        self._cols = tuple(cols)
        self.endResetModel()

    def set_data(self, df: pd.DataFrame) -> None:
        self.beginResetModel()
        self._nrows = len(df)
        self._text = []
        self._keys = []
        for col in self._cols:
            vals = (df[col.key].to_numpy() if col.key in df.columns
                    else np.full(self._nrows, None, dtype=object))
            texts: List[str] = []
            keys = np.empty(self._nrows, dtype=float)
            has_key = col.fmt != 's'
            for i, v in enumerate(vals):
                t, k = format_value(v, col)
                texts.append(t)
                if has_key:
                    keys[i] = k
            self._text.append(texts)
            self._keys.append(keys if has_key else None)
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else self._nrows

    def columnCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._cols)

    def data(self, index: QModelIndex, role=Qt.ItemDataRole.DisplayRole):
        if index.isValid() and role == Qt.ItemDataRole.DisplayRole:
            return self._text[index.column()][index.row()]
        return None

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if orientation != Qt.Orientation.Horizontal or section >= len(self._cols):
            return None
        col = self._cols[section]
        if role == Qt.ItemDataRole.DisplayRole:
            return col.header
        if role == Qt.ItemDataRole.ToolTipRole and col.tooltip:
            return col.tooltip
        return None

    def flags(self, index):
        return Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable

    def sort(self, column, order=Qt.SortOrder.AscendingOrder):
        if column < 0 or column >= len(self._cols) or self._nrows == 0:
            return
        keys = self._keys[column]
        if keys is None:  # string column: order by display text
            perm = np.argsort(np.asarray(self._text[column], dtype=object), kind='stable')
        else:
            perm = np.argsort(keys, kind='stable')
        if order == Qt.SortOrder.DescendingOrder:
            perm = perm[::-1]
        self.layoutAboutToBeChanged.emit()
        self._text = [[tcol[i] for i in perm] for tcol in self._text]
        self._keys = [None if k is None else k[perm] for k in self._keys]
        self.layoutChanged.emit()

    def value_at(self, row: int, col_index: int) -> str:
        return self._text[col_index][row]


class ResultsTableView(QTableView):
    """Schema + DataFrame driven sortable table, virtualized for large data.

    Same 'set_schema' / 'set_data' surface as 'ResultsTable'. Ctrl+C copies
    the selection (or the whole table) as TSV, matching the item table.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._model = _DataFrameModel(self)
        self.setModel(self._model)
        self.setSortingEnabled(True)
        self.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
        self.verticalHeader().setVisible(False)

    # Public API mirroring ResultsTable
    def set_schema(self, cols: Sequence[Column]) -> None:
        self._model.set_schema(cols)

    def set_data(self, df: pd.DataFrame) -> None:
        self._model.set_data(df)

    def rowCount(self) -> int:
        return self._model.rowCount()

    def select_first_by_text(self, col_index: int, text: str) -> bool:
        """Select + scroll to the first row whose column text equals 'text'."""
        for row in range(self._model.rowCount()):
            if self._model.value_at(row, col_index) == text:
                self.selectRow(row)
                self.scrollTo(self._model.index(row, col_index))
                return True
        return False

    # Copy (TSV, headers, whole table if nothing selected) -- matches CopyableTableWidget
    def keyPressEvent(self, event):
        if event.matches(QKeySequence.StandardKey.Copy):
            self._copy_selection()
        else:
            super().keyPressEvent(event)

    def _copy_selection(self) -> None:
        cols = list(range(self._model.columnCount()))
        sel = self.selectionModel().selectedRows()
        rows = ([i.row() for i in sel] if sel
                else list(range(self._model.rowCount())))
        header = [self._model.headerData(c, Qt.Orientation.Horizontal) for c in cols]
        lines = ['\t'.join(str(h or '') for h in header)]
        for r in rows:
            lines.append('\t'.join(self._model.value_at(r, c) for c in cols))
        QApplication.clipboard().setText('\n'.join(lines))


__all__ = ["ResultsTableView"]
