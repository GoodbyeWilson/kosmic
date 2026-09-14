# Shared GUI Widgets
# Reusable PyQt6 widget components used across tabs.

from PyQt6.QtWidgets import QTableWidget, QApplication
from PyQt6.QtGui import QKeySequence


class CopyableTableWidget(QTableWidget):
    """QTableWidget with Ctrl+C support.

    Copies selected cells as tab-separated text (pastes cleanly into
    Excel, Google Sheets, and plain-text editors).  If nothing is
    selected, copies the entire table including headers.
    """

    def keyPressEvent(self, event):
        if event.matches(QKeySequence.StandardKey.Copy):
            self.copy_selection()
        else:
            super().keyPressEvent(event)

    def copy_selection(self):
        selection = self.selectedRanges()
        if selection:
            rows: set = set()
            cols: set = set()
            for r in selection:
                rows.update(range(r.topRow(), r.bottomRow() + 1))
                cols.update(range(r.leftColumn(), r.rightColumn() + 1))
            rows = sorted(rows)
            cols = sorted(cols)
        else:
            rows = list(range(self.rowCount()))
            cols = list(range(self.columnCount()))

        header = []
        for c in cols:
            h = self.horizontalHeaderItem(c)
            header.append(h.text() if h else f"Col {c}")
        lines = ['\t'.join(header)]

        for r in rows:
            row_data = []
            for c in cols:
                item = self.item(r, c)
                row_data.append(item.text() if item else '')
            lines.append('\t'.join(row_data))

        QApplication.clipboard().setText('\n'.join(lines))


