# "Add donor metadata from table..." on the Inspect step.
#
# A paper's supplementary table joined onto obs by sample id. The dialog
# shows how the table's key column lines up with the study's samples,
# lists the table's columns with their distinct values so the user can
# choose which to bring in (sex- and age-like columns pre-selected,
# nothing else), and applies the chosen ones. The analysis lives in
# kosmic.scrna.inspect.donor_metadata; this is the form around it.
from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QDialog, QFileDialog, QHBoxLayout, QHeaderView, QLabel, QLineEdit,
    QTableWidgetItem, QVBoxLayout,
)

from kosmic.gui.shared.theme import NoScrollComboBox
from kosmic.gui.shared.widgets import (
    CaptionLabel, CopyableTableWidget, PrimaryButton, SecondaryButton,
    SecondaryLabel, StatusLabel,
)
from kosmic.scrna.inspect.donor_metadata import (
    column_summary, match_column, match_report, read_table, safe_obs_name,
    suggested_columns,
)


class DonorMetadataDialog(QDialog):
    """Pick a table, match its key column to the sample column, choose
    columns, apply. ``result()`` is Accepted when something was applied;
    ``chosen`` then holds ``{table column: obs name}``, ``table`` the
    DataFrame and ``key_column`` the matched column."""

    def __init__(self, adata, sample_col: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Add donor metadata from table")
        self.resize(760, 560)
        self._adata = adata
        self._sample_col = sample_col
        self._samples = [str(s) for s in adata.obs[sample_col].astype(str).unique()]
        self.table = None
        self.key_column = None
        self.chosen: dict = {}

        lay = QVBoxLayout(self)
        lay.setSpacing(8)
        lay.addWidget(CaptionLabel(
            "A CSV, TSV or Excel table with one row per donor or sample, "
            "such as a paper's supplementary table. One of its columns has "
            f"to hold the same ids as the study's sample column "
            f"('{sample_col}', {len(self._samples)} samples). Only the "
            "columns you tick are added; each becomes one value per cell."))

        file_row = QHBoxLayout()
        self._path_edit = QLineEdit()
        self._path_edit.setPlaceholderText("Choose a table file...")
        self._path_edit.setReadOnly(True)
        browse = SecondaryButton("Browse...")
        browse.clicked.connect(self._browse)
        file_row.addWidget(self._path_edit, 1)
        file_row.addWidget(browse)
        lay.addLayout(file_row)

        key_row = QHBoxLayout()
        key_row.addWidget(QLabel("Column holding the sample ids:"))
        self._key_combo = NoScrollComboBox()
        self._key_combo.currentTextChanged.connect(self._on_key_changed)
        key_row.addWidget(self._key_combo, 1)
        lay.addLayout(key_row)

        self._match_label = StatusLabel("")
        self._match_label.setWordWrap(True)
        lay.addWidget(self._match_label)

        lay.addWidget(SecondaryLabel(
            "Columns to add (tick). Sex- and age-like columns are ticked "
            "for you; the rest are off, because every column added is one "
            "more entry in every metadata list."))
        self._cols = CopyableTableWidget()
        self._cols.setColumnCount(4)
        self._cols.setHorizontalHeaderLabels(
            ["Table column", "Distinct", "Examples", "Name in KOSMIC"])
        hdr = self._cols.horizontalHeader()
        hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        hdr.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        hdr.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        hdr.setSectionResizeMode(3, QHeaderView.ResizeMode.Interactive)
        self._cols.setColumnWidth(3, 200)
        lay.addWidget(self._cols, 1)

        btn_row = QHBoxLayout()
        btn_row.addStretch()
        cancel = SecondaryButton("Cancel")
        cancel.clicked.connect(self.reject)
        self._apply_btn = PrimaryButton("Add columns")
        self._apply_btn.setEnabled(False)
        self._apply_btn.clicked.connect(self._apply)
        btn_row.addWidget(cancel)
        btn_row.addWidget(self._apply_btn)
        lay.addLayout(btn_row)

    # ------------------------------------------------------------------
    def _browse(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Donor metadata table", "",
            "Tables (*.csv *.tsv *.txt *.xlsx *.xls);;All files (*)")
        if path:
            self.load_table(path)

    def load_table(self, path):
        """Read the table and populate the key combo and column list."""
        try:
            self.table = read_table(path)
        except Exception as exc:  # noqa: BLE001 - shown to the user
            self._match_label.setText(f"Could not read {Path(path).name}: {exc}")
            self._match_label.set_state('warning')
            return
        self._path_edit.setText(str(path))
        guess = match_column(self.table, self._samples)
        self._key_combo.blockSignals(True)
        self._key_combo.clear()
        self._key_combo.addItems([str(c) for c in self.table.columns])
        if guess is not None:
            self._key_combo.setCurrentText(str(guess))
        self._key_combo.blockSignals(False)
        self._on_key_changed(self._key_combo.currentText())

    def _on_key_changed(self, key):
        if self.table is None or not key:
            return
        self.key_column = key
        rep = match_report(self.table, key, self._samples)
        text = (f"{rep.n_matched} of {rep.n_samples} samples matched on "
                f"'{key}'.")
        state = 'success'
        if rep.unmatched_samples:
            shown = ', '.join(rep.unmatched_samples[:6])
            more = f" (+{len(rep.unmatched_samples) - 6} more)" if len(rep.unmatched_samples) > 6 else ""
            text += f" Not in the table, will get blank values: {shown}{more}."
            state = 'warning'
        if rep.unmatched_rows:
            text += f" {len(rep.unmatched_rows)} table row(s) match no sample and are ignored."
        if rep.duplicate_keys:
            text += (f" Duplicate ids in the table (first row used): "
                     f"{', '.join(rep.duplicate_keys[:4])}.")
            state = 'warning'
        if rep.n_matched == 0:
            state = 'warning'
        self._match_label.setText(text)
        self._match_label.set_state(state)
        self._populate_columns()
        self._apply_btn.setEnabled(rep.n_matched > 0)

    def _populate_columns(self):
        self._cols.setRowCount(0)
        existing = list(self._adata.obs.columns)
        pre = set(suggested_columns(self.table, self.key_column, existing))
        for col in self.table.columns:
            if col == self.key_column:
                continue
            row = self._cols.rowCount()
            self._cols.insertRow(row)
            item = QTableWidgetItem(str(col))
            item.setFlags((item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                          & ~Qt.ItemFlag.ItemIsEditable)
            item.setCheckState(Qt.CheckState.Checked if col in pre
                               else Qt.CheckState.Unchecked)
            self._cols.setItem(row, 0, item)
            n, examples = column_summary(self.table, col)
            for c, text in ((1, str(n)), (2, examples)):
                cell = QTableWidgetItem(text)
                cell.setFlags(cell.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self._cols.setItem(row, c, cell)
            name = QTableWidgetItem(safe_obs_name(col, existing))
            name.setToolTip("The obs column this will be written to. Edit "
                            "to rename; an existing name is not overwritten.")
            self._cols.setItem(row, 3, name)

    def selected(self) -> dict:
        """``{table column: obs name}`` for the ticked rows."""
        out = {}
        for r in range(self._cols.rowCount()):
            item = self._cols.item(r, 0)
            if item is not None and item.checkState() == Qt.CheckState.Checked:
                name = self._cols.item(r, 3).text().strip() or safe_obs_name(item.text())
                out[item.text()] = name
        return out

    def _apply(self):
        chosen = self.selected()
        if not chosen:
            self._match_label.setText("Tick at least one column to add.")
            self._match_label.set_state('warning')
            return
        clash = [n for n in chosen.values() if n in self._adata.obs.columns]
        if clash:
            self._match_label.setText(
                "These names already exist in the dataset; rename them: "
                + ", ".join(clash))
            self._match_label.set_state('warning')
            return
        self.chosen = chosen
        self.accept()
