# The log panel at the bottom of the main window.
#
# Every workspace's 'log_message' signal ends up here, timestamped.
# It is the user-visible record of what the app is doing during a long
# run and the first place to look when something went wrong.
from __future__ import annotations

from datetime import datetime

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QCursor
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QTextEdit,
)


class OutputPanel(QWidget):
    """Always-visible output log with a header row (title + Clear) and a read-only text area."""

    def __init__(self, parent=None):
        super().__init__(parent)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 6, 10, 6)
        layout.setSpacing(4)

        header_row = QHBoxLayout()
        header_row.setSpacing(0)

        self._title = QLabel("OUTPUT")
        self._title.setProperty("role", "section_header")
        _tf = self._title.font()
        _tf.setLetterSpacing(QFont.SpacingType.AbsoluteSpacing, 1.5)
        self._title.setFont(_tf)
        self._title.setFixedHeight(16)
        header_row.addWidget(self._title)
        header_row.addStretch()

        self._clear_btn = QPushButton("Clear")
        self._clear_btn.setProperty("role", "flat_text")
        self._clear_btn.setFixedHeight(16)
        self._clear_btn.setFlat(True)
        self._clear_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._clear_btn.clicked.connect(self.clear)
        header_row.addWidget(self._clear_btn)

        layout.addLayout(header_row)

        self._log_text = QTextEdit()
        self._log_text.setProperty("role", "log_output")
        self._log_text.setReadOnly(True)
        self._log_text.setFont(QFont("Consolas", 9))
        self._log_text.setPlaceholderText("No output yet.")
        self._log_text.setFrameShape(QTextEdit.Shape.NoFrame)
        layout.addWidget(self._log_text)

    def log(self, message: str):
        """Append a timestamped message and scroll to the bottom."""
        ts = datetime.now().strftime("%H:%M:%S")
        self._log_text.append(f"[{ts}] {message}")
        sb = self._log_text.verticalScrollBar()
        sb.setValue(sb.maximum())

    def clear(self):
        self._log_text.clear()

    def refresh_theme(self):
        for w in (self._title, self._clear_btn, self._log_text):
            w.style().unpolish(w)
            w.style().polish(w)
