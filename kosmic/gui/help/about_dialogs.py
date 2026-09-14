# The About dialog (Help -> About KOSMIC): name, version and the text in
# 'about.md' beside this module.
from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QDialog, QPushButton, QTextBrowser, QVBoxLayout,
)


_HELP_DIR = Path(__file__).parent
_ABOUT_TEXT = (_HELP_DIR / "about.md").read_text(encoding="utf-8")


def show_about(parent) -> None:
    dlg = QDialog(parent)
    dlg.setWindowTitle("About KOSMIC")
    dlg.resize(480, 360)
    layout = QVBoxLayout(dlg)

    text = QTextBrowser()
    text.setMarkdown(_ABOUT_TEXT)
    text.setOpenExternalLinks(True)
    layout.addWidget(text)

    close_btn = QPushButton("Close")
    close_btn.clicked.connect(dlg.accept)
    layout.addWidget(close_btn, alignment=Qt.AlignmentFlag.AlignRight)
    dlg.exec()
