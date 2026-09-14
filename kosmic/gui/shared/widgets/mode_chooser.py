# Shared mode-chooser page for workspaces that branch into 2-3 sub-flows.
#
# Both the DE workspace ('Scoring' vs 'Discovery') and the Meta-Analysis
# workspace ('Hypothesis' / 'Discovery' / 'Methods Comparison') used to
# own their own near-identical 220-line 'Choose*Page' files. This module
# collapses them into:
#
# - 'ModeOption' -- describes one card: key, icon, title, subtitle, steps.
# - '_ModeCard' -- the styled clickable card (private widget).
# - 'ModeChooserPage' -- a title/subtitle + horizontal row of cards.
#
# Each workspace's "choose mode" page is now a thin wrapper that supplies a
# list of 'ModeOption'.
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QCursor, QFont
from PyQt6.QtWidgets import (
    QFrame, QHBoxLayout, QLabel, QVBoxLayout, QWidget,
)

from kosmic.gui.shared.icon_provider import make_icon


@dataclass
class ModeOption:
    """One card in the picker."""
    key: str
    icon_name: str
    title: str
    subtitle: str
    steps: List[str]
    accent_color: str


class _ModeCard(QFrame):
    """Clickable card used by 'ModeChooserPage'."""

    clicked = pyqtSignal(str)

    def __init__(self, option: ModeOption, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._option = option
        self._selected = False

        self.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.setFixedWidth(360)
        self.setMinimumHeight(320)

        self._build_ui()
        self._apply_style(hover=False)

    def _build_ui(self):
        opt = self._option
        lay = QVBoxLayout(self)
        lay.setContentsMargins(24, 24, 24, 24)
        lay.setSpacing(0)

        # Icon + title row.
        header = QHBoxLayout()
        header.setSpacing(12)
        icon_label = QLabel()
        icon_label.setPixmap(
            make_icon(opt.icon_name, opt.accent_color, 32).pixmap(32, 32))
        icon_label.setFixedSize(32, 32)
        header.addWidget(icon_label)

        title = QLabel(opt.title)
        title.setFont(QFont("Segoe UI", 15, QFont.Weight.Bold))
        header.addWidget(title, 1)
        lay.addLayout(header)
        lay.addSpacing(8)

        # Subtitle.
        sub = QLabel(opt.subtitle)
        sub.setWordWrap(True)
        sub.setFont(QFont("Segoe UI", 10))
        lay.addWidget(sub)
        lay.addSpacing(16)

        # Divider.
        div = QFrame()
        div.setFrameShape(QFrame.Shape.HLine)
        div.setFixedHeight(1)
        lay.addWidget(div)
        lay.addSpacing(12)

        # WORKFLOW heading.
        steps_title = QLabel("WORKFLOW")
        steps_title.setObjectName("mode_card_steps_title")
        steps_title.setFont(QFont("Segoe UI", 9, QFont.Weight.Bold))
        lay.addWidget(steps_title)
        lay.addSpacing(8)

        # Numbered steps.
        for i, step in enumerate(opt.steps, 1):
            row = QHBoxLayout()
            row.setSpacing(10)
            row.setContentsMargins(0, 2, 0, 2)

            num = QLabel(str(i))
            num.setObjectName("step_badge")
            num.setFixedSize(22, 22)
            num.setAlignment(Qt.AlignmentFlag.AlignCenter)
            num.setFont(QFont("Segoe UI", 9, QFont.Weight.Bold))
            row.addWidget(num)

            step_label = QLabel(step)
            step_label.setWordWrap(True)
            step_label.setFont(QFont("Segoe UI", 10))
            row.addWidget(step_label, 1)

            lay.addLayout(row)

        lay.addStretch()

    def _apply_style(self, hover: bool = False):
        state = ("selected" if self._selected
                 else ("hover" if hover else "normal"))
        self.setProperty("state", state)
        style = self.style()
        if style is not None:
            style.unpolish(self)
            style.polish(self)

    def set_selected(self, selected: bool):
        self._selected = selected
        self._apply_style()

    @property
    def key(self) -> str:
        return self._option.key

    def enterEvent(self, event):
        if not self._selected:
            self._apply_style(hover=True)
        super().enterEvent(event)

    def leaveEvent(self, event):
        self._apply_style()
        super().leaveEvent(event)

    def mousePressEvent(self, event):
        self.clicked.emit(self._option.key)


class ModeChooserPage(QWidget):
    """Title + subtitle + a horizontal row of clickable mode cards.

    Emits 'mode_selected(key)' when the user clicks a card. Tracks the
    most recently selected card so that re-activating the page after the
    user navigated away keeps the visual selection in sync.
    """

    mode_selected = pyqtSignal(str)

    def __init__(self, title: str, subtitle: str,
                 options: List[ModeOption],
                 parent: Optional[QWidget] = None):
        super().__init__(parent)
        if not options:
            raise ValueError("ModeChooserPage requires at least one ModeOption.")
        self._title = title
        self._subtitle = subtitle
        self._options = options
        self._cards: List[_ModeCard] = []
        self._current_mode: Optional[str] = None
        self._setup_ui()

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(40, 40, 40, 40)
        layout.setSpacing(0)
        layout.addStretch(2)

        title = QLabel(self._title)
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setFont(QFont("Segoe UI", 18, QFont.Weight.Bold))
        layout.addWidget(title)
        layout.addSpacing(6)

        subtitle = QLabel(self._subtitle)
        subtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        subtitle.setWordWrap(True)
        subtitle.setFont(QFont("Segoe UI", 11))
        subtitle.setProperty("role", "secondary")
        layout.addWidget(subtitle)
        layout.addSpacing(32)

        card_row = QHBoxLayout()
        card_row.setSpacing(24)
        card_row.addStretch()
        for opt in self._options:
            card = _ModeCard(opt)
            card.clicked.connect(self._on_card_clicked)
            self._cards.append(card)
            card_row.addWidget(card)
        card_row.addStretch()

        layout.addLayout(card_row)
        layout.addStretch(3)

    def _on_card_clicked(self, key: str):
        self._current_mode = key
        for card in self._cards:
            card.set_selected(card.key == key)
        self.mode_selected.emit(key)

    def set_current_mode(self, key: Optional[str]):
        """Programmatically select a card (without emitting the signal)."""
        self._current_mode = key
        for card in self._cards:
            card.set_selected(card.key == key)

    def on_activated(self):
        """Refresh selection state when the page becomes visible."""
        self.set_current_mode(self._current_mode)

    def refresh_theme(self):
        for card in self._cards:
            card._apply_style()
