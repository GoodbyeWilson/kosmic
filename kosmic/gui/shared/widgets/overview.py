# Reusable overview primitives: titled key/value panels with status
# icons, big-number stat blocks, and rounded status chips. Composed by
# dataset / results overview screens across workspaces. Styling lives
# in theme.py (objectName 'info_panel'; roles 'panel_title', 'kv_key',
# 'kv_value', 'stat_value', 'chip_<state>').
from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QFrame, QGridLayout, QLabel, QVBoxLayout, QWidget

from kosmic.gui.shared.icon_provider import make_icon
from kosmic.gui.shared.theme import get_color
from kosmic.gui.shared.widgets.semantic import CaptionLabel


class StatusChip(QLabel):
    """Rounded pill for page-level state (e.g. 'Ready to continue').

    Use 'set_state' with one of 'success', 'info', 'warning'.
    """

    _STATES = ('success', 'info', 'warning')

    def __init__(self, text: str = "", parent=None, state: str = 'success'):
        super().__init__(text, parent)
        self.set_state(state)

    def set_state(self, state: str) -> None:
        if state not in self._STATES:
            raise ValueError(
                f"StatusChip state must be one of {self._STATES}, "
                f"got {state!r}")
        self.setProperty("role", f"chip_{state}")
        style = self.style()
        if style is not None:
            style.unpolish(self)
            style.polish(self)


class IconTile(QFrame):
    """Rounded icon tile with an optional tiny caption ('RDS', 'h5ad').

    With an icon, the caption renders small beneath it; without an
    icon, the caption renders larger as the tile's main text ('GEO').
    """

    def __init__(self, icon: str | None = None, caption: str = "",
                 size: int = 48, icon_size: int = 20, parent=None):
        super().__init__(parent)
        self.setObjectName("icon_tile")
        self.setFixedSize(size, size)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 4, 2, 4)
        layout.setSpacing(0)
        layout.addStretch()
        if icon:
            icon_label = QLabel()
            icon_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            icon_label.setPixmap(
                make_icon(icon, get_color('accent_primary'),
                          icon_size).pixmap(icon_size, icon_size))
            layout.addWidget(icon_label)
        self._caption_label = QLabel(caption)
        self._caption_label.setProperty(
            "role", "icon_tile_caption" if icon else "icon_tile_text")
        self._caption_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._caption_label.setVisible(bool(caption))
        layout.addWidget(self._caption_label)
        layout.addStretch()

    def set_caption(self, text: str) -> None:
        self._caption_label.setText(text)
        self._caption_label.setVisible(bool(text))


class StatBlock(QWidget):
    """Large value over a small caption (Cells / Genes / File size)."""

    def __init__(self, caption: str, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self._value = QLabel("")
        self._value.setProperty("role", "stat_value")
        layout.addWidget(self._value)
        layout.addWidget(CaptionLabel(caption))

    def set_value(self, text: str) -> None:
        self._value.setText(text)


class InfoPanel(QFrame):
    """Titled key/value panel with an optional per-row status icon.

    'add_row(key, label, tip)' declares a row up front;
    'set_value(key, text, ok=, tooltip=)' fills it later. 'ok=True'
    draws a green check, 'ok=False' an amber dot, 'ok=None' leaves the
    icon slot empty (purely informational row).
    """

    def __init__(self, title: str, parent=None):
        super().__init__(parent)
        self.setObjectName("info_panel")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 12, 16, 14)
        layout.setSpacing(10)

        title_label = QLabel(title)
        title_label.setProperty("role", "panel_title")
        layout.addWidget(title_label)

        self._grid = QGridLayout()
        self._grid.setHorizontalSpacing(10)
        self._grid.setVerticalSpacing(7)
        # Values get the spare width and wrap when there is none;
        # keys are short and take what they need. With the stretch
        # on the key column instead, a long value was clipped on the
        # left ('dy for meta-analysis').
        self._grid.setColumnStretch(0, 0)
        self._grid.setColumnStretch(1, 1)
        layout.addLayout(self._grid)
        layout.addStretch()

        self._rows: dict[str, tuple[QLabel, QLabel]] = {}

    def add_row(self, key: str, label: str, tip: str | None = None) -> None:
        r = len(self._rows)
        key_label = QLabel(label)
        key_label.setProperty("role", "kv_key")
        if tip:
            key_label.setToolTip(tip)
        value_label = QLabel("")
        value_label.setProperty("role", "kv_value")
        value_label.setWordWrap(True)
        value_label.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        icon_label = QLabel()
        icon_label.setFixedWidth(16)
        self._grid.addWidget(key_label, r, 0)
        self._grid.addWidget(value_label, r, 1)
        self._grid.addWidget(icon_label, r, 2)
        self._rows[key] = (value_label, icon_label)

    def set_value(self, key: str, text: str, *, ok: bool | None = None,
                  tooltip: str | None = None) -> None:
        value_label, icon_label = self._rows[key]
        value_label.setText(text)
        value_label.setToolTip(tooltip or "")
        if ok is None:
            icon_label.clear()
        else:
            colour = get_color('success' if ok else 'warning')
            icon = 'circle-check' if ok else 'circle-dot'
            icon_label.setPixmap(make_icon(icon, colour, 14).pixmap(14, 14))
