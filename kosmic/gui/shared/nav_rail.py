# Vertical workspace navigation rail (left edge of the main window).
from __future__ import annotations

from PyQt6.QtCore import Qt, QRectF, pyqtSignal
from PyQt6.QtGui import QPainter, QColor, QFont, QCursor
from PyQt6.QtWidgets import QWidget, QVBoxLayout

from kosmic.gui.shared.theme import get_color, get_current_mode, ThemeMode
from kosmic.gui.shared.icon_provider import render_svg_to_painter


class _NavItem(QWidget):
    """Clickable rail item: centred SVG icon + label, with hover + active states."""

    clicked = pyqtSignal()
    INDICATOR_WIDTH = 3

    def __init__(self, icon_name: str, label: str, enabled: bool = True):
        super().__init__()
        self._enabled = enabled
        self._active = False
        self._hovered = False
        self._icon_name = icon_name
        self._label = label
        self.setFixedSize(64, 64)
        if enabled:
            self.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))

    def set_icon(self, name: str):
        self._icon_name = name
        self.update()

    def set_active(self, active: bool):
        self._active = active
        self.update()

    def refresh_theme(self):
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()

        if self._active:
            painter.fillRect(0, 0, w, h, QColor(get_color('bg_active')))
            painter.fillRect(0, 0, self.INDICATOR_WIDTH, h,
                             QColor(get_color('accent_primary')))
        elif self._hovered and self._enabled:
            painter.fillRect(0, 0, w, h, QColor(get_color('bg_hover')))

        fg = get_color(
            'fg_tertiary' if not self._enabled
            else 'fg_primary' if self._active
            else 'fg_secondary'
        )

        render_svg_to_painter(painter, self._icon_name, fg,
                              QRectF((w - 20) / 2, 8, 20, 20))
        painter.setPen(QColor(fg))
        painter.setFont(QFont("Segoe UI", 8))
        painter.drawText(
            QRectF(0, 30, w, h - 30),
            Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop,
            self._label,
        )
        painter.end()

    def mousePressEvent(self, event):
        if self._enabled and event.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit()
        super().mousePressEvent(event)

    def enterEvent(self, event):
        if self._enabled:
            self._hovered = True
            self.update()
        super().enterEvent(event)

    def leaveEvent(self, event):
        self._hovered = False
        self.update()
        super().leaveEvent(event)


class NavigationRail(QWidget):
    """Vertical rail with mode icons + tutorial/theme buttons."""

    mode_selected = pyqtSignal(str)      # "project" | "scrna" | "de" | "figures" | "meta"
    theme_toggled = pyqtSignal()
    tutorial_clicked = pyqtSignal()
    home_clicked = pyqtSignal()          # DNA logo at top of rail

    # (mode_key, icon, label, enabled). mode_key=None means not selectable.
    _MODE_ITEMS = (
        ("project", "folder-open",   "Project",                  True),
        ("scrna",   "umap",          "scRNA",                    True),
        ("de",      "volcano",       "Differential\nExpression", True),
        ("meta",    "forest-plot",   "Meta",                     True),
        ("figures", "image",         "Figure\nExport",           True),
    )

    def __init__(self):
        super().__init__()
        self.setFixedWidth(65)  # 64 content + 1px separator

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 48, 1, 4)  # top 48 = logo, right 1 = separator
        layout.setSpacing(0)

        self._items: dict[str, _NavItem] = {}
        self._all_items: list[_NavItem] = []
        for mode, icon, label, enabled in self._MODE_ITEMS:
            item = _NavItem(icon, label, enabled=enabled)
            if mode is not None:
                item.clicked.connect(lambda m=mode: self._select(m))
                self._items[mode] = item
            else:
                item.setToolTip("Combine Datasets — coming soon")
            layout.addWidget(item)
            self._all_items.append(item)

        layout.addStretch()

        self._tutorial_item = _NavItem("play-circle", "Tutorial")
        self._tutorial_item.clicked.connect(self.tutorial_clicked.emit)
        layout.addWidget(self._tutorial_item)
        self._all_items.append(self._tutorial_item)

        self._theme_item = _NavItem("moon", "Theme")
        self._theme_item.clicked.connect(self.theme_toggled.emit)
        layout.addWidget(self._theme_item)
        self._all_items.append(self._theme_item)

    def set_active(self, mode: str):
        for key, item in self._items.items():
            item.set_active(key == mode)

    def refresh_theme(self):
        self.update()
        self._theme_item.set_icon(
            "moon" if get_current_mode() == ThemeMode.DARK else "sun")
        for item in self._all_items:
            item.refresh_theme()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        painter.fillRect(0, 0, w, h, QColor(get_color('bg_secondary')))
        painter.fillRect(w - 1, 0, 1, h, QColor(get_color('border')))
        # Logo in the 64x48 top area, 24x24 centred
        render_svg_to_painter(
            painter, "dna", get_color('accent_primary'),
            QRectF((64 - 24) / 2, (48 - 24) / 2, 24, 24),
        )
        painter.end()

    def _select(self, mode: str):
        self.set_active(mode)
        self.mode_selected.emit(mode)

    def mousePressEvent(self, event):
        # Clicks in the 64x48 logo area at the top of the rail open Home.
        if (event.button() == Qt.MouseButton.LeftButton
                and event.position().y() < 48):
            self.home_clicked.emit()
            return
        super().mousePressEvent(event)

    def showEvent(self, event):
        super().showEvent(event)
        self.refresh_theme()
