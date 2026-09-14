# The Explorer pane's workflow list and status panel.
#
# 'WorkflowSteps' draws the numbered steps of the visible workspace
# with active / complete state and emits 'step_clicked';
# 'StepStatusPanel' is the STATUS block beneath the file tree, with a
# one-line message and a progress bar. Both are owned by 'ExplorerPane'
# and driven by 'AppWindow' on workspace activation; workspaces never
# touch them directly.
from __future__ import annotations

from PyQt6.QtCore import Qt, QRectF, pyqtSignal
from PyQt6.QtGui import QPainter, QColor, QFont, QCursor
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QLabel, QProgressBar

from kosmic.gui.shared.borderless import borderless
from kosmic.gui.shared.theme import get_color
from kosmic.gui.shared.icon_provider import render_svg_to_painter


class WorkflowStepItem(QWidget):
    """Clickable step row with circle / circle-dot / circle-check icon."""

    clicked = pyqtSignal()

    def __init__(self, index: int, label: str, parent=None):
        super().__init__(parent)
        self._index = index
        self._label = label
        self._active = False
        self._complete = False
        self._hovered = False
        self.setFixedHeight(34)
        self.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))

    def set_active(self, active: bool):
        self._active = active
        self.update()

    def set_complete(self, complete: bool):
        self._complete = complete
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()

        if self._active:
            painter.fillRect(0, 0, w, h, QColor(get_color('bg_active')))
        elif self._hovered:
            painter.fillRect(0, 0, w, h, QColor(get_color('bg_hover')))
        else:
            painter.fillRect(0, 0, w, h, QColor(get_color('bg_secondary')))

        if self._active:
            painter.fillRect(0, 0, 3, h, QColor(get_color('accent_primary')))

        if self._complete:
            fg_hex = get_color('success')
            icon_name = "circle-check"
        elif self._active:
            fg_hex = get_color('fg_primary')
            icon_name = "circle-dot"
        else:
            fg_hex = get_color('fg_secondary')
            icon_name = "circle"

        icon_size = 12
        render_svg_to_painter(
            painter, icon_name, fg_hex,
            QRectF(12, (h - icon_size) / 2, icon_size, icon_size),
        )

        painter.setPen(QColor(fg_hex))
        font = QFont("Segoe UI", 10)
        if self._active:
            font.setWeight(QFont.Weight.Medium)
        painter.setFont(font)
        text_x = 28
        painter.drawText(text_x, 0, w - text_x - 30, h,
                         Qt.AlignmentFlag.AlignVCenter, self._label)

        if self._complete:
            check_size = 12
            render_svg_to_painter(
                painter, "check", get_color('success'),
                QRectF(w - 24, (h - check_size) / 2, check_size, check_size),
            )

        painter.end()

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit()
        super().mousePressEvent(event)

    def enterEvent(self, event):
        self._hovered = True
        self.update()
        super().enterEvent(event)

    def leaveEvent(self, event):
        self._hovered = False
        self.update()
        super().leaveEvent(event)


class WorkflowSteps(QWidget):
    """Container for WorkflowStepItem rows with a "WORKFLOW" section header."""

    step_clicked = pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._items: list[WorkflowStepItem] = []

        layout = borderless(QVBoxLayout, self)

        self._header = QLabel("WORKFLOW")
        self._header.setProperty("role", "section_header")
        _header_font = self._header.font()
        _header_font.setLetterSpacing(QFont.SpacingType.AbsoluteSpacing, 1.5)
        self._header.setFont(_header_font)
        self._header.setFixedHeight(24)
        self._header.setContentsMargins(10, 4, 0, 0)
        layout.addWidget(self._header)

        self._steps_container = QWidget()
        self._steps_layout = borderless(QVBoxLayout, self._steps_container)
        layout.addWidget(self._steps_container)

    def set_steps(self, steps):
        """Replace all steps. 'steps' is list[str] or list[tuple[str, str]] (name, hint)."""
        self.clear()
        for i, entry in enumerate(steps):
            name, hint = entry if isinstance(entry, tuple) else (entry, "")
            item = WorkflowStepItem(i, name)
            if hint:
                item.setToolTip(hint)
            item.clicked.connect(lambda idx=i: self._on_item_clicked(idx))
            self._steps_layout.addWidget(item)
            self._items.append(item)

    def set_active(self, index: int):
        for i, item in enumerate(self._items):
            item.set_active(i == index)

    def set_complete(self, index: int, complete: bool = True):
        if 0 <= index < len(self._items):
            self._items[index].set_complete(complete)

    def clear(self):
        for item in self._items:
            self._steps_layout.removeWidget(item)
            item.deleteLater()
        self._items.clear()

    def _on_item_clicked(self, index: int):
        self.set_active(index)
        self.step_clicked.emit(index)

    def refresh_theme(self):
        self._header.style().unpolish(self._header)
        self._header.style().polish(self._header)
        for item in self._items:
            item.update()


class StepStatusPanel(QWidget):
    """Fixed-height bottom-of-sidebar panel: title + status line + progress bar + progress text."""

    def __init__(self, parent=None):
        super().__init__(parent)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 8, 10, 10)
        layout.setSpacing(6)

        self._title_label = QLabel("STATUS")
        self._title_label.setProperty("role", "section_header")
        _title_font = self._title_label.font()
        _title_font.setLetterSpacing(QFont.SpacingType.AbsoluteSpacing, 1.5)
        self._title_label.setFont(_title_font)
        self._title_label.setFixedHeight(16)
        layout.addWidget(self._title_label)

        self._status_label = QLabel("Open a project folder to begin")
        self._status_label.setProperty("role", "status_warning")
        self._status_label.setWordWrap(True)
        self._status_label.setFixedHeight(36)
        layout.addWidget(self._status_label)

        self.progress_bar = QProgressBar()
        self.progress_bar.setFixedHeight(6)
        self.progress_bar.setTextVisible(False)
        self.progress_bar.setValue(0)
        layout.addWidget(self.progress_bar)

        self.progress_label = QLabel("")
        self.progress_label.setProperty("role", "caption")
        self.progress_label.setWordWrap(True)
        self.progress_label.setFixedHeight(14)
        layout.addWidget(self.progress_label)

        self.setFixedHeight(120)

    def update_status(self, text: str, state: str = "info"):
        """Update status text. 'state': 'info' (default) | 'success' | 'warning' | 'error'."""
        self._status_label.setText(text)
        self._status_label.setProperty(
            "role", f"status_{state}" if state else "status_info")
        self._status_label.style().unpolish(self._status_label)
        self._status_label.style().polish(self._status_label)

    def refresh_theme(self):
        for w in (self._title_label, self._status_label, self.progress_label):
            w.style().unpolish(w)
            w.style().polish(w)
