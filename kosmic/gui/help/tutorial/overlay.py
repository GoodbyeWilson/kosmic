# Semi-transparent overlay with a spotlight cutout + floating instruction card.
from typing import Optional

from PyQt6.QtWidgets import (
    QWidget, QFrame, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QSizePolicy,
)
from PyQt6.QtCore import Qt, pyqtSignal, QRect, QPoint
from PyQt6.QtGui import QPainter, QColor, QRegion, QPen, QBrush, QPainterPath

from kosmic.gui.shared.theme import get_color, get_current_mode, ThemeMode


# Instruction card

class TutorialCard(QFrame):
    """Floating instruction card shown alongside the spotlight."""

    next_clicked = pyqtSignal()
    skip_clicked = pyqtSignal()
    auto_action_clicked = pyqtSignal()

    CARD_WIDTH = 380

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("tutorial_card")
        self.setFixedWidth(self.CARD_WIDTH)
        self._build_ui()
        self.hide()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 14, 16, 14)
        layout.setSpacing(8)

        # Title
        self._title = QLabel()
        self._title.setObjectName("tutorial_title")
        self._title.setWordWrap(True)
        self._title.setFixedWidth(self.CARD_WIDTH - 32)
        layout.addWidget(self._title)

        # Body
        self._body = QLabel()
        self._body.setObjectName("tutorial_body")
        self._body.setWordWrap(True)
        self._body.setFixedWidth(self.CARD_WIDTH - 32)  # account for margins
        self._body.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Minimum)
        layout.addWidget(self._body)

        layout.addSpacing(4)

        # Step counter
        self._counter = QLabel()
        self._counter.setObjectName("tutorial_step_counter")
        layout.addWidget(self._counter)

        # Buttons row
        btn_layout = QHBoxLayout()
        btn_layout.setSpacing(8)

        self._auto_btn = QPushButton("Do It For Me")
        self._auto_btn.setObjectName("tutorial_btn_auto")
        self._auto_btn.clicked.connect(self.auto_action_clicked.emit)
        self._auto_btn.hide()
        btn_layout.addWidget(self._auto_btn)

        btn_layout.addStretch()

        self._skip_btn = QPushButton("Skip Tutorial")
        self._skip_btn.setObjectName("tutorial_btn")
        self._skip_btn.clicked.connect(self.skip_clicked.emit)
        btn_layout.addWidget(self._skip_btn)

        self._next_btn = QPushButton("Next")
        self._next_btn.setObjectName("tutorial_btn_primary")
        self._next_btn.clicked.connect(self.next_clicked.emit)
        btn_layout.addWidget(self._next_btn)

        layout.addLayout(btn_layout)

    # -- Public API --

    def set_step(self, title: str, body: str, step: int, total: int,
                 next_text: str = "Next", has_auto_action: bool = False,
                 allow_next: bool = True):
        """Update the card contents for a new tutorial step."""
        self._title.setText(title)
        self._body.setText(body)
        self._counter.setText(f"Step {step} of {total}")
        self._next_btn.setText(next_text)
        self._next_btn.setEnabled(allow_next)
        self._auto_btn.setVisible(has_auto_action)
        # Force layout recalculation — word-wrapped labels need explicit sizing
        self._body.adjustSize()
        self._title.adjustSize()
        self.layout().activate()
        self.setFixedHeight(self.layout().sizeHint().height())

    def enable_next(self):
        """Enable the Next button (called after a signal wait completes)."""
        self._next_btn.setEnabled(True)


# Full-window overlay

class TutorialOverlay(QWidget):
    """Semi-transparent overlay that covers the MainWindow, with a
    spotlight cutout around the target widget."""

    SPOTLIGHT_PADDING = 8
    SPOTLIGHT_RADIUS = 8

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents, False)
        self.setWindowFlags(Qt.WindowType.Widget)
        self._spotlight_rect: Optional[QRect] = None
        self._card = TutorialCard(self)
        self.hide()

    @property
    def card(self) -> TutorialCard:
        return self._card

    def set_spotlight(self, widget: Optional[QWidget] = None):
        """Set the widget to spotlight, or None for no spotlight."""
        if widget is None:
            self._spotlight_rect = None
        else:
            # Map widget geometry to overlay coordinates
            top_left = widget.mapTo(self.parentWidget(), QPoint(0, 0))
            rect = QRect(top_left, widget.size())
            rect.adjust(
                -self.SPOTLIGHT_PADDING, -self.SPOTLIGHT_PADDING,
                self.SPOTLIGHT_PADDING, self.SPOTLIGHT_PADDING,
            )
            self._spotlight_rect = rect

        self._update_mask()
        self._position_card()
        self.update()

    def show_step(self, title: str, body: str, step: int, total: int,
                  target_widget=None, next_text: str = "Next",
                  has_auto_action: bool = False, allow_next: bool = True):
        """Configure and show the overlay for one tutorial step."""
        self.set_spotlight(target_widget)
        self._card.set_step(title, body, step, total, next_text, has_auto_action, allow_next)
        # Re-position and update mask after card content change (size may differ)
        self._position_card()
        self._update_mask()
        self._card.show()
        self._card.raise_()
        self.show()
        self.raise_()
        self._card.raise_()

    def hide_overlay(self):
        """Hide overlay and card."""
        self._card.hide()
        self.hide()

    # -- Internals --

    def _update_mask(self):
        """Set the widget mask so clicks pass through the spotlight area.

        The card region is always kept in the mask so its buttons remain
        clickable even when the spotlight cutout is very large.
        """
        if self._spotlight_rect is None:
            # No spotlight — overlay covers everything
            self.clearMask()
            return

        full = QRegion(self.rect())
        cutout = QRegion(self._spotlight_rect, QRegion.RegionType.Rectangle)
        mask = full.subtracted(cutout)

        # Ensure card stays clickable when spotlight overlaps its position
        if self._card.isVisible():
            mask = mask.united(QRegion(self._card.geometry()))

        self.setMask(mask)

    def _position_card(self):
        """Position the card near the spotlight, preferring below."""
        card = self._card
        sr = self._spotlight_rect
        parent_rect = self.rect()

        if sr is None:
            # Centre the card
            x = (parent_rect.width() - card.width()) // 2
            y = (parent_rect.height() - card.height()) // 2
            card.move(x, y)
            return

        margin = 16
        card_h = card.sizeHint().height()
        card_w = card.width()

        # Try below
        if sr.bottom() + margin + card_h < parent_rect.bottom():
            x = max(margin, min(sr.left(), parent_rect.right() - card_w - margin))
            y = sr.bottom() + margin
        # Try above
        elif sr.top() - margin - card_h > parent_rect.top():
            x = max(margin, min(sr.left(), parent_rect.right() - card_w - margin))
            y = sr.top() - margin - card_h
        # Try right
        elif sr.right() + margin + card_w < parent_rect.right():
            x = sr.right() + margin
            y = max(margin, min(sr.top(), parent_rect.bottom() - card_h - margin))
        # Try left
        else:
            x = max(margin, sr.left() - margin - card_w)
            y = max(margin, min(sr.top(), parent_rect.bottom() - card_h - margin))

        card.move(int(x), int(y))

    def paintEvent(self, event):
        """Draw the semi-transparent overlay with spotlight cutout."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        # Semi-transparent fill
        is_dark = get_current_mode() == ThemeMode.DARK
        overlay_color = QColor(0, 0, 0, 160 if is_dark else 120)
        painter.fillRect(self.rect(), overlay_color)

        # Draw spotlight cutout (clear area + accent border)
        if self._spotlight_rect is not None:
            sr = self._spotlight_rect

            # Clear the spotlight area
            painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Clear)
            path = QPainterPath()
            path.addRoundedRect(
                float(sr.x()), float(sr.y()),
                float(sr.width()), float(sr.height()),
                self.SPOTLIGHT_RADIUS, self.SPOTLIGHT_RADIUS,
            )
            painter.fillPath(path, QBrush(Qt.GlobalColor.transparent))

            # Draw accent border around spotlight
            painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceOver)
            accent = QColor(get_color('accent_primary'))
            accent.setAlpha(200)
            pen = QPen(accent, 2.0)
            painter.setPen(pen)
            painter.drawRoundedRect(
                sr, self.SPOTLIGHT_RADIUS, self.SPOTLIGHT_RADIUS,
            )

        painter.end()

    def resizeEvent(self, event):
        """Keep mask and card positioned when window resizes."""
        super().resizeEvent(event)
        self._update_mask()
        self._position_card()
