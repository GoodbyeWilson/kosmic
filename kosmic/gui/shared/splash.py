# Startup splash screen.
#
# Shown while 'AppWindow' builds every workspace; 'set_status' reports
# progress and 'dismiss_when_ready' closes it once the window can show,
# after a minimum display time so it does not flash.
from __future__ import annotations

import math
import os
import time
from typing import Callable

from PyQt6.QtCore import Qt, QTimer, QRectF, QPointF
from PyQt6.QtGui import QPainter, QColor, QFont, QLinearGradient, QPixmap
from PyQt6.QtWidgets import QApplication, QWidget


def _heart_pixmap(height: int = 250) -> QPixmap:
    """Load the bundled heart-UMAP render (assets/heart_umap.png).

    Returns a null pixmap when the asset is missing; the splash then
    falls back to its text-only layout.
    """
    # Beside this module in a source checkout and in a packaged build alike.
    base = os.path.dirname(os.path.abspath(__file__))
    pm = QPixmap(os.path.join(base, "assets", "heart_umap.png"))
    if not pm.isNull():
        pm = pm.scaledToHeight(
            height, Qt.TransformationMode.SmoothTransformation)
    return pm


class SplashScreen(QWidget):
    """Splash: heart-UMAP mark, title, subtitle, status, pulsing dots.

    Owns its own minimum-display window via 'dismiss_when_ready'.
    """

    # Minimum on-screen time so the particle animation gets a beat.
    MIN_DISPLAY_SECONDS = 2.5

    def __init__(self, parent=None):
        # Window flags must be set on the constructor so the native window
        # is created frameless from the start; setting them after the
        # widget is shown briefly flashes a default-framed window.
        super().__init__(
            parent,
            Qt.WindowType.FramelessWindowHint
            | Qt.WindowType.WindowStaysOnTopHint,
        )
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, False)
        self._heart = _heart_pixmap()
        self.setFixedSize(460, 470 if not self._heart.isNull() else 260)
        self._status_text = "Starting..."
        self._start_time = time.monotonic()
        self._timer = QTimer(self)
        self._timer.timeout.connect(self.update)
        self._timer.start(40)  # 25 fps -- enough for the pulsing dots
        self._center_on_primary_screen()

    def _center_on_primary_screen(self):
        screen = QApplication.primaryScreen()
        if screen is None:
            return
        geo = screen.availableGeometry()
        self.move(
            geo.center().x() - self.width() // 2,
            geo.center().y() - self.height() // 2,
        )

    def set_status(self, text: str):
        self._status_text = text
        self.update()

    def stop_animation(self):
        self._timer.stop()

    def dismiss_when_ready(self, on_done: Callable[[], None]) -> None:
        """Stop the animation, close the splash, and call 'on_done' --
        but not before the splash has been visible for at least
        'MIN_DISPLAY_SECONDS'."""
        elapsed = time.monotonic() - self._start_time
        remaining_ms = max(0, int((self.MIN_DISPLAY_SECONDS - elapsed) * 1000))

        def _finish():
            self.stop_animation()
            self.close()
            on_done()

        if remaining_ms > 0:
            QTimer.singleShot(remaining_ms, _finish)
        else:
            _finish()

    def paintEvent(self, event):
        w, h = self.width(), self.height()
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)

        bg = QLinearGradient(0, 0, 0, h)
        bg.setColorAt(0.0, QColor(15, 23, 42))
        bg.setColorAt(1.0, QColor(30, 30, 46))
        p.fillRect(0, 0, w, h, bg)

        if self._heart.isNull():
            title_y = h * 0.22
        else:
            p.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
            p.drawPixmap((w - self._heart.width()) // 2, 22, self._heart)
            title_y = 22 + self._heart.height() + 12

        p.setPen(QColor(255, 255, 255))
        p.setFont(QFont("Segoe UI", 30, QFont.Weight.Bold))
        p.drawText(QRectF(0, title_y, w, 46),
                   Qt.AlignmentFlag.AlignCenter, "KOSMIC")

        p.setPen(QColor(148, 163, 184))
        p.setFont(QFont("Segoe UI", 10))
        p.drawText(QRectF(0, title_y + 50, w, 22),
                   Qt.AlignmentFlag.AlignCenter,
                   "Open-source Single-cell Meta-analysis")

        p.drawText(QRectF(0, title_y + 88, w, 22),
                   Qt.AlignmentFlag.AlignCenter, self._status_text)

        elapsed = time.monotonic() - self._start_time
        for i in range(3):
            phase = (elapsed * 3 - i * 0.4) % 1.0
            alpha = int(60 + 160 * (0.5 + 0.5 * math.sin(phase * math.pi * 2)))
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor(0, 122, 204, alpha))
            p.drawEllipse(QPointF(w / 2 - 16 + i * 16, h - 28), 3, 3)

        p.end()
