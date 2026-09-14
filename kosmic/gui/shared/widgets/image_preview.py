# ImagePreviewWidget -- a QLabel that loads an image from disk and
# auto-scales it to the widget size while keeping aspect ratio.
from __future__ import annotations


from PyQt6.QtCore import Qt
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import QLabel


class ImagePreviewWidget(QLabel):
    """Widget for displaying plot/image previews loaded from disk."""

    def __init__(self):
        super().__init__()
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setMinimumSize(400, 400)
        self.setText("No plot generated yet")
        self._current_pixmap: QPixmap | None = None

    def _update_scaled_pixmap(self):
        if self._current_pixmap:
            scaled = self._current_pixmap.scaled(
                self.size(),
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation,
            )
            self.setPixmap(scaled)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._update_scaled_pixmap()
