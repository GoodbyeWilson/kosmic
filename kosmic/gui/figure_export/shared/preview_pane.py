# Right-side figure preview pane: scaled QLabel + status text.
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import QFrame, QLabel, QVBoxLayout, QWidget


class PreviewPane(QWidget):
    """
    A frame containing a centred QLabel that scales pixmaps to fit.

    Use 'set_pixmap(QPixmap)', 'set_png_bytes(bytes)', or
    'set_message(str)'. The pane keeps the most recent pixmap and rescales
    it whenever the widget is resized.
    """

    def __init__(self, placeholder: str = "No figure to display."):
        super().__init__()
        self._placeholder = placeholder
        self._pixmap: QPixmap | None = None

        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(0)

        frame = QFrame()
        frame.setObjectName("preview_frame")
        frame.setFrameShape(QFrame.Shape.StyledPanel)
        frame_lay = QVBoxLayout(frame)
        frame_lay.setContentsMargins(0, 0, 0, 0)

        self._label = QLabel(placeholder)
        self._label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._label.setMinimumSize(400, 300)
        self._label.setSizePolicy(self._label.sizePolicy().Policy.Expanding,
                                  self._label.sizePolicy().Policy.Expanding)
        frame_lay.addWidget(self._label)
        lay.addWidget(frame)

    def set_message(self, text: str):
        self._pixmap = None
        self._label.setPixmap(QPixmap())
        self._label.setText(text)

    def reset(self):
        self.set_message(self._placeholder)

    def set_pixmap(self, pixmap: QPixmap):
        if pixmap is None or pixmap.isNull():
            self.set_message(self._placeholder)
            return
        self._pixmap = pixmap
        self._rescale()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self._pixmap is not None:
            self._rescale()

    def _rescale(self):
        if self._pixmap is None:
            return
        scaled = self._pixmap.scaled(
            self._label.width() - 4,
            self._label.height() - 4,
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        )
        self._label.setText("")
        self._label.setPixmap(scaled)
