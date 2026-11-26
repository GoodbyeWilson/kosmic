# Width x height widget with an Auto checkbox. Emits 'changed' on any edit.
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QCheckBox, QDoubleSpinBox, QHBoxLayout, QLabel, QWidget


def _size_box(value: float) -> QDoubleSpinBox:
    box = QDoubleSpinBox()
    box.setRange(1, 30)
    box.setDecimals(1)
    box.setSingleStep(0.5)
    box.setValue(value)
    box.setSuffix('"')
    return box


class FigsizeRow(QWidget):
    """
    Inline width / height boxes (inches) with an Auto checkbox.

    'get_figsize()' returns 'None' while Auto is checked, so the figure
    takes the size its plot function chooses; otherwise '(width, height)'.
    While Auto is checked the boxes show the size the last preview was
    drawn at ('show_auto_size'), so unticking Auto keeps the figure as it
    is until a size is changed.
    """

    changed = pyqtSignal()

    def __init__(self, default_w: int = 0, default_h: int = 0):
        super().__init__()
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)

        self.auto_cb = QCheckBox("Auto")
        self.auto_cb.setChecked(default_w == 0 or default_h == 0)
        self.auto_cb.toggled.connect(self._on_auto_toggled)
        self.auto_cb.toggled.connect(self.changed)
        lay.addWidget(self.auto_cb)

        self.w_spin = _size_box(default_w if default_w > 0 else 10)
        self.w_spin.valueChanged.connect(self.changed)
        lay.addWidget(self.w_spin)

        lay.addWidget(QLabel("x"))

        self.h_spin = _size_box(default_h if default_h > 0 else 8)
        self.h_spin.valueChanged.connect(self.changed)
        lay.addWidget(self.h_spin)

        self._on_auto_toggled(self.auto_cb.isChecked())

    def _on_auto_toggled(self, auto: bool):
        self.w_spin.setEnabled(not auto)
        self.h_spin.setEnabled(not auto)

    def show_auto_size(self, size) -> None:
        """Show *size* ((width, height) in inches, or None) in the boxes
        while Auto is checked, without emitting 'changed'."""
        if size is None or not self.auto_cb.isChecked():
            return
        for box, value in zip((self.w_spin, self.h_spin), size):
            box.blockSignals(True)
            box.setValue(round(float(value), 1))
            box.blockSignals(False)

    def get_figsize(self):
        if self.auto_cb.isChecked():
            return None
        return (self.w_spin.value(), self.h_spin.value())
