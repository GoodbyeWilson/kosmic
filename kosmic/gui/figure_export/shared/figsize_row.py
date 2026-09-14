# Width x height widget with an Auto checkbox. Emits 'changed' on any edit.
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QCheckBox, QHBoxLayout, QLabel, QSpinBox, QWidget


class FigsizeRow(QWidget):
    """
    Inline width / height spinboxes (inches) with an Auto checkbox.

    'get_figsize()' returns 'None' while Auto is checked, otherwise
    '(width, height)'.
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

        self.w_spin = QSpinBox()
        self.w_spin.setRange(4, 30)
        self.w_spin.setValue(default_w if default_w > 0 else 10)
        self.w_spin.setSuffix('"')
        self.w_spin.valueChanged.connect(self.changed)
        lay.addWidget(self.w_spin)

        lay.addWidget(QLabel("x"))

        self.h_spin = QSpinBox()
        self.h_spin.setRange(4, 30)
        self.h_spin.setValue(default_h if default_h > 0 else 8)
        self.h_spin.setSuffix('"')
        self.h_spin.valueChanged.connect(self.changed)
        lay.addWidget(self.h_spin)

        self._on_auto_toggled(self.auto_cb.isChecked())

    def _on_auto_toggled(self, auto: bool):
        self.w_spin.setEnabled(not auto)
        self.h_spin.setEnabled(not auto)

    def get_figsize(self):
        if self.auto_cb.isChecked():
            return None
        return (self.w_spin.value(), self.h_spin.value())
