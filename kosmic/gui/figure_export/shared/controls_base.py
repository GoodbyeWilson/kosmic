# Base class for per-figure control panels.
#
# Subclasses add their own widgets in '__init__' and connect them to
# 'self.changed'. The 'FigsizeRow' is added by the base class.
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QFormLayout, QWidget

from kosmic.gui.figure_export.shared.figsize_row import FigsizeRow
from kosmic.gui.shared.widgets import SecondaryLabel


class FigureControls(QWidget):
    """
    A QFormLayout with a FigsizeRow at the top.

    Subclasses call 'self._form.addRow(label, widget)' to add controls and
    connect the widget's value-change signal to 'self.changed'.
    """

    changed = pyqtSignal()

    def __init__(self, default_w: int = 0, default_h: int = 0,
                 description: str = ""):
        super().__init__()
        self._form = QFormLayout(self)
        self._form.setContentsMargins(0, 0, 0, 0)

        self.figsize = FigsizeRow(default_w, default_h)
        self.figsize.changed.connect(self.changed)
        self._form.addRow("Size:", self.figsize)

        if description:
            info = SecondaryLabel(description)
            info.setWordWrap(True)
            self._form.addRow(info)

    def add_row(self, label: str, widget: QWidget):
        """Add a labelled control row. Caller still wires its own signals."""
        self._form.addRow(label, widget)
