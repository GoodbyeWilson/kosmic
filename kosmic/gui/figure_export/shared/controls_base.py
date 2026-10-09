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

    def add_lfc_gate(self, effect: str = "log2 fold change"):
        """Add 'Min |log2FC|' as 'self.lfc_gate'. Off (0) means a gene is
        significant on FDR < 0.05 alone; a value also requires the absolute
        *effect* to exceed it."""
        from PyQt6.QtWidgets import QDoubleSpinBox
        self.lfc_gate = QDoubleSpinBox()
        self.lfc_gate.setRange(0.0, 5.0)
        self.lfc_gate.setSingleStep(0.05)
        self.lfc_gate.setDecimals(2)
        self.lfc_gate.setSpecialValueText("Off")
        self.lfc_gate.setToolTip(
            "Off: significant means FDR < 0.05. A value also requires\n"
            f"|{effect}| above it.")
        self.lfc_gate.valueChanged.connect(self.changed)
        self.add_row("Min |log2FC|:", self.lfc_gate)
        return self.lfc_gate
