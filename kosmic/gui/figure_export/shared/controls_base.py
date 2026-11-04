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
        # The sidebar has a fixed width: a field that does not fit beside
        # the widest label goes below its own label instead of being
        # squeezed below its minimum size.
        self._form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)

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

    def add_x_limit(self, unit: str):
        """Add 'X limit (+/- *unit*)' as 'self.x_limit': the axis runs from
        minus to plus the value. Auto (0) leaves the range to the figure."""
        from PyQt6.QtWidgets import QDoubleSpinBox
        self.x_limit = QDoubleSpinBox()
        self.x_limit.setRange(0.0, 50.0)
        self.x_limit.setSingleStep(0.5)
        self.x_limit.setDecimals(1)
        self.x_limit.setValue(0.0)
        self.x_limit.setSpecialValueText("Auto")
        self.x_limit.valueChanged.connect(self.changed)
        self.add_row(f"X limit (+/- {unit}):", self.x_limit)
        return self.x_limit

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

    def add_named_labels(self):
        """Add 'Label named genes only' as 'self.named_labels', checked: the
        automatic gene labels pass over genes with no descriptive symbol."""
        from PyQt6.QtWidgets import QCheckBox
        self.named_labels = QCheckBox("Label named genes only")
        self.named_labels.setChecked(True)
        self.named_labels.setToolTip(
            "Skip genes with no descriptive symbol when choosing labels:\n"
            "Ensembl IDs, clone names such as AC004936.1, LINC and antisense\n"
            "(-AS1) transcripts. The genes are still plotted.")
        self.named_labels.toggled.connect(self.changed)
        self.add_row("", self.named_labels)
        return self.named_labels
