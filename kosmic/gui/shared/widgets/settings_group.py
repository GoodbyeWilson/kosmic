# SettingsGroup -- a labelled cluster of label/widget rows.
#
# Uses 'QGroupBox' internally so the rendered look matches Qt Fusion
# defaults the rest of the app already uses (titled border, themed title
# text). Replaces the same hand-rolled 'QGroupBox' + 'QVBoxLayout' +
# margin/spacing boilerplate that grew across ~65 sites.
#
# Optional collapsibility ('collapsible=True') maps to
# 'QGroupBox.setCheckable(True)': the title bar gets a checkbox that
# toggles content visibility, which is the same pattern cluster_tab
# already uses for its PCA / Clustering / Sweep groups.
#
# Usage::
#
#     group = SettingsGroup("PCA")
#     group.add_row("Components:", components_spinbox)
#     group.add_row("Solver:", solver_combo)
#     group.add_widget(run_pca_button)        # full-width, no label
#     parent_layout.addWidget(group)
#
#     # Collapsible (header checkbox toggles content visibility)
#     pca = SettingsGroup("PCA", collapsible=True)
#     pca.add_row(...)
from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QFormLayout, QGroupBox, QVBoxLayout, QWidget,
)


class SettingsGroup(QGroupBox):
    """Titled cluster of label/widget rows (QGroupBox-based).

    Subclasses 'QGroupBox' so it inherits the themed title, border,
    and (when 'collapsible=True') the built-in checkable title that
    toggles content visibility.

    Rows added via 'add_row' get a 'QFormLayout'-style label on
    the left and the widget on the right. Rows added via
    'add_widget' / 'add_layout' span the full width.

    Parameters
    ----------
    title
        Section title rendered in the QGroupBox title bar.
    collapsible
        If True, the title bar shows a checkbox; checking/unchecking
        toggles the content. Default False.
    expanded
        Initial expanded state when 'collapsible=True'. Default True.
    parent
        Optional parent widget.
    """

    toggled_expanded = pyqtSignal(bool)  # collapsible groups only

    def __init__(self, title: str, *,
                 collapsible: bool = False,
                 expanded: bool = True,
                 parent: Optional[QWidget] = None):
        super().__init__(title, parent)
        self._collapsible = collapsible
        self._build_ui(collapsible, expanded)

    def _build_ui(self, collapsible: bool, expanded: bool) -> None:
        # Content holder so collapse/expand can toggle one widget.
        self._content = QWidget()
        self._form = QFormLayout(self._content)
        self._form.setContentsMargins(0, 0, 0, 0)
        self._form.setSpacing(6)
        # Let field widgets (combos, spinboxes) fill the row width
        # instead of clinging to their sizeHint (Windows default).
        self._form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)

        # Outer layout on the QGroupBox itself; matches the cluster_tab
        # PCA-group margins so migrated groups visually replace 1:1.
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.setSpacing(0)
        outer.addWidget(self._content)

        if collapsible:
            self.setCheckable(True)
            self.setChecked(expanded)
            self._content.setVisible(expanded)
            self.toggled.connect(self._on_toggled)

    # -- public API -----------------------------------------------------

    def add_row(self, label, widget: QWidget) -> None:
        """Add a label-on-left, widget-on-right row.

        'label' may be either a string or a QWidget (e.g. a QCheckBox
        used as the row's left-hand selector).
        """
        self._form.addRow(label, widget)

    def add_widget(self, widget: QWidget) -> None:
        """Add a full-width widget (no label column)."""
        self._form.addRow(widget)

    def add_layout(self, layout) -> None:
        """Add a full-width layout (no label column)."""
        self._form.addRow(layout)

    # -- collapsible behaviour -----------------------------------------

    def is_expanded(self) -> bool:
        if not self._collapsible:
            return True
        return self.isChecked()

    def set_expanded(self, expanded: bool) -> None:
        if not self._collapsible:
            return
        self.setChecked(expanded)  # triggers toggled -> _on_toggled

    def _on_toggled(self, checked: bool):
        self._content.setVisible(checked)
        self.toggled_expanded.emit(checked)
