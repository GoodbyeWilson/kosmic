# Shared pyqtgraph base: themed plot widget with title, axis labels,
# and an empty-state overlay driven by class attributes.
from __future__ import annotations

from typing import Optional

import pyqtgraph as pg
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import QWidget

from kosmic.gui.shared.theme import (
    get_color, get_font_sizes, style_pg_plot,
)


def make_tooltip_label() -> pg.TextItem:
    """Hidden 'pg.TextItem' for hover tooltips on plot widgets."""
    label = pg.TextItem(
        text='', color=get_color('fg_primary'), anchor=(0, 1))
    label.setFont(QFont('Arial', int(get_font_sizes()['tick'])))
    label.hide()
    return label


class InteractivePlot(pg.PlotWidget):
    """Themed 'pg.PlotWidget' with empty-state overlay + reset hook."""

    # -- subclass overrides --------------------------------------------
    TITLE: str = ""
    LEFT_LABEL: str = ""
    BOTTOM_LABEL: str = ""
    UNAVAILABLE_MESSAGE: str = "No data."
    LEFT_AXIS_WIDTH: int = 50

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        *,
        title: Optional[str] = None,
        left_label: Optional[str] = None,
        bottom_label: Optional[str] = None,
        unavailable_message: Optional[str] = None,
    ):
        """Create a themed plot widget.

        Subclasses normally set the class-level attributes; the keyword
        arguments here are an escape hatch for cases where a plot is
        instantiated in a loop with per-instance text (e.g. one panel
        per pathway).
        """
        super().__init__(parent)

        self._title = title if title is not None else self.TITLE
        self._left_label = (
            left_label if left_label is not None else self.LEFT_LABEL)
        self._bottom_label = (
            bottom_label if bottom_label is not None else self.BOTTOM_LABEL)
        self._unavailable_message = (
            unavailable_message if unavailable_message is not None
            else self.UNAVAILABLE_MESSAGE)

        style_pg_plot(
            self,
            title=self._title,
            left_label=self._left_label,
            bottom_label=self._bottom_label,
        )

        # Standard left-axis width for cross-page alignment.
        self.getAxis('left').setWidth(self.LEFT_AXIS_WIDTH)

        # Empty-state overlay: a centred TextItem in view coordinates.
        # Anchored at (0.5, 0.5) so it stays centred when the view
        # auto-ranges to (0..1, 0..1) -- which it does when no real
        # data has been plotted yet.
        self._unavailable_item: Optional[pg.TextItem] = None
        self._show_initial_unavailable_message()

    # -- empty-state overlay -------------------------------------------

    def _show_initial_unavailable_message(self) -> None:
        if not self._unavailable_message:
            return
        msg = pg.TextItem(
            self._unavailable_message,
            color=get_color('fg_tertiary'),
            anchor=(0.5, 0.5),
        )
        font = QFont()
        font.setPointSizeF(get_font_sizes()['annotation'])
        msg.setFont(font)
        self.addItem(msg, ignoreBounds=True)
        msg.setPos(0.5, 0.5)
        self.setXRange(0, 1, padding=0)
        self.setYRange(0, 1, padding=0)
        self._unavailable_item = msg

    def show_unavailable_message(
        self, message: Optional[str] = None,
    ) -> None:
        """Render (or re-render) the empty-state overlay.

        Pass 'message' to override the configured text for this
        invocation (e.g. '"No genes match the current filter."').
        Resets the view to (0, 1) on both axes so the centred TextItem
        is on-screen regardless of any prior auto-range.
        """
        if message is not None:
            self._unavailable_message = message
        if self._unavailable_item is None:
            self._show_initial_unavailable_message()
        else:
            self._unavailable_item.setText(self._unavailable_message)
            self._unavailable_item.setVisible(True)
            self.setXRange(0, 1, padding=0)
            self.setYRange(0, 1, padding=0)

    def hide_unavailable_message(self) -> None:
        """Hide the empty-state overlay (call once real data is rendered)."""
        if self._unavailable_item is not None:
            self._unavailable_item.setVisible(False)
        self.enableAutoRange()

    # -- data lifecycle helpers ----------------------------------------

    def clear_plot_items(self) -> None:
        """Remove every item except the empty-state overlay."""
        items = list(self.getPlotItem().items)
        for item in items:
            if item is self._unavailable_item:
                continue
            self.removeItem(item)

    def reset_view(self) -> None:
        """Reset zoom/pan to fit all data. Wire to an external button if needed."""
        self.autoRange()

    # -- title / label updates -----------------------------------------
