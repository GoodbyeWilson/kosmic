# SimplePage -- archetype-1 base: vertical body, no sidebar, no tabs.
#
# For pages that are just a column of widgets (set-directory, load-data,
# methods overview, pathway explorer, etc.). Optional scroll wrapper for
# pages that may overflow vertically.
#
# Subclass usage:
#
#     class SetupPage(SimplePage):
#         SCROLLABLE = True
#         CONTENT_MARGINS = (40, 30, 40, 30)
#
#         def __init__(self, workspace):
#             super().__init__()
#             self.body_layout.addWidget(QLabel("Set project directory"))
#             self.body_layout.addWidget(self._dir_picker)
#             self.body_layout.addStretch()
#
# The body widget lives at 'self.body' and its 'QVBoxLayout' at
# 'self.body_layout'.
#
# Class-level overrides:
#
# * 'SCROLLABLE'           -- wrap body in a 'QScrollArea' (default False)
# * 'CONTENT_MARGINS'      -- (l, t, r, b) margins inside the body
# * 'CONTENT_SPACING'      -- spacing between body children
from __future__ import annotations

from typing import Optional, Tuple

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QFrame, QScrollArea, QVBoxLayout, QWidget,
)


class SimplePage(QWidget):
    """Archetype-1 base: vertical body, no sidebar, no tabs."""

    #: Help-content path (relative to 'kosmic/gui/help/content/',
    #: without the '.md' extension). Subclasses override.
    help_id: str = "index"

    SCROLLABLE: bool = False
    CONTENT_MARGINS: Tuple[int, int, int, int] = (16, 16, 16, 16)
    CONTENT_SPACING: int = 12

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._build_page_layout()

    def _build_page_layout(self) -> None:
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # Body widget + layout. Subclasses add widgets to body_layout.
        self.body = QWidget()
        self.body_layout = QVBoxLayout(self.body)
        self.body_layout.setContentsMargins(*self.CONTENT_MARGINS)
        self.body_layout.setSpacing(self.CONTENT_SPACING)

        if self.SCROLLABLE:
            self._scroll: Optional[QScrollArea] = QScrollArea()
            self._scroll.setWidgetResizable(True)
            self._scroll.setFrameShape(QFrame.Shape.NoFrame)
            self._scroll.setHorizontalScrollBarPolicy(
                Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            self._scroll.setWidget(self.body)
            outer.addWidget(self._scroll)
        else:
            self._scroll = None
            outer.addWidget(self.body)
