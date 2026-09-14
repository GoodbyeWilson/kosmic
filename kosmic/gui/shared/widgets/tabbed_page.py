# TabbedPage -- archetype-2 base: no sidebar, owns a 'QTabWidget'.
#
# For pages whose entire body is a tab widget (download/convert tabs,
# results-only pages). The tab widget lives at 'self.tabs'; subclasses
# call 'add_tab' to register tabs.
#
# Subclass usage:
#
#     class DownloadTab(TabbedPage):
#         def __init__(self, main_window):
#             super().__init__()
#             self.add_tab(self._build_geo_tab(),    "Download from GEO")
#             self.add_tab(self._build_local_tab(),  "Local Files")
#             self.add_tab(self._build_assemble_tab(), "Assemble from Parts")
#
# Class-level overrides:
#
# * 'CONTENT_MARGINS'      -- (l, t, r, b) margins around the tab widget
from __future__ import annotations

from typing import Optional, Tuple

from PyQt6.QtWidgets import QTabWidget, QVBoxLayout, QWidget


class TabbedPage(QWidget):
    """Archetype-2 base: no sidebar, 'QTabWidget' fills the body."""

    #: Help-content path (relative to 'kosmic/gui/help/content/',
    #: without the '.md' extension). Subclasses override.
    help_id: str = "index"

    CONTENT_MARGINS: Tuple[int, int, int, int] = (0, 0, 0, 0)

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._build_page_layout()

    def _build_page_layout(self) -> None:
        outer = QVBoxLayout(self)
        outer.setContentsMargins(*self.CONTENT_MARGINS)
        outer.setSpacing(0)

        self.tabs = QTabWidget()
        outer.addWidget(self.tabs)

    def add_tab(self, widget: QWidget, label: str) -> int:
        """Add a tab and return its index."""
        return self.tabs.addTab(widget, label)
