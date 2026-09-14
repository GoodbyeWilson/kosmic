# SidebarTabbedPage -- archetype-4 base: sidebar + tabbed content.
#
# Extends 'SidebarPage' by installing a 'QTabWidget' on the content
# side. Subclasses populate the sidebar and call 'add_tab' to
# register tabs.
#
# Subclass usage:
#
#     class GeneDEPage(SidebarTabbedPage):
#         def __init__(self, workspace):
#             super().__init__()
#             self.workspace = workspace
#
#             # Sidebar
#             settings = SettingsGroup("Settings", collapsible=True)
#             settings.add_row("DE method:", self._method_combo)
#             self.sidebar_layout.addWidget(settings)
#
#             # Tabs
#             self.add_tab(self._volcano_widget, "Volcano")
#             self.add_tab(self._results_table, "Results")
#
# The 'QTabWidget' lives at 'self.tabs' for direct access when
# needed (e.g. 'self.tabs.setCurrentIndex(0)' or connecting
# 'currentChanged').
from __future__ import annotations

from typing import Optional

from PyQt6.QtWidgets import QTabWidget, QWidget

from kosmic.gui.shared.widgets.sidebar_page import SidebarPage


class SidebarTabbedPage(SidebarPage):
    """Archetype-4 base: sidebar on the left, 'QTabWidget' on the right."""

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        # Tabs replace the bare content widget. The base's
        # 'content_layout' already exists with CONTENT_MARGINS so
        # the tab widget gets the same left gap from the splitter.
        self.tabs = QTabWidget()
        self.content_layout.addWidget(self.tabs)

    def add_tab(self, widget: QWidget, label: str) -> int:
        """Add a tab and return its index."""
        return self.tabs.addTab(widget, label)
