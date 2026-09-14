# SidebarPage -- archetype-3 base: sidebar + content body.
#
# Replaces the hand-rolled 'QSplitter' + 'QScrollArea' + max-width
# sidebar + content widget pattern that grew up across ~10 pages with
# drifting widths (280/320/340), inconsistent scroll wrapping, and
# non-uniform splitter settings.
#
# Subclass usage:
#
#     class FilterTab(SidebarPage):
#         def __init__(self, workspace):
#             super().__init__()
#             self.workspace = workspace
#
#             # Populate the sidebar (a QVBoxLayout)
#             cell_types = SettingsGroup("Cell Types", collapsible=True)
#             cell_types.add_widget(...)
#             self.sidebar_layout.addWidget(cell_types)
#             self.sidebar_layout.addStretch()
#
#             # Install whatever layout / widgets you want on self.content
#             content_lay = QVBoxLayout(self.content)
#             content_lay.addWidget(self._results_table)
#
# Class-level overrides (set as class attributes on the subclass):
#
# * 'SIDEBAR_WIDTH'        -- max sidebar width (default 320)
# * 'SIDEBAR_MIN_WIDTH'    -- min sidebar width (default 280; None to skip)
# * 'SIDEBAR_MARGINS'      -- (l, t, r, b) margins inside sidebar
# * 'SIDEBAR_SPACING'      -- spacing between sidebar children
# * 'SIDEBAR_SCROLLABLE'   -- wrap sidebar in QScrollArea (default True)
# * 'INITIAL_SIZES'        -- (sidebar_px, content_px) initial splitter
#                               sizes. Default starts the sidebar at its
#                               minimum so the content gets the
#                               remaining space; the user can drag wider
#                               up to 'SIDEBAR_WIDTH' if the sidebar
#                               has any drag-room
# * 'SIDEBAR_STRETCH'      -- splitter stretch factor for sidebar (default 0)
# * 'CONTENT_STRETCH'      -- splitter stretch factor for content (default 1)
# * 'HANDLE_WIDTH'         -- splitter handle width in px. Default 4
#                               gives Qt enough allocation to draw the
#                               theme's 1px coloured divider line even
#                               when both sides are width-constrained
#                               (Qt collapses the handle to 0 if children
#                               have 'min == max' and no explicit
#                               handle width is set)
# * 'CONTENT_MARGINS'      -- (l, t, r, b) margins on the content
#                               area. Default (8, 0, 0, 0) puts a small
#                               gap between the splitter divider and
#                               the content's first visible edge so the
#                               divider doesn't visually merge with a
#                               'QTabWidget' pane border or table edge.
#
# Subclass content layout: SidebarPage installs 'self.content_layout'
# (a 'QVBoxLayout' with 'CONTENT_MARGINS'). Subclasses add to it
# instead of installing their own layout on 'self.content'.
from __future__ import annotations

from typing import Optional, Tuple

from PyQt6.QtCore import Qt
from kosmic import UI_FONT_SCALE

from PyQt6.QtWidgets import (
    QFrame, QScrollArea, QSplitter, QVBoxLayout, QWidget,
)


class SidebarPage(QWidget):
    """Archetype-3 base: left sidebar + right content body."""

    #: Help-content path (relative to 'kosmic/gui/help/content/',
    #: without the '.md' extension). Subclasses override.
    help_id: str = "index"

    # Scaled with the UI font: a sidebar holds a fixed amount of *text*,
    # not a fixed number of pixels, so raising the font size without this
    # leaves settings panels demanding more width than the panel has and
    # their controls clipped at the card edge.
    SIDEBAR_WIDTH: int = round(320 * UI_FONT_SCALE)
    SIDEBAR_MIN_WIDTH: Optional[int] = round(320 * UI_FONT_SCALE)
    SIDEBAR_MARGINS: Tuple[int, int, int, int] = (10, 10, 10, 10)
    SIDEBAR_SPACING: int = 8
    SIDEBAR_SCROLLABLE: bool = True
    INITIAL_SIZES: Tuple[int, int] = (round(320 * UI_FONT_SCALE), 800)
    SIDEBAR_STRETCH: int = 0
    CONTENT_STRETCH: int = 1
    HANDLE_WIDTH: int = 4
    CONTENT_MARGINS: Tuple[int, int, int, int] = (8, 0, 0, 0)

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._build_page_layout()

    # -- internals ------------------------------------------------------

    def _build_page_layout(self) -> None:
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        self._splitter = QSplitter(Qt.Orientation.Horizontal)
        self._splitter.setChildrenCollapsible(False)
        self._splitter.setHandleWidth(self.HANDLE_WIDTH)

        # -- sidebar widget + layout (always a QVBoxLayout) --
        self.sidebar = QWidget()
        self.sidebar.setMaximumWidth(self.SIDEBAR_WIDTH)
        if self.SIDEBAR_MIN_WIDTH is not None:
            self.sidebar.setMinimumWidth(self.SIDEBAR_MIN_WIDTH)
        self.sidebar_layout = QVBoxLayout(self.sidebar)
        self.sidebar_layout.setContentsMargins(*self.SIDEBAR_MARGINS)
        self.sidebar_layout.setSpacing(self.SIDEBAR_SPACING)

        # Optional scroll-area wrapper around the sidebar.
        if self.SIDEBAR_SCROLLABLE:
            self._sidebar_scroll = QScrollArea()
            self._sidebar_scroll.setWidgetResizable(True)
            self._sidebar_scroll.setFrameShape(QFrame.Shape.NoFrame)
            self._sidebar_scroll.setHorizontalScrollBarPolicy(
                Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            self._sidebar_scroll.setMaximumWidth(self.SIDEBAR_WIDTH)
            if self.SIDEBAR_MIN_WIDTH is not None:
                self._sidebar_scroll.setMinimumWidth(self.SIDEBAR_MIN_WIDTH)
            # Theme hook: theme.py styles 'QScrollArea[role="thin_scroll"]'
            # with a slim scrollbar; sidebars universally use this look.
            self._sidebar_scroll.setProperty("role", "thin_scroll")
            self._sidebar_scroll.setWidget(self.sidebar)
            sidebar_container: QWidget = self._sidebar_scroll
        else:
            self._sidebar_scroll = None
            sidebar_container = self.sidebar

        # -- content widget + default vertical layout with margins --
        # Subclasses add to 'self.content_layout' rather than installing
        # their own layout on 'self.content'. The CONTENT_MARGINS left
        # gap exists so the splitter divider stays visually distinct
        # from any QTabWidget pane border or table edge inside.
        self.content = QWidget()
        self.content_layout = QVBoxLayout(self.content)
        self.content_layout.setContentsMargins(*self.CONTENT_MARGINS)
        self.content_layout.setSpacing(0)

        # -- assemble splitter --
        self._splitter.addWidget(sidebar_container)
        self._splitter.addWidget(self.content)
        self._splitter.setSizes(list(self.INITIAL_SIZES))
        self._splitter.setStretchFactor(0, self.SIDEBAR_STRETCH)
        self._splitter.setStretchFactor(1, self.CONTENT_STRETCH)

        outer.addWidget(self._splitter)

    # -- public helpers -------------------------------------------------

    @property
    def splitter(self) -> QSplitter:
        """The main horizontal splitter (sidebar | content)."""
        return self._splitter
