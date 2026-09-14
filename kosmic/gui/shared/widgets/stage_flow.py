# Stage-flow primitive for the "summary first, parameters on demand"
# pattern: a sidebar of StageSummaryCards, each recording what its
# analysis stage ran with and expanding in place (accordion-style) to
# show the stage's settings. The owning page coordinates the accordion
# so opening one card closes the others. Styling lives in theme.py
# (objectName 'stage_card'; roles 'stage_card_title', 'link_accent').
from __future__ import annotations

from PyQt6.QtCore import QSize, Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QFrame, QHBoxLayout, QLabel, QPushButton, QScrollArea, QSizePolicy,
    QVBoxLayout, QWidget,
)

from kosmic.gui.shared.icon_provider import make_icon
from kosmic.gui.shared.theme import get_color
from kosmic.gui.shared.widgets.semantic import SecondaryLabel


class _SettingsScroll(QScrollArea):
    """Scroll container for a card's settings.

    Keeps a small, fixed size hint so an expanded card never inflates
    the sidebar's preferred size (which would overflow the sidebar's
    own scroll area on first show); the accordion's layout stretch, not
    the hint, decides the real height.
    """

    def sizeHint(self):
        return QSize(super().sizeHint().width(), 160)


class StageSummaryCard(QFrame):
    """Sidebar card for one analysis stage: icon, title, completed
    check, short summary lines of what was run, and an inline settings
    area that expands accordion-style via the Edit link.
    """

    edit_clicked = pyqtSignal()

    #: Uniform collapsed height (fits header + two summary lines), so
    #: cards match across stages and pages regardless of line count.
    COLLAPSED_HEIGHT = 88

    def __init__(self, title: str, icon: str | None = None, parent=None):
        super().__init__(parent)
        self.setObjectName("stage_card")
        self._expanded = False
        self._read_only = False
        self._expandable = True
        self._settings_scroll: QScrollArea | None = None
        self.setFixedHeight(self.COLLAPSED_HEIGHT)

        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(12, 10, 12, 12)
        self._layout.setSpacing(6)

        header = QHBoxLayout()
        header.setSpacing(8)
        if icon:
            icon_label = QLabel()
            icon_label.setPixmap(
                make_icon(icon, get_color('accent_primary'), 18).pixmap(18, 18))
            header.addWidget(icon_label)
        title_label = QLabel(title)
        title_label.setProperty("role", "stage_card_title")
        header.addWidget(title_label)
        self._check = QLabel()
        self._check.setPixmap(
            make_icon('circle-check', get_color('success'), 14).pixmap(14, 14))
        self._check.hide()
        header.addWidget(self._check)
        header.addStretch()
        self._edit_btn = QPushButton("Edit  ›")
        self._edit_btn.setProperty("role", "link_accent")
        self._edit_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._edit_btn.clicked.connect(self.edit_clicked.emit)
        # Shown by 'set_settings_widget'. A card with nothing behind it
        # must not offer a link that does nothing when clicked.
        self._edit_btn.hide()
        header.addWidget(self._edit_btn)
        self._layout.addLayout(header)

        self._summary_box = QVBoxLayout()
        self._summary_box.setSpacing(2)
        self._layout.addLayout(self._summary_box)

    # -- summary -----------------------------------------------------

    def set_expandable(self, expandable: bool = True) -> None:
        """Allow or forbid opening this card.

        For a card whose settings are all inapplicable -- nothing
        estimable to adjust for, say -- the summary is the whole story.
        Opening onto an empty list tells the user less than the closed
        card already did, so the link goes away and clicks do nothing.
        """
        self._expandable = bool(expandable)
        if not self._expandable and self._expanded:
            self.set_expanded(False)
        self._edit_btn.setVisible(
            self._expandable and self._settings_scroll is not None)

    def set_read_only(self, read_only: bool = True) -> None:
        """Mark the card as informational: expanding shows detail, not settings.

        Only the affordance changes -- 'Details' rather than 'Edit' -- so a
        card with nothing to change does not invite you to change it.
        """
        self._read_only = bool(read_only)
        self._sync_edit_label()

    def _sync_edit_label(self) -> None:
        if self._expanded:
            self._edit_btn.setText("Close")
        else:
            self._edit_btn.setText(
                "Details  \u203a" if self._read_only else "Edit  \u203a")

    def set_completed(self, completed: bool) -> None:
        self._check.setVisible(completed)

    def set_summary(self, lines: list[str]) -> None:
        """Replace the summary lines under the header."""
        while self._summary_box.count():
            item = self._summary_box.takeAt(0)
            w = item.widget()
            if w:
                # Hide first: deleteLater waits for the event loop, and
                # until then a still-visible label paints over its
                # replacement (two summaries on top of each other).
                w.hide()
                w.deleteLater()
        for line in lines:
            label = SecondaryLabel(line)
            label.setWordWrap(True)
            # Add to the layout BEFORE showing it. A widget with no parent
            # is a top-level window, so setVisible(True) here flashes an
            # empty frame on screen for the instant before addWidget
            # reparents it -- visible as a popup whenever summaries are
            # refreshed on every keystroke.
            self._summary_box.addWidget(label)
            label.setVisible(not self._expanded)

    # -- inline settings (accordion) ----------------------------------

    def set_settings_widget(self, widget: QWidget) -> None:
        """Install the stage's settings page, shown when expanded."""
        scroll = _SettingsScroll()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setProperty("role", "thin_scroll")
        # 'Ignored' horizontally, not 'Preferred': a settings panel with a
        # long combo entry or a wide spin row would otherwise demand its
        # content width from the layout and push the card -- and the
        # controls inside it -- past the 320px sidebar. A scroll area
        # exists so content can exceed the viewport; asking to be as wide
        # as that content defeats it.
        scroll.setSizePolicy(QSizePolicy.Policy.Ignored,
                             QSizePolicy.Policy.Expanding)
        scroll.setWidget(widget)
        scroll.hide()
        self._settings_scroll = scroll
        self._layout.addWidget(scroll, 1)
        self._edit_btn.setVisible(self._expandable)
        self._sync_edit_label()

    def is_expanded(self) -> bool:
        return self._expanded

    def set_expanded(self, expanded: bool) -> None:
        """Show or hide the inline settings; summary lines swap out."""
        if self._settings_scroll is None or not self._expandable:
            expanded = False
        self._expanded = expanded
        if self._settings_scroll is not None:
            self._settings_scroll.setVisible(expanded)
        for i in range(self._summary_box.count()):
            w = self._summary_box.itemAt(i).widget()
            if w:
                w.setVisible(not expanded)
        self._sync_edit_label()
        # Collapsed cards are all exactly COLLAPSED_HEIGHT; the expanded
        # card is freed to take the accordion's spare space.
        if expanded:
            self.setMinimumHeight(self.COLLAPSED_HEIGHT)
            self.setMaximumHeight(16777215)  # QWIDGETSIZE_MAX
        else:
            self.setFixedHeight(self.COLLAPSED_HEIGHT)
        policy = self.sizePolicy()
        policy.setVerticalPolicy(
            QSizePolicy.Policy.Expanding if expanded
            else QSizePolicy.Policy.Preferred)
        self.setSizePolicy(policy)

    def mousePressEvent(self, event):
        # A click anywhere on the collapsed card opens it; collapsing
        # is explicit (the Close link) so stray clicks while editing
        # never fold the card up. A card with no settings has nothing to
        # open, and must not collapse whichever card is open by asking.
        if (not self._expanded and self._settings_scroll is not None
                and self._expandable):
            self.edit_clicked.emit()
        super().mousePressEvent(event)


class StageAccordion:
    """Coordinates StageSummaryCards in one sidebar layout: opening a
    card closes the others, and the open card takes the layout's spare
    height (earlier cards pin to the top, later ones to the bottom).

    Usage: construct with the sidebar QVBoxLayout, 'add_card' each card
    after adding it to the layout, then call 'finalize()' immediately
    after the layout's trailing 'addStretch()'.
    """

    def __init__(self, layout):
        self._layout = layout
        self._cards: dict[str, StageSummaryCard] = {}
        self._stretch_idx: int | None = None

    def add_card(self, key: str, card: StageSummaryCard) -> None:
        self._cards[key] = card
        card.edit_clicked.connect(lambda k=key: self.toggle(k))

    def finalize(self) -> None:
        """Record the trailing stretch's index (call right after addStretch)
        and normalise the initial stretch state so the first layout pass
        matches the post-interaction one."""
        self._stretch_idx = self._layout.count() - 1
        self._apply_stretch()

    def toggle(self, key: str) -> None:
        target = self._cards[key]
        expanding = not target.is_expanded()
        for card in self._cards.values():
            card.set_expanded(False)
        if expanding:
            target.set_expanded(True)
        self._apply_stretch()

    def collapse_all(self) -> None:
        for card in self._cards.values():
            card.set_expanded(False)
        self._apply_stretch()

    def _apply_stretch(self) -> None:
        for card in self._cards.values():
            self._layout.setStretch(
                self._layout.indexOf(card), 1 if card.is_expanded() else 0)
        if self._stretch_idx is not None:
            any_open = any(c.is_expanded() for c in self._cards.values())
            self._layout.setStretch(self._stretch_idx, 0 if any_open else 1)


def show_methods_report(parent, h5ad_path) -> None:
    """Show the study's recorded methods (provenance methods.md) for the
    working h5ad, or an info dialog when nothing is recorded yet."""
    from pathlib import Path

    from kosmic.gui.shared import dialogs

    methods = None
    if h5ad_path:
        candidate = Path(h5ad_path).parent / "methods.md"
        if candidate.exists():
            methods = candidate.read_text(encoding="utf-8")
    if not methods:
        dialogs.info(
            parent, "No Report Yet",
            "No methods have been recorded for this study yet.\n"
            "Run an analysis stage first.")
        return
    from PyQt6.QtWidgets import QDialog, QTextBrowser
    dlg = QDialog(parent)
    dlg.setWindowTitle("Analysis Report")
    dlg.resize(560, 620)
    dlg_layout = QVBoxLayout(dlg)
    browser = QTextBrowser()
    browser.setMarkdown(methods)
    dlg_layout.addWidget(browser)
    dlg.exec()
