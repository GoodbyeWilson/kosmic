# Semantic widget classes for KOSMIC.
#
# Thin subclasses that set the right 'objectName' / dynamic 'role'
# property so the global stylesheet in 'gui.shared.theme' picks
# them up automatically.  Using these instead of inline
# 'setStyleSheet' calls keeps styling centralised -- a theme change
# in 'theme.py' flows to every instance without touching call sites.
#
# Each class takes optional 'text' + 'parent' mirroring the Qt
# constructors they wrap.  No other behaviour is added; callers keep
# using the underlying 'QPushButton' / 'QLabel' / 'QFrame' API
# for everything else.
#
# Companion selectors live in 'gui.shared.theme' (see
# 'QPushButton#primary_button', 'QLabel[role="secondary"]', etc.).
from __future__ import annotations

from PyQt6.QtWidgets import QFrame, QLabel, QPushButton


# Buttons

class PrimaryButton(QPushButton):
    """Accent-coloured action button (Run / Apply / Generate / etc.)."""

    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setObjectName("primary_button")


class SecondaryButton(QPushButton):
    """Neutral-bg button with subtle border (Cancel / Close / Skip)."""

    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setObjectName("secondary_button")


class DangerButton(QPushButton):
    """Destructive action button (Delete / Reset / Clear)."""

    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setObjectName("danger_button")


class TabButton(QPushButton):
    """In-dialog tab-strip button.  Use 'set_active' to flip state.

    Stylesheet reads the dynamic 'active' property, so toggling
    forces a re-polish.
    """

    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setObjectName("tab_button")
        self.setProperty("active", "false")

    def set_active(self, active: bool) -> None:
        self.setProperty("active", "true" if active else "false")
        # Qt doesn't re-evaluate stylesheet selectors on a dynamic
        # property change unless we poke the style engine.
        style = self.style()
        if style is not None:
            style.unpolish(self)
            style.polish(self)


# Labels

class SectionHeader(QLabel):
    """Uppercase small-caps section header for sidebars + group titles."""

    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setProperty("role", "section_header")


class HeaderLabel(QLabel):
    """Primary heading inside a page body (slightly larger, bold)."""

    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setProperty("role", "header")


class SecondaryLabel(QLabel):
    """Grey secondary text (subtitles, metadata, footnotes)."""

    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setProperty("role", "secondary")


class HintLabel(QLabel):
    """Small italic grey hint text for inline help.

    Wraps by default. A QLabel that does not wrap reports its whole
    string as its minimum width, so one sentence of help demanded 1,560
    pixels inside a 274-pixel card and dragged every control in that
    card out of the container with it. Half the callers were already
    setting word wrap by hand; inline help that cannot wrap is not
    useful anyway.
    """

    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setProperty("role", "hint")
        self.setWordWrap(True)


class CaptionLabel(QLabel):
    """11px grey label for 'Caption: value' grid rows."""

    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setProperty("role", "caption")


class ValueLabel(QLabel):
    """Bold primary-colour label for the value cell of a KPI grid."""

    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setProperty("role", "value")


class StatusLabel(QLabel):
    """Inline status line with state-dependent colour.

    Use 'set_state' with one of 'info', 'success',
    'warning', 'error'.  Default state is 'info'.
    """

    _STATES = ('info', 'success', 'warning', 'error')

    def __init__(self, text: str = "", parent=None, state: str = 'info'):
        super().__init__(text, parent)
        self.set_state(state)

    def set_state(self, state: str) -> None:
        if state not in self._STATES:
            raise ValueError(
                f"StatusLabel state must be one of {self._STATES}, "
                f"got {state!r}")
        self.setProperty("role", f"status_{state}")
        style = self.style()
        if style is not None:
            style.unpolish(self)
            style.polish(self)


# Frames

class CardFrame(QFrame):
    """Bordered content panel (bg_secondary + rounded corners)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setProperty("role", "card")
