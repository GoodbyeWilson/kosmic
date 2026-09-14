# Build a Qt layout with zero margins and zero spacing in one call.

from __future__ import annotations

from typing import Optional, Type, TypeVar

from PyQt6.QtWidgets import QLayout, QWidget


_LayoutT = TypeVar("_LayoutT", bound=QLayout)


def borderless(
    layout_cls: Type[_LayoutT],
    parent: Optional[QWidget] = None,
) -> _LayoutT:
    """Construct 'layout_cls' with zero margins and zero spacing.

    Parameters
    ----------
    layout_cls
        'QVBoxLayout' / 'QHBoxLayout' / 'QGridLayout' / etc.
    parent
        Optional parent widget. If given, the layout is installed on
        the parent (same as calling 'layout_cls(parent)' directly).
    """
    layout: _LayoutT = layout_cls(parent) if parent is not None else layout_cls()
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(0)
    return layout
