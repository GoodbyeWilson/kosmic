# Render Lucide SVG icons recoloured at runtime; SVGs in kosmic/gui/icons/
# carry stroke="currentColor" which is substituted with the requested hex.
from __future__ import annotations

import os
import sys
from functools import lru_cache

from PyQt6.QtCore import Qt, QRectF, QByteArray, QFileInfo
from PyQt6.QtGui import QPainter, QPixmap, QIcon, QAbstractFileIconProvider
from PyQt6.QtSvg import QSvgRenderer


def _icons_dir() -> str:
    # sys._MEIPASS is set when running from a PyInstaller bundle.
    if getattr(sys, "frozen", False):
        base = sys._MEIPASS  # type: ignore[attr-defined]
    else:
        base = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base, "icons")


def _read_svg(name: str) -> str:
    """Read the raw SVG text for *name* (without .svg extension)."""
    path = os.path.join(_icons_dir(), f"{name}.svg")
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def _recolour(svg_text: str, color_hex: str) -> str:
    """Substitute both stroke="currentColor" and fill="currentColor"."""
    svg_text = svg_text.replace('stroke="currentColor"', f'stroke="{color_hex}"')
    svg_text = svg_text.replace('fill="currentColor"', f'fill="{color_hex}"')
    return svg_text


def render_svg_to_painter(
    painter: QPainter,
    name: str,
    color_hex: str,
    rect: QRectF,
) -> None:
    """Render the named SVG icon into *rect* on *painter*.

    'currentColor' in the SVG is replaced with *color_hex* (e.g. '"#cccccc"').
    """
    svg_text = _recolour(_read_svg(name), color_hex)
    data = QByteArray(svg_text.encode("utf-8"))
    renderer = QSvgRenderer(data)
    if renderer.isValid():
        renderer.render(painter, rect)


@lru_cache(maxsize=128)
def make_icon(name: str, color_hex: str, size: int = 16) -> QIcon:
    """Return a 'QIcon' for the named SVG, coloured with *color_hex*.

    Results are cached by (name, color_hex, size) so repeated calls are free.
    """
    svg_text = _recolour(_read_svg(name), color_hex)
    data = QByteArray(svg_text.encode("utf-8"))
    renderer = QSvgRenderer(data)

    pixmap = QPixmap(size, size)
    pixmap.fill(Qt.GlobalColor.transparent)

    p = QPainter(pixmap)
    if renderer.isValid():
        renderer.render(p, QRectF(0, 0, size, size))
    p.end()

    return QIcon(pixmap)


class LucideFileIconProvider(QAbstractFileIconProvider):
    """QFileSystemModel icon provider using the Lucide folder/file SVGs."""

    def icon(self, info_or_type):
        from kosmic.gui.shared.theme import get_color
        fg = get_color('fg_secondary')
        is_dir = (
            (isinstance(info_or_type, QFileInfo) and info_or_type.isDir())
            or info_or_type == QAbstractFileIconProvider.IconType.Folder
        )
        return make_icon("folder" if is_dir else "file", fg, 16)
