"""Long plain-text tooltips are converted to rich text, which Qt wraps.

A plain-text tooltip is shown as one line per paragraph; 52 of KOSMIC's
tooltips had a line over 90 characters and ran across the screen.
"""
from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication, QCheckBox  # noqa: E402

from kosmic.gui.shared.theme import ThemeMode, apply_theme, wrapped_tooltip  # noqa: E402

LONG = ("Hartung-Knapp-Sidik-Jonkman adjustment: t-test on k-1 df plus an "
        "adjusted SE. Applied to DL and REML only.")


@pytest.fixture(scope="module")
def app():
    app = QApplication.instance() or QApplication([])
    apply_theme(app, ThemeMode.DARK)
    return app


def test_a_long_tooltip_becomes_wrapped_rich_text(app):
    cb = QCheckBox("HKSJ")
    cb.setToolTip(LONG + "\n\nSecond paragraph & more.")
    assert cb.toolTip() == ("<p>" + LONG + "</p>"
                            "<p>Second paragraph &amp; more.</p>")


@pytest.mark.parametrize("text", ["Short tooltip.", "<p>Already rich text.</p>", ""])
def test_short_and_rich_tooltips_are_left_alone(app, text):
    cb = QCheckBox("x")
    cb.setToolTip(text)
    assert cb.toolTip() == text


def test_single_line_breaks_are_kept():
    assert wrapped_tooltip(LONG + "\nnext line") == f"<p>{LONG}<br>next line</p>"
