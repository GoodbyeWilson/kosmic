"""Refreshing a stage card's summary must not flash a top-level window.

A QWidget with no parent IS a window. Building the summary labels and
calling setVisible() on them before adding them to the layout showed an
empty frame on screen for an instant -- harmless once, but the DE page
refreshes its card summaries on every settings keystroke, so it read as
something being launched over and over.
"""
from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication, QLabel, QVBoxLayout, QWidget  # noqa: E402

from kosmic.gui.shared.widgets import StageSummaryCard  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def summary_labels(card):
    box = card._summary_box
    return [box.itemAt(i).widget() for i in range(box.count())
            if box.itemAt(i).widget() is not None]


def test_summary_labels_are_parented_before_being_shown(app):
    card = StageSummaryCard("Filters", icon=None)
    card.set_summary(["first line", "second line"])
    labels = summary_labels(card)
    assert len(labels) == 2
    for label in labels:
        # A parented widget is not a window; an unparented one is.
        assert label.parentWidget() is not None
        assert not label.isWindow()


def test_repeated_refresh_leaves_exactly_the_new_lines(app):
    card = StageSummaryCard("Filters", icon=None)
    for i in range(5):
        card.set_summary([f"line {i}a", f"line {i}b"])
        app.processEvents()
    labels = summary_labels(card)
    assert [w.text() for w in labels] == ["line 4a", "line 4b"]


def test_no_stray_top_level_windows_appear(app):
    """The symptom itself: a window appearing during a refresh."""
    card = StageSummaryCard("Filters", icon=None)
    card.set_summary(["one", "two"])
    app.processEvents()
    before = {w for w in app.topLevelWidgets() if w.isVisible()}
    for i in range(10):
        card.set_summary([f"{i}", f"{i}"])
    after = {w for w in app.topLevelWidgets() if w.isVisible()}
    assert after == before


def test_summary_visibility_follows_expansion(app):
    card = StageSummaryCard("Filters", icon=None)
    page = QWidget()
    QVBoxLayout(page).addWidget(QLabel("settings"))
    card.set_settings_widget(page)
    card.set_summary(["a", "b"])

    card.set_expanded(True)
    assert all(not w.isVisibleTo(card) for w in summary_labels(card))
    card.set_summary(["c", "d"])
    assert all(not w.isVisibleTo(card) for w in summary_labels(card))

    card.set_expanded(False)
    assert all(w.isVisibleTo(card) for w in summary_labels(card))
