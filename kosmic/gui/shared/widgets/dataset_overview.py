# DatasetOverview -- the one layout a "here is your dataset" screen uses.
#
# scRNA's Load Data step and DE's Dataset step both show the same thing:
# a page header with a readiness chip, a hero card (icon, name, caption,
# status line, stat blocks, action buttons), a row of InfoPanels, and an
# optional note. They were hand-built separately and drifted -- different
# widths, spacing, centring, overlaps. This widget owns the geometry;
# the pages supply names, numbers, panels and buttons.
#
#     ov = DatasetOverview(
#         "Dataset", "Review the study's dataset before running DE.",
#         stats=(("cells", "Cells"), ("genes", "Genes"), ("samples", "Samples")))
#     ov.add_action(SecondaryButton("Change study..."))
#     ov.add_panel(InfoPanel("DATASET"))
#     ov.set_name("GSE292067"); ov.stats["cells"].set_value("154,207")
#     ov.set_status("Loaded", 'success'); ov.set_ready(True)
#
# Anything a page needs beneath the panels goes in 'ov.extra_layout'.
from __future__ import annotations

from typing import Sequence

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QFrame, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget,
)

from kosmic.gui.shared.widgets.overview import IconTile, InfoPanel, StatBlock, StatusChip
from kosmic.gui.shared.widgets.semantic import CaptionLabel, StatusLabel


class DatasetOverview(QWidget):
    """Header + hero card + info panels + note, laid out once."""

    #: Width of the column every dataset page sits in.
    COLUMN_MAX_WIDTH = 1380
    ACTION_WIDTH = 200

    def __init__(self, title: str, caption: str, *,
                 stats: Sequence[tuple[str, str]],
                 icon: str = "dataset", icon_label: str = "h5ad",
                 parent: QWidget | None = None):
        super().__init__(parent)
        self.setMaximumWidth(self.COLUMN_MAX_WIDTH)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(12)

        # -- header: title + caption left, readiness chip right --
        header_row = QHBoxLayout()
        header_col = QVBoxLayout()
        header_col.setSpacing(2)
        self.title = QLabel(title)
        self.title.setProperty("role", "page_header")
        header_col.addWidget(self.title)
        self.caption = CaptionLabel(caption)
        header_col.addWidget(self.caption)
        header_row.addLayout(header_col)
        header_row.addStretch()
        self.ready_chip = StatusChip("✓ Ready to continue")
        self.ready_chip.hide()
        header_row.addWidget(self.ready_chip, alignment=Qt.AlignmentFlag.AlignTop)
        outer.addLayout(header_row)

        # -- hero card --
        self.hero = QFrame()
        self.hero.setObjectName("dataset_hero_card")
        hero_row = QHBoxLayout(self.hero)
        hero_row.setContentsMargins(24, 18, 24, 18)
        hero_row.setSpacing(20)
        hero_row.addWidget(IconTile(icon, icon_label, size=64, icon_size=28),
                           alignment=Qt.AlignmentFlag.AlignTop)

        identity = QVBoxLayout()
        identity.setSpacing(2)
        name_row = QHBoxLayout()
        name_row.setSpacing(10)
        self.name = QLabel("")
        self.name.setObjectName("hero_dataset_name")
        self.name.setWordWrap(True)
        name_row.addWidget(self.name)
        self.badge = QLabel("Active")
        self.badge.setProperty("role", "active_badge")
        self.badge.hide()
        name_row.addWidget(self.badge, alignment=Qt.AlignmentFlag.AlignVCenter)
        name_row.addStretch()
        identity.addLayout(name_row)
        self.subtitle = CaptionLabel("")
        identity.addWidget(self.subtitle)
        self.status = StatusLabel("", state='info')
        self.status.setWordWrap(True)
        self.status.hide()
        identity.addWidget(self.status)
        identity.addSpacing(12)

        stats_row = QHBoxLayout()
        stats_row.setSpacing(28)
        self.stats: dict[str, StatBlock] = {}
        for i, (key, label) in enumerate(stats):
            if i:
                divider = QFrame()
                divider.setProperty("role", "stat_divider")
                divider.setFixedWidth(1)
                divider.setFixedHeight(38)
                stats_row.addWidget(divider)
            block = StatBlock(label)
            stats_row.addWidget(block)
            self.stats[key] = block
        stats_row.addStretch()
        identity.addLayout(stats_row)
        hero_row.addLayout(identity, stretch=1)

        self._actions = QVBoxLayout()
        self._actions.setSpacing(6)
        self._actions_title = QLabel("ACTIONS")
        self._actions_title.setProperty("role", "panel_title")
        self._actions_title.hide()
        self._actions.addWidget(self._actions_title)
        self._actions.addStretch()
        hero_row.addLayout(self._actions)
        outer.addWidget(self.hero)

        # -- info panels --
        self.panels_widget = QWidget()
        self._panels = QHBoxLayout(self.panels_widget)
        self._panels.setContentsMargins(0, 0, 0, 0)
        self._panels.setSpacing(12)
        outer.addWidget(self.panels_widget)

        # -- note under the panels --
        self.note = QLabel("")
        self.note.setProperty("role", "info_banner")
        self.note.setWordWrap(True)
        self.note.hide()
        outer.addWidget(self.note)

        # -- page-specific extras --
        self.extra_layout = QVBoxLayout()
        self.extra_layout.setContentsMargins(0, 0, 0, 0)
        self.extra_layout.setSpacing(12)
        outer.addLayout(self.extra_layout)

    # -- building --------------------------------------------------------

    def add_action(self, button: QPushButton) -> QPushButton:
        """Put a button in the hero's action column (fixed width, stacked)."""
        button.setFixedWidth(self.ACTION_WIDTH)
        self._actions.insertWidget(self._actions.count() - 1, button)
        self._actions_title.show()
        return button

    def add_panel(self, panel: InfoPanel) -> InfoPanel:
        self._panels.addWidget(panel, 1)
        return panel

    # -- state ------------------------------------------------------------

    def set_name(self, text: str, tooltip: str = "") -> None:
        self.name.setText(text)
        self.name.setToolTip(tooltip)

    def set_subtitle(self, text: str) -> None:
        self.subtitle.setText(text)

    def set_status(self, text: str, state: str = 'info', tooltip: str = "") -> None:
        """Status line under the subtitle; empty text hides it."""
        self.status.setText(text)
        self.status.set_state(state)
        self.status.setToolTip(tooltip)
        self.status.setVisible(bool(text))

    def set_ready(self, ready: bool) -> None:
        self.ready_chip.setVisible(ready)

    def set_note(self, text: str) -> None:
        self.note.setText(text)
        self.note.setVisible(bool(text))
