# Left-pane file explorer: workflow steps, directory indicator, file tree, status panel.
from __future__ import annotations

import os
import subprocess
import sys

from PyQt6.QtCore import Qt, QDir, pyqtSignal
from PyQt6.QtGui import QFileSystemModel, QFont, QPainter, QColor
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QTreeView,
    QMenu,
)

from kosmic.gui.shared.borderless import borderless
from kosmic.gui.shared.theme import get_color
from kosmic.gui.shared.icon_provider import make_icon, LucideFileIconProvider
from kosmic.gui.shared.workflow import WorkflowSteps, StepStatusPanel
from kosmic.gui.shared import dialogs


class ExplorerPane(QWidget):
    """Left pane: workflow steps | FILES header | directory | tree | status."""

    file_activated = pyqtSignal(str)  # double-click on a file
    folder_opened = pyqtSignal(str)   # user picked a folder via Open Folder
    step_clicked = pyqtSignal(int)    # workflow step clicked

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumWidth(220)
        self._root_path = ""
        self._root_set = False

        layout = borderless(QVBoxLayout, self)

        self._workflow_steps = WorkflowSteps()
        self._workflow_steps.step_clicked.connect(self.step_clicked.emit)
        layout.addWidget(self._workflow_steps)

        self._divider = QWidget()
        self._divider.setProperty("role", "divider")
        self._divider.setFixedHeight(1)
        layout.addWidget(self._divider)

        self._files_header = QLabel("FILES")
        self._files_header.setProperty("role", "section_header")
        _hf = self._files_header.font()
        _hf.setLetterSpacing(QFont.SpacingType.AbsoluteSpacing, 1.5)
        self._files_header.setFont(_hf)
        self._files_header.setFixedHeight(24)
        self._files_header.setContentsMargins(10, 4, 0, 0)
        layout.addWidget(self._files_header)

        # Directory indicator (folder icon + name)
        self._dir_indicator = QWidget()
        self._dir_indicator.setFixedHeight(36)
        dir_layout = QHBoxLayout(self._dir_indicator)
        dir_layout.setContentsMargins(10, 4, 10, 4)
        dir_layout.setSpacing(6)
        self._dir_icon_label = QLabel()
        self._dir_icon_label.setFixedSize(16, 16)
        dir_layout.addWidget(self._dir_icon_label)
        self._dir_name_label = QLabel("No directory set")
        _df = QFont("Segoe UI", 10)
        _df.setWeight(QFont.Weight.Medium)
        self._dir_name_label.setFont(_df)
        dir_layout.addWidget(self._dir_name_label, 1)
        layout.addWidget(self._dir_indicator)

        self._model = QFileSystemModel()
        self._model.setReadOnly(True)
        self._model.setIconProvider(LucideFileIconProvider())
        self._model.setFilter(
            QDir.Filter.AllDirs | QDir.Filter.Files | QDir.Filter.NoDotAndDotDot
        )

        self._tree = QTreeView()
        self._tree.setModel(self._model)
        self._tree.setHeaderHidden(True)
        for col in (1, 2, 3):  # hide Size, Type, Date Modified — show only Name
            self._tree.hideColumn(col)
        self._tree.setAnimated(True)
        self._tree.setIndentation(16)
        self._tree.doubleClicked.connect(self._on_double_click)
        self._tree.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self._tree.customContextMenuRequested.connect(self._on_tree_context_menu)
        layout.addWidget(self._tree, 1)

        self._status_divider = QWidget()
        self._status_divider.setProperty("role", "divider")
        self._status_divider.setFixedHeight(1)
        layout.addWidget(self._status_divider)

        self._data_status = StepStatusPanel()
        layout.addWidget(self._data_status)

    # ---- workflow + status pass-throughs ---------------------------------
    def set_workflow_steps(self, steps):
        self._workflow_steps.set_steps(steps)

    def set_active_step(self, index: int):
        self._workflow_steps.set_active(index)

    def set_step_complete(self, index: int, complete: bool = True):
        self._workflow_steps.set_complete(index, complete)

    def update_data_status(self, text: str, state: str = "info"):
        """state: 'info' | 'success' | 'warning' | 'error'."""
        self._data_status.update_status(text, state)

    @property
    def progress_panel(self) -> StepStatusPanel:
        return self._data_status

    # ---- root / navigation -----------------------------------------------
    def set_root(self, path: str):
        self._root_path = path
        idx = self._model.setRootPath(path)
        self._tree.setRootIndex(idx)
        self._root_set = True
        self._dir_name_label.setText(os.path.basename(path) or path)
        self._dir_name_label.setToolTip(path)
        self._update_dir_indicator_style()

    def _on_double_click(self, index):
        path = self._model.filePath(index)
        if path and not self._model.isDir(index):
            self.file_activated.emit(path)

    # ---- context menu ----------------------------------------------------
    def _on_tree_context_menu(self, pos):
        index = self._tree.indexAt(pos)
        if not index.isValid():
            return
        path = self._model.filePath(index)
        is_dir = self._model.isDir(index)
        name = os.path.basename(path)

        menu = QMenu(self)
        reveal = menu.addAction("Reveal in Finder")
        reveal.triggered.connect(lambda: self._reveal_in_finder(path))
        menu.addSeparator()
        delete = menu.addAction(f"Delete {'Folder' if is_dir else 'File'}")
        delete.triggered.connect(lambda: self._delete_path(path, name, is_dir))
        menu.exec(self._tree.viewport().mapToGlobal(pos))

    def _reveal_in_finder(self, path):
        if sys.platform == 'darwin':
            subprocess.Popen(['open', '-R', path])
        elif sys.platform == 'win32':
            subprocess.Popen(['explorer', '/select,', os.path.normpath(path)])
        else:
            subprocess.Popen(['xdg-open', os.path.dirname(path)])

    def _delete_path(self, path, name, is_dir):
        item_type = "folder and all its contents" if is_dir else "file"
        if not dialogs.confirm(
            self,
            f"Delete {name}",
            f"Are you sure you want to delete this {item_type}?\n\n{path}\n\nThis cannot be undone.",
        ):
            return
        try:
            if is_dir:
                import shutil
                shutil.rmtree(path)
            else:
                os.remove(path)
        except Exception as e:
            dialogs.warning(self, "Delete Failed", f"Could not delete:\n{e}")

    # ---- painting + theme ------------------------------------------------
    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(0, 0, self.width(), self.height(),
                         QColor(get_color('bg_secondary')))
        painter.end()
        super().paintEvent(event)

    def _update_dir_indicator_style(self):
        fg = get_color('fg_primary') if self._root_set else get_color('fg_tertiary')
        self._dir_name_label.setStyleSheet(f"color: {fg}; background: transparent;")
        icon_color = (get_color('accent_primary') if self._root_set
                      else get_color('fg_tertiary'))
        icon_name = "folder-open" if self._root_set else "folder"
        self._dir_icon_label.setPixmap(
            make_icon(icon_name, icon_color, 16).pixmap(16, 16))

    def refresh_theme(self):
        self._update_dir_indicator_style()
        for w in (self._divider, self._status_divider, self._files_header):
            w.style().unpolish(w)
            w.style().polish(w)
        self._workflow_steps.refresh_theme()
        self._data_status.refresh_theme()
        self.update()
