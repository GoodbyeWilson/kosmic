# Intent-named wrappers around QMessageBox: info / warning / error / choose /
# confirm.
from __future__ import annotations

from typing import Optional, Sequence

from PyQt6.QtWidgets import QMessageBox, QWidget


def info(parent: Optional[QWidget], title: str, text: str) -> None:
    """Show a non-blocking informational message."""
    QMessageBox.information(parent, title, text)


def warning(parent: Optional[QWidget], title: str, text: str) -> None:
    """Show a warning. Use for recoverable issues (bad input, partial result).

    For genuine errors (operation failed, exception caught) use 'error'.
    """
    QMessageBox.warning(parent, title, text)


def error(parent: Optional[QWidget], title: str, text: str) -> None:
    """Show an error. Use for failed operations / caught exceptions."""
    QMessageBox.critical(parent, title, text)


def choose(parent: Optional[QWidget], title: str, text: str,
           options: Sequence[str], *, escape: Optional[str] = None) -> Optional[str]:
    """Ask the user to pick one of 'options', shown as buttons in order.

    Returns the label of the button clicked. Esc or closing the box
    returns 'escape' when it is one of the options, else None.
    """
    box = QMessageBox(QMessageBox.Icon.Question, title, text, parent=parent)
    buttons = {}
    for label in options:
        buttons[label] = box.addButton(label, QMessageBox.ButtonRole.ActionRole)
    if options:
        box.setDefaultButton(buttons[options[0]])
    if escape in buttons:
        box.setEscapeButton(buttons[escape])
    box.exec()
    clicked = box.clickedButton()
    for label, button in buttons.items():
        if button is clicked:
            return label
    return escape if escape in buttons else None


def confirm(parent: Optional[QWidget], title: str, text: str, *,
            default_yes: bool = False) -> bool:
    """Yes / No confirmation. Returns True for Yes, False for No or Esc.

    Set 'default_yes=True' to make Yes the keyboard-focused button.
    Default is No-focused (safer for destructive prompts).
    """
    default = (QMessageBox.StandardButton.Yes if default_yes
               else QMessageBox.StandardButton.No)
    reply = QMessageBox.question(
        parent, title, text,
        QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        default,
    )
    return reply == QMessageBox.StandardButton.Yes
