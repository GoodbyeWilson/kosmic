# Tutorial Controller -- state machine driving the guided in-app tutorial.
# Step copy lives in 'steps.toml'; this file owns navigation, signal
# waiting, and the auto-action implementations the steps reference by id.

import contextlib
import tomllib
from dataclasses import dataclass
from importlib.resources import files
from typing import Callable, Optional

from PyQt6.QtCore import QObject, QTimer, pyqtSignal

from kosmic.gui.help.tutorial.overlay import TutorialOverlay
from kosmic.gui.shared import dialogs


@dataclass
class TutorialStep:
    tab_index: int                          # Which tab (0-5), -1 for project page
    title: str                              # Card title
    body: str                               # Instruction text
    target_widget: Optional[str] = None     # Dot-path to widget for spotlight
    wait_for_signal: Optional[str] = None   # Signal dot-path to wait for before enabling Next
    auto_action: Optional[Callable] = None  # "Do It For Me" callback
    allow_manual_next: bool = True          # Can user click Next without completing action?
    next_text: str = "Next"                 # Button label
    on_advance: Optional[Callable] = None   # Runs when Next is clicked (before advancing)
    on_show: Optional[Callable] = None      # Runs when step is displayed


def _load_steps(ctrl: 'TutorialController') -> list[TutorialStep]:
    """Parse steps.toml and resolve *_id fields to bound controller methods."""
    raw = files('kosmic.gui.help.tutorial').joinpath('steps.toml').read_bytes()
    doc = tomllib.loads(raw.decode('utf-8'))

    def _resolve(name: Optional[str]) -> Optional[Callable]:
        if not name:
            return None
        method = getattr(ctrl, f'_{name}', None)
        if method is None:
            raise ValueError(f"steps.toml refers to unknown callback id: {name!r}")
        return method

    steps = []
    for entry in doc['steps']:
        steps.append(TutorialStep(
            tab_index=entry['tab_index'],
            title=entry['title'],
            body=entry['body'].strip(),
            target_widget=entry.get('target_widget'),
            wait_for_signal=entry.get('wait_for_signal'),
            allow_manual_next=entry.get('allow_manual_next', True),
            next_text=entry.get('next_text', 'Next'),
            auto_action=_resolve(entry.get('auto_action_id')),
            on_advance=_resolve(entry.get('on_advance_id')),
            on_show=_resolve(entry.get('on_show_id')),
        ))
    return steps


class TutorialController(QObject):
    """Drives the tutorial overlay through its step sequence."""

    tutorial_started = pyqtSignal()
    tutorial_ended = pyqtSignal()

    def __init__(self, main_window):
        super().__init__(main_window)
        self.main_window = main_window
        self.overlay = TutorialOverlay(main_window.centralWidget())
        self._steps: list = []
        self._current_step = 0
        self._tutorial_active = False
        self._temp_dir: Optional[str] = None
        self._waiting_signal = None
        self._waiting_connection = None

        self.overlay.card.next_clicked.connect(self._on_next)
        self.overlay.card.skip_clicked.connect(self.end_tutorial)
        self.overlay.card.auto_action_clicked.connect(self._on_auto_action)

    @property
    def is_active(self) -> bool:
        return self._tutorial_active

    def start_tutorial(self):
        """Begin the tutorial from step 0."""
        if self._tutorial_active:
            if not dialogs.confirm(self.main_window, "Restart Tutorial",
                                   "A tutorial is already in progress. Restart from the beginning?"):
                return
            self._disconnect_waiting_signal()

        self._steps = _load_steps(self)
        self._current_step = 0
        self._tutorial_active = True

        cw = self.main_window.centralWidget()
        self.overlay.setGeometry(cw.rect())

        self.tutorial_started.emit()
        self._show_current_step()

    def end_tutorial(self):
        """End the tutorial and clean up."""
        self._tutorial_active = False
        self._disconnect_waiting_signal()
        self.overlay.hide_overlay()
        self.tutorial_ended.emit()
        self.main_window.status_bar.showMessage("Tutorial ended. Your data is still loaded.")

    def resize_overlay(self):
        """Called by MainWindow on resize to keep overlay sized correctly."""
        if self._tutorial_active:
            cw = self.main_window.centralWidget()
            self.overlay.setGeometry(cw.rect())
            self._show_current_step()

    # Step navigation

    def _on_next(self):
        if not self._tutorial_active:
            return

        step = self._steps[self._current_step]
        if step.on_advance:
            step.on_advance()

        self._disconnect_waiting_signal()

        self._current_step += 1
        if self._current_step >= len(self._steps):
            self.end_tutorial()
            return

        step = self._steps[self._current_step]
        if step.tab_index >= 0:
            self.main_window.switch_tab(step.tab_index)

        # Defer until Qt finishes layout from switch_tab().
        QTimer.singleShot(0, self._show_current_step)

    def _on_auto_action(self):
        if not self._tutorial_active:
            return
        step = self._steps[self._current_step]
        if step.auto_action:
            step.auto_action()

    def _show_current_step(self):
        if not self._tutorial_active or self._current_step >= len(self._steps):
            return

        step = self._steps[self._current_step]

        if step.on_show:
            step.on_show()

        target = self._resolve_widget(step.target_widget) if step.target_widget else None

        self.overlay.show_step(
            title=step.title,
            body=step.body,
            step=self._current_step + 1,
            total=len(self._steps),
            target_widget=target,
            next_text=step.next_text,
            has_auto_action=step.auto_action is not None,
            allow_next=step.allow_manual_next,
        )

        if step.wait_for_signal:
            self._connect_waiting_signal(step.wait_for_signal)

    # Widget and signal resolution

    def _resolve_widget(self, dot_path: str):
        """Resolve a dot-separated path from main_window to a widget.

        E.g. "cluster_tab.run_cluster_btn" -> main_window.cluster_tab.run_cluster_btn
        """
        obj = self.main_window
        for attr in dot_path.split('.'):
            obj = getattr(obj, attr, None)
            if obj is None:
                return None
        return obj

    def _connect_waiting_signal(self, signal_path: str):
        """Connect to a signal and auto-advance when it fires."""
        self._disconnect_waiting_signal()

        signal = self._resolve_widget(signal_path)
        if signal is None:
            return

        def _on_signal_fired(*args):
            if self._tutorial_active:
                self.overlay.card.enable_next()
                self.overlay.card._next_btn.setText("Next")
                # Pacing -- let the user see the "complete" state.
                QTimer.singleShot(500, self._on_next)

        try:
            signal.connect(_on_signal_fired)
            self._waiting_signal = signal
            self._waiting_connection = _on_signal_fired
        except Exception:
            pass

    def _disconnect_waiting_signal(self):
        if self._waiting_signal is not None and self._waiting_connection is not None:
            with contextlib.suppress(Exception):
                self._waiting_signal.disconnect(self._waiting_connection)
            self._waiting_signal = None
            self._waiting_connection = None

    # Auto-action implementations (referenced by *_id fields in steps.toml)
