# Uniform base class for every QThread worker in KOSMIC.
#
# The signal contract is split: 'finished_ok = pyqtSignal(object)' for
# success (payload only) and 'failed = pyqtSignal(str)' for the error
# path (traceback string). Subclasses override '_run()'; this base
# handles the try/except wrapping and signal emission.
#
# Subclass contract
#
# - Override '_run()'. Return the payload the caller needs (a DataFrame,
#   dict, list, or 'None'). Raising an exception is caught and emitted
#   as 'failed' with the traceback.
# - Emit 'progress' for human-readable status messages.
# - Emit 'progress_pct' for 0-100 integer progress (e.g. CC permutation loop).
# - Subclasses may declare *additional* signals (e.g. 'fold_progress' on
#   LOO workers) -- these coexist with BaseWorker's declared signals.
#
# Migration notes
#
# **From Pattern A** ('finished.emit(True, msg, data)' / 'finished.emit(False, err, None)'):
#
# Before::
#
#     class FooWorker(QThread):
#         progress = pyqtSignal(str)
#         finished = pyqtSignal(bool, str, object)
#
#         def run(self):
#             try:
#                 ...
#                 self.finished.emit(True, "Done", result)
#             except Exception as e:
#                 self.finished.emit(False, str(e), None)
#
#     # Caller:
#     def _on_finished(self, success, msg, data):
#         if success: ...
#         else: dialogs.warning(...)
#
# After::
#
#     class FooWorker(BaseWorker):
#         # progress inherited
#
#         def _run(self):
#             ...
#             self.progress.emit("Done")
#             return result   # <- BaseWorker.run emits finished_ok(result)
#
#     # Caller:
#     worker.finished_ok.connect(self._on_result)
#     worker.failed.connect(self._on_error)
from __future__ import annotations

import traceback

from PyQt6.QtCore import QThread, pyqtSignal


class BaseWorker(QThread):
    """QThread with uniform success/error signal shape + exception safety.

    Signals
    -------
    progress : str
        Human-readable status text. Typically routed to OutputPanel.
    progress_pct : int
        0-100 integer progress, for progress bars on long loops.
    finished_ok : object
        Emitted once on success with the subclass's result payload.
    failed : str
        Emitted once on failure with the exception message + traceback.
    """

    progress = pyqtSignal(str)
    progress_pct = pyqtSignal(int)
    finished_ok = pyqtSignal(object)
    failed = pyqtSignal(str)

    # Subclasses override this.
    def _run(self):
        """Perform the actual work. Return the payload the caller
        should receive via 'finished_ok'. Any exception raised here is
        caught and routed through 'failed' with a full traceback.
        """
        raise NotImplementedError(
            f"{type(self).__name__} must override _run()")

    def run(self) -> None:  # final -- subclasses override _run() instead
        try:
            result = self._run()
        except Exception as exc:  # noqa: BLE001 -- by design
            self.failed.emit(f"{exc}\n{traceback.format_exc()}")
        else:
            self.finished_ok.emit(result)
