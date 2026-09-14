# Wire a BaseWorker's signals to caller-supplied callbacks and start it.
# Caller must hold a reference to the worker (e.g. self._worker = ...) so
# Qt doesn't garbage-collect it mid-run.
from __future__ import annotations

from typing import Any, Callable, Optional

from kosmic.gui.shared.widgets.base_worker import BaseWorker


def run_worker(
    worker: BaseWorker,
    *,
    on_finished: Callable[[Any], None],
    on_failed: Optional[Callable[[str], None]] = None,
    on_progress: Optional[Callable[[str], None]] = None,
    on_progress_pct: Optional[Callable[[int], None]] = None,
) -> BaseWorker:
    """Wire BaseWorker signals to the supplied callbacks and start the
    worker.

    Parameters
    ----------
    worker
        The 'BaseWorker' instance to start. The caller must keep a
        reference (e.g. 'self._foo_worker = worker') so it is not
        garbage-collected.
    on_finished
        Called with the worker's success payload ('finished_ok').
    on_failed
        Called with the failure message + traceback ('failed').
        If omitted, errors are silently dropped; pass an explicit
        handler unless you have a specific reason not to.
    on_progress
        Called with each human-readable progress string.
    on_progress_pct
        Called with each 0-100 progress integer (e.g. for a progress
        bar).

    Returns
    -------
    The same worker, so call sites can chain assignments inline:

        self._worker = run_worker(MyWorker(...), on_finished=...)
    """
    if on_progress is not None:
        worker.progress.connect(on_progress)
    if on_progress_pct is not None:
        worker.progress_pct.connect(on_progress_pct)
    worker.finished_ok.connect(on_finished)
    if on_failed is not None:
        worker.failed.connect(on_failed)
    worker.start()
    return worker
