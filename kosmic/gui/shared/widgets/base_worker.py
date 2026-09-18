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


def _memory_gb():
    """(resident, peak resident) of this process in GB, or None.

    psutil when present; otherwise the operating system directly
    (/proc on Linux, GetProcessMemoryInfo on Windows, rusage on macOS),
    so the readout does not depend on an optional package.
    """
    try:
        import psutil
        info = psutil.Process().memory_info()
        peak = getattr(info, 'peak_wset', None)
        return info.rss / 1e9, (peak / 1e9 if peak else None)
    except Exception:
        pass
    try:
        import sys
        if sys.platform.startswith('linux'):
            rss = peak = None
            with open('/proc/self/status') as f:
                for line in f:
                    if line.startswith('VmRSS:'):
                        rss = int(line.split()[1]) * 1024
                    elif line.startswith('VmHWM:'):
                        peak = int(line.split()[1]) * 1024
            if rss is not None:
                return rss / 1e9, (peak / 1e9 if peak else None)
        elif sys.platform == 'win32':
            import ctypes
            from ctypes import wintypes

            class _PMC(ctypes.Structure):
                _fields_ = [('cb', wintypes.DWORD), ('PageFaultCount', wintypes.DWORD),
                            ('PeakWorkingSetSize', ctypes.c_size_t),
                            ('WorkingSetSize', ctypes.c_size_t),
                            ('QuotaPeakPagedPoolUsage', ctypes.c_size_t),
                            ('QuotaPagedPoolUsage', ctypes.c_size_t),
                            ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t),
                            ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                            ('PagefileUsage', ctypes.c_size_t),
                            ('PeakPagefileUsage', ctypes.c_size_t)]
            pmc = _PMC()
            pmc.cb = ctypes.sizeof(_PMC)
            kernel32 = ctypes.windll.kernel32
            kernel32.GetCurrentProcess.restype = wintypes.HANDLE
            psapi = ctypes.windll.psapi
            psapi.GetProcessMemoryInfo.argtypes = [
                wintypes.HANDLE, ctypes.POINTER(_PMC), wintypes.DWORD]
            psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
            if psapi.GetProcessMemoryInfo(kernel32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb):
                return pmc.WorkingSetSize / 1e9, pmc.PeakWorkingSetSize / 1e9
        else:
            import resource
            peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss  # bytes on macOS
            return peak / 1e9, peak / 1e9
    except Exception:
        pass
    return None


def memory_report(step: str, before, after=None) -> str:
    """One line saying what a step cost: resident before, after, and the
    process peak. The peak is per process, not per step, so it only moves
    when a step sets a new high; the step's own cost is the difference
    between after and before plus whatever transient the peak shows.

    Every worker records these, and run_worker writes the line to the
    output panel when the step completes, so the panel always says what
    is in memory. Nothing else in KOSMIC did.
    """
    if after is None:
        after = _memory_gb()
    if before is None or after is None:
        return f"{step}: memory not available"
    text = f"{step}: memory {before[0]:.1f} GB before, {after[0]:.1f} GB after"
    if after[1] is not None:
        text += f", process peak {after[1]:.1f} GB"
    return text


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
        # Memory before and after the step is recorded on the worker; the
        # owner reads it on the GUI thread once the step has completed
        # (run_worker does so). Nothing is emitted from this thread that
        # the completion handlers did not already expect.
        self.memory_before = _memory_gb()
        try:
            result = self._run()
        except Exception as exc:  # noqa: BLE001 -- by design
            self.memory_after = _memory_gb()
            self.failed.emit(f"{exc}\n{traceback.format_exc()}")
        else:
            self.memory_after = _memory_gb()
            self.finished_ok.emit(result)
        finally:
            self._release_inputs()

    def memory_line(self) -> str:
        """What this step cost, for the output panel; see memory_report."""
        return memory_report(type(self).__name__,
                             getattr(self, 'memory_before', None),
                             getattr(self, 'memory_after', None))

    def _release_inputs(self) -> None:
        """Drop references to any AnnData the worker was given.

        Owners keep their last worker as an attribute so the QThread
        stays alive, and the worker keeps the dataset it was handed. Once
        the result has been emitted that dataset is dead weight: after a
        subset the parent (16 GB on Reichart) stayed in memory because the
        finished FilterWorker still pointed at it. The result object is
        the caller's; only the inputs are released here.
        """
        try:
            import anndata as ad
        except ImportError:  # pragma: no cover
            return
        for name, value in list(vars(self).items()):
            if isinstance(value, ad.AnnData):
                setattr(self, name, None)
