# Base class for Figure Export pages. Subclasses override
# 'dependencies_met()' (does the workspace have what this figure needs?)
# and '_make_render_func()' (closure returning a Figure or None).

from __future__ import annotations

import io
from pathlib import Path
from typing import Callable, Optional

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import QFileDialog, QLabel

from kosmic.gui.figure_export.shared.preview_pane import PreviewPane
from kosmic.gui.shared import run_worker
from kosmic.gui.shared.theme import get_export_dpi, get_font_sizes
from kosmic.gui.shared.widgets import (
    BaseWorker, SecondaryButton, SecondaryLabel, SidebarPage,
)


def _fig_to_png_bytes(fig, dpi: int = 150) -> bytes:
    """Serialise a matplotlib Figure to PNG bytes and close it."""
    import matplotlib.pyplot as plt

    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight',
                facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close(fig)
    buf.seek(0)
    return buf.getvalue()


def _save_png_pdf(fig, png_path: Path) -> Path:
    """Save fig to png_path and a sibling .pdf at the user-configured DPI.

    facecolor is preserved so dark-mode previews export accurately.
    Closes the figure.
    """
    import matplotlib.pyplot as plt

    dpi = get_export_dpi()
    fig.savefig(png_path, dpi=dpi, bbox_inches='tight',
                facecolor=fig.get_facecolor(), edgecolor='none')
    if png_path.suffix.lower() == '.png':
        # Vector PDF with editable text/axes; any rasterized artist (e.g. a
        # dense UMAP scatter) is rendered crisply, not at the ~100 dpi default.
        fig.savefig(png_path.with_suffix('.pdf'), dpi=max(dpi, 400),
                    bbox_inches='tight',
                    facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close(fig)
    return png_path


class RenderWorker(BaseWorker):
    """
    Run 'func()' (returning a Figure or None) on a worker thread;
    emit 'finished_ok' with PNG bytes (b'' when 'func' returns None).
    """

    def __init__(self, func, parent=None):
        super().__init__(parent)
        self._func = func

    def _run(self):
        fig = self._func()
        if fig is None:
            return b''
        return _fig_to_png_bytes(fig)


class FigurePage(SidebarPage):
    """
    Abstract base for one figure type.

    Subclasses must override 'TITLE', '_build_controls()',
    'dependencies_met()', and '_make_render_func()'.
    """

    # SidebarPage geometry overrides for figure pages: figure controls
    # are short, so no scroll; margins slightly larger to suit the
    # title-on-top layout.
    SIDEBAR_SCROLLABLE = False
    SIDEBAR_MARGINS = (12, 12, 12, 12)

    TITLE = "Figure"
    UNAVAILABLE_MESSAGE = "This figure is not available in the current state."

    # DE-driven pages set this True to show a banner when the loaded DE / pathway
    # results predate the current data (provenance upstream-token mismatch).
    CHECKS_DE_STALENESS = False

    log_message = pyqtSignal(str)

    def __init__(self, workspace):
        super().__init__()
        self.workspace = workspace
        self._render_worker: Optional[RenderWorker] = None
        self._cached_pixmap: Optional[QPixmap] = None
        self._cached_data_version: int = -1

        self._setup_ui()

    # -- UI scaffolding -------------------------------------------------

    def _setup_ui(self):
        # Sidebar: title + controls + export + status (uses SidebarPage's
        # self.sidebar_layout).
        title_lbl = QLabel(self.TITLE)
        title_lbl.setProperty("role", "section_title")
        self.sidebar_layout.addWidget(title_lbl)

        self._controls = self._build_controls()
        if self._controls is not None:
            self._controls.changed.connect(self._invalidate_and_render)
            self.sidebar_layout.addWidget(self._controls)

        self.sidebar_layout.addStretch()

        self._export_btn = SecondaryButton("Export PNG + PDF")
        self._export_btn.clicked.connect(self._export_current)
        self.sidebar_layout.addWidget(self._export_btn)

        self._status_lbl = SecondaryLabel("")
        self._status_lbl.setWordWrap(True)
        self.sidebar_layout.addWidget(self._status_lbl)

        # Content: an optional staleness banner above the preview pane.
        self._staleness_banner = SecondaryLabel("")
        self._staleness_banner.setWordWrap(True)
        self._staleness_banner.setProperty("role", "status_warning")
        self._staleness_banner.setVisible(False)
        self.content_layout.addWidget(self._staleness_banner)

        self._preview = PreviewPane(self.UNAVAILABLE_MESSAGE)
        self.content_layout.addWidget(self._preview)

    def _build_controls(self):
        """
        Override to return a 'FigureControls' subclass instance.

        Default: no controls.
        """
        return None

    # -- Lifecycle ------------------------------------------------------

    def on_activated(self):
        """Called by the workspace when this page becomes the visible one."""
        if not self.dependencies_met():
            self._staleness_banner.setVisible(False)
            self._preview.set_message(self.UNAVAILABLE_MESSAGE)
            self._export_btn.setEnabled(False)
            return
        self._update_staleness()
        self._export_btn.setEnabled(True)
        version = self.workspace.data_version()
        if (self._cached_pixmap is not None
                and version == self._cached_data_version):
            self._preview.set_pixmap(self._cached_pixmap)
            return
        self._render()

    def refresh_theme(self):
        """Theme changed -- drop the cached pixmap and re-render if visible."""
        self._cached_pixmap = None
        if self.isVisible() and self.dependencies_met():
            self._render()

    def invalidate(self):
        """Drop the cached pixmap (call when underlying data changes)."""
        self._cached_pixmap = None
        self._cached_data_version = -1

    # -- Rendering ------------------------------------------------------

    def _update_staleness(self):
        """Show/hide the 'results out of date' banner for DE-driven pages."""
        note = None
        if self.CHECKS_DE_STALENESS:
            from kosmic.gui.figure_export.shared.staleness import (
                de_results_stale_note,
            )
            note = de_results_stale_note(getattr(self.workspace, 'de_ws', None))
        self._staleness_banner.setText(note or "")
        self._staleness_banner.setVisible(bool(note))

    def _invalidate_and_render(self):
        self._cached_pixmap = None
        self._render()

    def _render(self):
        func = self._make_render_func()
        if func is None:
            self._preview.set_message(self.UNAVAILABLE_MESSAGE)
            self._export_btn.setEnabled(False)
            return

        if self._render_worker is not None and self._render_worker.isRunning():
            self._render_worker.requestInterruption()
            self._render_worker.wait(50)

        self._status_lbl.setText("Rendering...")
        self._render_worker = RenderWorker(func, parent=self)
        run_worker(
            self._render_worker,
            on_finished=self._on_render_finished,
            on_failed=self._on_render_failed,
        )

    def _on_render_finished(self, png_bytes: bytes):
        self._render_worker = None
        if not png_bytes:
            self._preview.set_message("Nothing to render.")
            self._status_lbl.setText("Nothing to render.")
            return
        from PyQt6.QtGui import QImage
        img = QImage()
        img.loadFromData(png_bytes)
        pm = QPixmap.fromImage(img)
        self._cached_pixmap = pm
        self._cached_data_version = self.workspace.data_version()
        self._preview.set_pixmap(pm)
        self._status_lbl.setText("Ready.")

    def _on_render_failed(self, message: str):
        self._render_worker = None
        first_line = message.splitlines()[0] if message else "Render failed."
        self._preview.set_message(f"Render failed: {first_line}")
        self._status_lbl.setText(f"Render failed: {first_line}")
        log = getattr(self.workspace, 'log_message', None)
        if log is not None:
            log.emit(f"Figure render failed:\n{message}")

    # -- Export ---------------------------------------------------------

    def _default_export_dir(self) -> Optional[Path]:
        """
        Override to provide a default export folder for this page.

        Returns None to fall back to a folder picker.
        """
        return None

    def _default_export_basename(self) -> str:
        return self.TITLE.lower().replace(' ', '_').replace('/', '_')

    def _export_current(self):
        func = self._make_render_func()
        if func is None:
            self._status_lbl.setText("Nothing to export.")
            return

        default_dir = self._default_export_dir()
        if default_dir:
            default_dir.mkdir(parents=True, exist_ok=True)
            default_path = str(default_dir / f"{self._default_export_basename()}.png")
        else:
            default_path = f"{self._default_export_basename()}.png"

        path, _ = QFileDialog.getSaveFileName(
            self, f"Export {self.TITLE}", default_path,
            "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)",
        )
        if not path:
            return

        try:
            fig = func()
            if fig is None:
                self._status_lbl.setText("Render returned no figure; export aborted.")
                return
            saved = _save_png_pdf(fig, Path(path))
            self._status_lbl.setText(f"Exported to {saved}")
            self.log_message.emit(f"Figure exported: {saved}")
        except Exception as e:
            self._status_lbl.setText(f"Export failed: {e}")

    # -- Hooks ----------------------------------------------------------

    def dependencies_met(self) -> bool:
        """Override: return True iff the workspace state can drive a render."""
        return False

    def _make_render_func(self) -> Optional[Callable]:
        """Override: return a closure that renders a Figure, or None."""
        return None

    # -- Helpers --------------------------------------------------------

    @staticmethod
    def _font_sizes():
        return get_font_sizes()
