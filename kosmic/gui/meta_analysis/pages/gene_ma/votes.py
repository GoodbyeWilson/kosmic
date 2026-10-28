# Votes tab: genes by the number of studies agreeing (0/k ... k/k, and
# conflict), split up / down, with the run's other method's significant
# genes overlaid. Drawn by kosmic.visualisation.meta.vote_levels, the
# same figure as Figure Export's "Study Agreement" page.
#
# Parent drives via 'populate' after each run.
from __future__ import annotations

from typing import Optional

import pandas as pd
from PyQt6.QtWidgets import QVBoxLayout, QWidget

from kosmic.gui.shared.widgets import SecondaryLabel


class VotesTab(QWidget):
    """Study-agreement bar chart for the current pooling run."""

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._lay = QVBoxLayout(self)
        self._lay.setContentsMargins(4, 4, 4, 4)
        self.status_label = SecondaryLabel(
            "Run a meta-analysis to see how many studies agree on each gene.")
        self.status_label.setWordWrap(True)
        self._lay.addWidget(self.status_label)
        self._canvas = None
        self.levels: Optional[pd.DataFrame] = None

    def populate(self, meta_df: Optional[pd.DataFrame], method_keys,
                 n_studies: int) -> None:
        from kosmic.meta_analysis.pooling.vote_count import run_vote_levels

        self._clear_canvas()
        self.levels = None
        if (meta_df is None or meta_df.empty or n_studies < 1
                or not {'n_up', 'n_down'} <= set(meta_df.columns)):
            self.status_label.setText(
                "No study-level calls to count: the studies have no adjusted "
                "p-values, or no run yet.")
            return
        levels, sig_label = run_vote_levels(meta_df, method_keys, n_studies)
        self.levels = levels
        self.status_label.setText(
            "Genes by the number of studies that call them significant "
            "(study FDR < 0.05) in the same direction."
            + (f" Hatched: genes with {sig_label}." if sig_label else ""))

        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

        from kosmic.visualisation.meta.vote_levels import create_vote_levels_plot
        fig = create_vote_levels_plot(levels, title='Study agreement',
                                      significant_label=sig_label)
        if fig is None:
            return
        self._canvas = FigureCanvasQTAgg(fig)
        self._lay.addWidget(self._canvas, 1)

    def clear(self) -> None:
        self._clear_canvas()
        self.levels = None

    def _clear_canvas(self) -> None:
        if self._canvas is not None:
            import matplotlib.pyplot as plt
            plt.close(self._canvas.figure)
            self._lay.removeWidget(self._canvas)
            self._canvas.deleteLater()
            self._canvas = None


__all__ = ["VotesTab"]
