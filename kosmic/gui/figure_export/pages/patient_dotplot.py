# Patient pathway dotplot page.
#
# A matplotlib mirror of the DE window's per-donor pathway score plot. It never
# recomputes scores: it renders the in-memory scores from the last scoring run
# ('de_ws.pathway_score_data') and falls back to the per-donor scores saved to
# disk by that run, so the figure always reflects the latest actual computation
# (whatever scoring method was used) rather than a recomputed approximation.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import QCheckBox

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class _PatientDotplotControls(FigureControls):
    def __init__(self):
        # Default to Auto sizing -- the grid height scales with the pathway
        # count, so a fixed size crops a large panel.
        super().__init__()
        self.standardize = QCheckBox("Standardise per pathway (z-score)")
        self.standardize.setChecked(True)
        self.standardize.setToolTip(
            "Z-score each pathway's per-donor scores (mean 0, SD 1) so all "
            "pathways share a comparable scale. Off shows the raw per-donor "
            "scores in the scoring method's own units.")
        self.standardize.toggled.connect(self.changed)
        self.add_row("", self.standardize)


class PatientDotplotPage(FigurePage):
    TITLE = "Patient Dotplot"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run pathway scoring in the DE window first -- this figure mirrors those "
        "per-donor scores."
    )

    def _build_controls(self) -> FigureControls:
        return _PatientDotplotControls()

    def _scores_path(self) -> Optional[Path]:
        """On-disk per-donor scores saved by the last scoring run, or None."""
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_pathway_scoring_dir
        from kosmic.de.pathway_scoring import PER_DONOR_SCORES_FILE
        return de_pathway_scoring_dir(de.project_dir) / PER_DONOR_SCORES_FILE

    def _score_data(self):
        """Per-donor scores to render: the in-memory scores from this session's
        scoring run, else the ones saved to disk by the last run. None if
        neither exists (scoring never run)."""
        de = self.workspace.de_ws
        if de is None:
            return None
        cache = getattr(de, 'pathway_score_data', None)
        if cache and cache.get('names'):
            return cache
        from kosmic.de.pathway_scoring import load_per_donor_scores
        path = self._scores_path()
        return load_per_donor_scores(path) if path else None

    def dependencies_met(self) -> bool:
        de = self.workspace.de_ws
        if de is None:
            return False
        if not (getattr(de, 'disease_label', None)
                and getattr(de, 'control_label', None)):
            return False
        return self._score_data() is not None

    def _make_render_func(self):
        data = self._score_data()
        if data is None:
            return None
        de = self.workspace.de_ws
        d, c = de.disease_label, de.control_label
        standardize = self._controls.standardize.isChecked()
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()
        from kosmic.gui.shared.theme import get_color
        theme_colors = {'plot_control': get_color('plot_control'),
                        'plot_disease': get_color('plot_disease')}

        def _render(data=data, d=d, c=c, std=standardize,
                    f=font_sizes, sz=figsize, tc=theme_colors):
            from kosmic.visualisation.de.dotplots import (
                pathway_scores_to_sample_df, create_patient_dotplot,
            )
            sdf, pn = pathway_scores_to_sample_df(data, standardize=std)
            ylabel = ("Pathway score (z-scored)" if std
                      else "Pathway score (per donor)")
            fig, _ = create_patient_dotplot(
                sdf, pn, d, c, font_sizes=f, figsize=sz, ylabel=ylabel,
                shared_y=std, theme_colors=tc)
            return fig
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        de = self.workspace.de_ws
        if de is None or not getattr(de, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(de.project_dir)

    def _default_export_basename(self) -> str:
        return "patient_dotplot"
