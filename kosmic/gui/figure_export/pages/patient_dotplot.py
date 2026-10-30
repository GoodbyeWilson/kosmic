# Patient pathway dotplot page.
#
# A matplotlib mirror of the DE window's per-donor pathway score plot. It never
# recomputes scores: it renders the per-donor scores of the last scoring run,
# from the DE workspace or as saved in the study's folder (FigurePage._de),
# so the figure always reflects the latest actual computation (whatever
# scoring method was used) rather than a recomputed approximation.

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

    STALENESS_ATTR = 'pathway_score_data'

    def _score_data(self):
        """Per-donor scores to render, or None if scoring was never run."""
        de = self._de()
        data = de.pathway_score_data if de is not None else None
        return data if data and data.get('names') else None

    def dependencies_met(self) -> bool:
        de = self._de()
        return (de is not None and de.available('disease_label')
                and de.available('control_label')
                and de.available('pathway_score_data'))

    def _make_render_func(self):
        data = self._score_data()
        if data is None:
            return None
        de = self._de()
        d, c = de.disease_label, de.control_label
        standardize = self._controls.standardize.isChecked()
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()
        # The Pathway DE FDR for each pathway, shown above its panel.
        fdr = {}
        pdr = getattr(de, 'pathway_de_results', None)
        if pdr is not None and {'names', 'pvals_adj'} <= set(pdr.columns):
            fdr = dict(zip(pdr['names'], pdr['pvals_adj'].astype(float)))

        def _render(data=data, d=d, c=c, std=standardize,
                    f=font_sizes, sz=figsize, fdr=fdr):
            from kosmic.visualisation.de.dotplots import (
                pathway_scores_to_sample_df, create_patient_dotplot,
            )
            sdf, pn = pathway_scores_to_sample_df(data, standardize=std)
            ylabel = ("Pathway score (z-scored)" if std
                      else "Pathway score (per donor)")
            return create_patient_dotplot(
                sdf, pn, d, c, font_sizes=f, figsize=sz, ylabel=ylabel,
                shared_y=std, fdr_by_pathway=fdr)
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        return self._de_export_dir()

    def _default_export_basename(self) -> str:
        return "patient_dotplot"
