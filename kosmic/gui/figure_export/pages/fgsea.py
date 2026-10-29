# Preranked GSEA (fgsea) enrichment page.
#
# A matplotlib mirror of the DE window's NES bar chart. It never recomputes:
# it renders the fgsea summary of the last run, from the DE workspace or as
# saved in the study's folder (FigurePage._de), so the figure always
# reflects the latest actual GSEA (against the genome-wide gene ranking)
# rather than a recomputed approximation.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls


class FgseaPage(FigurePage):
    TITLE = "Pathway Enrichment (GSEA)"
    CHECKS_DE_STALENESS = True
    UNAVAILABLE_MESSAGE = (
        "Run preranked GSEA (fgsea) on the Pathway DE page first -- this figure "
        "mirrors those enrichment results."
    )

    def _build_controls(self) -> FigureControls:
        return FigureControls()

    STALENESS_ATTR = 'fgsea_results'

    def _results(self):
        """fgsea summary to render, or None if fgsea was never run."""
        de = self._de()
        res = de.fgsea_results if de is not None else None
        return res if res is not None and len(res) > 0 else None

    def dependencies_met(self) -> bool:
        de = self._de()
        return de is not None and de.available('fgsea_results')

    def _make_render_func(self):
        res = self._results()
        if res is None:
            return None
        de = self._de()
        d = getattr(de, 'disease_label', None)
        c = getattr(de, 'control_label', None)
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()

        def _render(res=res, d=d, c=c, f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.fgsea_plot import create_fgsea_nes_plot
            return create_fgsea_nes_plot(
                res, disease_label=d, control_label=c,
                figsize=sz, font_sizes=f)
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        return self._de_export_dir()

    def _default_export_basename(self) -> str:
        return "pathway_enrichment_gsea"
