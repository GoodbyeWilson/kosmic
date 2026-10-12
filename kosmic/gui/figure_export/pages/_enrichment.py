# Base for figure pages drawn from an enrichment run, which can come from
# the DE workspace (one study's DE) or the Meta-Analysis workspace's
# Enrichment step (the pooled gene list, via figure_inputs(); ADR-010).
#
# Each source is presented with the attribute names the DE workspace uses
# ('enrichment_results', 'gsea_library_results', 'gsea_library_details',
# 'gsea_library_name', 'de_results', 'disease_label', 'control_label'), so
# a page reads 'self._source()' where it used to read 'workspace.de_ws'.
# A Source selector appears when both workspaces hold a run; Meta-analysis
# is the default.

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Optional

from PyQt6.QtWidgets import QComboBox

from kosmic.gui.figure_export.pages._base import FigurePage

SOURCE_META, SOURCE_STUDY = "Meta-analysis", "Study DE"


class EnrichmentFigurePage(FigurePage):
    """A figure page whose data is one enrichment run."""

    #: Which kind of run the page draws: 'ora' (GO / Enrichr
    #: over-representation) or 'gsea' (preranked GSEA of a library).
    KIND = 'ora'

    def _add_source_row(self, ctrl):
        """Add the Source selector to *ctrl*; hidden until both sources hold
        a run."""
        ctrl.source = QComboBox()
        ctrl.source.currentIndexChanged.connect(ctrl.changed)
        ctrl.add_row("Source:", ctrl.source)
        self._set_source_row_visible(ctrl, False)
        return ctrl

    @staticmethod
    def _set_source_row_visible(ctrl, visible):
        ctrl.source.setVisible(visible)
        label = ctrl._form.labelForField(ctrl.source)
        if label is not None:
            label.setVisible(visible)

    # -- sources ------------------------------------------------------------

    def _meta_source(self):
        meta = getattr(self.workspace, 'meta_ws', None)
        inputs = meta.figure_inputs() if meta is not None else None
        if not inputs:
            return None
        sel = str(inputs.get('selection') or '').replace('_', ' ')
        return SimpleNamespace(
            enrichment_results=inputs.get('enrichment_results'),
            gsea_library_results=inputs.get('gsea_library_results'),
            gsea_library_details=inputs.get('gsea_library_details'),
            gsea_library_name=inputs.get('gsea_library_name'),
            de_results=inputs.get('enrichment_ranking'),
            disease_label='Disease', control_label='Control',
            contrast=f"{sel} (pooled)" if sel else 'Pooled',
            project_folder=inputs.get('project_folder'),
            selection=inputs.get('selection'),
        )

    def _has_run(self, src, is_meta: bool) -> bool:
        """Whether *src* holds a run of this page's KIND."""
        if src is None:
            return False
        if self.KIND == 'gsea':
            res = getattr(src, 'gsea_library_results', None)
        else:
            res = getattr(src, 'enrichment_results', None)
        return res is not None and len(res) > 0

    def _sources(self) -> dict:
        """{name: source} for the sources holding a run, Meta-analysis first."""
        out = {}
        meta = self._meta_source()
        if self._has_run(meta, True):
            out[SOURCE_META] = meta
        de = self.workspace.de_ws
        if self._has_run(de, False):
            out[SOURCE_STUDY] = de
        return out

    def _source_name(self) -> Optional[str]:
        sources = self._sources()
        chosen = self._controls.source.currentText()
        if chosen in sources:
            return chosen
        return next(iter(sources), None)

    def _source(self):
        """The selected source, or None when neither holds a run."""
        name = self._source_name()
        return self._sources().get(name) if name else None

    def _is_meta(self) -> bool:
        return self._source_name() == SOURCE_META

    def _contrast(self, src) -> Optional[str]:
        """'DCM vs NF' for a study, '<selection> (pooled)' for the
        meta-analysis; None when the study's labels are unknown."""
        if getattr(src, 'contrast', None):
            return src.contrast
        d, c = getattr(src, 'disease_label', None), getattr(src, 'control_label', None)
        return f"{d} vs {c}" if d and c else None

    def _sync_sources(self):
        """Refill the Source selector, keeping the choice while it is
        available, and show it only when there is a choice to make."""
        combo = self._controls.source
        names = list(self._sources())
        current = combo.currentText()
        combo.blockSignals(True)
        combo.clear()
        combo.addItems(names)
        if current in names:
            combo.setCurrentText(current)
        combo.blockSignals(False)
        self._set_source_row_visible(self._controls, len(names) > 1)

    # -- lifecycle -----------------------------------------------------------

    def on_activated(self):
        self._sync_sources()
        super().on_activated()

    def dependencies_met(self) -> bool:
        return self._source() is not None

    def _update_staleness(self):
        # The staleness check compares DE results with the current data;
        # it says nothing about a meta-analysis run.
        if self._is_meta():
            self._staleness_banner.setVisible(False)
            return
        super()._update_staleness()

    def _default_export_dir(self) -> Optional[Path]:
        src = self._source()
        if self._is_meta():
            if not getattr(src, 'project_folder', None):
                return None
            from kosmic.paths import meta_output_dir
            return meta_output_dir(src.project_folder, src.selection)
        if src is None or not getattr(src, 'project_dir', None):
            return None
        from kosmic.paths import de_figures_dir
        return de_figures_dir(src.project_dir)
