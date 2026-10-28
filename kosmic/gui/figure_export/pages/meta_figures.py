# Meta-analysis figure pages: pooled volcano, per-study forest and
# leave-one-out sensitivity, for gene- or pathway-level results, and the
# study agreement (vote levels) of a gene-level run.
#
# Read the Meta-Analysis workspace through 'workspace.meta_ws.figure_inputs()'
# (ADR-010): the pooled tables of the latest runs and the studies' own DE
# tables for the current selection.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import QComboBox, QDoubleSpinBox, QSpinBox

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls

LEVEL_GENE, LEVEL_PATHWAY = "Gene", "Pathway"


class MetaFigurePage(FigurePage):
    """Base for pages drawn from a meta-analysis result at one LEVEL."""

    LEVEL = 'gene'   # 'gene' or 'pathway'

    @property
    def UNAVAILABLE_MESSAGE(self):
        return (f"Run a {self._level()}-level meta-analysis in the "
                "Meta-Analysis workspace to view this figure.")

    def _level(self) -> str:
        return self.LEVEL

    def _inputs(self) -> Optional[dict]:
        meta = getattr(self.workspace, 'meta_ws', None)
        return meta.figure_inputs() if meta is not None else None

    def _results(self):
        inputs = self._inputs()
        df = inputs.get(f'{self._level()}_results') if inputs else None
        return df if df is not None and len(df) else None

    def _studies(self) -> list:
        inputs = self._inputs()
        return inputs.get(f'{self._level()}_studies', []) if inputs else []

    def _selection_label(self) -> str:
        inputs = self._inputs() or {}
        sel = str(inputs.get('selection') or '').replace('_', ' ')
        n = len(self._studies())
        return f"{sel} ({n} studies)" if sel else f"{n} studies"

    def dependencies_met(self) -> bool:
        return self._results() is not None

    def _default_export_dir(self) -> Optional[Path]:
        inputs = self._inputs() or {}
        if not inputs.get('project_folder'):
            return None
        from kosmic.paths import meta_output_dir
        return meta_output_dir(inputs['project_folder'], inputs.get('selection'))


# -- Volcano -----------------------------------------------------------------

class _MetaVolcanoControls(FigureControls):
    def __init__(self, max_labels: int, point_size: int):
        super().__init__(8, 8)
        self.max_labels = QSpinBox()
        self.max_labels.setRange(0, 60)
        self.max_labels.setValue(max_labels)
        self.max_labels.valueChanged.connect(self.changed)
        self.add_row("Max labels:", self.max_labels)

        self.point_size = QSpinBox()
        self.point_size.setRange(1, 30)
        self.point_size.setValue(point_size)
        self.point_size.valueChanged.connect(self.changed)
        self.add_row("Point size:", self.point_size)

        # Off: significant on FDR alone, as the Meta-Analysis workspace counts.
        self.add_lfc_gate("pooled log2 fold change")


class MetaVolcanoPage(MetaFigurePage):
    TITLE = "Gene Meta-Analysis Volcano"
    HEADLINE = "Meta-Analysed Differential Gene Expression"

    def _build_controls(self):
        return _MetaVolcanoControls(max_labels=10, point_size=3)

    def _make_render_func(self):
        df = self._results()
        if df is None:
            return None
        dr = df[['names', 'logfoldchanges', 'fdr']].rename(columns={'fdr': 'pvals_adj'})
        ctrl = self._controls
        ml, ps = ctrl.max_labels.value(), ctrl.point_size.value() ** 2
        gate = ctrl.lfc_gate.value()
        sub, headline = self._selection_label(), self.HEADLINE
        figsize = ctrl.figsize.get_figsize() or (8, 8)
        font_sizes = self._font_sizes()

        def _render(dr=dr, ml=ml, ps=ps, gate=gate, sub=sub, headline=headline,
                    sz=figsize, f=font_sizes):
            from kosmic.visualisation.de.volcano import create_volcano_plot
            return create_volcano_plot(
                dr, '', '', title=headline, subtitle=sub,
                xlabel='Pooled Log$_2$ fold change', max_labels=ml,
                logfc_threshold=gate,
                point_size=ps, figsize=sz, font_sizes=f)
        return _render

    def _default_export_basename(self) -> str:
        return f"{self._level()}_meta_volcano"


class PathwayMetaVolcanoPage(MetaVolcanoPage):
    TITLE = "Pathway Meta-Analysis Volcano"
    HEADLINE = "Meta-Analysed Pathway Activity"
    LEVEL = 'pathway'

    def _build_controls(self):
        return _MetaVolcanoControls(max_labels=20, point_size=6)


# -- Forest and leave-one-out --------------------------------------------------

class _ItemControls(FigureControls):
    """Pick one gene / pathway (most significant first) and the x range."""

    def __init__(self, with_level: bool = False):
        super().__init__()
        self.level = None
        if with_level:
            self.level = QComboBox()
            self.level.addItems([LEVEL_GENE, LEVEL_PATHWAY])
            self.add_row("Level:", self.level)

        self.item = QComboBox()
        self.item.setEditable(True)
        self.item.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self.item.currentIndexChanged.connect(self.changed)
        self.add_row("Gene / pathway:", self.item)

        # 0 = symmetric range fitted to the data.
        self.x_limit = QDoubleSpinBox()
        self.x_limit.setRange(0.0, 20.0)
        self.x_limit.setSingleStep(0.5)
        self.x_limit.setDecimals(1)
        self.x_limit.setSpecialValueText("Auto")
        self.x_limit.valueChanged.connect(self.changed)
        self.add_row("X limit (+/- log2FC):", self.x_limit)

    def populate_items(self, names: list[str]):
        prev = self.item.currentText()
        self.item.blockSignals(True)
        self.item.clear()
        self.item.addItems(names)
        idx = self.item.findText(prev)
        self.item.setCurrentIndex(idx if idx >= 0 else 0)
        self.item.blockSignals(False)


class _ItemPage(MetaFigurePage):
    """Shared item handling for the forest and leave-one-out pages."""

    def on_activated(self):
        df = self._results()
        if df is not None:
            order = df.sort_values('fdr')['names'].astype(str).tolist()
            self._controls.populate_items(order)
        super().on_activated()

    def _item_row(self):
        df = self._results()
        item = self._controls.item.currentText()
        if df is None or not item:
            return None, None
        hit = df[df['names'].astype(str) == item]
        return item, (hit.iloc[0] if len(hit) else None)

    def _display_name(self, item: str) -> str:
        return item.replace('_', ' ')

    def _write_report(self, fig_path: Path, text: str) -> list[Path]:
        out = Path(fig_path).with_suffix('.txt')
        out.write_text(text, encoding='utf-8')
        return [out]


class StudyForestPage(_ItemPage):
    TITLE = "Gene Forest Plot"

    def _build_controls(self):
        return _ItemControls()

    def _make_render_func(self):
        item, row = self._item_row()
        if row is None:
            return None
        studies = self._studies()
        name = self._display_name(item)
        xl = self._controls.x_limit.value() or None
        figsize = self._controls.figsize.get_figsize()
        font_sizes = self._font_sizes()

        def _render(item=item, row=row, studies=studies, name=name, xl=xl,
                    sz=figsize, f=font_sizes):
            from kosmic.visualisation.meta.study_forest import (
                create_study_forest, study_effects,
            )
            return create_study_forest(name, study_effects(item, studies), row,
                                       figsize=sz, font_sizes=f, x_limit=xl)
        return _render

    def _export_companions(self, fig_path: Path) -> list[Path]:
        item, row = self._item_row()
        if row is None:
            return []
        from kosmic.visualisation.meta.study_forest import forest_report, study_effects
        text = forest_report(self._display_name(item),
                             study_effects(item, self._studies()), row)
        return self._write_report(fig_path, text)

    def _default_export_basename(self) -> str:
        return f"forest_{self._controls.item.currentText()}"


class PathwayForestPage(StudyForestPage):
    TITLE = "Pathway Forest Plot"
    LEVEL = 'pathway'


class LeaveOneOutPage(_ItemPage):
    TITLE = "Leave-One-Out Sensitivity"

    def _build_controls(self):
        ctrl = _ItemControls(with_level=True)
        ctrl.level.currentIndexChanged.connect(self._on_level_changed)
        return ctrl

    def _level(self) -> str:
        ctrl = getattr(self, '_controls', None)   # unset while the page is built
        level = ctrl.level.currentText() if ctrl is not None else LEVEL_GENE
        return 'pathway' if level == LEVEL_PATHWAY else 'gene'

    def _on_level_changed(self):
        self.on_activated()
        self._invalidate_and_render()

    def _make_render_func(self):
        item, row = self._item_row()
        if row is None:
            return None
        studies = self._studies()
        name = self._display_name(item)
        xl = self._controls.x_limit.value() or None
        figsize = self._controls.figsize.get_figsize()
        font_sizes = self._font_sizes()

        def _render(item=item, row=row, studies=studies, name=name, xl=xl,
                    sz=figsize, f=font_sizes):
            from kosmic.visualisation.meta.study_forest import (
                create_loo_plot, loo_estimates,
            )
            return create_loo_plot(name, loo_estimates(item, studies, row), row,
                                   figsize=sz, font_sizes=f, x_limit=xl)
        return _render

    def _export_companions(self, fig_path: Path) -> list[Path]:
        item, row = self._item_row()
        if row is None:
            return []
        from kosmic.visualisation.meta.study_forest import loo_estimates, loo_report
        text = loo_report(self._display_name(item),
                          loo_estimates(item, self._studies(), row), row)
        return self._write_report(fig_path, text)

    def _default_export_basename(self) -> str:
        return f"leave_one_out_{self._controls.item.currentText()}"


# -- Study agreement (vote levels) ---------------------------------------------

_NO_OVERLAY = "None"


class _VoteLevelsControls(FigureControls):
    def __init__(self):
        super().__init__(7, 5)
        self.overlay = QComboBox()
        self.overlay.setToolTip(
            "Hatch, at each level, the genes this method of the run calls "
            "significant (FDR < 0.05).")
        self.overlay.currentIndexChanged.connect(self.changed)
        self.add_row("Overlay:", self.overlay)

    def populate_methods(self, keys: list[str]):
        prev = self.overlay.currentText()
        self.overlay.blockSignals(True)
        self.overlay.clear()
        self.overlay.addItems(keys + [_NO_OVERLAY])
        idx = self.overlay.findText(prev)
        self.overlay.setCurrentIndex(idx if idx >= 0 else 0)
        self.overlay.blockSignals(False)


class VoteLevelsPage(MetaFigurePage):
    """Genes by the number of studies agreeing (0/k ... k/k, conflict)."""

    TITLE = "Study Agreement (Vote Levels)"

    def _build_controls(self):
        return _VoteLevelsControls()

    def dependencies_met(self) -> bool:
        df = self._results()
        return df is not None and {'n_up', 'n_down'} <= set(df.columns)

    def _methods(self) -> list[str]:
        """The run's methods with an FDR column, other than vote counting."""
        df = self._results()
        if df is None:
            return []
        return [c[len('fdr_'):] for c in df.columns
                if c.startswith('fdr_') and c != 'fdr_vote_count']

    def on_activated(self):
        self._controls.populate_methods(self._methods())
        super().on_activated()

    def _levels(self):
        from kosmic.meta_analysis.pooling.vote_count import run_vote_levels
        df = self._results()
        method = self._controls.overlay.currentText()
        keys = [] if method in ('', _NO_OVERLAY) else [method]
        return run_vote_levels(df, keys, len(self._studies()))

    def _make_render_func(self):
        if not self.dependencies_met():
            return None
        levels, sig_label = self._levels()
        title = f"Study Agreement: {self._selection_label()}"
        figsize = self._controls.figsize.get_figsize()
        font_sizes = self._font_sizes()

        def _render(levels=levels, sig_label=sig_label, title=title,
                    sz=figsize, f=font_sizes):
            from kosmic.visualisation.meta.vote_levels import create_vote_levels_plot
            return create_vote_levels_plot(levels, title=title,
                                           significant_label=sig_label,
                                           figsize=sz, font_sizes=f)
        return _render

    def _export_companions(self, fig_path: Path) -> list[Path]:
        """The counts behind the bars."""
        if not self.dependencies_met():
            return []
        out = Path(fig_path).with_suffix('.csv')
        self._levels()[0].to_csv(out)
        return [out]

    def _default_export_basename(self) -> str:
        return "study_agreement_vote_levels"
