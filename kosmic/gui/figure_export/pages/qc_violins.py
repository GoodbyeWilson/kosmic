# Single-cell QC pages: QC Violins and Cell Filtering.
#
# QC Violins reads the obs of the study's scRNA dataset (FigurePage._adata):
# genes and counts per cell and mitochondrial %, optionally DecontX
# contamination and the doublet score, pooled or by condition or sample.
# Both pages read the QC limits and the cells removed at each step from the
# study's provenance record, which sits beside its h5ad.

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd
from PyQt6.QtWidgets import QCheckBox, QComboBox, QHBoxLayout, QLabel, QLineEdit, QWidget

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls

GROUP_POOLED, GROUP_CONDITION, GROUP_SAMPLE = "Pooled", "By condition", "By sample"
START_LOADED, START_FILTERED = "Loaded cells", "After the pre-analysis filter"
# Default figure height (inches) per grouping; the width is 12 for all.
DEFAULT_HEIGHT = {GROUP_POOLED: 6, GROUP_CONDITION: 6, GROUP_SAMPLE: 10}
# Role keys in obs['_role'] and their names on the figure, controls first.
ROLES = (('control', 'Control'), ('disease', 'Disease'))


def _provenance_record(page: FigurePage) -> Optional[dict]:
    """The study's provenance record (beside its h5ad), or None."""
    src = page._study_source()
    path = src.h5ad_path() if src is not None else None
    if path is None:
        return None
    from kosmic import provenance
    return provenance.load(path.parent)


class _QCViolinControls(FigureControls):
    def __init__(self):
        super().__init__(default_w=12, default_h=DEFAULT_HEIGHT[GROUP_POOLED])
        self.group = QComboBox()
        self.group.addItems([GROUP_POOLED, GROUP_CONDITION, GROUP_SAMPLE])
        self.group.setToolTip("Pooled: one violin per metric for all cells.\n"
                              "By condition: control and disease side by side.\n"
                              "By sample: one violin per sample, controls first,\n"
                              "with the metrics stacked in rows.")
        self._group_shown = GROUP_POOLED
        self.group.currentIndexChanged.connect(self._group_changed)
        self.group.currentIndexChanged.connect(self.changed)
        self.add_row("Violins:", self.group)

        self.log_counts = QCheckBox("Log scale for genes and counts")
        self.log_counts.setChecked(True)
        self.log_counts.toggled.connect(self.changed)
        self.add_row("", self.log_counts)

        self.limits = QCheckBox("Show the QC limits")
        self.limits.setChecked(True)
        self.limits.setToolTip("Dotted lines at the limits of the study's last\n"
                               "recorded QC step.")
        self.limits.toggled.connect(self.changed)
        self.add_row("", self.limits)

        self.extra = QCheckBox("Add DecontX and doublet score")
        self.extra.setToolTip("DecontX contamination (% ambient RNA) and Scrublet's\n"
                              "doublet score, drawn when the dataset has them.")
        self.extra.toggled.connect(self.changed)
        self.add_row("", self.extra)

        # Y-axis range per metric, so the same range can be typed for every
        # study and figures placed side by side share their scales.
        from kosmic.visualisation.scrna.qc_plots import EXTRA_QC_METRICS, QC_METRICS
        short = {'n_genes_by_counts': "Genes:", 'total_counts': "Counts:",
                 'pct_counts_mt': "MT %:", 'decontX_contamination': "DecontX %:",
                 'doublet_score': "Doublet:"}
        self._form.addRow(QLabel("Y-axis range:"))
        self.y_range: dict[str, tuple[QLineEdit, QLineEdit]] = {}
        for metric in (*QC_METRICS, *EXTRA_QC_METRICS):
            row = QWidget()
            lay = QHBoxLayout(row)
            lay.setContentsMargins(0, 0, 0, 0)
            pair = []
            for i, tip in enumerate(("Lowest value on the axis", "Highest value on the axis")):
                # Empty is Auto; a number is a limit (0 included). Redraws on
                # Enter or leaving the box, not on every keystroke.
                box = QLineEdit()
                box.setPlaceholderText("Auto")
                box.setToolTip(f"{tip}. Leave empty for Auto: set from the data.")
                box.editingFinished.connect(self._range_edited)
                box.setMinimumWidth(0)
                lay.addWidget(box, 1)
                if i == 0:
                    lay.addWidget(QLabel("to"))
                pair.append(box)
            self.y_range[metric] = tuple(pair)
            self.add_row(short[metric], row)
        self.extra.toggled.connect(self._enable_extra_ranges)
        self._enable_extra_ranges(False)

    def _group_changed(self) -> None:
        """Switch to the new grouping's default height, unless the height
        was changed by hand. By sample stacks the metrics in rows, which
        need more height than panels side by side."""
        new = self.group.currentText()
        h = self.figsize.h_spin
        if h.value() == DEFAULT_HEIGHT[self._group_shown]:
            h.blockSignals(True)     # the group change redraws once already
            h.setValue(DEFAULT_HEIGHT[new])
            h.blockSignals(False)
        self._group_shown = new

    def _enable_extra_ranges(self, on: bool) -> None:
        from kosmic.visualisation.scrna.qc_plots import EXTRA_QC_METRICS
        for metric in EXTRA_QC_METRICS:
            for box in self.y_range[metric]:
                box.setEnabled(on)

    @staticmethod
    def _number(box: QLineEdit) -> Optional[float]:
        """The box's number (commas allowed), or None when empty or not a number."""
        text = box.text().replace(',', '').strip()
        try:
            return float(text) if text else None
        except ValueError:
            return None

    def _range_edited(self) -> None:
        # Text that is not a number is cleared, so the box shows Auto as it
        # is drawn; editingFinished fires on Enter and again on leaving, so
        # redraw only when the limits changed.
        box = self.sender()
        if box is not None and box.text().strip() and self._number(box) is None:
            box.clear()
        limits = self.y_limits()
        if limits != getattr(self, '_drawn_limits', {}):
            self._drawn_limits = limits
            self.changed.emit()

    def y_limits(self) -> dict:
        """{metric: (low, high)} for the metrics with a range set; None for
        an end left on Auto."""
        out = {}
        for metric, (lo, hi) in self.y_range.items():
            low, high = self._number(lo), self._number(hi)
            if low is not None or high is not None:
                out[metric] = (low, high)
        return out


class QCViolinsPage(FigurePage):
    TITLE = "QC Violins"
    UNAVAILABLE_MESSAGE = (
        "Compute QC metrics in the scRNA workflow to view this figure."
    )

    def _build_controls(self) -> FigureControls:
        return _QCViolinControls()

    def dependencies_met(self) -> bool:
        adata = self._adata()
        if adata is None or 'n_genes_by_counts' not in adata.obs.columns:
            return False
        group = self._controls.group.currentText()
        if group == GROUP_POOLED:
            return True
        cols = adata.obs.columns
        return '_role' in cols and (group == GROUP_CONDITION or 'sample' in cols)

    def _groups(self, obs: pd.DataFrame):
        """(group per cell, order, colour per group), or (None, None, None)
        for the pooled figure."""
        from kosmic.visualisation.style import CONTROL_COLOR, DISEASE_COLOR
        colour_of = {'control': CONTROL_COLOR, 'disease': DISEASE_COLOR}
        group = self._controls.group.currentText()
        if group == GROUP_POOLED:
            return None, None, None
        role = obs['_role'].astype(str)
        if group == GROUP_CONDITION:
            names = dict(ROLES)
            return (role.map(names).fillna(role), [n for _, n in ROLES],
                    {names[r]: colour_of[r] for r, _ in ROLES})
        sample = obs['sample'].astype(str)
        # Each sample under its majority role, controls first.
        sample_role = role.groupby(sample).agg(lambda r: r.value_counts().index[0])
        order = [s for r, _ in ROLES for s in sorted(sample_role.index[sample_role == r])]
        order += sorted(set(sample_role.index) - set(order))
        colours = {s: colour_of.get(sample_role[s]) for s in order
                   if sample_role[s] in colour_of}
        return sample, order, colours

    def _metrics(self) -> list:
        from kosmic.visualisation.scrna.qc_plots import EXTRA_QC_METRICS, QC_METRICS
        return [*QC_METRICS, *(EXTRA_QC_METRICS if self._controls.extra.isChecked() else ())]

    def _make_render_func(self):
        if not self.dependencies_met():
            return None
        obs = self._adata().obs
        groups, order, colours = self._groups(obs)
        limits = None
        if self._controls.limits.isChecked():
            from kosmic.provenance import qc_limits
            limits = qc_limits(_provenance_record(self))
        metrics = self._metrics()
        log_counts = self._controls.log_counts.isChecked()
        y_limits = self._controls.y_limits()
        figsize = self._controls.figsize.get_figsize()
        font_sizes = self._font_sizes()

        def _render(obs=obs, groups=groups, order=order, colours=colours, limits=limits,
                    metrics=metrics, log_counts=log_counts, y_limits=y_limits, sz=figsize,
                    f=font_sizes):
            from kosmic.visualisation.scrna.qc_plots import create_qc_violin_plots
            return create_qc_violin_plots(
                obs, metrics=metrics, thresholds=limits, groups=groups, group_order=order,
                group_colours=colours, log_counts=log_counts, y_limits=y_limits,
                figsize=sz, font_sizes=f)
        return _render

    def _export_companions(self, fig_path: Path) -> list[Path]:
        """Cells and the median of each metric per violin."""
        obs = self._adata().obs
        groups, order, _ = self._groups(obs)
        metrics = [m for m in self._metrics() if m in obs.columns]
        values = obs[metrics].apply(pd.to_numeric, errors='coerce')
        if 'decontX_contamination' in values:
            values['decontX_contamination'] *= 100
        key = pd.Series('All cells' if groups is None else groups.to_numpy(), index=obs.index)
        table = values.groupby(key).median().add_prefix('median_')
        table.insert(0, 'cells', key.value_counts())
        if order:
            table = table.reindex([g for g in order if g in table.index])
        table.index.name = 'group'
        out = Path(fig_path).with_suffix('.csv')
        table.round(4).to_csv(out)
        return [out]

    def _default_export_dir(self) -> Optional[Path]:
        return self._scrna_export_dir()

    def _default_export_basename(self) -> str:
        group = self._controls.group.currentText()
        return {GROUP_POOLED: "qc_violins", GROUP_CONDITION: "qc_violins_by_condition",
                GROUP_SAMPLE: "qc_violins_by_sample"}[group]


class _CellFilteringControls(FigureControls):
    def __init__(self):
        super().__init__(default_w=6, default_h=6)
        self.start = QComboBox()
        self.start.addItems([START_LOADED, START_FILTERED])
        self.start.setToolTip("The pre-analysis filter removes whole samples or\n"
                              "conditions; starting after it shows the QC steps\n"
                              "at a readable scale.")
        self.start.currentIndexChanged.connect(self.changed)
        self.add_row("Start at:", self.start)


class CellFilteringPage(FigurePage):
    TITLE = "Cell Filtering"
    UNAVAILABLE_MESSAGE = (
        "Run QC in the scRNA workflow to view this figure; it reads the cells "
        "removed at each step from the study's provenance record.")

    def _build_controls(self) -> FigureControls:
        return _CellFilteringControls()

    def _steps(self) -> list:
        from kosmic.provenance import cell_filtering_steps
        steps = cell_filtering_steps(_provenance_record(self))
        if self._controls.start.currentText() == START_FILTERED:
            last = max((i for i, s in enumerate(steps) if s['step'] == 'Pre-analysis filter'),
                       default=None)
            if last is not None:
                steps = [{**steps[last], 'step': 'After pre-analysis filter', 'n_removed': 0},
                         *steps[last + 1:]]
        return steps

    def dependencies_met(self) -> bool:
        return self._study_dir() is not None and len(self._steps()) > 1

    def _make_render_func(self):
        steps = self._steps()
        if len(steps) < 2:
            return None
        figsize = self._controls.figsize.get_figsize()
        font_sizes = self._font_sizes()

        def _render(steps=steps, sz=figsize, f=font_sizes):
            from kosmic.visualisation.scrna.qc_plots import create_cell_filtering_plot
            return create_cell_filtering_plot(steps, figsize=sz, font_sizes=f)
        return _render

    def _export_companions(self, fig_path: Path) -> list[Path]:
        """The steps as a table: step, cells remaining, cells removed."""
        out = Path(fig_path).with_suffix('.csv')
        pd.DataFrame(self._steps()).to_csv(out, index=False)
        return [out]

    def _default_export_dir(self) -> Optional[Path]:
        return self._scrna_export_dir()

    def _default_export_basename(self) -> str:
        return "cell_filtering"
