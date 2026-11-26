# Methods Comparison figure pages (ADR-017): the Venn diagram of the genes
# two or three methods of a comparison call.
#
# Read the comparison through the workspace's comparison_inputs(): the
# Methods Comparison workspace's comparison (ComparisonWorkspace.
# figure_inputs()), or a comparison saved in the project's
# methods_comparison/ folder.

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import QComboBox

from kosmic.comparison.standard import METHOD_LABELS
from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls

_NONE = "None"
#: The methods the diagram opens on, when the comparison has them.
_DEFAULT_METHODS = ('pb_merged', 'ma_reml', 'ma_vote')


class _VennControls(FigureControls):
    def __init__(self):
        super().__init__(5, 5)
        self.method_boxes = []
        for i, label in enumerate(("First method:", "Second method:", "Third method:")):
            box = QComboBox()
            box.currentIndexChanged.connect(self.changed)
            self.add_row(label, box)
            self.method_boxes.append(box)
        self.method_boxes[2].setToolTip("None draws a diagram of two methods.")

    def populate(self, methods: list):
        """Offer ``methods``; keep the choices made where still offered."""
        for i, box in enumerate(self.method_boxes):
            # A box filled before keeps its choice ('None' included).
            chose = box.count() > 0
            prev = box.currentData()
            box.blockSignals(True)
            box.clear()
            if i == 2:
                box.addItem(_NONE, None)
            for m in methods:
                box.addItem(METHOD_LABELS.get(m, m), m)
            index = box.findData(prev) if chose else -1
            if index < 0:
                wanted = [m for m in _DEFAULT_METHODS if m in methods] + methods
                index = box.findData(wanted[i] if i < len(wanted) else None)
            box.setCurrentIndex(max(index, 0))
            box.blockSignals(False)

    def chosen(self) -> list:
        """The methods chosen, each once, in the order of the boxes."""
        out = []
        for box in self.method_boxes:
            m = box.currentData()
            if m is not None and m not in out:
                out.append(m)
        return out


class MethodsVennPage(FigurePage):
    """Venn diagram of the genes two or three methods call."""

    TITLE = "Genes Called (Venn)"
    UNAVAILABLE_MESSAGE = ("Run a comparison in the Methods Comparison workspace to view "
                           "this figure.")

    def _build_controls(self):
        return _VennControls()

    def _bundle(self):
        inputs = getattr(self.workspace, 'comparison_inputs', lambda: None)()
        return (inputs or {}).get('bundle')

    def dependencies_met(self) -> bool:
        bundle = self._bundle()
        return bundle is not None and len(bundle.methods) >= 2

    def on_activated(self):
        bundle = self._bundle()
        self._controls.populate(bundle.methods if bundle is not None else [])
        super().on_activated()

    def _sets(self) -> Optional[dict]:
        bundle = self._bundle()
        chosen = [m for m in self._controls.chosen() if bundle is not None
                  and m in bundle.methods]
        if bundle is None or len(chosen) < 2:
            return None
        sets = bundle.sets
        return {m: sets[m] for m in chosen}

    def _title(self, bundle) -> str:
        cell_type = bundle.params.get('cell_type') or 'Genes called'
        return f"{cell_type}: genes called ({bundle.rule.text})"

    def _make_render_func(self):
        bundle, sets = self._bundle(), self._sets()
        if bundle is None or sets is None:
            return None
        labels = {m: METHOD_LABELS.get(m, m) for m in sets}
        n, title = bundle.n_universe, self._title(bundle)
        figsize = self._controls.figsize.get_figsize()
        font_sizes = self._font_sizes()

        def _render(sets=sets, labels=labels, n=n, title=title, sz=figsize, f=font_sizes):
            from kosmic.visualisation.comparison.venn import create_venn_plot
            return create_venn_plot(sets, labels, n, title=title, figsize=sz, font_sizes=f)
        return _render

    def _export_companions(self, fig_path: Path) -> list[Path]:
        """Each region's genes, and the comparison's parameters, beside the figure."""
        from kosmic.comparison.agreement import venn_regions
        bundle, sets = self._bundle(), self._sets()
        if bundle is None or sets is None:
            return []
        regions = venn_regions(sets).rename(columns={m: METHOD_LABELS.get(m, m)
                                                     for m in sets})
        regions['region'] = regions['region'].map(
            lambda r: ' + '.join(METHOD_LABELS.get(m, m) for m in r.split(' + ')))
        csv = Path(fig_path).with_suffix('.csv')
        regions.to_csv(csv, index=False)
        params = Path(fig_path).with_name(Path(fig_path).stem + '_params.json')
        params.write_text(json.dumps({**bundle.params, 'methods_shown': list(sets),
                                      'rule': bundle.rule.as_dict(),
                                      'rule_text': bundle.rule.text}, indent=2, default=str),
                          encoding='utf-8')
        return [csv, params]

    def _default_export_dir(self) -> Optional[Path]:
        bundle = self._bundle()
        return bundle.folder if bundle is not None else None

    def _default_export_basename(self) -> str:
        return "methods_venn"
