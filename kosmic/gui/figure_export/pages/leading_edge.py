# Leading-edge heatmap page.
#
# The strongest leading-edge genes of one or more gene sets from the last
# library GSEA run, and each gene's log2 fold change in every dataset. From
# the Meta-Analysis Enrichment step (pooled GSEA) the columns are the
# studies of the selection and the pooled estimate; from the DE workspace
# the single column is that study. The leading edge always comes from the
# run being drawn: the pooled one for a meta-analysis.

from __future__ import annotations

from pathlib import Path

import pandas as pd
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QDoubleSpinBox, QListWidget, QListWidgetItem, QSpinBox

from kosmic.gui.figure_export.pages._enrichment import EnrichmentFigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls
from kosmic.reference.pathways.libraries import term_label

#: Gene sets ticked when the page first shows a run.
_DEFAULT_TICKED = 2


class _LeadingEdgeControls(FigureControls):
    """Gene sets to draw, genes per set and the colour range."""

    def __init__(self):
        super().__init__()
        self.sets = QListWidget()
        self.sets.setMinimumHeight(140)
        self.sets.setToolTip("Gene sets from the GSEA run, most significant "
                             "first. Each ticked set is one block of rows.")
        self.sets.itemChanged.connect(self.changed)
        self.add_row("Gene sets:", self.sets)

        self.n_genes = QSpinBox()
        self.n_genes.setRange(3, 50)
        self.n_genes.setValue(15)
        self.n_genes.valueChanged.connect(self.changed)
        self.add_row("Genes per set:", self.n_genes)

        self.vmax = QDoubleSpinBox()
        self.vmax.setRange(0.5, 10.0)
        self.vmax.setSingleStep(0.5)
        self.vmax.setDecimals(1)
        self.vmax.setValue(2.0)
        self.vmax.valueChanged.connect(self.changed)
        self.add_row("Colour range (+/- log2FC):", self.vmax)

    def ticked(self) -> list[str]:
        return [self.sets.item(i).data(Qt.ItemDataRole.UserRole)
                for i in range(self.sets.count())
                if self.sets.item(i).checkState() == Qt.CheckState.Checked]


class LeadingEdgeHeatmapPage(EnrichmentFigurePage):
    TITLE = "Leading-Edge Heatmap"
    CHECKS_DE_STALENESS = True
    KIND = 'gsea'
    UNAVAILABLE_MESSAGE = (
        "Run preranked GSEA against a library on the DE or Meta-Analysis "
        "Enrichment step first -- this figure shows the leading-edge genes "
        "of that run."
    )

    def _build_controls(self) -> FigureControls:
        ctrl = self._add_source_row(_LeadingEdgeControls())
        ctrl.source.currentIndexChanged.connect(self._sync_sets)
        return ctrl

    # -- data -------------------------------------------------------------

    def _meta_source(self):
        src = super()._meta_source()
        if src is not None:
            inputs = self.workspace.meta_ws.figure_inputs()
            src.gene_results = inputs.get('gene_results')
            src.gene_studies = inputs.get('gene_studies') or []
        return src

    def _summary(self):
        src = self._source()
        res = getattr(src, 'gsea_library_results', None)
        if res is None or not len(res) or 'leading_edge' not in res.columns:
            return None
        return res

    def _columns(self, src) -> list[tuple]:
        """(label, DE table, adjusted-p column) for each heatmap column."""
        if self._is_meta():
            cols = [(name, df, 'pvals_adj') for name, df in src.gene_studies]
            pooled = getattr(src, 'gene_results', None)
            if pooled is None:
                pooled = getattr(src, 'de_results', None)
            if pooled is not None:
                cols.append(('Pooled', pooled, 'fdr'))
            return cols
        de = getattr(src, 'de_results', None)
        if de is None:
            return []
        label = Path(src.project_dir).name if getattr(src, 'project_dir', None) else 'This study'
        return [(label, de, 'pvals_adj')]

    def _ranking(self, src) -> pd.Series:
        """The statistic the GSEA ranked genes by, indexed by gene."""
        if self._is_meta():
            table = getattr(src, 'de_results', None)
            if table is None:
                table = getattr(src, 'gene_results', None)
            col = 'z_stat'
        else:
            table, col = getattr(src, 'de_results', None), 'logfoldchanges'
        if table is None or col not in table.columns:
            col = 'logfoldchanges'
        if table is None:
            return pd.Series(dtype=float)
        return pd.to_numeric(table.drop_duplicates('names').set_index('names')[col],
                             errors='coerce')

    # -- lifecycle --------------------------------------------------------

    def on_activated(self):
        self._sync_sources()
        self._sync_sets()
        self.invalidate()        # the run lives in memory and may have changed
        super().on_activated()

    def _sync_sets(self):
        """List the run's gene sets, most significant first, keeping ticks."""
        res = self._summary()
        lst = self._controls.sets
        kept = set(self._controls.ticked())
        lst.blockSignals(True)
        lst.clear()
        if res is not None:
            order = res.assign(_p=pd.to_numeric(res['pvals_adj'], errors='coerce'),
                               _n=-pd.to_numeric(res['nes'], errors='coerce').abs()
                               ).sort_values(['_p', '_n'])
            for k, (name, nes, fdr) in enumerate(zip(order['names'], order['nes'], order['_p'])):
                item = QListWidgetItem(f"{term_label(name)}  (NES {nes:+.2f}, FDR {fdr:.2g})")
                item.setData(Qt.ItemDataRole.UserRole, str(name))
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                tick = (str(name) in kept) if kept else k < _DEFAULT_TICKED
                item.setCheckState(Qt.CheckState.Checked if tick else Qt.CheckState.Unchecked)
                lst.addItem(item)
        lst.blockSignals(False)

    # -- render -----------------------------------------------------------

    def _blocks(self):
        """The heatmap blocks for the ticked sets, or None."""
        from kosmic.visualisation.de.leading_edge import (
            leading_edge_genes, leading_edge_table,
        )
        src, res = self._source(), self._summary()
        if src is None or res is None:
            return None
        columns, ranking = self._columns(src), self._ranking(src)
        if not columns:
            return None
        n = self._controls.n_genes.value()
        by_name = res.drop_duplicates('names').set_index(res.drop_duplicates('names')['names'].astype(str))
        blocks = []
        for name in self._controls.ticked():
            if name not in by_name.index:
                continue
            row = by_name.loc[name]
            nes = float(row['nes'])
            genes = leading_edge_genes(row['leading_edge'], ranking, nes, n=n)
            lfc, padj = leading_edge_table(genes, columns)
            blocks.append({
                'name': name,
                'title': f"{term_label(name)}  (NES {nes:+.2f}, "
                         f"FDR {float(row['pvals_adj']):.2g})",
                'lfc': lfc, 'padj': padj})
        return blocks or None

    def _make_render_func(self):
        blocks = self._blocks()
        if not blocks:
            return None
        src = self._source()
        lib = getattr(src, 'gsea_library_name', None)
        contrast = self._contrast(src)
        title = f"{lib}: {contrast}" if lib and contrast else lib
        gap = 'Pooled' if self._is_meta() else None
        vmax = self._controls.vmax.value()
        font_sizes = self._font_sizes()
        figsize = self._controls.figsize.get_figsize()

        def _render(blocks=blocks, vmax=vmax, gap=gap, title=title,
                    f=font_sizes, sz=figsize):
            from kosmic.visualisation.de.leading_edge import create_leading_edge_heatmap
            return create_leading_edge_heatmap(
                blocks, vmax=vmax, gap_before=gap, title=title,
                figsize=sz, font_sizes=f)
        return _render

    def _export_companions(self, fig_path: Path) -> list[Path]:
        """The numbers behind the cells, one row per gene set, gene and column."""
        blocks = self._blocks()
        if not blocks:
            return []
        rows = []
        for b in blocks:
            long = (b['lfc'].stack(future_stack=True).rename('log2fc').to_frame()
                    .join(b['padj'].stack(future_stack=True).rename('padj')))
            long.index.names = ['gene', 'column']
            long.insert(0, 'gene_set', b['name'])       # the library's own name
            rows.append(long.reset_index())
        out = Path(fig_path).with_suffix('.csv')
        pd.concat(rows, ignore_index=True).round(6).to_csv(out, index=False)
        return [out]

    def _default_export_basename(self) -> str:
        prefix = 'meta_' if self._is_meta() else ''
        return f"{prefix}leading_edge_heatmap"
