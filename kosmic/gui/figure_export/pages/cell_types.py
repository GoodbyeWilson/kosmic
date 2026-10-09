# Cell-type figure pages: a marker dot plot, one embedding per marker panel
# coloured by its module score, and cell-type composition by condition or
# sample.
#
# Read 'workspace.scrna_ws.current_adata'. Disease and control come from
# 'adata.obs["_role"]', so studies that name their conditions differently
# are grouped the same way.

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
from PyQt6.QtWidgets import QComboBox, QDoubleSpinBox, QFileDialog, QSpinBox

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.pages.umap import label_columns
from kosmic.gui.figure_export.shared.controls_base import FigureControls
from kosmic.gui.shared import dialogs

#: Labels that mean "no cell type"; their cells are left out of the figures.
UNLABELLED = {'', 'nan', 'none', 'unknown', 'unassigned'}

SOURCE_GENERAL = "General (built-in)"
SOURCE_FILE = "From file..."
BARS_CONDITION, BARS_SAMPLE = "By condition", "By sample"
ROLES = (('control', 'Control'), ('disease', 'Disease'))


def _marker_sources() -> dict:
    """{source name: {cell type: [genes]}}: the bundled curated panels, then
    the general built-in markers."""
    from kosmic.reference.marker_panels import builtin_panels
    from kosmic.scrna.annotate.score import DEFAULT_MARKERS
    return {**builtin_panels(), SOURCE_GENERAL: dict(DEFAULT_MARKERS)}


def _labelled(values) -> np.ndarray:
    """Boolean mask of cells that carry a cell-type label."""
    v = np.char.lower(np.char.strip(np.asarray(values, dtype=str)))
    return ~np.isin(v, list(UNLABELLED))


def _expression_label(adata) -> str:
    """Colour bar label naming the unit of 'adata.X' when it is recognisably
    log1p of counts per 10,000 (each cell's expm1 sums to 10,000)."""
    import scipy.sparse as sp
    x = adata.X[:50]
    x = x.toarray() if sp.issparse(x) else np.asarray(x)
    sums = np.expm1(x.astype(float)).sum(axis=1)
    if len(sums) and np.allclose(sums, 1e4, rtol=0.01):
        return 'Mean expression\nlog1p(CP10K)'
    return 'Mean expression'


class _ColumnControls(FigureControls):
    """Controls with a cell-type column selector."""

    def __init__(self, default_w: int = 0, default_h: int = 0):
        super().__init__(default_w, default_h)
        self.column = QComboBox()
        self.column.currentIndexChanged.connect(self.changed)
        self.add_row("Cell types:", self.column)

    def populate_columns(self, options: list[str]):
        prev = self.column.currentText()
        self.column.blockSignals(True)
        self.column.clear()
        self.column.addItems(options)
        target = prev if prev in options else next(
            (c for c in ('cell_type', 'Names', 'predicted_labels') if c in options),
            options[0] if options else '')
        self.column.setCurrentIndex(max(self.column.findText(target), 0))
        self.column.blockSignals(False)


class _MarkerSource(QComboBox):
    """Marker source: a bundled panel, the general built-in markers, or a
    gene-set file chosen when 'From file...' is picked."""

    def __init__(self, controls: FigureControls):
        super().__init__()
        self._controls = controls
        self.file_markers: dict = {}
        self.addItems(list(_marker_sources()) + [SOURCE_FILE])
        self.activated.connect(self._on_activated)
        self.currentIndexChanged.connect(controls.changed)
        controls.add_row("Markers:", self)

    def _on_activated(self):
        if self.currentText() != SOURCE_FILE:
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Open marker set", "",
            "GMT gene sets (*.gmt);;All files (*)")
        markers = {}
        if path:
            from kosmic.reference.pathways.formats import load_gmt
            try:
                markers = {k: v for k, v in load_gmt(path).items() if v}
            except Exception as exc:
                dialogs.warning(self, "Could not read file", str(exc))
        if markers:
            self.file_markers = markers
            self._controls.changed.emit()
        else:
            self.setCurrentIndex(0)

    def markers(self) -> dict:
        if self.currentText() == SOURCE_FILE:
            return self.file_markers
        return _marker_sources().get(self.currentText(), {})


class _ScRNAPage(FigurePage):
    """Base: reads the scRNA workspace's current dataset."""

    def _adata(self):
        scrna = self.workspace.scrna_ws
        return getattr(scrna, 'current_adata', None) if scrna is not None else None

    def _default_export_dir(self) -> Optional[Path]:
        scrna = self.workspace.scrna_ws
        if scrna is None or not getattr(scrna, 'current_project_dir', None):
            return None
        from kosmic.paths import scrna_figures_dir
        return scrna_figures_dir(scrna.current_project_dir)


class _CellTypePage(_ScRNAPage):
    """Base for pages grouped by a cell-type column."""

    def _column(self) -> str:
        return self._controls.column.currentText()

    def dependencies_met(self) -> bool:
        adata = self._adata()
        return adata is not None and self._column() in adata.obs.columns

    def on_activated(self):
        adata = self._adata()
        if adata is not None:
            self._controls.populate_columns(label_columns(adata))
        super().on_activated()


# -- Marker dot plot -------------------------------------------------------------

class MarkerDotPlotPage(_CellTypePage):
    TITLE = "Marker Dot Plot"
    UNAVAILABLE_MESSAGE = "Load an annotated dataset in the scRNA workflow to view this figure."

    def __init__(self, workspace):
        self._cache: dict = {}
        super().__init__(workspace)

    def _build_controls(self):
        ctrl = _ColumnControls()
        ctrl.source = _MarkerSource(ctrl)
        return ctrl

    def _summary(self, adata, col, source, markers):
        """(mean, pct, owner) for the panel genes present, cached per dataset,
        column and marker source."""
        key = (id(adata), self.workspace.data_version(), col, source,
               tuple(markers))
        if key not in self._cache:
            from kosmic.scrna.inspect.gene_group import marker_panel, summarise_gene_group
            genes, owner = marker_panel(markers, cell_types=list(markers))
            keep = _labelled(adata.obs[col])
            s = summarise_gene_group(adata[keep], genes, group_col=col)
            owner = {g: owner[g] for g in s.present}
            self._cache = {key: (s.mean_matrix, s.pct_matrix, owner)}
        return self._cache[key]

    def _make_render_func(self):
        adata = self._adata()
        col = self._column()
        markers = self._controls.source.markers()
        if adata is None or col not in adata.obs.columns or not markers:
            return None
        source = self._controls.source.currentText()
        figsize = self._controls.figsize.get_figsize()
        font_sizes = self._font_sizes()

        def _render(adata=adata, col=col, source=source, markers=markers,
                    sz=figsize, f=font_sizes):
            from kosmic.visualisation.scrna.markers import create_marker_dotplot
            mean, pct, owner = self._summary(adata, col, source, markers)
            if mean.empty:
                return None
            return create_marker_dotplot(mean, pct, owner, figsize=sz, font_sizes=f,
                                         colour_label=_expression_label(adata))
        return _render

    def _export_companions(self, fig_path: Path) -> list[Path]:
        """The numbers behind the dots, one row per gene and cell type."""
        adata = self._adata()
        mean, pct, owner = self._summary(adata, self._column(),
                                         self._controls.source.currentText(),
                                         self._controls.source.markers())
        long = (mean.stack().rename('mean_expression').to_frame()
                .join(pct.mul(100).stack().rename('pct_expressing')))
        long.index.names = ['gene', 'cell_type']
        long.insert(0, 'marks', [owner[g] for g, _ in long.index])
        out = Path(fig_path).with_suffix('.csv')
        long.round(4).to_csv(out)
        return [out]

    def _default_export_basename(self) -> str:
        return "marker_dotplot"


# -- Marker score embeddings -----------------------------------------------------

class _ScoreControls(FigureControls):
    def __init__(self):
        super().__init__()
        self.source = _MarkerSource(self)

        self.embedding = QComboBox()
        self.embedding.addItems(["UMAP", "t-SNE"])
        self.embedding.currentIndexChanged.connect(self.changed)
        self.add_row("Embedding:", self.embedding)

        self.n_cells = QSpinBox()
        self.n_cells.setRange(10_000, 2_000_000)
        self.n_cells.setSingleStep(10_000)
        self.n_cells.setValue(150_000)
        self.n_cells.setGroupSeparatorShown(True)
        self.n_cells.setToolTip("Cells drawn, sampled at random. Scoring and\n"
                                "drawing every cell of a large atlas is slow and\n"
                                "looks the same.")
        self.n_cells.valueChanged.connect(self.changed)
        self.add_row("Cells shown:", self.n_cells)

        self.point_size = QDoubleSpinBox()
        self.point_size.setRange(0.1, 10.0)
        self.point_size.setDecimals(2)
        self.point_size.setSingleStep(0.25)
        self.point_size.setValue(1.0)
        self.point_size.valueChanged.connect(self.changed)
        self.add_row("Point size:", self.point_size)

        self.ncols = QSpinBox()
        self.ncols.setRange(1, 8)
        self.ncols.setValue(4)
        self.ncols.valueChanged.connect(self.changed)
        self.add_row("Columns:", self.ncols)


class MarkerScoreUMAPPage(_ScRNAPage):
    TITLE = "Marker Score UMAPs"
    UNAVAILABLE_MESSAGE = (
        "Run dimensionality reduction in the scRNA workflow to view this figure.")

    def __init__(self, workspace):
        self._cache: dict = {}
        super().__init__(workspace)

    def _build_controls(self):
        return _ScoreControls()

    def _embed_key(self) -> str:
        return 'X_tsne' if self._controls.embedding.currentText() == 't-SNE' else 'X_umap'

    def dependencies_met(self) -> bool:
        adata = self._adata()
        return adata is not None and self._embed_key() in adata.obsm

    def _make_render_func(self):
        adata = self._adata()
        markers = self._controls.source.markers()
        key = self._embed_key()
        if adata is None or key not in adata.obsm or not markers:
            return None
        c = self._controls
        source, n_cells = c.source.currentText(), c.n_cells.value()
        point_size, ncols = c.point_size.value(), c.ncols.value()
        figsize = c.figsize.get_figsize()
        font_sizes = self._font_sizes()

        def _render(adata=adata, markers=markers, key=key, source=source,
                    n_cells=n_cells, point_size=point_size, ncols=ncols,
                    sz=figsize, f=font_sizes):
            from kosmic.visualisation.scrna.markers import (
                create_marker_score_grid, panel_scores,
            )
            ck = (id(adata), self.workspace.data_version(), source, tuple(markers), n_cells)
            if ck not in self._cache:
                n = adata.n_obs
                idx = (np.arange(n) if n <= n_cells else
                       np.sort(np.random.default_rng(0).choice(n, n_cells, replace=False)))
                self._cache = {ck: (idx, panel_scores(adata[idx], markers))}
            idx, scores = self._cache[ck]
            if not scores:
                return None
            return create_marker_score_grid(
                np.asarray(adata.obsm[key])[idx], scores, ncols=ncols,
                point_size=point_size, title='Marker module scores',
                figsize=sz, font_sizes=f)
        return _render

    def _default_export_basename(self) -> str:
        return "marker_score_umaps"


# -- Composition -----------------------------------------------------------------

class _CompositionControls(_ColumnControls):
    def __init__(self):
        super().__init__()
        self.bars = QComboBox()
        self.bars.addItems([BARS_CONDITION, BARS_SAMPLE])
        self.bars.currentIndexChanged.connect(self.changed)
        self.add_row("Bars:", self.bars)

        from kosmic.visualisation.scrna.umap import PALETTE_NAMES
        self.colours = QComboBox()
        self.colours.addItems(PALETTE_NAMES)
        self.colours.setToolTip("The same scheme as the UMAP page gives each\n"
                                "cell type the same colour in both figures.")
        self.colours.currentIndexChanged.connect(self.changed)
        self.add_row("Colours:", self.colours)


class CompositionPage(_CellTypePage):
    TITLE = "Cell-Type Composition"
    UNAVAILABLE_MESSAGE = (
        "Load an annotated dataset with disease and control roles set (Inspect "
        "tab) to view this figure.")

    def _build_controls(self):
        return _CompositionControls()

    def dependencies_met(self) -> bool:
        adata = self._adata()
        return (super().dependencies_met() and '_role' in adata.obs.columns
                and (self._controls.bars.currentText() == BARS_CONDITION
                     or 'sample' in adata.obs.columns))

    def _table(self):
        """(percentages, bracket blocks or None, colours) for the current
        controls."""
        from kosmic.visualisation.scrna.composition import composition_table
        from kosmic.visualisation.scrna.umap import category_colours
        obs = self._adata().obs
        col = self._column()
        colours = category_colours(obs[col].astype(str), self._controls.colours.currentText())
        role = obs['_role'].astype(str)
        keep = _labelled(obs[col]) & role.isin([r for r, _ in ROLES]).to_numpy()
        obs, role = obs[keep], role[keep]
        if self._controls.bars.currentText() == BARS_CONDITION:
            names = dict(ROLES)
            return (composition_table(obs[col].astype(str), role.map(names),
                                      [n for _, n in ROLES]), None, colours)
        sample = obs['sample'].astype(str)
        # Each sample under its majority role, controls first.
        sample_role = role.groupby(sample).agg(lambda r: r.value_counts().index[0])
        order, blocks = [], []
        for r, name in ROLES:
            members = sorted(sample_role.index[sample_role == r])
            if members:
                order += members
                blocks.append((name, len(members)))
        return composition_table(obs[col].astype(str), sample, order), blocks, colours

    def _make_render_func(self):
        if not self.dependencies_met():
            return None
        pct, blocks, colours = self._table()
        if pct.empty:
            return None
        figsize = self._controls.figsize.get_figsize()
        font_sizes = self._font_sizes()

        def _render(pct=pct, blocks=blocks, colours=colours, sz=figsize, f=font_sizes):
            from kosmic.visualisation.scrna.composition import create_composition_plot
            return create_composition_plot(pct, colours, blocks=blocks, figsize=sz,
                                           font_sizes=f)
        return _render

    def _export_companions(self, fig_path: Path) -> list[Path]:
        """The percentages, one column per bar, cell types in stack order."""
        pct, _, _ = self._table()
        out = Path(fig_path).with_suffix('.csv')
        t = pct.T.round(4)
        t.index.name = 'cell_type'
        t.to_csv(out)
        return [out]

    def _default_export_basename(self) -> str:
        by = 'sample' if self._controls.bars.currentText() == BARS_SAMPLE else 'condition'
        return f"composition_by_{by}"
