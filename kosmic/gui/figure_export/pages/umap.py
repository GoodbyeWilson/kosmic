# UMAP / t-SNE embedding plot page.
#
# Single page that switches between an Overview (2x2) layout and
# per-column plots via a 'Color by' combo. Reads from 'workspace.scrna_ws'.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import QComboBox, QDoubleSpinBox

from kosmic.gui.figure_export.pages._base import FigurePage
from kosmic.gui.figure_export.shared.controls_base import FigureControls

# 'Overview' = 2x2 grid of top categorical columns.
COLOR_OVERVIEW = "Overview (2x2 grid)"
# Highlight selector default: colour every category (no single-category focus).
HL_NONE = "(show all categories)"

_COL_PRIORITY = ['cell_type', 'leiden', 'sample', 'condition',
                 'clusters', 'seurat_clusters', 'batch']


def _pick_overview_columns(adata, max_cols: int = 4):
    chosen = []
    for col in _COL_PRIORITY:
        if col in adata.obs.columns and adata.obs[col].nunique() <= 50:
            chosen.append(col)
            if len(chosen) >= max_cols:
                return chosen
    for col in adata.obs.columns:
        if col in chosen:
            continue
        dtype = adata.obs[col].dtype
        if str(dtype) == 'category' or dtype == 'object':
            if 1 < adata.obs[col].nunique() <= 30:
                chosen.append(col)
                if len(chosen) >= max_cols:
                    break
    return chosen


class _UMAPControls(FigureControls):

    def __init__(self, default_w: int = 8, default_h: int = 8):
        super().__init__(default_w, default_h)

        self.color_by = QComboBox()
        self.color_by.currentIndexChanged.connect(self.changed)
        self.add_row("Color by:", self.color_by)

        self.embedding = QComboBox()
        self.embedding.addItems(["UMAP", "t-SNE"])
        self.embedding.currentIndexChanged.connect(self.changed)
        self.add_row("Embedding:", self.embedding)

        self.point_size = QDoubleSpinBox()
        self.point_size.setRange(0.1, 30.0)
        self.point_size.setDecimals(2)
        self.point_size.setSingleStep(0.25)
        self.point_size.setValue(2.0)
        self.point_size.valueChanged.connect(self.changed)
        self.add_row("Point size:", self.point_size)

        from kosmic.visualisation.scrna.umap import PALETTE_NAMES
        self.colours = QComboBox()
        self.colours.addItems(PALETTE_NAMES)
        self.colours.currentIndexChanged.connect(self.changed)
        self.add_row("Colours:", self.colours)

        # Colour only one category of the 'Color by' column; grey the rest.
        self.highlight = QComboBox()
        self.highlight.currentIndexChanged.connect(self.changed)
        self.add_row("Highlight:", self.highlight)

    def populate_highlight(self, values: list[str]):
        """Fill the highlight combo with '(show all)' plus category values."""
        prev = self.highlight.currentText()
        self.highlight.blockSignals(True)
        self.highlight.clear()
        self.highlight.addItems([HL_NONE] + list(values))
        idx = self.highlight.findText(prev)
        self.highlight.setCurrentIndex(idx if idx >= 0 else 0)
        self.highlight.blockSignals(False)

    def populate_columns(self, options: list[str]):
        """Populate the colour-by combo: single columns first (a biological
        label preferred as default), with the 2x2 overview last."""
        prev = self.color_by.currentText()
        self.color_by.blockSignals(True)
        self.color_by.clear()
        self.color_by.addItems(options + [COLOR_OVERVIEW])
        target = prev
        if not target:
            target = next(
                (c for c in ('cell_type', 'Names', 'predicted_labels', 'leiden')
                 if c in options),
                options[0] if options else COLOR_OVERVIEW)
        idx = self.color_by.findText(target)
        self.color_by.setCurrentIndex(idx if idx >= 0 else 0)
        self.color_by.blockSignals(False)


class UMAPPage(FigurePage):
    TITLE = "UMAP / t-SNE"
    UNAVAILABLE_MESSAGE = (
        "Run dimensionality reduction in the scRNA workflow to view this figure."
    )

    def _build_controls(self) -> _UMAPControls:
        ctrl = _UMAPControls()
        # Refresh the highlight values whenever the colour-by column changes.
        ctrl.color_by.currentIndexChanged.connect(self._on_color_by_changed)
        # Initial column population happens at on_activated time when
        # we know which adata is loaded.
        return ctrl

    def _sync_highlight(self, adata):
        """Populate the highlight combo with the current colour-by column's
        categories (or empty for the overview / no column)."""
        col = self._controls.color_by.currentText()
        if col and col != COLOR_OVERVIEW and col in adata.obs.columns:
            vals = sorted(adata.obs[col].astype(str).unique())
            self._controls.populate_highlight(vals)
        else:
            self._controls.populate_highlight([])

    def _on_color_by_changed(self):
        adata = self._adata()
        if adata is not None:
            self._sync_highlight(adata)

    def _embed_key(self) -> str:
        if self._controls.embedding.currentText() == 't-SNE':
            return 'X_tsne'
        return 'X_umap'

    def _adata(self):
        scrna = self.workspace.scrna_ws
        if scrna is None:
            return None
        return getattr(scrna, 'current_adata', None)

    def dependencies_met(self) -> bool:
        adata = self._adata()
        if adata is None:
            return False
        # Need at least one embedding key.
        return any(k in adata.obsm for k in ('X_umap', 'X_tsne'))

    def on_activated(self):
        adata = self._adata()
        if adata is not None:
            cols = self._available_columns(adata)
            self._controls.populate_columns(cols)
            self._sync_highlight(adata)
        super().on_activated()

    def _available_columns(self, adata) -> list[str]:
        """Categorical columns suitable for colouring, coarse-to-fine label
        columns first. The cap is high enough to include fine subtype columns
        (e.g. 'predicted_labels', ~70 categories)."""
        out: list[str] = []
        for c in ('cell_type', 'Names', 'predicted_labels', 'leiden',
                  'clusters', 'seurat_clusters', 'condition', 'sample', 'batch'):
            if c in adata.obs.columns and c not in out:
                out.append(c)
        for c in adata.obs.columns:
            if c in out:
                continue
            dtype = adata.obs[c].dtype
            if str(dtype) == 'category' or dtype == 'object':
                n = adata.obs[c].nunique()
                if 1 < n <= 100:
                    out.append(c)
        return out

    def _make_render_func(self):
        adata = self._adata()
        if adata is None:
            return None

        embed_key = self._embed_key()
        if embed_key not in adata.obsm:
            return None
        embed_name = 'UMAP' if 'umap' in embed_key else 't-SNE'
        ctrl = self._controls
        figsize = ctrl.figsize.get_figsize()
        point_size = ctrl.point_size.value()
        scheme = ctrl.colours.currentText()
        font_sizes = self._font_sizes()
        coords = adata.obsm[embed_key]

        choice = ctrl.color_by.currentText()

        if choice == COLOR_OVERVIEW or not choice:
            cols = _pick_overview_columns(adata)
            if not cols:
                return None
            obs = adata.obs

            def _render(coords=coords, obs=obs, cols=cols, embed_name=embed_name,
                        figsize=figsize, font_sizes=font_sizes, scheme=scheme):
                from kosmic.visualisation.scrna.umap import create_overview_figure
                kw = dict(embedding_name=embed_name, font_sizes=font_sizes, scheme=scheme)
                if figsize:
                    kw['figsize_per_panel'] = (figsize[0] // 2, figsize[1] // 2)
                return create_overview_figure(coords, obs, cols, **kw)
            return _render

        # Per-column plot
        col = choice
        if col not in adata.obs.columns:
            return None
        labels = adata.obs[col].astype(str).values

        # Highlight mode: colour only one category, grey the rest.
        hl = ctrl.highlight.currentText()
        if hl and hl != HL_NONE and hl in set(labels):
            title = f"{embed_name}: {hl}"

            def _render(coords=coords, labels=labels, hl=hl, title=title,
                        embed_name=embed_name, figsize=figsize,
                        point_size=point_size, font_sizes=font_sizes, scheme=scheme):
                from kosmic.visualisation.scrna.umap import create_highlight_plot
                return create_highlight_plot(
                    coords, labels, [hl], title=title,
                    embedding_name=embed_name, figsize=figsize,
                    point_size=point_size, font_sizes=font_sizes, scheme=scheme,
                )
            return _render

        title = f"{embed_name} coloured by {col}"

        def _render(coords=coords, labels=labels, title=title,
                    embed_name=embed_name, figsize=figsize,
                    point_size=point_size, font_sizes=font_sizes, scheme=scheme):
            from kosmic.visualisation.scrna.umap import create_embedding_plot
            return create_embedding_plot(
                coords, labels, title=title, embedding_name=embed_name,
                figsize=figsize, point_size=point_size, font_sizes=font_sizes,
                scheme=scheme,
            )
        return _render

    def _default_export_dir(self) -> Optional[Path]:
        scrna = self.workspace.scrna_ws
        if scrna is None or not getattr(scrna, 'current_project_dir', None):
            return None
        from kosmic.paths import scrna_figures_dir
        return scrna_figures_dir(scrna.current_project_dir)

    def _default_export_basename(self) -> str:
        choice = self._controls.color_by.currentText()
        if choice == COLOR_OVERVIEW:
            return "embedding_overview"
        hl = self._controls.highlight.currentText()
        if hl and hl != HL_NONE:
            safe = ''.join(ch if ch.isalnum() else '_' for ch in hl)
            return f"embedding_{choice}_{safe}"
        return f"embedding_{choice}"
