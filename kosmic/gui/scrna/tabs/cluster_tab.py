# Cluster Tab
# Preview-first clustering layout with integrated PCA + batch correction.
# Large interactive UMAP fills most of the screen; slim parameter panel
# on the left (PCA + clustering controls) and a resolution-sweep
# thumbnail strip below. "Elbow Plot" tab sits alongside the UMAP tab.

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QLineEdit, QCompleter,
    QCheckBox, QFrame,
    QScrollArea, QTabWidget,
)
from PyQt6.QtCore import Qt, pyqtSignal, QTimer
from PyQt6.QtGui import QPixmap, QFont, QColor
from pathlib import Path
from typing import Optional
import numpy as np
import pyqtgraph as pg

from kosmic import UI_EMBEDDING_MAX_POINTS

from kosmic.gui.shared.theme import get_color, style_pg_plot, NoScrollComboBox, NoScrollSpinBox, NoScrollDoubleSpinBox
from kosmic.gui.shared.plots import InteractivePlot
from kosmic.gui.shared.widgets import (
    BaseWorker, SidebarPage, SecondaryLabel, SectionHeader, PrimaryButton,
    SecondaryButton, StageAccordion, StageSummaryCard, show_methods_report,
)
from kosmic.paths import processed_data_dir, processed_h5ad_path
from kosmic.gui.shared import borderless, dialogs, run_worker


def _refresh_style(widget) -> None:
    """
    Force Qt to re-evaluate stylesheet selectors after a dynamic
    property change (Qt doesn't re-polish on setProperty by default).
    """
    style = widget.style()
    if style is not None:
        style.unpolish(widget)
        style.polish(widget)


# ---------------------------------------------------------------------------
# Workers: thin wrappers around the analysis package
# ---------------------------------------------------------------------------

class PCAHarmonyWorker(BaseWorker):
    """
    Worker: normalise -> HVG selection -> PCA -> optional Harmony.

    Runs the full pre-clustering pipeline in a background thread so the
    GUI stays responsive. Emits 'finished_ok' with '(adata, message)'.
    """

    def __init__(self, adata, n_hvg: int, n_pcs: int, *,
                 harmony: bool = False,
                 scale: bool = False,
                 honour_roles: bool = True,
                 harmony_key: str = "",
                 harmony_condition_key: str = "",
                 output_path: Optional[str] = None):
        super().__init__()
        self.adata = adata
        self.n_hvg = n_hvg
        self.n_pcs = n_pcs
        self.harmony = harmony
        self.scale = scale
        self.honour_roles = honour_roles
        self.harmony_key = harmony_key
        self.harmony_condition_key = harmony_condition_key
        self.output_path = output_path

    def _run(self):
        import scanpy as sc
        from scipy import sparse as sp
        from kosmic.scrna.cluster.hvg import find_hvg
        from kosmic.scrna.cluster.harmony import run_harmony
        from kosmic.scrna.cluster.workset import apply_results, working_subset
        from kosmic.scrna.inspect.roles import included_mask

        # Cells marked '_role == exclude' are left out of the analysis but
        # not out of the file. Without this an excluded arm still drives
        # HVG selection, PCA and the graph -- on GSE292067 that is 44% of
        # the cells, from a cohort the DE contrast never touches.
        self._full = self.adata
        self._mask = included_mask(self.adata) if self.honour_roles else None

        # --- Is this counts or already log-normalised? ------------------
        X = self._full.X
        if sp.issparse(X):
            max_val = float(np.max(X.data)) if len(X.data) > 0 else 0
        else:
            max_val = float(np.max(X))
        already_normalized = max_val < 20

        # The step works on one matrix for the included cells and its own
        # obs/var, never on a copy of the whole study (see workset.py). The
        # matrix is shared when every cell is included and it is only read;
        # counts that still need normalising get a private copy.
        if self._mask is not None:
            n_out = int((~self._mask).sum())
            self.progress.emit(
                f"Excluding {n_out:,} cell(s) marked 'exclude' from the "
                f"embedding ({int(self._mask.sum()):,} remain)...")
        else:
            self.progress.emit("Preparing dataset...")
        self.adata = working_subset(
            self._full, self._mask, copy_matrix=not already_normalized)
        self.progress_pct.emit(2)

        # --- HVG ---------------------------------------------------------
        # Before normalisation, not after: 'seurat_v3' fits its
        # mean-variance curve to raw counts, and it is the selector that
        # copes with the depth differences between studies. Running it here
        # avoids stashing a second copy of the counts on a dataset that can
        # be tens of GB.
        #
        # Selection is per batch when there is one. Pooling every cell
        # first lets genes that differ between studies score as "variable",
        # so the gene set carries the batch effect into PCA -- before
        # Harmony ever gets a chance to remove it.
        hvg_batch = self.harmony_key if self.harmony else None
        self.progress.emit(
            f"Selecting up to {self.n_hvg:,} highly variable genes"
            + (f" within '{hvg_batch}'..." if hvg_batch else "..."))
        self.progress_pct.emit(10)
        find_hvg(self.adata, n_hvg=self.n_hvg, subset=False,
                 batch_key=hvg_batch)
        self._hvg_flavor = str(
            self.adata.uns.get('hvg', {}).get('flavor', 'unknown'))
        self.progress.emit(f"HVG method: {self._hvg_flavor}")

        # --- Normalise (skip if already log-transformed) ----------------
        if not already_normalized:
            # Only the step's private matrix is normalised; the study keeps
            # its counts in X. The QC tab's Normalise is the route that
            # normalises the study itself.
            self.progress.emit(
                "X holds counts: normalising a private copy for the PCA "
                "(the study is unchanged; run Normalise on the QC tab)...")
            self.progress_pct.emit(20)
            sc.pp.normalize_total(self.adata, target_sum=10000)
            sc.pp.log1p(self.adata)
        else:
            self.progress.emit("Data already normalised -- skipping")
            self.progress_pct.emit(20)

        # --- Scale (optional) --------------------------------------------
        # Z-scores each gene so PCA weighs them equally. Off by default:
        # it turns the sparse matrix dense (every zero becomes a negative
        # number), which on a large master costs far more memory than the
        # rest of the pipeline combined.
        if self.scale:
            self.progress.emit("Scaling genes (z-score, clipped at 10)...")
            self.progress_pct.emit(30)
            sc.pp.scale(self.adata, max_value=10, zero_center=True)

        # --- PCA ---------------------------------------------------------
        n_pcs = min(self.n_pcs, self.adata.n_obs - 1, self.adata.n_vars - 1)
        # Default centred PCA densifies the matrix (n_cells x n_genes
        # floats). At ~500K+ cells that runs out of RAM; switch to the
        # sparse-native zero_center=False path used by atlas pipelines.
        large = self.adata.n_obs > 500_000
        self.progress.emit(
            f"Running PCA ({n_pcs} components"
            f"{', sparse mode' if large else ''})...")
        self.progress_pct.emit(40)
        sc.tl.pca(self.adata, n_comps=n_pcs, zero_center=not large)
        self.progress_pct.emit(70)

        # --- Harmony (optional) ------------------------------------------
        if self.harmony and self.harmony_key:
            self.progress.emit(f"Running Harmony on '{self.harmony_key}'...")
            self.progress_pct.emit(75)
            try:
                run_harmony(
                    self.adata,
                    batch_key=self.harmony_key,
                    condition_key=self.harmony_condition_key or None,
                )
            except ImportError as e:
                raise RuntimeError(
                    "harmonypy is not installed. Install it with: pip install harmonypy") from e
            self.progress_pct.emit(95)
            self.progress.emit("Harmony correction complete.")
        else:
            # Drop stale harmony embedding if user toggled it off
            if 'X_pca_harmony' in self.adata.obsm:
                del self.adata.obsm['X_pca_harmony']

        n_hvg_found = int(self.adata.var['highly_variable'].sum())
        excluded_msg = ""

        # Put the analysis back on the full object. Excluded cells survive
        # in the file with their metadata and NaN where the analysis never
        # saw them. A new embedding invalidates any earlier graph.
        self.progress.emit("Writing results back onto the dataset...")
        self.adata.var['highly_variable'] = self.adata.var[
            'highly_variable'].astype(bool)
        self._full.var = self.adata.var.reindex(self._full.var_names)
        apply_results(
            self._full, self.adata, self._mask,
            obsm_keys=('X_pca', 'X_pca_harmony'),
            uns_keys=('hvg', 'pca'),
            obsp_keys=('connectivities', 'distances'),
            drop_uns=('neighbors',))
        if 'X_pca_harmony' not in self.adata.obsm:
            self._full.obsm.pop('X_pca_harmony', None)
        if self._mask is None:
            self._full.obsp.pop('connectivities', None)
            self._full.obsp.pop('distances', None)
        else:
            excluded_msg = f", {int((~self._mask).sum()):,} cells excluded"
        self.adata = self._full

        # Save inside the worker so the slot doesn't have to write_h5ad
        # on the main thread.
        if self.output_path:
            from pathlib import Path
            self.progress.emit(f"Saving to {Path(self.output_path).name}...")
            self.progress_pct.emit(98)
            Path(self.output_path).parent.mkdir(parents=True, exist_ok=True)
            self.adata.write_h5ad(self.output_path)

        self.progress_pct.emit(100)
        harmony_msg = " + Harmony" if self.harmony and self.harmony_key else ""
        message = (
            f"PCA complete{harmony_msg} -- {n_pcs} components, "
            f"{n_hvg_found:,} HVGs{excluded_msg}")
        return (self.adata, message)

class ClusteringWorker(BaseWorker):
    """
    Worker for running the clustering pipeline.

    Emits 'finished_ok' with '(adata, message)'.
    """

    def __init__(self, adata, params: dict, output_path: str):
        super().__init__()
        self.adata = adata
        self.params = params
        self.output_path = output_path

    def _run(self):
        import scanpy as sc
        from pathlib import Path
        from kosmic.scrna.cluster.neighbors import compute_neighbors
        from kosmic.scrna.cluster.leiden import cluster_leiden
        from kosmic.scrna.annotate.marker_genes import compute_cluster_marker_genes

        self.progress.emit("Preparing dataset...")
        self.progress_pct.emit(2)

        n_pcs = self.params.get('n_pcs', 30)
        n_neighbors = self.params.get('n_neighbors', 15)
        resolution = self.params.get('resolution', 1.0)
        preserve_embeddings = self.params.get('preserve_embeddings', False)

        from kosmic.scrna.cluster.workset import apply_results, working_subset
        from kosmic.scrna.inspect.roles import included_mask

        # Same rule as the PCA step: an excluded arm must not shape the
        # graph or the clusters. The PCA step already left NaN rows in
        # X_pca for those cells, so they could not be clustered anyway.
        # The step reads X (marker genes) and never writes it, so the
        # matrix is shared rather than copied when every cell is included.
        full = self.adata
        mask = (included_mask(self.adata)
                if self.params.get('honour_roles', True) else None)
        self.adata = working_subset(
            full, mask,
            obsm_keys=('X_pca', 'X_pca_harmony', 'X_umap'),
            uns_keys=('neighbors', 'pca', 'umap', 'log1p'),
            obsp_keys=('connectivities', 'distances'))
        if mask is not None:
            self.progress.emit(
                f"Clustering {self.adata.n_obs:,} cells "
                f"({int((~mask).sum()):,} excluded)...")
        else:
            self.progress.emit(f"Clustering {self.adata.n_obs:,} cells...")
        self.progress_pct.emit(5)

        # The graph does not depend on resolution, and rebuilding it is by
        # far the most expensive step (~35s per 40k cells). Reuse it when
        # the parameters have not moved, so trying another resolution costs
        # a Leiden run rather than the whole pipeline.
        from kosmic.scrna.cluster.neighbors import neighbors_are_current
        use_rep = ('X_pca_harmony' if 'X_pca_harmony' in self.adata.obsm
                   else 'X_pca')
        reused = neighbors_are_current(
            self.adata, n_neighbors,
            min(n_pcs, self.adata.obsm[use_rep].shape[1]), use_rep)
        self.progress.emit(
            f"Reusing existing neighbour graph ({n_neighbors} neighbors, "
            f"{n_pcs} PCs)" if reused
            else f"Computing {n_neighbors} neighbors ({n_pcs} PCs)...")
        compute_neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        self.progress_pct.emit(30)

        self.progress.emit(f"Leiden clustering (resolution={resolution})...")
        cluster_leiden(self.adata, resolution=resolution)
        n_clusters = self.adata.obs['leiden'].nunique()
        self.progress.emit(f"Found {n_clusters} clusters")
        self.progress_pct.emit(55)

        # Per-cluster top differential genes -- biologists need these
        # for manual annotation, and computing them now (while we're
        # already iterating over clusters) means the Annotate tab can
        # render them instantly on click instead of triggering a
        # one-off worker on first selection.
        self.progress.emit("Computing top genes per cluster (Wilcoxon)...")
        try:
            compute_cluster_marker_genes(self.adata, cluster_col='leiden')
        except Exception as exc:
            # Don't fail the whole clustering if rank_genes_groups errors;
            # Annotate tab has a lazy on-click fallback.
            self.progress.emit(f"  warning: rank_genes_groups failed ({exc})")
        self.progress_pct.emit(75)

        # UMAP is a function of the graph, not of the resolution: if the
        # graph was reused there is nothing for it to recompute, and it
        # costs ~16s per 40k cells. Changing resolution recolours the same
        # embedding, which is what you want to compare resolutions on
        # anyway.
        if 'X_umap' in self.adata.obsm and (preserve_embeddings or reused):
            self.progress.emit("Reusing existing UMAP")
        else:
            self.progress.emit("Computing UMAP...")
            sc.tl.umap(self.adata)
        self.progress_pct.emit(90)

        self.progress.emit("Writing results back onto the dataset...")
        apply_results(
            full, self.adata, mask,
            obsm_keys=('X_umap',), obs_cols=('leiden',),
            uns_keys=('neighbors', 'rank_genes_groups', 'leiden', 'umap'),
            obsp_keys=('connectivities', 'distances'))
        self.adata = full

        self.progress.emit(f"Saving to {self.output_path}...")
        Path(self.output_path).parent.mkdir(parents=True, exist_ok=True)
        self.adata.write_h5ad(self.output_path)

        self.progress_pct.emit(100)
        message = (
            f"Clustering complete: {n_clusters} clusters, "
            f"{self.adata.n_obs:,} cells")
        return (self.adata, message)


class ResolutionSweepWorker(BaseWorker):
    """
    Worker for sweeping multiple Leiden resolutions and generating UMAP previews.

    Emits 'finished_ok' with '(results_list, message)'. Cancellation
    raises so the caller's 'failed' slot surfaces it uniformly.
    """

    def __init__(self, adata, resolutions: list, output_dir: str):
        super().__init__()
        self.adata = adata
        self.resolutions = sorted(resolutions)
        self.output_dir = output_dir
        self._cancelled = False

    def cancel(self):
        self._cancelled = True

    def _run(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from pathlib import Path
        from kosmic.scrna.cluster.leiden import sweep_leiden_resolutions

        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        if self._cancelled:
            raise RuntimeError("Sweep cancelled")

        self.progress.emit(f"Sweeping {len(self.resolutions)} resolutions...")
        sweep_results = sweep_leiden_resolutions(self.adata, self.resolutions)

        umap_coords = self.adata.obsm['X_umap']
        total = len(sweep_results)
        results = []

        for step, sr in enumerate(sweep_results):
            if self._cancelled:
                raise RuntimeError("Sweep cancelled")

            r = sr['resolution']
            labels = sr['labels']
            n_clusters = sr['n_clusters']
            min_size = sr['min_size']

            self.progress.emit(
                f"Rendering thumbnail for resolution {r:.2f} ({step + 1}/{total})...")

            umap_path = str(Path(self.output_dir) / f"sweep_r{r:.2f}.png")
            try:
                fig, ax = plt.subplots(figsize=(3.2, 2.6), dpi=100)
                unique_labels = sorted(np.unique(labels), key=lambda x: int(x) if str(x).isdigit() else 0)
                cmap = plt.get_cmap('tab20', max(len(unique_labels), 1))
                label_to_int = {lb: i for i, lb in enumerate(unique_labels)}
                colors = [cmap(label_to_int[lb]) for lb in labels]
                ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
                           c=colors, s=0.5, alpha=0.6, linewidths=0, rasterized=True)
                ax.set_title(f"res={r:.2f} ({n_clusters} clusters)", fontsize=7, pad=2)
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xlabel("")
                ax.set_ylabel("")
                for spine in ax.spines.values():
                    spine.set_visible(False)
                fig.tight_layout(pad=0.3)
                fig.savefig(umap_path, dpi=100, bbox_inches='tight')
                plt.close(fig)
            except Exception:
                umap_path = ""
                plt.close('all')

            results.append({
                'resolution': r,
                'n_clusters': n_clusters,
                'min_size': min_size,
                'labels': labels,
                'umap_path': umap_path,
            })

            self.progress_pct.emit(int((step + 1) / total * 100))

        message = f"Sweep complete: {total} resolutions evaluated"
        return (results, message)


class ClustreWorker(BaseWorker):
    """
    Worker for Leiden stability analysis across resolutions.

    Emits 'finished_ok' with '(output_path, message)'.
    """

    def __init__(self, adata, resolutions: list, output_dir: str):
        super().__init__()
        self.adata = adata
        self.resolutions = sorted(resolutions)
        self.output_dir = output_dir

    def _run(self):
        import matplotlib.pyplot as plt
        from pathlib import Path
        from kosmic.scrna.cluster.clustree import compute_clustree_data
        from kosmic.visualisation.scrna.clustree import create_clustree_figure

        output_dir = Path(self.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.progress.emit(
            f"Computing stability data for {len(self.resolutions)} resolutions...")
        self.progress_pct.emit(10)

        data = compute_clustree_data(self.adata, self.resolutions)
        self.progress_pct.emit(65)

        res_map = "  |  ".join(
            f"L{i+1}=res{r:.2f} ({len(set(lbl))}cl)"
            for i, (r, lbl) in enumerate(data['leiden_per_level'])
        )
        self.progress.emit(f"Cluster counts: {res_map}")
        self.progress.emit("Rendering stability graph...")
        self.progress_pct.emit(70)

        fig = create_clustree_figure(data['leiden_per_level'], data['level_positions'])

        output_path = str(output_dir / "stability_analysis.png")
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)

        self.progress_pct.emit(100)
        cluster_summary = " \u2192 ".join(
            str(len(set(lbl))) for _, lbl in data['leiden_per_level']
        )
        message = (
            f"Stability analysis complete  \u2014  clusters: {cluster_summary}")
        return (output_path, message)


class TSNEWorker(BaseWorker):
    """
    Background worker for computing t-SNE embedding using openTSNE (FFT-accelerated).

    Emits 'finished_ok' with the adata.
    """

    def __init__(self, adata, perplexity=100):
        super().__init__()
        self.adata = adata
        self.perplexity = perplexity

    def _run(self):
        from kosmic.scrna.cluster.tsne import compute_tsne

        compute_tsne(self.adata, perplexity=self.perplexity)
        return self.adata


# ---------------------------------------------------------------------------
# Interactive UMAP Widget
# ---------------------------------------------------------------------------

# 20 distinct colors matching scanpy/matplotlib tab20
_TAB20 = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
    '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
]


def _viridis_color(t: float) -> str:
    """Simple viridis-like gradient: dark purple → teal → yellow."""
    t = max(0.0, min(1.0, t))
    if t < 0.5:
        s = t / 0.5
        r = int(68 + s * (32 - 68))
        g = int(1 + s * (144 - 1))
        b = int(84 + s * (140 - 84))
    else:
        s = (t - 0.5) / 0.5
        r = int(32 + s * (253 - 32))
        g = int(144 + s * (231 - 144))
        b = int(140 + s * (37 - 140))
    return f"#{r:02x}{g:02x}{b:02x}"


def _label_sort_key(label):
    """Numbered clusters first in numeric order, then names; 'nan' last.

    A study whose excluded donors were never clustered carries NaN in
    obs['leiden'], which reaches here as the string 'nan'. Mixing int and
    str keys made sorted() raise and the plot never drew.
    """
    text = str(label)
    if text.isdigit():
        return (0, int(text), '')
    if text.lower() in ('nan', 'none', ''):
        return (2, 0, text)
    return (1, 0, text)


class _UMAPWidget(QWidget):
    """
    Interactive pyqtgraph scatter plot for UMAP/t-SNE embeddings.

    Supports coloring by categorical obs columns (discrete tab20 palette),
    continuous obs columns (viridis gradient), or gene expression.
    """

    color_by_changed = pyqtSignal(str)  # emitted when user changes selection
    tsne_computed = pyqtSignal(object)  # emitted when background t-SNE finishes

    def __init__(self, parent=None):
        super().__init__(parent)
        self._coords = None
        self._values = None       # current color values (str for categorical, float for continuous)
        self._is_continuous = False
        self._color_key = ""      # current color-by key
        self._adata = None        # reference for gene lookups
        self._embedding = 'umap'  # 'umap' or 'tsne'
        self._tsne_worker = None
        self._computing_tsne = False

        layout = borderless(QVBoxLayout, self)

        # Controls are created here but added to the ClusterTab sidebar
        # via setup_sidebar_controls()

        # UMAP / t-SNE toggle
        self._umap_btn = QPushButton("UMAP")
        self._tsne_btn = QPushButton("t-SNE")
        for btn in (self._umap_btn, self._tsne_btn):
            btn.setFixedWidth(60)
            btn.setCheckable(True)
            btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._umap_btn.setChecked(True)
        self._umap_btn.clicked.connect(lambda: self._switch_embedding('umap'))
        self._tsne_btn.clicked.connect(lambda: self._switch_embedding('tsne'))
        self._update_toggle_style()

        # Recompute t-SNE button (visible only in t-SNE mode)
        self._recompute_tsne_btn = QPushButton("↻")
        self._recompute_tsne_btn.setFixedWidth(28)
        self._recompute_tsne_btn.setToolTip("Recompute t-SNE with optimized parameters")
        self._recompute_tsne_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._recompute_tsne_btn.clicked.connect(self._recompute_tsne)
        self._recompute_tsne_btn.setVisible(False)

        self._color_combo = NoScrollComboBox()
        self._color_combo.currentTextChanged.connect(self._on_color_changed)

        self._gene_edit = QLineEdit()
        self._gene_edit.setPlaceholderText("Search gene (e.g. PTPRC, CD3E)")
        self._gene_edit.setToolTip(
            "Colour the embedding by a gene's expression.\n"
            "Type a gene symbol and press Enter.")
        self._gene_edit.returnPressed.connect(self._on_gene_search)
        self._gene_completer = QCompleter()
        self._gene_completer.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        self._gene_completer.setMaxVisibleItems(10)
        self._gene_edit.setCompleter(self._gene_completer)

        # Info label
        self._info_label = QLabel()

        # Key button — toggles legend overlay
        self._key_btn = QPushButton("Hide Legend")
        self._key_btn.setProperty("role", "toggle")
        self._key_btn.setProperty("state", "inactive")
        self._key_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._key_btn.clicked.connect(self._toggle_legend)
        self._key_btn.setVisible(False)

        # Legend state
        self._legend_entries = []  # list of (color_hex, label_str)
        self._legend_overlay = None
        self._legend_auto_show = True  # auto-show legend after categorical plot

        # ── pyqtgraph plot ──
        self._plot = InteractivePlot(
            title='UMAP',
            left_label='UMAP 2',
            bottom_label='UMAP 1',
            unavailable_message='Run clustering to see the UMAP embedding.',
        )
        self._plot.setMenuEnabled(False)
        self._plot.setMinimumHeight(300)
        self._plot.setAspectLocked(True)

        self._scatter = pg.ScatterPlotItem(size=2, pen=pg.mkPen(None))
        self._plot.addItem(self._scatter)

        # Hover tooltip
        self._tooltip = pg.TextItem(anchor=(0, 1))
        self._tooltip.setVisible(False)
        self._plot.addItem(self._tooltip)
        self._scatter.sigHovered.connect(self._on_hover)

        layout.addWidget(self._plot, 1)

    # ── Public API ──

    def set_adata(self, adata):
        """Store reference and populate color-by combo with obs columns."""
        self._adata = adata
        self._populate_color_combo()
        self._populate_gene_completer()

    def build_toolbar(self):
        """Build a horizontal toolbar widget with visualization controls."""
        toolbar = QWidget()
        row = QHBoxLayout(toolbar)
        row.setContentsMargins(0, 2, 0, 6)
        row.setSpacing(8)

        row.addWidget(self._umap_btn)
        row.addWidget(self._tsne_btn)
        row.addWidget(self._recompute_tsne_btn)

        row.addSpacing(8)
        row.addWidget(QLabel("Color by"))
        self._color_combo.setMinimumWidth(130)
        row.addWidget(self._color_combo)

        row.addSpacing(8)
        self._gene_edit.setFixedWidth(220)
        row.addWidget(self._gene_edit)

        row.addStretch()
        row.addWidget(self._info_label)
        row.addWidget(self._key_btn)

        return toolbar

    def set_data_categorical(self, coords: np.ndarray, labels: np.ndarray,
                              info_text: str = ""):
        """Plot embedding colored by categorical labels (clusters, sample, etc)."""
        # Clear legend overlay if showing
        if self._legend_overlay is not None:
            self._legend_overlay.deleteLater()
            self._legend_overlay = None
        self._coords = coords
        self._values = labels
        self._is_continuous = False

        label_arr = np.asarray(labels)
        unique = sorted(set(label_arr.tolist()), key=_label_sort_key)
        n_categories = len(unique)
        label_to_idx = {lb: i for i, lb in enumerate(unique)}

        # Legend counts are over every cell, whatever is drawn below.
        self._legend_entries = [
            (_TAB20[i % len(_TAB20)],
             f"{lb}   ({int((label_arr == lb).sum()):,})")
            for i, lb in enumerate(unique)
        ]
        self._key_btn.setVisible(True)

        n_total = len(coords)
        # Shuffle draw order so no single label (e.g. whichever was concatenated
        # last when the AnnData was built) systematically paints over the others
        # in overlap/boundary regions. The same permutation also picks the
        # displayed subset when the dataset is over the display cap.
        order = self._draw_order(n_total)
        n_shown = len(order)

        # Alpha scales down with cell count: opaque for small datasets, translucent
        # for large ones so overlapping points blend instead of fully occluding
        # whichever happens to be drawn last (which otherwise makes densely
        # overplotted regions look like arbitrary noise even when the underlying
        # clusters separate cleanly -- most visible above ~100k cells).
        alpha = int(max(60, min(255, 3_000_000 / n_shown)))

        # One brush per category, shared by every point of that category:
        # building a dict and a QBrush per cell was the slow part for
        # atlas-sized studies.
        brushes = []
        for i in range(n_categories):
            color = pg.mkColor(_TAB20[i % len(_TAB20)])
            color.setAlpha(alpha)
            brushes.append(pg.mkBrush(color))
        shown_labels = label_arr[order].tolist()
        idx = [label_to_idx[lb] for lb in shown_labels]

        self._scatter.setData(
            x=coords[order, 0].astype(float), y=coords[order, 1].astype(float),
            brush=[brushes[i] for i in idx],
            data=shown_labels,
        )
        self._plot.hide_unavailable_message()
        # Auto-size points for cell count
        auto_size = max(1, min(4, 5000 / n_shown))
        self._point_size = int(auto_size)
        self._scatter.setSize(self._point_size)

        if not info_text:
            info_text = f"{self._shown_text(n_shown, n_total)}  |  {n_categories} groups"
        self._info_label.setText(info_text)
        self._plot.update()

        # Auto-show legend on first categorical plot
        if self._legend_auto_show and self._legend_entries:
            QTimer.singleShot(100, self._ensure_legend_shown)

    def set_data_continuous(self, coords: np.ndarray, values: np.ndarray,
                             info_text: str = ""):
        """Plot embedding colored by continuous values (gene expr, QC metric)."""
        self._coords = coords
        self._values = values
        self._is_continuous = True
        self._legend_entries = []
        self._key_btn.setVisible(False)

        values = np.asarray(values, dtype=float)
        vmin = np.nanmin(values)
        vmax = np.nanmax(values)
        rng = vmax - vmin if vmax > vmin else 1.0

        n_total = len(coords)
        shown = self._draw_order(n_total)
        # Sort so high values render on top
        shown = shown[np.argsort(values[shown], kind='stable')]
        n_shown = len(shown)
        # See set_data_categorical for why alpha scales down with cell count.
        alpha = int(max(60, min(255, 3_000_000 / n_shown)))

        # Quantise the colour scale to 64 shared brushes instead of one per cell.
        n_levels = 64
        levels = (values[shown] - vmin) / rng * (n_levels - 1)
        levels = np.clip(np.nan_to_num(levels, nan=0.0), 0, n_levels - 1).astype(np.int64)
        brushes = []
        for k in range(n_levels):
            color = pg.mkColor(_viridis_color(k / (n_levels - 1)))
            color.setAlpha(alpha)
            brushes.append(pg.mkBrush(color))

        self._scatter.setData(
            x=coords[shown, 0].astype(float), y=coords[shown, 1].astype(float),
            brush=[brushes[k] for k in levels.tolist()],
            data=values[shown].tolist(),
        )
        self._plot.hide_unavailable_message()
        auto_size = max(1, min(4, 5000 / n_shown))
        self._point_size = int(auto_size)
        self._scatter.setSize(self._point_size)

        if not info_text:
            info_text = (f"{self._shown_text(n_shown, n_total)}  |  "
                         f"range {vmin:.2f} \u2013 {vmax:.2f}")
        self._info_label.setText(info_text)

    @staticmethod
    def _draw_order(n_total: int) -> np.ndarray:
        """Indices to draw, shuffled; a uniform random subset above the cap.

        Every cell is still counted in the legend and the hover tooltip;
        only the drawn points are capped, because pyqtgraph repaints every
        point on each pan or zoom and an atlas-sized study made that
        unusable. The seed is fixed so the same cells are drawn each time.
        """
        order = np.random.default_rng(0).permutation(n_total)
        cap = int(UI_EMBEDDING_MAX_POINTS)
        if cap > 0 and n_total > cap:
            order = order[:cap]
        return order

    @staticmethod
    def _shown_text(n_shown: int, n_total: int) -> str:
        if n_shown < n_total:
            return f"{n_shown:,} of {n_total:,} cells drawn"
        return f"{n_total:,} cells"

    def clear_plot(self):
        self._scatter.clear()
        self._info_label.setText("No data \u2014 run clustering to see preview")
        self._coords = None
        self._values = None
        self._legend_entries = []
        self._key_btn.setVisible(False)
        if self._legend_overlay is not None:
            self._legend_overlay.deleteLater()
            self._legend_overlay = None

    # ── Color-by selector ──

    def _populate_color_combo(self):
        self._color_combo.blockSignals(True)
        prev = self._color_combo.currentText()
        self._color_combo.clear()

        if self._adata is None:
            self._color_combo.blockSignals(False)
            return

        # Priority categorical columns first
        priority_cat = ['leiden', 'clusters', 'seurat_clusters', 'cell_type',
                        'sample', 'Sample', 'donor', 'Donor', 'batch', 'Batch',
                        'condition', 'Condition', 'orig.ident', 'patient', 'Patient']
        added = []
        for col in priority_cat:
            if col in self._adata.obs.columns:
                self._color_combo.addItem(col)
                added.append(col)

        # Other categorical obs columns (2-200 unique values)
        for col in self._adata.obs.columns:
            if col not in added:
                try:
                    n = self._adata.obs[col].nunique()
                    if 2 <= n <= 200:
                        self._color_combo.addItem(col)
                        added.append(col)
                except Exception:
                    pass

        # Separator + continuous QC metrics
        if added:
            self._color_combo.insertSeparator(self._color_combo.count())
        for col in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt',
                     'n_genes', 'n_counts', 'doublet_score']:
            if col in self._adata.obs.columns and col not in added:
                self._color_combo.addItem(col)
                added.append(col)

        # Restore previous selection if still available
        idx = self._color_combo.findText(prev)
        if idx >= 0:
            self._color_combo.setCurrentIndex(idx)
        elif self._color_combo.count() > 0:
            self._color_combo.setCurrentIndex(0)

        self._color_combo.blockSignals(False)

    def _populate_gene_completer(self):
        if self._adata is None:
            return
        try:
            from PyQt6.QtCore import QStringListModel
            genes = list(self._adata.var_names[:5000])  # cap for performance
            model = QStringListModel(genes)
            self._gene_completer.setModel(model)
        except Exception:
            pass

    def _recompute_tsne(self):
        """Force recompute t-SNE even if it already exists."""
        if self._adata is None or self._computing_tsne:
            return
        # Remove cached t-SNE so _switch_embedding triggers recompute
        if 'X_tsne' in self._adata.obsm:
            del self._adata.obsm['X_tsne']
        self._embedding = 'umap'  # reset so _switch_embedding doesn't early-return
        self._switch_embedding('tsne')

    def _switch_embedding(self, embedding: str):
        """Switch between UMAP and t-SNE views."""
        if embedding == self._embedding:
            return

        self._embedding = embedding
        self._umap_btn.setChecked(embedding == 'umap')
        self._tsne_btn.setChecked(embedding == 'tsne')
        self._recompute_tsne_btn.setVisible(embedding == 'tsne')
        self._update_toggle_style()

        obsm_key = 'X_umap' if embedding == 'umap' else 'X_tsne'
        label = 'UMAP' if embedding == 'umap' else 't-SNE'
        self._plot.setLabel('bottom', f'{label} 1')
        self._plot.setLabel('left', f'{label} 2')

        if self._adata is None:
            return

        if obsm_key in self._adata.obsm:
            self._recolor()
        elif embedding == 'tsne' and not self._computing_tsne:
            # Compute t-SNE lazily in background
            self._computing_tsne = True
            self._scatter.setData([], [])  # clear plot
            self._info_label.setText("Computing t-SNE — this may take a minute...")
            self._tsne_btn.setText("t-SNE ...")
            self._tsne_btn.setProperty("role", "toggle")
            self._tsne_btn.setProperty("state", "warning")
            _refresh_style(self._tsne_btn)
            # Higher perplexity for large datasets gives tighter clusters
            perplexity = min(150, max(30, self._adata.n_obs // 200))
            self._tsne_worker = TSNEWorker(self._adata, perplexity=perplexity)
            run_worker(
                self._tsne_worker,
                on_finished=self._on_tsne_finished,
                on_failed=self._on_tsne_failed,
            )
        elif embedding == 'tsne':
            self._info_label.setText("Computing t-SNE — this may take a minute...")

    def _on_tsne_finished(self, adata):
        self._computing_tsne = False
        self._tsne_btn.setText("t-SNE")
        self._update_toggle_style()
        if adata is not None:
            self._adata = adata
            self.tsne_computed.emit(adata)
            self._info_label.setText("")
            if self._embedding == 'tsne':
                self._recolor()

    def _on_tsne_failed(self, _message: str):
        self._computing_tsne = False
        self._tsne_btn.setText("t-SNE")
        self._update_toggle_style()
        self._info_label.setText("t-SNE failed")

    def _update_toggle_style(self):
        for btn, active in [(self._umap_btn, self._embedding == 'umap'),
                            (self._tsne_btn, self._embedding == 'tsne')]:
            if active:
                btn.setProperty("role", "toggle")
                btn.setProperty("state", "active")
                _refresh_style(btn)
            else:
                btn.setProperty("role", "toggle")
                btn.setProperty("state", "inactive")
                _refresh_style(btn)

    @property
    def _obsm_key(self):
        return 'X_umap' if self._embedding == 'umap' else 'X_tsne'

    def _on_color_changed(self, text: str):
        if not text or self._adata is None or self._obsm_key not in self._adata.obsm:
            return
        self._color_key = text
        self._gene_edit.clear()
        self._recolor()

    def _on_gene_search(self):
        gene = self._gene_edit.text().strip()
        if not gene or self._adata is None or self._obsm_key not in self._adata.obsm:
            return

        if gene not in self._adata.var_names:
            # Try case-insensitive match
            matches = [g for g in self._adata.var_names if g.upper() == gene.upper()]
            if matches:
                gene = matches[0]
                self._gene_edit.setText(gene)
            else:
                self._info_label.setText(f"Gene '{gene}' not found")
                return

        self._color_key = f"gene:{gene}"
        coords = self._adata.obsm[self._obsm_key]

        # Get expression values
        gene_idx = list(self._adata.var_names).index(gene)
        from scipy import sparse
        if sparse.issparse(self._adata.X):
            expr = np.asarray(self._adata.X[:, gene_idx].todense()).flatten()
        else:
            expr = self._adata.X[:, gene_idx].flatten()

        self.set_data_continuous(
            coords, expr,
            f"{gene}  |  {len(coords):,} cells  |  "
            f"expr {expr.min():.2f} \u2013 {expr.max():.2f}"
        )

    def _recolor(self):
        """Re-render the embedding with the current color-by selection."""
        if self._adata is None or self._obsm_key not in self._adata.obsm:
            return

        coords = self._adata.obsm[self._obsm_key]
        col = self._color_key

        # Handle gene expression (set via _on_gene_search)
        if col.startswith('gene:'):
            gene = col[5:]
            if gene in self._adata.var_names:
                gene_idx = list(self._adata.var_names).index(gene)
                from scipy import sparse
                if sparse.issparse(self._adata.X):
                    expr = np.asarray(self._adata.X[:, gene_idx].todense()).flatten()
                else:
                    expr = self._adata.X[:, gene_idx].flatten()
                self.set_data_continuous(
                    coords, expr,
                    f"{gene}  |  {len(coords):,} cells  |  "
                    f"expr {expr.min():.2f} – {expr.max():.2f}"
                )
            return

        if col not in self._adata.obs.columns:
            # Fall back to first available: leiden, combo selection, or first obs column
            for fallback in ['leiden', 'clusters', 'seurat_clusters']:
                if fallback in self._adata.obs.columns:
                    col = fallback
                    break
            else:
                combo_text = self._color_combo.currentText()
                if combo_text and combo_text in self._adata.obs.columns:
                    col = combo_text
                else:
                    return
            self._color_key = col

        series = self._adata.obs[col]

        # Decide categorical vs continuous
        if series.dtype.name == 'category' or series.dtype == object or series.nunique() <= 50:
            labels = series.astype(str).values
            n_groups = len(set(labels))
            self.set_data_categorical(
                coords, labels,
                f"{col}  |  {len(coords):,} cells  |  {n_groups} groups"
            )
        else:
            values = series.values.astype(float)
            self.set_data_continuous(
                coords, values,
                f"{col}  |  {len(coords):,} cells  |  "
                f"range {np.nanmin(values):.1f} \u2013 {np.nanmax(values):.1f}"
            )

    # ── Hover ──

    def _on_hover(self, _plot, points, _ev):
        if points:
            pt = points[0]
            val = pt.data()
            if self._is_continuous:
                self._tooltip.setText(f"{float(val):.3f}")
            else:
                if self._values is not None:
                    count = np.sum(np.array(self._values) == val)
                    self._tooltip.setText(f"{val}\n{count:,} cells")
                else:
                    self._tooltip.setText(str(val))
            self._tooltip.setPos(pt.pos())
            self._tooltip.setVisible(True)
        else:
            self._tooltip.setVisible(False)

    def _toggle_legend(self):
        """Toggle the legend overlay on/off."""
        if self._legend_overlay is not None:
            self._legend_overlay.deleteLater()
            self._legend_overlay = None
            self._key_btn.setText("Show Legend")
            self._legend_auto_show = False  # user explicitly hid it
        else:
            self._build_legend()
            self._legend_auto_show = True

    def _ensure_legend_shown(self):
        """Show legend if not already visible."""
        if self._legend_overlay is None and self._legend_entries:
            self._build_legend()

    def _build_legend(self):
        """Build and display the legend overlay on top of the plot."""
        if not self._legend_entries:
            return
        # Remove existing if any
        if self._legend_overlay is not None:
            self._legend_overlay.deleteLater()
            self._legend_overlay = None

        overlay = QWidget(self._plot)
        overlay.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents, False)
        lay = QVBoxLayout(overlay)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(2)

        overlay.setObjectName("legend_overlay")

        for color, label in self._legend_entries:
            row = QHBoxLayout()
            row.setContentsMargins(0, 0, 0, 0)
            row.setSpacing(6)

            # Runtime-variable swatch colour comes from the plot's
            # category palette, so we paint a solid-colour pixmap on
            # the swatch QLabel rather than using a stylesheet.  This
            # keeps the "no inline setStyleSheet" rule intact while
            # still allowing per-instance colouring.
            swatch = QLabel()
            swatch.setFixedSize(10, 10)
            swatch_pix = QPixmap(10, 10)
            swatch_pix.fill(QColor(color))
            swatch.setPixmap(swatch_pix)
            row.addWidget(swatch)

            text = QLabel(label)
            text.setProperty("role", "legend_text")
            row.addWidget(text, 1)

            container = QWidget()
            container.setLayout(row)
            lay.addWidget(container)

        overlay.adjustSize()
        plot_rect = self._plot.rect()
        overlay.move(plot_rect.width() - overlay.width() - 10, 10)
        overlay.show()
        overlay.raise_()
        self._legend_overlay = overlay
        self._key_btn.setText("Hide Legend")

    def refresh_theme(self):
        style_pg_plot(self._plot, title='UMAP',
                      left_label='UMAP 2', bottom_label='UMAP 1')
        # Re-apply role-based styles (Qt re-polishes after setProperty).
        self._update_toggle_style()
        _refresh_style(self._info_label)
        _refresh_style(self._key_btn)


# ---------------------------------------------------------------------------
# Sweep Thumbnail Strip
# ---------------------------------------------------------------------------

class _SweepThumbnail(QLabel):
    """Clickable thumbnail for a single resolution in the sweep strip."""
    clicked = pyqtSignal(int)

    def __init__(self, index: int, parent=None):
        super().__init__(parent)
        self._index = index
        self._selected = False
        self.setFixedSize(120, 96)
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setScaledContents(False)
        self._update_style()

    def set_selected(self, selected: bool):
        self._selected = selected
        self._update_style()

    def _update_style(self):
        self.setProperty("selected", "true" if self._selected else "false")
        _refresh_style(self)

    def mousePressEvent(self, ev):
        self.clicked.emit(self._index)


class _SweepStrip(QWidget):
    """Horizontal scrollable strip of resolution sweep thumbnails."""
    resolution_selected = pyqtSignal(int)  # index

    def __init__(self, parent=None):
        super().__init__(parent)
        self._thumbnails: list[_SweepThumbnail] = []
        self._current = -1

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)

        # Header
        self._header = QLabel("Resolution Sweep")
        font = QFont()
        font.setBold(True)
        font.setPointSize(10)
        self._header.setFont(font)
        layout.addWidget(self._header)

        # Scroll area for thumbnails
        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._scroll.setFixedHeight(120)
        self._scroll.setFrameShape(QFrame.Shape.NoFrame)

        self._strip_widget = QWidget()
        self._strip_layout = QHBoxLayout(self._strip_widget)
        self._strip_layout.setContentsMargins(0, 0, 0, 0)
        self._strip_layout.setSpacing(6)
        self._strip_layout.addStretch()

        self._scroll.setWidget(self._strip_widget)
        layout.addWidget(self._scroll)

        self.setVisible(False)

    def clear(self):
        """Remove all thumbnails and hide the strip."""
        for t in self._thumbnails:
            self._strip_layout.removeWidget(t)
            t.deleteLater()
        self._thumbnails.clear()
        self._current = -1
        self.setVisible(False)

    def set_results(self, results: list):
        """Populate strip with sweep result thumbnails."""
        # Clear existing
        for t in self._thumbnails:
            self._strip_layout.removeWidget(t)
            t.deleteLater()
        self._thumbnails.clear()
        self._current = -1

        # Remove trailing stretch
        self._strip_layout.takeAt(self._strip_layout.count() - 1)

        for i, r in enumerate(results):
            thumb = _SweepThumbnail(i)
            path = r.get('umap_path', '')
            if path:
                pixmap = QPixmap(path)
                if not pixmap.isNull():
                    thumb.setPixmap(pixmap.scaled(
                        116, 92,
                        Qt.AspectRatioMode.KeepAspectRatio,
                        Qt.TransformationMode.SmoothTransformation
                    ))
                else:
                    thumb.setText(f"res={r['resolution']:.2f}")
            else:
                thumb.setText(f"res={r['resolution']:.2f}")

            thumb.setToolTip(f"Resolution {r['resolution']:.2f}\n{r['n_clusters']} clusters")
            thumb.clicked.connect(self._on_thumb_clicked)
            self._strip_layout.addWidget(thumb)
            self._thumbnails.append(thumb)

        self._strip_layout.addStretch()
        self.setVisible(True)

        if results:
            self.select(0)

    def select(self, index: int):
        if index < 0 or index >= len(self._thumbnails):
            return
        if self._current >= 0 and self._current < len(self._thumbnails):
            self._thumbnails[self._current].set_selected(False)
        self._current = index
        self._thumbnails[index].set_selected(True)
        self._scroll.ensureWidgetVisible(self._thumbnails[index])

    def _on_thumb_clicked(self, index: int):
        self.select(index)
        self.resolution_selected.emit(index)

    def refresh_theme(self):
        for t in self._thumbnails:
            t._update_style()


# ---------------------------------------------------------------------------
# ClusterTab
# ---------------------------------------------------------------------------

class ClusterTab(SidebarPage):
    """Preview-first clustering tab with integrated PCA, interactive UMAP, and sweep strip."""

    help_id = "scrna/cluster"

    clustering_complete = pyqtSignal(str)
    pca_complete = pyqtSignal(object)   # emits adata after PCA (+Harmony)
    log_message = pyqtSignal(str)

    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.project_dir = None
        self.adata = None
        self.file_path = None
        self.h5ad_path = None
        self.clustering_worker = None
        self.sweep_worker = None
        self.clustree_worker = None
        self._pca_worker = None
        self._sweep_results = []
        self._suggested_pcs = None

        self._setup_ui()

    def _setup_ui(self):
        # Sidebar lives on self.sidebar_layout (provided by SidebarPage).
        left_layout = self.sidebar_layout

        # ═══ Analysis summary: one card per stage; each card expands in
        # place (accordion) to show that stage's settings. ═══
        left_layout.addWidget(SectionHeader("ANALYSIS SUMMARY"))
        self._accordion = StageAccordion(left_layout)

        self._pca_card = StageSummaryCard("PCA", icon="bar-chart-2")
        self._pca_card.set_summary(["Not run yet"])
        left_layout.addWidget(self._pca_card)
        self._accordion.add_card('pca', self._pca_card)

        self._cluster_card = StageSummaryCard("Clustering", icon="umap")
        self._cluster_card.set_summary(["Not run yet"])
        left_layout.addWidget(self._cluster_card)
        self._accordion.add_card('clustering', self._cluster_card)

        self._annot_card = StageSummaryCard("Annotate", icon="cells")
        self._annot_card.set_summary(["Not run yet"])
        left_layout.addWidget(self._annot_card)
        self._accordion.add_card('annotate', self._annot_card)

        left_layout.addStretch()
        self._accordion.finalize()
        report_btn = SecondaryButton("View full report")
        report_btn.setToolTip(
            "Show the methods recorded for this study (provenance).")
        report_btn.clicked.connect(
            lambda: show_methods_report(self, self.h5ad_path))
        left_layout.addWidget(report_btn)

        # ═══ Drawer page: PCA settings ═══
        pca_page = QWidget()
        pca_page_layout = QVBoxLayout(pca_page)
        pca_page_layout.setContentsMargins(0, 4, 0, 0)
        pca_page_layout.setSpacing(8)

        pca_grid = QGridLayout()
        pca_grid.setSpacing(4)
        pca_grid.addWidget(QLabel("HVG count"), 0, 0)
        self.hvg_spin = NoScrollSpinBox()
        self.hvg_spin.setRange(500, 5000)
        self.hvg_spin.setValue(2000)
        self.hvg_spin.setSingleStep(500)
        self.hvg_spin.setToolTip("Number of highly variable genes. 2,000 is standard.")
        pca_grid.addWidget(self.hvg_spin, 0, 1)

        pca_grid.addWidget(QLabel("Compute PCs"), 1, 0)
        self.pca_compute_spin = NoScrollSpinBox()
        self.pca_compute_spin.setRange(10, 100)
        self.pca_compute_spin.setValue(50)
        self.pca_compute_spin.setToolTip("PCs to compute. Check elbow plot to choose how many to use.")
        pca_grid.addWidget(self.pca_compute_spin, 1, 1)
        pca_page_layout.addLayout(pca_grid)

        self.run_pca_btn = PrimaryButton("Run PCA")
        self.run_pca_btn.setMinimumHeight(32)
        self.run_pca_btn.setEnabled(False)
        self.run_pca_btn.clicked.connect(self._run_pca)
        pca_page_layout.addWidget(self.run_pca_btn)

        self.pca_status = SecondaryLabel("")
        self.pca_status.setWordWrap(True)
        pca_page_layout.addWidget(self.pca_status)

        # ═══ Drawer page: Clustering settings ═══
        cluster_page = QWidget()
        cluster_page_layout = QVBoxLayout(cluster_page)
        cluster_page_layout.setContentsMargins(0, 4, 0, 0)
        cluster_page_layout.setSpacing(8)

        graph_grid = QGridLayout()
        graph_grid.setSpacing(4)

        graph_grid.addWidget(QLabel("PCs to use"), 0, 0)
        self.n_pcs_spin = NoScrollSpinBox()
        self.n_pcs_spin.setRange(5, 100)
        self.n_pcs_spin.setValue(30)
        self.n_pcs_spin.setToolTip("PCA components for neighbor graph. Auto-set from elbow plot.")
        graph_grid.addWidget(self.n_pcs_spin, 0, 1)

        graph_grid.addWidget(QLabel("Neighbors (k)"), 1, 0)
        self.neighbors_spin = NoScrollSpinBox()
        self.neighbors_spin.setRange(5, 100)
        self.neighbors_spin.setValue(15)
        self.neighbors_spin.setToolTip("KNN neighbors. 10-15 = finer, 30-50 = broader.")
        graph_grid.addWidget(self.neighbors_spin, 1, 1)

        graph_grid.addWidget(QLabel("Resolution"), 2, 0)
        self.resolution_spin = NoScrollDoubleSpinBox()
        self.resolution_spin.setRange(0.1, 5.0)
        self.resolution_spin.setValue(0.5)
        self.resolution_spin.setSingleStep(0.1)
        self.resolution_spin.setToolTip(
            "How finely cells are split into clusters. Start at 0.5; raise\n"
            "it if distinct cell types share one cluster, lower it if one\n"
            "cell type splits across several. Check the choice with marker\n"
            "genes (Annotate) or the resolution sweep below.")
        graph_grid.addWidget(self.resolution_spin, 2, 1)

        cluster_page_layout.addLayout(graph_grid)

        self.preserve_embed_check = QCheckBox("Keep existing embeddings")
        self.preserve_embed_check.setChecked(False)
        self.preserve_embed_check.setToolTip("Keep existing UMAP/t-SNE instead of recomputing.")
        cluster_page_layout.addWidget(self.preserve_embed_check)

        self.run_cluster_btn = PrimaryButton("Run Clustering")
        self.run_cluster_btn.setMinimumHeight(32)
        self.run_cluster_btn.clicked.connect(self._run_clustering)
        self.run_cluster_btn.setEnabled(False)
        cluster_page_layout.addWidget(self.run_cluster_btn)

        cluster_elbow_link = QPushButton("View elbow plot")
        cluster_elbow_link.setProperty("role", "link_accent")
        cluster_elbow_link.setCursor(Qt.CursorShape.PointingHandCursor)
        cluster_elbow_link.clicked.connect(
            lambda: self._right_tabs.setCurrentIndex(0))
        cluster_page_layout.addWidget(
            cluster_elbow_link, alignment=Qt.AlignmentFlag.AlignLeft)

        # Interactive UMAP / t-SNE
        self._umap_widget = _UMAPWidget()

        # ── Annotate (controls from embedded AnnotateTab) ──
        from kosmic.gui.scrna.tabs.annotate_tab import AnnotateTab
        self.annotate_tab = AnnotateTab(self.main_window, embedded=True)
        self.annotate_tab.annotation_complete.connect(self._on_annotation_complete)

        # ── PCA drawer page: batch correction (below the run button) ──
        harmony_label = SectionHeader("BATCH CORRECTION")
        pca_page_layout.addWidget(harmony_label)

        self.honour_roles_check = QCheckBox("Ignore cells marked 'exclude'")
        self.honour_roles_check.setChecked(True)
        self.honour_roles_check.setToolTip(
            "Leave cells whose condition you set to Exclude in Inspect\n"
            "out of HVG selection, PCA, the neighbour graph, clustering\n"
            "and the UMAP. They stay in the file with their metadata;\n"
            "they simply get no coordinates and no cluster." + "\n\n"
            "Untick to embed every cell. Exclude then applies to the\n"
            "differential-expression contrast only, which is what it\n"
            "used to mean everywhere.")
        pca_page_layout.addWidget(self.honour_roles_check)

        self.scale_check = QCheckBox("Scale genes before PCA")
        self.scale_check.setChecked(False)
        self.scale_check.setToolTip(
            "Z-score each gene to mean 0, variance 1 (clipped at 10 SD)\n"
            "so PCA weighs a lowly-expressed transcription factor as\n"
            "heavily as a ribosomal gene. Without it the first PCs\n"
            "track the most abundant genes, because log-scale variance\n"
            "still rises with the mean." + "\n\n"
            "Off by default: scaling makes the matrix dense, which on a\n"
            "large combined dataset costs more memory than everything\n"
            "else put together. HVG selection has already removed most\n"
            "of the low-variance genes it would help with.")
        pca_page_layout.addWidget(self.scale_check)

        self.harmony_check = QCheckBox("Harmony batch correction")
        self.harmony_check.setChecked(False)   # turned on once a batch column is known
        self._harmony_default_applied = False
        self.harmony_check.setToolTip(
            "Correct the PCA space for the batch key (the donor / sample) so "
            "cells of one type from different donors overlap. On by default "
            "for multi-donor data; it changes clustering only, never DE.")
        self.harmony_check.stateChanged.connect(self._on_harmony_toggled)
        pca_page_layout.addWidget(self.harmony_check)

        batch_grid = QGridLayout()
        batch_grid.setSpacing(4)
        batch_grid.addWidget(QLabel("Batch key"), 0, 0)
        self.batch_combo = NoScrollComboBox()
        self.batch_combo.setEditable(True)
        self.batch_combo.setEnabled(False)
        self.batch_combo.setToolTip(
            "The grouping to correct for. On a combined master this is\n"
            "'study' -- the cohort effect is what has to come out. Within\n"
            "one study it is the sample or donor column.")
        batch_grid.addWidget(self.batch_combo, 0, 1)
        batch_grid.addWidget(QLabel("Merge across"), 1, 0)
        self.condition_combo = NoScrollComboBox()
        self.condition_combo.setEditable(True)
        self.condition_combo.setEnabled(False)
        self.condition_combo.addItem("")
        self.condition_combo.setToolTip(
            "Optional second variable to correct, so its groups land in\n"
            "the same clusters. Set it to the condition column when you\n"
            "want one cardiomyocyte cluster rather than a diseased and a\n"
            "healthy one." + "\n\n"
            "Leave blank to keep disease-specific cell states separate.\n"
            "This affects clustering and the UMAP only -- differential\n"
            "expression runs on raw counts and never sees the corrected\n"
            "embedding.")
        batch_grid.addWidget(self.condition_combo, 1, 1)
        pca_page_layout.addLayout(batch_grid)

        self._batch_hint = SecondaryLabel(
            "Correcting on 'Merge across' pulls its groups together, it does not "
            "shield them. Use it to group by cell type; leave it blank to keep "
            "disease-specific states apart. DE is unaffected either way."
        )
        self._batch_hint.setWordWrap(True)
        self._batch_hint.setVisible(False)
        pca_page_layout.addWidget(self._batch_hint)

        self.batch_status_label = QLabel("")
        self.batch_status_label.setWordWrap(True)
        self.batch_status_label.setVisible(False)
        pca_page_layout.addWidget(self.batch_status_label)

        pca_elbow_link = QPushButton("View elbow plot")
        pca_elbow_link.setProperty("role", "link_accent")
        pca_elbow_link.setCursor(Qt.CursorShape.PointingHandCursor)
        pca_elbow_link.clicked.connect(
            lambda: self._right_tabs.setCurrentIndex(0))
        pca_page_layout.addWidget(
            pca_elbow_link, alignment=Qt.AlignmentFlag.AlignLeft)
        pca_page_layout.addStretch()

        # ── Clustering drawer page: resolution tools ──
        sweep_label = SectionHeader("RESOLUTION SWEEP")
        cluster_page_layout.addWidget(sweep_label)

        sweep_grid = QGridLayout()
        sweep_grid.setSpacing(4)

        sweep_grid.addWidget(QLabel("Min"), 0, 0)
        self.sweep_min_spin = NoScrollDoubleSpinBox()
        self.sweep_min_spin.setRange(0.05, 4.9)
        self.sweep_min_spin.setValue(0.1)
        self.sweep_min_spin.setSingleStep(0.1)
        sweep_grid.addWidget(self.sweep_min_spin, 0, 1)

        sweep_grid.addWidget(QLabel("Max"), 1, 0)
        self.sweep_max_spin = NoScrollDoubleSpinBox()
        self.sweep_max_spin.setRange(0.1, 5.0)
        self.sweep_max_spin.setValue(0.8)
        self.sweep_max_spin.setSingleStep(0.1)
        sweep_grid.addWidget(self.sweep_max_spin, 1, 1)

        sweep_grid.addWidget(QLabel("Step"), 2, 0)
        self.sweep_step_spin = NoScrollDoubleSpinBox()
        self.sweep_step_spin.setRange(0.01, 1.0)
        self.sweep_step_spin.setValue(0.1)
        self.sweep_step_spin.setSingleStep(0.05)
        self.sweep_step_spin.setDecimals(2)
        sweep_grid.addWidget(self.sweep_step_spin, 2, 1)

        cluster_page_layout.addLayout(sweep_grid)

        self.run_sweep_btn = QPushButton("Run Sweep")
        self.run_sweep_btn.setEnabled(False)
        self.run_sweep_btn.clicked.connect(self._run_resolution_sweep)
        cluster_page_layout.addWidget(self.run_sweep_btn)

        self.run_clustree_btn = QPushButton("Stability Analysis")
        self.run_clustree_btn.setEnabled(False)
        self.run_clustree_btn.setToolTip(
            "Shows how clusters split and merge across resolutions.\n"
            "Blue edges = stable (>50% of cells), Orange = unstable.\n"
            "Pick the resolution where your clusters stabilise."
        )
        self.run_clustree_btn.clicked.connect(self._run_clustree)
        cluster_page_layout.addWidget(self.run_clustree_btn)
        cluster_page_layout.addStretch()

        # ═══ Install the settings pages into their summary cards ═══
        self._pca_card.set_settings_widget(pca_page)
        self._cluster_card.set_settings_widget(cluster_page)
        self._annot_card.set_settings_widget(self.annotate_tab.controls_widget)

        # ---- Right: preview area (tabbed: UMAP | Elbow Plot) ----
        right_layout = self.content_layout

        self._right_tabs = QTabWidget()

        # Tab 1: Elbow Plot (shown first — PCA is the first step)
        elbow_container = QWidget()
        elbow_layout = borderless(QVBoxLayout, elbow_container)

        self._elbow_plot = InteractivePlot(
            title='Elbow Plot',
            left_label='Variance Explained',
            bottom_label='Principal Component',
            unavailable_message='Run PCA to see the elbow plot.',
        )
        elbow_layout.addWidget(self._elbow_plot, 1)

        self._elbow_info_label = QLabel("")
        self._elbow_info_label.setWordWrap(True)
        self._elbow_info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._elbow_info_label.setProperty("role", "success_banner")
        self._elbow_info_label.setVisible(False)
        elbow_layout.addWidget(self._elbow_info_label)

        self._right_tabs.addTab(elbow_container, "Elbow Plot")

        # Tab 2: UMAP + sweep strip
        umap_container = QWidget()
        umap_container_layout = borderless(QVBoxLayout, umap_container)
        self._umap_widget.tsne_computed.connect(self._on_tsne_ready)

        # Visualization toolbar sits above the plot
        umap_container_layout.addWidget(self._umap_widget.build_toolbar())
        umap_container_layout.addWidget(self._umap_widget, 1)

        # Sweep thumbnail strip (hidden until sweep is run)
        self._sweep_strip = _SweepStrip()
        self._sweep_strip.resolution_selected.connect(self._on_sweep_thumb_selected)
        umap_container_layout.addWidget(self._sweep_strip)

        # Stats footer: clusters / cell types / cells
        self._umap_stats_bar = QFrame()
        self._umap_stats_bar.setObjectName("info_panel")
        self._umap_stats_layout = QHBoxLayout(self._umap_stats_bar)
        self._umap_stats_layout.setContentsMargins(16, 8, 16, 8)
        self._umap_stats_layout.setSpacing(16)
        self._umap_stats_bar.setVisible(False)
        umap_container_layout.addWidget(self._umap_stats_bar)

        self._right_tabs.addTab(umap_container, "UMAP")

        # Tab 3: Annotate (table only — controls are in the drawer)
        self._right_tabs.addTab(self.annotate_tab, "Annotate")

        right_layout.addWidget(self._right_tabs, 1)

        # Side panel progress bar — assigned by main.py
        self.progress_bar = None

    # ------------------------------------------------------------------
    # Data flow
    # ------------------------------------------------------------------

    def reset_state(self):
        """Clear all cached state (called on raw data reload)."""
        # Kill running workers
        for attr in ('clustering_worker', 'sweep_worker', 'clustree_worker', '_pca_worker'):
            w = getattr(self, attr, None)
            if w is not None and w.isRunning():
                w.quit()
                w.wait(2000)
        self.clustering_worker = None
        self.sweep_worker = None
        self.clustree_worker = None
        self._pca_worker = None
        self.adata = None
        self.file_path = None
        self.h5ad_path = None
        self._sweep_results = []
        self._suggested_pcs = None
        self._last_adata_version = -1
        # Reset UI
        self._umap_widget.clear_plot()
        self._elbow_plot.clear_plot_items()
        self._elbow_plot.show_unavailable_message()
        self._sweep_strip.clear()
        for card in (self._pca_card, self._cluster_card, self._annot_card):
            card.set_completed(False)
            card.set_summary(["Not run yet"])
        self._accordion.collapse_all()
        self._update_state()
        # Reset annotate tab too
        if hasattr(self, 'annotate_tab'):
            self.annotate_tab.reset_state()

    def set_project_directory(self, directory: str):
        self.project_dir = Path(directory)

    def set_data(self, adata, file_path=None):
        """Receive adata from the workspace / upstream tab (PCA tab compat)."""
        prev_adata = self.adata
        self.adata = adata
        self.file_path = file_path

        # Only repopulate combos if adata object actually changed
        if adata is not prev_adata:
            self._umap_widget.set_adata(adata)
            self._populate_harmony_keys()
            # Forward to embedded annotate tab
            if hasattr(self, 'annotate_tab'):
                self.annotate_tab.set_adata(adata)

        self._update_state()

        # Detect pre-existing analysis and populate plots
        if adata is not prev_adata and adata is not None:
            if 'X_pca' in adata.obsm:
                self._show_elbow_plot(adata)
                if self._suggested_pcs:
                    self.n_pcs_spin.setValue(self._suggested_pcs)

        # Defer UMAP rendering so the tab appears immediately
        self._umap_widget._info_label.setText("Rendering...")
        QTimer.singleShot(50, self._update_preview)

    def _on_tsne_ready(self, adata):
        """Update our adata reference when background t-SNE completes."""
        self.adata = adata

    def _update_state(self):
        """Enable/disable buttons based on data state and show detected status."""
        has_data = self.adata is not None
        self.run_pca_btn.setEnabled(has_data)
        self.run_cluster_btn.setEnabled(has_data)

        has_pca = has_data and 'X_pca' in self.adata.obsm
        has_neighbors = has_data and 'neighbors' in self.adata.uns
        has_umap = has_data and 'X_umap' in self.adata.obsm
        has_leiden = has_data and 'leiden' in self.adata.obs.columns

        self.run_sweep_btn.setEnabled(has_pca and has_neighbors and has_umap)
        self.run_clustree_btn.setEnabled(has_pca and has_neighbors)

        # Show status of pre-existing analysis
        if has_data:
            if has_pca:
                n_pcs = self.adata.obsm['X_pca'].shape[1]
                n_hvg = self.adata.var['highly_variable'].sum() if 'highly_variable' in self.adata.var.columns else '?'
                harmony = 'X_pca_harmony' in self.adata.obsm
                msg = f"PCA present: {n_pcs} PCs, {n_hvg} HVGs"
                if harmony:
                    msg += " + Harmony"
                self.pca_status.setText(msg)
                self.pca_status.setProperty("role", "status_success")
                _refresh_style(self.pca_status)
                # Auto-fill spinbox from existing PCA
                self.pca_compute_spin.setValue(n_pcs)
                self._pca_card.set_completed(True)
                self._pca_card.set_summary([f"{n_pcs} PCs computed"])
            else:
                self.pca_status.setText("")
                self.pca_status.setProperty("role", "secondary")
                _refresh_style(self.pca_status)
                self._pca_card.set_completed(False)
                self._pca_card.set_summary(["Not run yet"])

            # Auto-detect pre-existing clustering and mark step complete
            if has_umap and has_leiden:
                self._right_tabs.setCurrentIndex(1)  # UMAP tab
                self.main_window.mark_step_complete(4)  # Cluster step
                self._cluster_card.set_completed(True)
                n_cl = self.adata.obs['leiden'].nunique()
                self._cluster_card.set_summary([f"{n_cl} Leiden clusters"])
            else:
                self._cluster_card.set_completed(False)
                self._cluster_card.set_summary(["Not run yet"])

            # Auto-detect pre-existing annotation
            if 'cell_type' in self.adata.obs.columns:
                n_types = self.adata.obs['cell_type'].nunique()
                self._annot_card.set_completed(True)
                self._annot_card.set_summary([f"{n_types} cell types"])

    def _update_preview(self):
        """Update the interactive embedding with current data."""
        if self.adata is None:
            self._umap_widget.clear_plot()
            self._umap_stats_bar.setVisible(False)
            return

        obsm_key = self._umap_widget._obsm_key
        if obsm_key not in self.adata.obsm:
            self._umap_widget.clear_plot()
            self._umap_widget._info_label.setText("No embedding — run clustering first")
            self._umap_stats_bar.setVisible(False)
            return

        # Determine what to color by:
        # Prefer cell_type if annotation has been done, otherwise leiden
        selected = self._umap_widget._color_combo.currentText()
        has_cell_type = 'cell_type' in self.adata.obs.columns

        if selected and selected in self.adata.obs.columns:
            self._umap_widget._color_key = selected
            self._umap_widget._recolor()
        elif has_cell_type:
            # Auto-switch to cell_type after annotation
            idx = self._umap_widget._color_combo.findText('cell_type')
            if idx >= 0:
                self._umap_widget._color_combo.setCurrentIndex(idx)
            else:
                coords = self.adata.obsm[obsm_key]
                labels = self.adata.obs['cell_type'].astype(str).values
                n = len(set(labels))
                self._umap_widget.set_data_categorical(coords, labels)
        else:
            # Fallback: color by cluster labels
            coords = self.adata.obsm[obsm_key]
            for col in ['leiden', 'clusters', 'seurat_clusters']:
                if col in self.adata.obs.columns:
                    labels = self.adata.obs[col].astype(str).values
                    n = len(set(labels))
                    self._umap_widget.set_data_categorical(
                        coords, labels,
                        f"{len(coords):,} cells  |  {n} clusters  |  "
                        f"res = {self.resolution_spin.value():.2f}"
                    )
                    break
            else:
                self._umap_widget.clear_plot()
                self._umap_widget._info_label.setText("UMAP available but no cluster labels")
                self._umap_stats_bar.setVisible(False)
                return

        # Update the stats footer below the UMAP
        self._update_umap_info_label()

    def _on_annotation_complete(self, message: str):
        """Handle annotation completing — refresh UMAP to show cell types."""
        # Pull updated adata from annotate tab
        if self.annotate_tab.adata is not None:
            self.adata = self.annotate_tab.adata
            self._umap_widget.set_adata(self.adata)
            self.main_window.current_adata = self.adata
            self.main_window.notify_adata_changed()
            self._last_adata_version = self.main_window._adata_version
        # Switch color-by to cell_type and refresh
        self._update_preview()
        # Update the Annotate summary card
        if self.adata is not None and 'cell_type' in self.adata.obs.columns:
            n_types = self.adata.obs['cell_type'].nunique()
            method = self.annotate_tab.method_combo.currentText()
            self._annot_card.set_completed(True)
            self._annot_card.set_summary([method, f"{n_types} cell types"])
        # Switch to UMAP tab to show annotated results
        self._right_tabs.setCurrentIndex(1)

    def _update_umap_info_label(self):
        """Update the stats footer below the UMAP with cluster/annotation summary."""
        if self.adata is None:
            self._umap_stats_bar.setVisible(False)
            return

        parts = []
        # Cluster count
        if 'leiden' in self.adata.obs.columns:
            n_clusters = self.adata.obs['leiden'].nunique()
            parts.append(f"{n_clusters} Leiden clusters")
        elif 'clusters' in self.adata.obs.columns:
            n_clusters = self.adata.obs['clusters'].nunique()
            parts.append(f"{n_clusters} clusters")

        # Cell type count
        if 'cell_type' in self.adata.obs.columns:
            n_types = self.adata.obs['cell_type'].nunique()
            parts.append(f"{n_types} cell types annotated")

        parts.append(f"{self.adata.n_obs:,} cells")

        # Rebuild the footer: evenly spaced stat cells with dividers
        while self._umap_stats_layout.count():
            item = self._umap_stats_layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        self._umap_stats_layout.addStretch()
        for i, part in enumerate(parts):
            if i:
                divider = QFrame()
                divider.setProperty("role", "stat_divider")
                divider.setFixedWidth(1)
                divider.setFixedHeight(18)
                self._umap_stats_layout.addWidget(divider)
                self._umap_stats_layout.addStretch()
            stat = QLabel(part)
            stat.setProperty("role", "kv_value")
            self._umap_stats_layout.addWidget(stat)
            self._umap_stats_layout.addStretch()
        self._umap_stats_bar.setVisible(bool(parts))

    # ------------------------------------------------------------------
    # Clustering
    # ------------------------------------------------------------------

    def _run_clustering(self):
        if self.adata is None:
            return

        if hasattr(self.main_window, 'current_h5ad_path') and self.main_window.current_h5ad_path:
            self.h5ad_path = self.main_window.current_h5ad_path
        elif hasattr(self.main_window, 'inspect_tab') and self.main_window.inspect_tab.h5ad_path:
            self.h5ad_path = self.main_window.inspect_tab.h5ad_path
        elif self.project_dir:
            self.h5ad_path = processed_h5ad_path(self.project_dir, "clustered.h5ad")
        else:
            dialogs.warning(self, "No Path", "Set a project directory first.")
            return

        # Use the spinner value for PCs (capped by what's available)
        has_pca = 'X_pca' in self.adata.obsm
        max_pcs = self.adata.obsm['X_pca'].shape[1] if has_pca else 50
        n_pcs = min(self.n_pcs_spin.value(), max_pcs)

        params = {
            'n_hvg': self.adata.var['highly_variable'].sum() if 'highly_variable' in self.adata.var.columns else 2000,
            'n_pcs': n_pcs,
            'n_neighbors': self.neighbors_spin.value(),
            'resolution': self.resolution_spin.value(),
            'compute_tsne': False,
            'preserve_embeddings': self.preserve_embed_check.isChecked(),
            'honour_roles': self.honour_roles_check.isChecked(),
            # Harmony is now handled in the PCA tab; use corrected PCA if available
            'harmony': 'X_pca_harmony' in self.adata.obsm,
            'harmony_key': '',
        }

        self._last_clustering_params = params
        self.run_cluster_btn.setEnabled(False)
        self._reset_progress()

        # Pass adata directly; ClusteringWorker copies it on the worker
        # thread so the main thread doesn't freeze for several seconds
        # on big datasets.
        self.clustering_worker = ClusteringWorker(self.adata, params, str(self.h5ad_path))
        run_worker(
            self.clustering_worker,
            on_finished=self._on_clustering_finished,
            on_failed=self._on_clustering_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_progress,
        )

    def _on_clustering_finished(self, payload):
        adata, message = payload
        self._reset_progress()
        self.run_cluster_btn.setEnabled(True)

        self.adata = adata
        self._umap_widget.set_adata(adata)
        self.main_window.set_adata(adata, str(self.h5ad_path))
        if hasattr(self.main_window, 'record_provenance'):
            n_clusters = (int(adata.obs['leiden'].nunique())
                          if 'leiden' in adata.obs.columns else None)
            self.main_window.record_provenance('cluster', {
                'resolution': round(float(self.resolution_spin.value()), 3),
                'n_clusters': n_clusters,
            })
        self._last_adata_version = self.main_window._adata_version
        # Forward to annotate tab (now has leiden clusters)
        if hasattr(self, 'annotate_tab'):
            self.annotate_tab.set_adata(adata)
        self._update_state()
        self._update_preview()
        self._log_resolution_guidance(adata)
        self._cluster_card.set_completed(True)
        self._cluster_card.set_summary([
            f"{self.n_pcs_spin.value()} PCs · k={self.neighbors_spin.value()}",
            f"resolution {self.resolution_spin.value():.2f}",
        ])
        # Switch to UMAP tab to show results
        self._right_tabs.setCurrentIndex(1)  # UMAP is tab index 1
        dialogs.info(self, "Clustering Complete", message)
        self.clustering_complete.emit(str(self.h5ad_path))

    def _on_clustering_failed(self, message: str):
        self._reset_progress()
        self.run_cluster_btn.setEnabled(True)
        dialogs.warning(self, "Clustering Failed", message)

    def _log_resolution_guidance(self, adata):
        if 'leiden' not in adata.obs.columns:
            return
        cluster_counts = adata.obs['leiden'].value_counts()
        n = len(cluster_counts)
        res = self.resolution_spin.value()
        self._on_status(
            f"{n} clusters at resolution {res:.2f}  "
            f"(min={cluster_counts.min():,}, max={cluster_counts.max():,}, "
            f"median={int(np.median(cluster_counts)):,})"
        )
        if n < 10:
            self._on_status(f"Few clusters — consider increasing resolution above {res:.1f}")
        elif n > 30:
            self._on_status(f"Many clusters — consider decreasing resolution below {res:.1f}")

    # ------------------------------------------------------------------
    # Resolution Sweep
    # ------------------------------------------------------------------

    def _run_resolution_sweep(self):
        if self.adata is None:
            return
        if 'X_pca' not in self.adata.obsm:
            dialogs.warning(self, "PCA Needed",
                            "Run PCA first (step 1 in the sidebar).")
            return
        if 'X_umap' not in self.adata.obsm:
            dialogs.warning(self, "Clustering Needed",
                            "Run Clustering once first so a UMAP exists "
                            "to preview resolutions on.")
            return
        if 'neighbors' not in self.adata.uns:
            dialogs.warning(self, "Clustering Needed",
                            "Run Clustering once first to build the "
                            "neighbour graph the sweep reuses.")
            return

        min_r = self.sweep_min_spin.value()
        max_r = self.sweep_max_spin.value()
        step = self.sweep_step_spin.value()
        if min_r >= max_r:
            dialogs.warning(self, "Invalid Range", "Min must be less than max.")
            return

        resolutions = list(np.round(np.arange(min_r, max_r + step * 0.5, step), 3))

        if self.project_dir:
            output_dir = str(processed_data_dir(self.project_dir))
        else:
            import tempfile
            output_dir = tempfile.mkdtemp(prefix="sweep_plots_")

        self.run_sweep_btn.setEnabled(False)
        self._reset_progress()

        self.sweep_worker = ResolutionSweepWorker(self.adata, resolutions, output_dir)
        run_worker(
            self.sweep_worker,
            on_finished=self._on_sweep_finished,
            on_failed=self._on_sweep_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_progress,
        )

    def _on_sweep_finished(self, payload):
        results, message = payload
        self._reset_progress()
        self.run_sweep_btn.setEnabled(True)
        self._on_status(message)

        self._sweep_results = results
        self._sweep_strip.set_results(results)

        # Show first result in main preview
        if results:
            self._on_sweep_thumb_selected(0)

    def _on_sweep_failed(self, message: str):
        self._reset_progress()
        self.run_sweep_btn.setEnabled(True)
        self._on_status(message)
        dialogs.warning(self, "Sweep Failed", message)

    def _on_sweep_thumb_selected(self, index: int):
        """User clicked a sweep thumbnail — update main UMAP and resolution spinner."""
        if index < 0 or index >= len(self._sweep_results):
            return

        r = self._sweep_results[index]
        self._sweep_strip.select(index)

        # Update main preview with this resolution's labels
        if 'labels' in r and 'X_umap' in self.adata.obsm:
            coords = self.adata.obsm['X_umap']
            labels = np.array(r['labels']).astype(str)
            n = len(set(labels))
            self._umap_widget.set_data_categorical(
                coords, labels,
                f"{len(coords):,} cells  |  {n} clusters  |  "
                f"res = {r['resolution']:.2f}"
            )

        # Update resolution spinner
        self.resolution_spin.setValue(r['resolution'])
        self._on_status(
            f"Resolution {r['resolution']:.2f} — {r['n_clusters']} clusters  "
            f"(click Run Clustering to apply)"
        )

        # Brief green-border flash on resolution spinner; theme.py
        # defines QDoubleSpinBox[flash="true"] with the border style.
        self.resolution_spin.setProperty("flash", "true")
        _refresh_style(self.resolution_spin)

        def _clear_flash():
            self.resolution_spin.setProperty("flash", "false")
            _refresh_style(self.resolution_spin)
        QTimer.singleShot(1500, _clear_flash)

    # ------------------------------------------------------------------
    # Clustree
    # ------------------------------------------------------------------

    def _run_clustree(self):
        if self.adata is None:
            return
        if 'neighbors' not in self.adata.uns:
            dialogs.warning(self, "Clustering Needed",
                            "Run Clustering once first to build the "
                            "neighbour graph the analysis reuses.")
            return

        min_r = self.sweep_min_spin.value()
        max_r = self.sweep_max_spin.value()
        step = self.sweep_step_spin.value()
        if min_r >= max_r:
            dialogs.warning(self, "Invalid Range", "Min must be less than max.")
            return

        resolutions = list(np.round(np.arange(min_r, max_r + step * 0.5, step), 3))

        if self.project_dir:
            output_dir = str(processed_data_dir(self.project_dir))
        else:
            import tempfile
            output_dir = tempfile.mkdtemp(prefix="clustree_output_")

        self.run_clustree_btn.setEnabled(False)
        self.run_sweep_btn.setEnabled(False)
        self._reset_progress()

        self.clustree_worker = ClustreWorker(self.adata, resolutions, output_dir)
        run_worker(
            self.clustree_worker,
            on_finished=self._on_clustree_finished,
            on_failed=self._on_clustree_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_progress,
        )

    def _on_clustree_finished(self, payload):
        image_path, message = payload
        self._reset_progress()
        self.run_clustree_btn.setEnabled(True)
        self.run_sweep_btn.setEnabled('neighbors' in (self.adata.uns if self.adata else {}))
        self._on_status(message)
        self._show_clustree_dialog(image_path)

    def _on_clustree_failed(self, message: str):
        self._reset_progress()
        self.run_clustree_btn.setEnabled(True)
        self.run_sweep_btn.setEnabled('neighbors' in (self.adata.uns if self.adata else {}))
        self._on_status(message)
        dialogs.warning(self, "Clustree Failed", message)

    def _show_clustree_dialog(self, image_path: str):
        from PyQt6.QtWidgets import QDialog, QVBoxLayout, QLabel, QDialogButtonBox
        from PyQt6.QtGui import QPixmap
        from PyQt6.QtCore import QSize

        screen = self.screen()
        if screen:
            avail = screen.availableGeometry()
            dialog_w = int(avail.width() * 0.88)
            dialog_h = int(avail.height() * 0.88)
        else:
            dialog_w, dialog_h = 1100, 880

        dialog = QDialog(self)
        dialog.setWindowTitle("Resolution Stability Analysis")
        dialog.resize(dialog_w, dialog_h)

        dlg_layout = QVBoxLayout(dialog)
        dlg_layout.setContentsMargins(8, 8, 8, 8)
        dlg_layout.setSpacing(4)

        img_label = QLabel()
        img_label.setAlignment(Qt.AlignmentFlag.AlignCenter)

        pixmap = QPixmap(image_path)
        if not pixmap.isNull():
            scaled = pixmap.scaled(
                QSize(dialog_w - 16, dialog_h - 60),
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation
            )
            img_label.setPixmap(scaled)
        else:
            img_label.setText(f"Could not load image:\n{image_path}")

        dlg_layout.addWidget(img_label, stretch=1)

        path_label = SecondaryLabel(f"Saved to: {image_path}")
        dlg_layout.addWidget(path_label)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(dialog.close)
        dlg_layout.addWidget(buttons)

        dialog.exec()


    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _on_status(self, text: str):
        self.log_message.emit(text)

    def _on_progress(self, value: int):
        if self.progress_bar:
            self.progress_bar.setValue(value)

    def _reset_progress(self):
        if self.progress_bar:
            self.progress_bar.setValue(0)

    def on_tab_activated(self):
        """Called when tab becomes visible — pull latest adata from workspace."""
        ws = self.main_window
        if ws.current_adata is None:
            return
        version = getattr(ws, '_adata_version', 0)
        if version != getattr(self, '_last_adata_version', -1):
            self._last_adata_version = version
            self.set_data(ws.current_adata, file_path=ws.current_h5ad_path)

    def refresh_theme(self):
        self._umap_widget.refresh_theme()
        self._sweep_strip.refresh_theme()
        # Re-theme elbow plot
        style_pg_plot(self._elbow_plot, title='Elbow Plot',
                      left_label='Variance Explained',
                      bottom_label='Principal Component')

    # ------------------------------------------------------------------
    # PCA: Harmony toggle, populate keys, run, finished, elbow plot
    # ------------------------------------------------------------------

    def _on_harmony_toggled(self, state):
        enabled = Qt.CheckState(state) == Qt.CheckState.Checked
        self.batch_combo.setEnabled(enabled)
        self.condition_combo.setEnabled(enabled)
        self._batch_hint.setVisible(enabled)

    def _populate_harmony_keys(self):
        """Detect and populate batch/condition columns from adata.obs."""
        self.batch_combo.clear()
        if self.adata is None:
            return

        # Batch key combo — prioritise common batch column names
        priority = [
            'study', 'Study',          # a combined master's cohort column
            'sample', 'Sample', 'donor', 'Donor', 'patient', 'Patient',
            'batch', 'Batch', 'orig.ident', 'library_id',
        ]
        added: list[str] = []
        for col in priority:
            if col in self.adata.obs.columns:
                self.batch_combo.addItem(col)
                added.append(col)
        for col in self.adata.obs.columns:
            if col not in added:
                try:
                    n = self.adata.obs[col].nunique()
                    if 2 <= n <= 200:
                        self.batch_combo.addItem(col)
                except Exception:
                    pass
        if self.batch_combo.count() == 0:
            self.batch_combo.addItem('sample')

        # Auto-select detected batch column
        try:
            from kosmic.scrna.inspect.batch import detect_batch_column
            detected_batch = detect_batch_column(self.adata)
            if detected_batch:
                idx = self.batch_combo.findText(detected_batch)
                if idx >= 0:
                    self.batch_combo.setCurrentIndex(idx)
        except Exception:
            pass

        # Merge-across combo. Offered, never pre-selected: correcting on the
        # condition pulls disease and control cells into shared clusters,
        # which hides disease-specific states in the embedding. The user
        # opts in; the default is to correct on the sample only.
        self.condition_combo.clear()
        self.condition_combo.addItem("")  # no second variable
        condition_priority = [
            'condition', 'Condition', 'disease', 'Disease',
            'group', 'Group', 'treatment', 'Treatment',
            'status', 'Status', 'phenotype', 'Phenotype',
        ]
        for col in condition_priority:
            if col in self.adata.obs.columns:
                self.condition_combo.addItem(col)
        for col in self.adata.obs.columns:
            if col not in condition_priority and col not in added:
                try:
                    n = self.adata.obs[col].nunique()
                    if 2 <= n <= 10:
                        self.condition_combo.addItem(col)
                except Exception:
                    pass

        self.condition_combo.setCurrentIndex(0)

        # Harmony on by default when there is a sample column to correct
        # on: every multi-donor study has donor as its dominant technical
        # batch, and the authors of the deposited datasets corrected on it.
        if self.batch_combo.currentText().strip() and not self._harmony_default_applied:
            self.harmony_check.setChecked(True)
            self._harmony_default_applied = True

    def _get_save_path(self) -> str | None:
        """Resolve the h5ad path to save to."""
        if self.file_path:
            return str(self.file_path)
        ws = self.main_window
        if ws and hasattr(ws, 'current_h5ad_path') and ws.current_h5ad_path:
            return ws.current_h5ad_path
        if self.project_dir:
            return str(processed_h5ad_path(self.project_dir, "pca.h5ad"))
        return None

    def _run_pca(self):
        """Launch the PCA pipeline worker."""
        if self.adata is None:
            dialogs.warning(self, "No Data", "Load a dataset first.")
            return

        self.run_pca_btn.setEnabled(False)
        self.pca_status.setText("Running PCA...")
        self._reset_progress()

        # Worker copies adata internally and saves to disk; main thread
        # stays free of the (slow) deep copy and h5ad write.
        save_path = self._get_save_path()
        self._pca_worker = PCAHarmonyWorker(
            self.adata,
            self.hvg_spin.value(),
            self.pca_compute_spin.value(),
            harmony=self.harmony_check.isChecked(),
            scale=self.scale_check.isChecked(),
            honour_roles=self.honour_roles_check.isChecked(),
            harmony_key=self.batch_combo.currentText().strip(),
            harmony_condition_key=self.condition_combo.currentText().strip(),
            output_path=str(save_path) if save_path else None,
        )
        run_worker(
            self._pca_worker,
            on_finished=self._on_pca_finished,
            on_failed=self._on_pca_failed,
            on_progress=self._on_pca_status,
            on_progress_pct=self._on_progress,
        )

    def _on_pca_status(self, text: str):
        self.pca_status.setText(text)
        self.log_message.emit(text)

    def _on_pca_failed(self, message: str):
        self.run_pca_btn.setEnabled(True)
        self._reset_progress()
        self.pca_status.setText(f"PCA failed: {message}")
        dialogs.warning(self, "PCA Failed", message)

    def _on_pca_finished(self, payload):
        adata, message = payload
        self.run_pca_btn.setEnabled(True)
        self._reset_progress()
        self._on_pca_status(message)

        # Store the processed adata. Save was already done inside the
        # worker (output_path passed at launch) so no main-thread h5ad
        # write here.
        self.adata = adata
        self._umap_widget.set_adata(adata)

        # Every one of these changes the embedding, and the 'cluster'
        # stage records only the resolution -- so without this the methods
        # could not say how the space the clusters live in was built.
        if hasattr(self.main_window, 'record_provenance'):
            harmony_on = self.harmony_check.isChecked()
            self.main_window.record_provenance('embedding', {
                'n_hvg': int(self.hvg_spin.value()),
                'hvg_flavor': str(adata.uns.get('hvg', {}).get('flavor', '')),
                'hvg_batch_key': (self.batch_combo.currentText().strip()
                                  if harmony_on else None),
                'scaled': bool(self.scale_check.isChecked()),
                'n_pcs_computed': int(self.pca_compute_spin.value()),
                'harmony': harmony_on,
                'harmony_batch_key': (self.batch_combo.currentText().strip()
                                      if harmony_on else None),
                'excluded_cells_omitted': bool(
                    self.honour_roles_check.isChecked()),
                'harmony_merge_across': (
                    self.condition_combo.currentText().strip() or None
                    if harmony_on else None),
            })

        # Show elbow plot and auto-set PCA spinner
        self._show_elbow_plot(adata)

        # Auto-set the "PCs to use" spinner from elbow suggestion
        if self._suggested_pcs:
            self.n_pcs_spin.setValue(self._suggested_pcs)

        # Switch to Elbow Plot tab to show the result
        self._right_tabs.setCurrentIndex(0)
        self._pca_card.set_completed(True)
        self._pca_card.set_summary([
            f"{self.hvg_spin.value():,} HVGs · {self.pca_compute_spin.value()} PCs"])

        # Update Harmony status label
        if 'X_pca_harmony' in adata.obsm:
            batch_key = self.batch_combo.currentText().strip()
            self.batch_status_label.setText(
                f"Harmony applied on '{batch_key}'. "
                f"Corrected PCA stored in X_pca_harmony."
            )
            self.batch_status_label.setVisible(True)
        else:
            self.batch_status_label.setVisible(False)

        self._update_state()

        # Emit for downstream (analysis_tab wiring, workspace)
        self.pca_complete.emit(adata)

    def _show_elbow_plot(self, adata):
        """Generate and display the interactive elbow plot using pyqtgraph."""
        if 'pca' not in adata.uns or 'variance_ratio' not in adata.uns['pca']:
            self._on_pca_status("PCA variance data not available.")
            return

        from kosmic.visualisation.scrna.pca_plots import elbow_report

        variance_ratio = np.asarray(adata.uns['pca']['variance_ratio'], dtype=float)
        report = elbow_report(variance_ratio)
        suggested = report['suggested']
        cumulative = report['cumulative']
        n_pcs = len(variance_ratio)
        cum_at_suggested = cumulative[min(suggested - 1, n_pcs - 1)] * 100
        verdict = (f"Elbow at PC {report['elbow']}" if report['found']
                   else f"No clear elbow within {n_pcs} PCs; using {suggested}")

        # Store suggestion for auto-setting the PCs spinner
        self._suggested_pcs = suggested
        self._on_pca_status(f"{verdict} ({cum_at_suggested:.0f}% variance)")

        # Clear and draw
        self._elbow_plot.clear_plot_items()
        self._elbow_plot.hide_unavailable_message()
        accent = get_color('accent_primary')
        warning = get_color('warning')

        x = np.arange(1, n_pcs + 1).astype(float)

        # Variance ratio curve
        self._elbow_plot.plot(
            x, variance_ratio,
            pen=pg.mkPen(accent, width=2),
            symbol='o', symbolSize=5, symbolBrush=accent,
        )

        # Elbow vertical line
        elbow_line = pg.InfiniteLine(
            pos=suggested, angle=90,
            pen=pg.mkPen(warning, width=2, style=Qt.PenStyle.DashLine),
        )
        self._elbow_plot.addItem(elbow_line)

        # Elbow label
        elbow_text = pg.TextItem(
            f"Elbow: PC {suggested}" if report['found'] else f"No clear elbow (using {suggested})",
            color=warning, anchor=(0, 1),
        )
        elbow_text.setPos(suggested + 0.5, variance_ratio[0] * 0.9)
        self._elbow_plot.addItem(elbow_text)

        # Auto-range to fit data
        self._elbow_plot.enableAutoRange()

        # Show info below the plot
        total_var = cumulative[-1] * 100 if len(cumulative) > 0 else 0
        self._elbow_info_label.setText(
            f"{verdict}  |  "
            f"Cumulative variance at PC {suggested}: {cum_at_suggested:.1f}%  |  "
            f"Total in {n_pcs} PCs: {total_var:.1f}%"
        )
        self._elbow_info_label.setVisible(True)

        # Mirror the verdict under the Run PCA button, where the user
        # picks "PCs to use" next.
        self.pca_status.setText(f"{verdict} ({cum_at_suggested:.0f}% variance)")
        self.pca_status.setProperty("role", "status_success")
        _refresh_style(self.pca_status)
        self._pca_card.set_summary([
            f"{self.hvg_spin.value():,} HVGs · {self.pca_compute_spin.value()} PCs",
            f"elbow at PC {suggested}",
        ])
