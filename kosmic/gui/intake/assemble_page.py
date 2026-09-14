# Assemble-from-parts intake tool: combine matrix + metadata +
# coordinate files into one h5ad in the active study's processed_data/.
# Hosted by the Project workspace's Add Data dialog (ADR-001 Phase E).
from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QFileDialog, QFrame, QHBoxLayout, QLabel, QLineEdit, QProgressBar,
    QPushButton, QScrollArea, QVBoxLayout, QWidget,
)

from kosmic.gui.intake.workers import MTXConvertWorker
from kosmic.gui.shared import dialogs, run_worker
from kosmic.gui.shared.widgets import (
    CaptionLabel, HeaderLabel, IconTile, InfoPanel, PrimaryButton,
    SecondaryButton,
)
from kosmic.paths import processed_data_dir


class _DropZone(QFrame):
    """Dashed drop target; hands a dropped local folder to a callback."""

    def __init__(self, on_folder):
        super().__init__()
        self.setProperty("role", "drop_zone")
        self.setAcceptDrops(True)
        self._on_folder = on_folder

    def dragEnterEvent(self, event):
        md = event.mimeData()
        if md.hasUrls() and any(u.isLocalFile() for u in md.urls()):
            event.acceptProposedAction()

    def dropEvent(self, event):
        for u in event.mimeData().urls():
            if u.isLocalFile():
                p = Path(u.toLocalFile())
                self._on_folder(p if p.is_dir() else p.parent)
                event.acceptProposedAction()
                return


class AssemblePage(QWidget):
    """Assemble-from-parts intake tool, decoupled from the scRNA workspace.

    Context: 'study_dir_provider' returns the target study Path (or
    None), 'log_cb' receives log lines, 'on_h5ad_ready(path)' fires when
    the assembled h5ad has been written.
    """

    def __init__(self, study_dir_provider, log_cb, on_h5ad_ready,
                 parent=None):
        super().__init__(parent)
        self._study_dir_provider = study_dir_provider
        self._log_cb = log_cb
        self._on_h5ad_ready = on_h5ad_ready
        self.convert_worker = None

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(6)
        self._build_ui(outer)

        footer = QHBoxLayout()
        self.status_label = CaptionLabel("")
        footer.addWidget(self.status_label, 1)
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setFixedWidth(220)
        footer.addWidget(self.progress_bar)
        outer.addLayout(footer)
        self._tab_status_label = self.status_label

    @property
    def project_dir(self):
        return self._study_dir_provider()

    def _log(self, msg: str):
        self._log_cb(msg)

    def _on_status(self, status: str):
        self.status_label.setText(status)

    def _build_ui(self, outer):
        # --- Tab 3: Assemble from Parts ---
        assemble_scroll = QScrollArea()
        assemble_scroll.setWidgetResizable(True)
        assemble_scroll.setFrameShape(QFrame.Shape.NoFrame)
        assemble_page = QWidget()
        sl = QVBoxLayout(assemble_page)
        sl.setContentsMargins(20, 16, 20, 16)
        sl.setSpacing(12)

        asm_title = QLabel("Assemble from Parts")
        asm_title.setProperty("role", "page_header")
        sl.addWidget(asm_title)
        sl.addWidget(CaptionLabel(
            "Combine separate matrix, metadata, and embedding files into a "
            "single dataset for import."))

        asm_row = QHBoxLayout()
        asm_row.setSpacing(12)

        # Left card: the three file slots
        select_card = QFrame()
        select_card.setObjectName("info_panel")
        sc_layout = QVBoxLayout(select_card)
        sc_layout.setContentsMargins(16, 14, 16, 16)
        sc_layout.setSpacing(10)
        sc_layout.addWidget(HeaderLabel("Select files"))

        self.mtx_path_edit = QLineEdit()
        self.mtx_path_edit.setPlaceholderText("matrix.mtx, .mtx.gz, or .h5")
        self.mtx_browse_btn = QPushButton("Browse...")
        self.mtx_browse_btn.clicked.connect(self._browse_mtx)
        sc_layout.addWidget(self._make_part_row(
            "dataset", "Expression matrix", "Required",
            self.mtx_path_edit, self.mtx_browse_btn, required=True))

        self.meta_path_edit = QLineEdit()
        self.meta_path_edit.setPlaceholderText(
            "cell annotations, cluster labels, sample information")
        self.meta_browse_btn = QPushButton("Browse...")
        self.meta_browse_btn.clicked.connect(self._browse_meta)
        sc_layout.addWidget(self._make_part_row(
            "file", "Metadata", "Optional",
            self.meta_path_edit, self.meta_browse_btn))

        self.clusters_path_edit = QLineEdit()
        self.clusters_path_edit.setPlaceholderText("UMAP / t-SNE coordinates")
        self.clusters_browse_btn = QPushButton("Browse...")
        self.clusters_browse_btn.clicked.connect(self._browse_clusters)
        sc_layout.addWidget(self._make_part_row(
            "umap", "Coordinates", "Optional",
            self.clusters_path_edit, self.clusters_browse_btn))

        strip_row = QHBoxLayout()
        strip_row.setSpacing(10)
        formats_note = QLabel(
            "ⓘ  Supports Matrix Market (.mtx + genes/features + barcodes), "
            "10x output folders, .h5, and related metadata files.")
        formats_note.setProperty("role", "info_banner")
        formats_note.setWordWrap(True)
        strip_row.addWidget(formats_note, stretch=1)
        self.auto_detect_btn = SecondaryButton("Auto-detect files in folder")
        self.auto_detect_btn.clicked.connect(self._auto_detect_files)
        strip_row.addWidget(self.auto_detect_btn)
        self.mtx_convert_btn = PrimaryButton("Assemble & Import")
        self.mtx_convert_btn.clicked.connect(self._convert_mtx)
        strip_row.addWidget(self.mtx_convert_btn)
        sc_layout.addLayout(strip_row)

        asm_row.addWidget(select_card, stretch=2)

        # Right panel: live preview of what will be assembled
        self._assembly_panel = InfoPanel("ASSEMBLY PREVIEW")
        self._assembly_panel.add_row(
            'matrix', "Matrix", "The expression matrix. Required.")
        self._assembly_panel.add_row(
            'metadata', "Metadata",
            "Per-cell annotations merged into obs. Optional.")
        self._assembly_panel.add_row(
            'coords', "Coordinates",
            "Precomputed UMAP / t-SNE coordinates. Optional.")
        self._assembly_panel.add_row(
            'output', "Output", "Everything is assembled into one h5ad file.")
        asm_row.addWidget(self._assembly_panel, stretch=1)

        sl.addLayout(asm_row)

        for edit in (self.mtx_path_edit, self.meta_path_edit,
                     self.clusters_path_edit):
            edit.textChanged.connect(self._update_assembly_preview)

        # Drag & drop target: drop a folder to auto-detect its parts
        drop_zone = _DropZone(self._detect_in_folder)
        dz_row = QHBoxLayout(drop_zone)
        dz_row.setContentsMargins(16, 14, 16, 14)
        dz_row.setSpacing(16)
        dz_row.addWidget(IconTile("folder-open", size=48, icon_size=22))
        dz_col = QVBoxLayout()
        dz_col.setSpacing(2)
        dz_title = QLabel("Prefer drag & drop?")
        dz_title.setProperty("role", "file_row_name")
        dz_col.addWidget(dz_title)
        dz_col.addWidget(CaptionLabel(
            "Drop a folder containing 10x output or related files to "
            "auto-detect and populate the fields above."))
        dz_row.addLayout(dz_col)
        dz_row.addStretch()
        sl.addWidget(drop_zone)

        sl.addStretch()
        self._update_assembly_preview()
        assemble_scroll.setWidget(assemble_page)

        outer.addWidget(assemble_scroll, 1)

    def _make_part_row(self, icon, label_text, tag_text, edit, browse_btn,
                       *, required=False):
        """One file slot on the Assemble page: icon tile, label + tag, path field, browse."""
        row = QFrame()
        row.setProperty("role", "file_row")
        if required:
            row.setProperty("selected", "true")  # accent border marks the required slot
        h = QHBoxLayout(row)
        h.setContentsMargins(10, 8, 10, 8)
        h.setSpacing(12)
        h.addWidget(IconTile(icon, size=40, icon_size=18))
        label_col = QVBoxLayout()
        label_col.setSpacing(2)
        name = QLabel(label_text)
        name.setProperty("role", "file_row_name")
        name.setMinimumWidth(140)
        label_col.addWidget(name)
        tag = QLabel(tag_text)
        tag.setProperty("role", "tag")
        label_col.addWidget(tag, alignment=Qt.AlignmentFlag.AlignLeft)
        h.addLayout(label_col)
        h.addWidget(edit, stretch=1)
        h.addWidget(browse_btn)
        return row

    def _update_assembly_preview(self):
        """Mirror the three path fields into the Assembly preview panel."""
        def fill(key, text, required=False):
            if text:
                self._assembly_panel.set_value(
                    key, Path(text).name, ok=True, tooltip=text)
            elif required:
                self._assembly_panel.set_value(key, "not selected", ok=False)
            else:
                self._assembly_panel.set_value(key, "optional")

        fill('matrix', self.mtx_path_edit.text().strip(), required=True)
        fill('metadata', self.meta_path_edit.text().strip())
        fill('coords', self.clusters_path_edit.text().strip())
        self._assembly_panel.set_value('output', "AnnData (.h5ad)")

    def _auto_detect_files(self):
        """Pick a folder, then auto-detect MTX, metadata, and cluster files in it."""
        folder = QFileDialog.getExistingDirectory(
            self, "Select Folder with Downloaded Files",
            str(self.project_dir) if self.project_dir else ""
        )
        if not folder:
            return
        self._detect_in_folder(Path(folder))

    def _detect_in_folder(self, folder: Path):
        """Auto-detect MTX, metadata, and cluster files in 'folder'."""
        # Look for expression matrix
        for pattern in ['matrix.mtx*', '*.mtx', '*.mtx.gz', '*expression*.mtx*',
                        'filtered_feature_bc_matrix.h5', '*.h5']:
            matches = list(folder.glob(pattern))
            if matches:
                self.mtx_path_edit.setText(str(matches[0]))
                self._log(f"Found matrix: {matches[0].name}")
                break

        # Look for metadata
        for pattern in ['*metadata*.tsv', '*metadata*.csv', '*meta*.txt',
                        'scp_meta.txt', '*annotation*.tsv']:
            matches = list(folder.glob(pattern))
            if matches:
                self.meta_path_edit.setText(str(matches[0]))
                self._log(f"Found metadata: {matches[0].name}")
                break

        # Look for cluster/UMAP coordinates
        for pattern in ['*cluster*.tsv', '*cluster*.csv', '*umap*.tsv',
                        '*tsne*.tsv', '*coordinates*.tsv']:
            matches = list(folder.glob(pattern))
            if matches:
                self.clusters_path_edit.setText(str(matches[0]))
                self._log(f"Found coordinates: {matches[0].name}")
                break

    def _browse_mtx(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Expression Matrix",
            str(self.project_dir) if self.project_dir else "",
            "Matrix Files (*.mtx *.mtx.gz *.h5);;All Files (*)"
        )
        if path:
            self.mtx_path_edit.setText(path)

    def _browse_meta(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Metadata File",
            str(self.project_dir) if self.project_dir else "",
            "Metadata Files (*.tsv *.csv *.txt);;All Files (*)"
        )
        if path:
            self.meta_path_edit.setText(path)

    def _browse_clusters(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Cluster/Coordinates File",
            str(self.project_dir) if self.project_dir else "",
            "Cluster Files (*.tsv *.csv *.txt);;All Files (*)"
        )
        if path:
            self.clusters_path_edit.setText(path)

    def _convert_mtx(self):
        """Convert MTX + metadata to h5ad."""
        mtx_path = self.mtx_path_edit.text().strip()
        if not mtx_path:
            dialogs.warning(self, "Error", "Select an expression matrix file")
            return

        if not Path(mtx_path).exists():
            dialogs.warning(self, "Error", f"Matrix file not found: {mtx_path}")
            return

        # Determine output path
        if self.project_dir:
            output_dir = processed_data_dir(self.project_dir)
        else:
            output_dir = Path(mtx_path).parent

        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate output filename
        mtx_name = Path(mtx_path).stem.replace('.mtx', '').replace('.gz', '')
        output_path = output_dir / f"{mtx_name}.h5ad"

        # Get optional files
        meta_path = self.meta_path_edit.text().strip() or None
        clusters_path = self.clusters_path_edit.text().strip() or None

        self._log("Converting to h5ad...")
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.mtx_convert_btn.setEnabled(False)

        self.mtx_worker = MTXConvertWorker(
            mtx_path=mtx_path,
            metadata_path=meta_path,
            output_path=str(output_path),
            clusters_path=clusters_path
        )
        run_worker(
            self.mtx_worker,
            on_finished=self._on_mtx_convert_done,
            on_failed=self._on_mtx_convert_failed,
            on_progress=self._on_status,
            on_progress_pct=lambda v: self.progress_bar.setValue(v),
        )

    def _on_mtx_convert_done(self, payload):
        output_path, message = payload
        self.mtx_convert_btn.setEnabled(True)
        self._log(message)
        self.progress_bar.setValue(100)
        self.status_label.setText(message)
        if output_path:
            self._on_h5ad_ready(str(output_path))

    def _on_mtx_convert_failed(self, message: str):
        self.mtx_convert_btn.setEnabled(True)
        self._log(message)
        self.progress_bar.setValue(0)
        dialogs.warning(self, "Conversion Failed", message)

