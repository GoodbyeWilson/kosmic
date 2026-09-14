# Download Tab - Download scRNA-seq data from GEO with auto-detection.

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QMessageBox, QFrame, QScrollArea,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QCursor
from pathlib import Path
from typing import TYPE_CHECKING
import shutil

if TYPE_CHECKING:
    pass

from kosmic.gui.shared.widgets import (
    DatasetOverview,
    TabbedPage,
    SecondaryLabel, SectionHeader, CaptionLabel, PrimaryButton,
    SecondaryButton, InfoPanel,
)
from kosmic.paths import processed_data_dir, raw_data_dir
from kosmic.gui.shared import borderless, dialogs, run_worker
from kosmic.gui.intake.workers import (
    BatchConvertWorker, CSVConvertWorker, ExtractAndConvert10xWorker,
    ExtractWorker, H5adImportWorker, RDSConvertWorker, TenXFolderConvertWorker,
    _H5adLoadWorker, _convert_h5, _convert_matrix_txt, _convert_tsv,
)




class DownloadTab(TabbedPage):
    """Tab for downloading GEO datasets with auto-detection."""

    help_id = "scrna/load_data"

    data_status_changed = pyqtSignal(str)  # status_text for sidebar panel
    log_message = pyqtSignal(str)  # routed to OutputPanel

    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.project_dir = None
        self.extract_worker = None
        self._load_worker = None
        self.convert_worker = None
        self.extract_convert_worker = None
        self._batch_total_files = 0
        self._batch_completed_files = 0
        self._setup_ui()

    def _setup_ui(self):
        # Single Dataset page (intake tools moved to the Project-hosted
        # Add Data dialog, ADR-001 Phase E); the tab bar is hidden.
        outer_layout = self.layout()
        self._import_tabs = self.tabs
        self.tabs.tabBar().setVisible(False)

        # --- Tab 1: Dataset ---
        local_scroll = QScrollArea()
        local_scroll.setWidgetResizable(True)
        local_scroll.setFrameShape(QFrame.Shape.NoFrame)
        local_page = QWidget()
        page_layout = QVBoxLayout(local_page)
        page_layout.setContentsMargins(20, 20, 20, 20)
        page_layout.setSpacing(0)

        # Content sits in a width-capped centered column so the dataset
        # overview reads as a focused panel rather than a full-width strip.
        column = QWidget()
        column.setMaximumWidth(DatasetOverview.COLUMN_MAX_WIDTH)
        ll = QVBoxLayout(column)
        ll.setContentsMargins(0, 0, 0, 0)
        ll.setSpacing(0)

        # === State A: Data loaded (shown when h5ad is loaded) ===
        # The dataset screen is the shared scaffold; this page only
        # supplies its stats, actions and panels.
        ov = DatasetOverview(
            "Dataset", "Review your loaded dataset and its basic properties.",
            stats=(("cells", "Cells"), ("genes", "Genes"),
                   ("size", "File size")))
        self._overview = ov
        self._loaded_frame = ov
        self._hero_name = ov.name
        self._hero_stats = ov.stats
        self._ready_chip = ov.ready_chip
        ov.set_subtitle("AnnData  •  h5ad format")
        ov.badge.show()
        ov.set_ready(True)

        # Data management (adding / replacing a study's data) lives on
        # the Project page; this screen only shows what is registered.
        self._reload_btn = ov.add_action(
            SecondaryButton("Restore imported data..."))
        self._reload_btn.setToolTip(
            "Copy the originally imported file back as the working dataset,\n"
            "discarding QC, clustering, and annotation."
        )
        self._reload_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._reload_btn.clicked.connect(self._on_reload_raw)
        self._reload_btn.hide()

        self._overview_panel = ov.add_panel(InfoPanel("DATA OVERVIEW"))
        self._overview_panel.add_row(
            'cell_meta', "Cell metadata",
            "Per-cell annotations (obs): sample, condition, cell type.")
        self._overview_panel.add_row(
            'gene_meta', "Gene metadata",
            "Per-gene annotations (var): IDs, symbols, QC flags.")
        self._overview_panel.add_row(
            'layers', "Layers",
            "Alternative count matrices stored alongside the main one.")
        self._overview_panel.add_row(
            'raw_counts', "Raw counts (.raw)",
            "Original counts snapshot used by DE and pathway scoring.")
        self._overview_panel.add_row(
            'sparse', "Sparse matrix",
            "Sparse storage keeps memory usage low.")
        self._overview_panel.add_row(
            'normalised', "Normalisation",
            "Whether the main matrix holds raw integer counts or already\n"
            "log-normalised values. The QC Normalize step expects raw counts.")

        self._source_panel = ov.add_panel(InfoPanel("DATA SOURCE"))
        self._source_panel.add_row('file_path', "File path")
        self._source_panel.add_row('file_type', "File type")
        self._source_panel.add_row(
            'imported', "Imported",
            "When the originally imported file was added to this study.")
        self._source_panel.add_row(
            'anndata_ver', "anndata version",
            "Version of the anndata library reading this file.")

        self._quality_panel = ov.add_panel(InfoPanel("DATA QUALITY CHECKS"))
        self._quality_panel.add_row('dims', "Matrix dimensions")
        self._quality_panel.add_row(
            'total_counts', "Total counts",
            "Sum of all counts in the matrix.")
        self._quality_panel.add_row(
            'detected', "Detected genes (≥1 count)",
            "Genes with at least one count in at least one cell. A very low\n"
            "fraction can indicate a wrong species or a heavily filtered matrix.")
        self._quality_panel.add_row(
            'median_genes', "Median genes per cell",
            "Below ~200 often indicates empty droplets or low-quality\n"
            "cells; refine at the Quality Control step.")
        self._quality_panel.add_row(
            'mito', "Mitochondrial genes (%)",
            "Share of counts from MT- genes. High values (>20%) suggest\n"
            "stressed or dying cells; filter at the Quality Control step.")

        ov.set_note(
            "ⓘ  These are summary statistics for the imported data. Detailed "
            "exploration and quality control are available in the next steps.")

        # Other h5ads in processed_data/ -- switch targets, shown only when
        # more than one dataset exists.
        self._others_section = QWidget()
        others_layout = borderless(QVBoxLayout, self._others_section)
        others_layout.setSpacing(4)
        others_layout.addSpacing(8)
        others_layout.addWidget(SectionHeader("OTHER DATASETS"))
        self._dataset_list = QVBoxLayout()
        self._dataset_list.setSpacing(2)
        others_layout.addLayout(self._dataset_list)
        self._others_section.hide()
        ov.extra_layout.addWidget(self._others_section)

        self._loaded_frame.hide()
        ll.addWidget(self._loaded_frame)

        # === State B: Convertible data found ===
        self._detected_frame = QFrame()
        self._detected_frame.setProperty("role", "panel")
        det_layout = QVBoxLayout(self._detected_frame)
        det_layout.setContentsMargins(12, 12, 12, 12)
        det_layout.setSpacing(8)

        self._detected_label = QLabel("")
        self._detected_label.setWordWrap(True)
        self._detected_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._detected_label.setProperty("role", "selection_header")
        det_layout.addWidget(self._detected_label)

        det_btn_row = QHBoxLayout()
        det_btn_row.addStretch()
        self._convert_btn = PrimaryButton("  Convert to h5ad")
        self._convert_btn.setMinimumHeight(32)
        self._convert_btn.setFixedWidth(200)
        self._convert_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._convert_btn.clicked.connect(self._auto_convert)
        det_btn_row.addWidget(self._convert_btn)
        det_btn_row.addStretch()
        det_layout.addLayout(det_btn_row)

        self._detected_frame.hide()
        ll.addWidget(self._detected_frame)

        ll.addSpacing(16)

        # === State C: nothing registered for this study. A message, not
        # a form: intake starts from the Project page ("Import
        # External..."), which opens the tools hosted on the sub-tabs
        # above. ===
        self._import_section = QWidget()
        imp_layout = borderless(QVBoxLayout, self._import_section)

        self._import_title = QLabel("")
        self._import_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._import_title.setProperty("role", "page_header")
        imp_layout.addWidget(self._import_title)

        self._import_subtitle = CaptionLabel("")
        self._import_subtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._import_subtitle.setWordWrap(True)
        imp_layout.addWidget(self._import_subtitle)

        imp_layout.addSpacing(4)
        ll.addWidget(self._import_section)

        center_row = QHBoxLayout()
        center_row.addStretch(1)
        center_row.addWidget(column)
        center_row.addStretch(1)
        page_layout.addStretch(2)
        page_layout.addLayout(center_row)
        page_layout.addStretch(3)

        # Forward navigation: the next pipeline step, bottom-right.
        nav_row = QHBoxLayout()
        nav_row.addStretch()
        self._continue_btn = PrimaryButton("Continue to Gene Names  →")
        self._continue_btn.setMinimumHeight(32)
        self._continue_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._continue_btn.setToolTip("Go to the next step: harmonise gene symbols.")
        self._continue_btn.clicked.connect(lambda: self.main_window.switch_tab(1))
        self._continue_btn.hide()
        nav_row.addWidget(self._continue_btn)
        page_layout.addLayout(nav_row)

        local_scroll.setWidget(local_page)
        self._import_tabs.addTab(local_scroll, "Dataset")



        # --- In-tab status label (progress bar is in sidebar) ---
        self._tab_status_label = SecondaryLabel("No data loaded")
        self._tab_status_label.setObjectName("tab_status_strip")
        outer_layout.addWidget(self._tab_status_label)

        self._skip_confirm = False

        # ==================================================================
        # Shared widgets (progress bar / status label wired from sidebar)
        # ==================================================================
        self.progress_bar = None  # wired to sidebar by main.py
        self.status_label = None  # wired to sidebar by main.py via set_status_label

    def set_status_label(self, label) -> None:
        """Wire a foreign QLabel (sidebar progress label) for status updates."""
        self.status_label = label

    def set_project_directory(self, directory: str):
        self.project_dir = Path(directory)
        self._log(f"Project: {directory}")
        self._refresh_local()
        self._run_auto_detect()

    def _detect_convertible_format(self, *, skip_h5ad_check=False):
        """
        Scan raw_data/ and return info about the best convertible format found.

        Returns a dict with keys: format_name, format_key, files, description
        or None if nothing convertible is found. Pass 'skip_h5ad_check=True'
        to skip the "already converted" short-circuit (used after extraction
        when fresh raw files have just appeared next to an old h5ad).
        """
        if not self.project_dir:
            return None

        # 1. Check if h5ad already exists (before requiring raw_data/)
        if not skip_h5ad_check:
            processed_dir = processed_data_dir(self.project_dir)
            if processed_dir.exists():
                h5ad_files = list(processed_dir.glob("*.h5ad"))
                if h5ad_files:
                    return {
                        'format_name': 'h5ad already converted',
                        'format_key': 'h5ad_exists',
                        'files': h5ad_files,
                        'description': f"Found {len(h5ad_files)} existing h5ad file(s) in processed_data/"
                    }

        raw_dir = raw_data_dir(self.project_dir)
        if not raw_dir.exists():
            return None

        # 1b. An h5ad dropped straight into raw_data/. Nothing to convert, but
        # it still needs importing: depositors disagree about which slot holds
        # raw counts (CellxGene normalises X and puts counts in .raw), so the
        # import step picks the counts matrix rather than assuming X.
        raw_h5ads = sorted(raw_dir.glob("*.h5ad"))
        if raw_h5ads:
            return {
                'format_name': 'h5ad (raw_data)',
                'format_key': 'h5ad_raw',
                'files': raw_h5ads,
                'description': f"{len(raw_h5ads)} h5ad file(s) in raw_data/ — "
                               "will select the raw-counts matrix on import"
            }

        # 2. Sample tarballs (e.g., GSM*_filtered_feature_bc_matrix.tar.gz)
        sample_tarballs = self._find_sample_tarballs()
        if sample_tarballs:
            return {
                'format_name': '10X Sample Archives (tar.gz)',
                'format_key': 'sample_tarballs',
                'files': sample_tarballs,
                'description': f"{len(sample_tarballs)} sample tar.gz archive(s)"
            }

        # 3. 10X MTX folders
        mtx_folders = self._find_10x_folders()
        # Filter to only folders that actually have MTX files (not just H5)
        real_mtx_folders = []
        for folder in mtx_folders:
            if list(folder.glob("*.mtx.gz")) or list(folder.glob("*.mtx")) or \
               list(folder.glob("*_matrix.mtx.gz")) or list(folder.glob("*_matrix.mtx")) or \
               list(folder.glob("*.matrix.mtx.gz")) or list(folder.glob("*.matrix.mtx")):
                real_mtx_folders.append(folder)
        if real_mtx_folders:
            return {
                'format_name': '10X Genomics (MTX)',
                'format_key': '10x_mtx',
                'files': real_mtx_folders,
                'description': f"{len(real_mtx_folders)} folder(s) with 10X matrix files"
            }

        # 4. 10X H5 files
        h5_files = [f for f in raw_dir.glob("*.h5") if not f.name.endswith('.h5ad')]
        if h5_files:
            return {
                'format_name': '10X Genomics (H5)',
                'format_key': '10x_h5',
                'files': h5_files,
                'description': f"{len(h5_files)} H5 file(s)"
            }

        # 5. Dense matrix.txt.gz files
        matrix_txt_files = self._find_matrix_txt_files()
        if matrix_txt_files:
            return {
                'format_name': 'Dense Matrix (TXT)',
                'format_key': 'matrix_txt',
                'files': matrix_txt_files,
                'description': f"{len(matrix_txt_files)} matrix.txt file(s)"
            }

        # 6. TSV/CSV expression matrices
        tsv_files = self._find_tsv_matrices()
        if tsv_files:
            return {
                'format_name': 'TSV/CSV Expression Matrix',
                'format_key': 'tsv_csv',
                'files': tsv_files,
                'description': f"{len(tsv_files)} TSV/CSV matrix file(s)"
            }

        # 7. RDS/Robj files
        rds_files = (list(raw_dir.glob("*.rds")) + list(raw_dir.glob("*.RDS")) +
                     list(raw_dir.glob("*.rds.gz")) + list(raw_dir.glob("*.RDS.gz")) +
                     list(raw_dir.glob("*.Robj.gz")) + list(raw_dir.glob("*.robj.gz")) +
                     list(raw_dir.glob("*.Robj")) + list(raw_dir.glob("*.robj")))
        if rds_files:
            return {
                'format_name': 'Seurat RDS/Robj',
                'format_key': 'rds',
                'files': rds_files,
                'description': f"{len(rds_files)} RDS/Robj file(s)"
            }

        # 8. CSV count matrices
        csv_files = []
        for pattern in ['*count*.csv.gz', '*count*.csv', '*matrix*.csv.gz', '*matrix*.csv',
                        '*expression*.csv.gz', '*expression*.csv', '*Count*.csv.gz', '*Count*.csv']:
            csv_files.extend(raw_dir.glob(pattern))
        csv_files = list(set(csv_files))
        if csv_files:
            return {
                'format_name': 'CSV Count Matrix',
                'format_key': 'csv_matrix',
                'files': csv_files,
                'description': f"{len(csv_files)} CSV count matrix file(s)"
            }

        # 9. Unextracted TAR archives
        tar_files = list(raw_dir.glob("*.tar")) + list(raw_dir.glob("*.tar.gz")) + list(raw_dir.glob("*.tgz"))
        if tar_files:
            return {
                'format_name': 'Unextracted TAR Archive',
                'format_key': 'tar_archive',
                'files': tar_files,
                'description': f"{len(tar_files)} TAR archive(s) — will extract first, then convert"
            }

        return None

    def _auto_convert(self):
        """Auto-detect and convert files in raw_data/ without user prompts."""
        detected = self._detect_convertible_format()
        if not detected or detected['format_key'] == 'h5ad_exists':
            return

        fmt_key = detected['format_key']
        self._log(f"Converting: {detected['format_name']}")
        self._detected_frame.hide()
        self._tab_status_label.setText(f"Converting {detected['format_name']}...")

        self._skip_confirm = True
        try:
            if fmt_key == 'tar_archive':
                self._auto_extract_then_convert(detected['files'])
            elif fmt_key == 'sample_tarballs':
                self._extract_and_convert_tarballs()
            elif fmt_key in ('10x_mtx', '10x_h5', 'tsv_csv'):
                self._convert_10x_to_h5ad()
            elif fmt_key == 'matrix_txt':
                self._convert_matrix_txt_to_h5ad()
            elif fmt_key == 'h5ad_raw':
                self._start_h5ad_import(detected['files'])
            elif fmt_key in ('rds', 'csv_matrix'):
                self._convert_to_h5ad()
        finally:
            self._skip_confirm = False

    def _start_h5ad_import(self, h5ad_files):
        """Import an h5ad from raw_data/, selecting its raw-counts matrix."""
        h5ad_file = max(h5ad_files, key=lambda f: f.stat().st_mtime)
        processed_dir = processed_data_dir(self.project_dir)

        self.convert_worker = H5adImportWorker(
            str(h5ad_file), str(processed_dir))
        run_worker(
            self.convert_worker,
            on_finished=self._on_convert_done,
            on_failed=self._on_convert_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_convert_progress,
        )

        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)

        self._tab_status_label.setText(f"Importing {h5ad_file.name}...")
        self._log(f"Importing h5ad: {h5ad_file.name}")

    def _auto_extract_then_convert(self, tar_files):
        """Extract TAR archive(s) and then auto-detect and convert the extracted files."""
        if not tar_files:
            return

        raw_dir = raw_data_dir(self.project_dir)
        tar_file = max(tar_files, key=lambda f: f.stat().st_mtime)

        self._log(f"Auto-extracting: {tar_file.name} before conversion...")

        self.extract_worker = ExtractWorker(str(tar_file), str(raw_dir))
        run_worker(
            self.extract_worker,
            on_finished=self._on_auto_extract_done,
            on_failed=self._on_auto_extract_failed,
            on_progress=self._on_status,
        )

    def _on_auto_extract_failed(self, message: str):
        self._log(message)
        dialogs.warning(self, "Extraction Failed", message)

    def _on_auto_extract_done(self, payload):
        """Callback after auto TAR extraction. Re-detects format and dispatches conversion."""
        files, message = payload
        self._log(message)

        self._refresh_local()
        self._log("Extraction complete. Detecting convertible format...")

        # Re-detect now that files are extracted (skip h5ad check)
        detected = self._detect_convertible_format(skip_h5ad_check=True)
        if detected is None:
            dialogs.warning(self, "No Convertible Data",
                "TAR was extracted but no convertible data files were found inside.")
            return

        self._log(f"Found: {detected['format_name']} — converting...")
        fmt_key = detected['format_key']

        self._skip_confirm = True
        try:
            if fmt_key == 'sample_tarballs':
                self._extract_and_convert_tarballs()
            elif fmt_key in ('10x_mtx', '10x_h5', 'tsv_csv'):
                self._convert_10x_to_h5ad()
            elif fmt_key == 'matrix_txt':
                self._convert_matrix_txt_to_h5ad()
            elif fmt_key == 'h5ad_raw':
                self._start_h5ad_import(detected['files'])
            elif fmt_key in ('rds', 'csv_matrix'):
                self._convert_to_h5ad()
            else:
                self._log(f"No converter for '{detected['format_name']}'")
        finally:
            self._skip_confirm = False

    def _on_extract_only_done(self, payload):
        """Callback after extract-only TAR extraction."""
        files, message = payload
        self._log(message)
        self._refresh_local()
        n_files = len(files) if files else 0
        self._log(
            f"Extracted {n_files} files. Review the samples in raw_data/, "
            "delete any you don't want, then re-import to auto-convert."
        )
        dialogs.info(
            self, "Extraction Complete",
            f"Extracted {n_files} files to raw_data/.\n\n"
            "Review the extracted samples and delete any unwanted ones, "
            "then re-import to auto-convert the remaining samples."
        )

    def _on_extract_only_failed(self, message: str):
        self._log(message)
        dialogs.warning(self, "Extraction Failed", message)

    def _convert_to_h5ad(self):
        raw_dir = raw_data_dir(self.project_dir)

        # Find RDS/Robj files
        rds_files = (list(raw_dir.glob("*.rds")) + list(raw_dir.glob("*.RDS")) +
                     list(raw_dir.glob("*.rds.gz")) + list(raw_dir.glob("*.RDS.gz")) +
                     list(raw_dir.glob("*.Robj.gz")) + list(raw_dir.glob("*.robj.gz")) +
                     list(raw_dir.glob("*.Robj")) + list(raw_dir.glob("*.robj")))

        # Find CSV count matrix files
        csv_files = []
        for pattern in ['*count*.csv.gz', '*count*.csv', '*matrix*.csv.gz', '*matrix*.csv',
                        '*expression*.csv.gz', '*expression*.csv', '*Count*.csv.gz', '*Count*.csv']:
            csv_files.extend(raw_dir.glob(pattern))
        csv_files = list(set(csv_files))  # Remove duplicates

        processed_dir = processed_data_dir(self.project_dir)
        processed_dir.mkdir(exist_ok=True)

        if rds_files and csv_files:
            # Both available - ask user
            msg = QMessageBox(self)
            msg.setWindowTitle("Select Format")
            msg.setText("Found both RDS/Robj and CSV files. Which would you like to convert?")
            rds_btn = msg.addButton("RDS/Robj (needs R)", QMessageBox.ButtonRole.ActionRole)
            csv_btn = msg.addButton("CSV (Python only)", QMessageBox.ButtonRole.ActionRole)
            msg.addButton(QMessageBox.StandardButton.Cancel)
            msg.exec()

            if msg.clickedButton() == rds_btn:
                self._start_rds_conversion(rds_files, processed_dir)
            elif msg.clickedButton() == csv_btn:
                self._start_csv_conversion(csv_files, processed_dir)
        elif rds_files:
            self._start_rds_conversion(rds_files, processed_dir)
        elif csv_files:
            self._start_csv_conversion(csv_files, processed_dir)
        else:
            dialogs.warning(self, "Error", "No RDS/Robj or CSV count matrix files found")

    def _start_rds_conversion(self, rds_files, processed_dir):
        """Start RDS/Robj to h5ad conversion."""
        rds_file = max(rds_files, key=lambda f: f.stat().st_mtime)

        self.convert_worker = RDSConvertWorker(str(rds_file), str(processed_dir))
        run_worker(
            self.convert_worker,
            on_finished=self._on_convert_done,
            on_failed=self._on_convert_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_convert_progress,
        )

        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)

        self._tab_status_label.setText(f"Converting RDS: {rds_file.name}...")
        self._log(f"Converting RDS: {rds_file.name}...")

    def _start_csv_conversion(self, csv_files, processed_dir):
        """Start CSV count matrix to h5ad conversion."""
        csv_file = max(csv_files, key=lambda f: f.stat().st_mtime)

        self.convert_worker = CSVConvertWorker(str(csv_file), str(processed_dir))
        run_worker(
            self.convert_worker,
            on_finished=self._on_convert_done,
            on_failed=self._on_convert_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_convert_progress,
        )

        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)

        self._tab_status_label.setText(f"Converting CSV: {csv_file.name}...")
        self._log(f"Converting CSV: {csv_file.name}...")

    def _on_convert_progress(self, value: int):
        """Handle conversion progress update."""
        if self.progress_bar:
            self.progress_bar.setValue(value)

    def _on_convert_done(self, payload):
        output_path, message = payload
        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
        self._log(message)
        if self.progress_bar:
            self.progress_bar.setValue(100)
        self._tab_status_label.setText(f"Converted: {Path(output_path).name}")
        self._refresh_local()
        self._run_auto_detect()
        self.main_window.mark_step_complete(0)
        self._load_and_set_adata(output_path)

    def _on_convert_failed(self, message: str):
        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
        self._log(message)
        if self.progress_bar:
            self.progress_bar.setValue(0)
        self._tab_status_label.setText(f"Conversion failed: {message[:80]}")
        dialogs.warning(self, "Conversion Failed", message)

    def _ensure_raw_copy(self, h5ad_path: str):
        """
        Ensure a raw (immutable) copy exists in raw_data/.

        If the file is already in raw_data/, just record the path.
        If it's in processed_data/ or elsewhere, copy it to raw_data/.
        Sets main_window.raw_h5ad_path to the raw copy.
        """
        from pathlib import Path

        h5ad_path = Path(h5ad_path)
        if not h5ad_path.exists():
            return

        raw_dir = raw_data_dir(self.project_dir)
        raw_dir.mkdir(exist_ok=True)
        raw_copy = raw_dir / h5ad_path.name

        if h5ad_path.parent.resolve() == raw_dir.resolve():
            # Already in raw_data/
            self.main_window.raw_h5ad_path = str(raw_copy)
        elif not raw_copy.exists():
            # Copy to raw_data/ as the immutable reference
            shutil.copy2(h5ad_path, raw_copy)
            self._log(f"Raw copy saved: raw_data/{h5ad_path.name}")
            self.main_window.raw_h5ad_path = str(raw_copy)
        else:
            # Raw copy already exists
            self.main_window.raw_h5ad_path = str(raw_copy)

    def _load_and_set_adata(self, h5ad_path: str):
        """Load the h5ad file in a background thread so the UI stays responsive."""
        from pathlib import Path

        # If a load is already in progress, skip to avoid destroying the running thread
        if self._load_worker is not None and self._load_worker.isRunning():
            return

        h5ad_path = str(h5ad_path)
        file_name = Path(h5ad_path).name

        # Show loading status in sidebar + in-tab
        if self.status_label:
            self.status_label.setText(f"Loading {file_name}...")
        if self.progress_bar:
            self.progress_bar.setValue(0)
            self.progress_bar.setMaximum(0)  # indeterminate / animated bar
        self._tab_status_label.setText(f"Loading {file_name}...")
        self.data_status_changed.emit(f"Loading {file_name}...")
        # The hero area must say what is happening now, not repeat
        # whatever state it showed last ('No data registered' from the
        # previous study would sit there for the whole load).
        self._loaded_frame.hide()
        self._detected_frame.hide()
        self._import_title.setText(Path(self.project_dir).name
                                   if self.project_dir else file_name)
        self._import_subtitle.setText(f"Loading {file_name}...")
        self._import_section.show()

        # Load in background thread
        self._load_worker = _H5adLoadWorker(h5ad_path)
        run_worker(
            self._load_worker,
            on_finished=self._on_h5ad_loaded,
            on_failed=self._on_h5ad_load_failed,
        )

    def _on_h5ad_loaded(self, payload):
        """Handle h5ad load completion from background thread."""
        from pathlib import Path
        adata, h5ad_path, message = payload
        file_name = Path(h5ad_path).name

        if self.progress_bar:
            self.progress_bar.setMaximum(100)
            self.progress_bar.setValue(100)

        self._log(f"Loaded: {message}")
        if self.status_label:
            self.status_label.setText(f"Loaded {message}")
        self._tab_status_label.setText("Ready")
        self.data_status_changed.emit(f"<b>Loaded:</b> {file_name} ({message})")
        # Dataset overview becomes the hero; data management lives on
        # the Project page.
        self._loaded_frame.show()
        self._detected_frame.hide()
        self._import_section.hide()
        self._ensure_raw_copy(h5ad_path)
        self.main_window.set_adata(adata, h5ad_path)
        self.main_window.mark_step_complete(0)
        self._update_reload_btn_visibility()
        self._refresh_dataset_combo()
        self._update_forward_nav()

    def _on_h5ad_load_failed(self, message: str):
        if self.progress_bar:
            self.progress_bar.setMaximum(100)
            self.progress_bar.setValue(0)
        self._log(f"Load error: {message}")
        if self.status_label:
            self.status_label.setText(f"Load failed: {message}")
        self._tab_status_label.setText(f"Load failed: {message}")
        self.data_status_changed.emit(f"<b>Error:</b> {message}")
        self._import_subtitle.setText(f"Could not load the dataset: {message}")

    def _find_sample_tarballs(self):
        """Find sample-level tar.gz files (e.g., GSM*_filtered_feature_bc_matrix.tar.gz)."""
        raw_dir = raw_data_dir(self.project_dir)
        if not raw_dir.exists():
            return []

        sample_tarballs = []
        for pattern in ['*_filtered_feature_bc_matrix.tar.gz', '*_raw_feature_bc_matrix.tar.gz',
                        'GSM*.tar.gz', '*_matrix.tar.gz']:
            sample_tarballs.extend(raw_dir.glob(pattern))

        # Filter out the main RAW.tar file
        sample_tarballs = [t for t in sample_tarballs if '_RAW.tar' not in t.name and t.name != 'RAW.tar']
        return sample_tarballs

    def _extract_and_convert_tarballs(self):
        """Extract sample tar.gz files and convert to h5ad."""
        if not self.project_dir:
            dialogs.warning(self, "No Project", "Please open a project first.")
            return

        raw_dir = raw_data_dir(self.project_dir)
        if not raw_dir.exists():
            dialogs.warning(self, "No Data", "No raw_data folder found in project.")
            return

        # Check for sample tar.gz files
        sample_tarballs = self._find_sample_tarballs()

        if not sample_tarballs:
            dialogs.warning(self, "No Sample Archives",
                "No sample tar.gz files found in raw_data.\n\n"
                "Expected files like:\n"
                "  - GSM*_filtered_feature_bc_matrix.tar.gz\n"
                "  - GSM*_raw_feature_bc_matrix.tar.gz")
            return

        # Determine output path
        processed_dir = processed_data_dir(self.project_dir)
        processed_dir.mkdir(exist_ok=True)
        output_name = f"{self.project_dir.name}_combined.h5ad"
        output_path = processed_dir / output_name

        # Show confirmation dialog (skip if already confirmed by unified button)
        if not self._skip_confirm:
            sample_list = "\n".join([f"  - {t.name}" for t in sample_tarballs[:10]])
            if len(sample_tarballs) > 10:
                sample_list += f"\n  ... and {len(sample_tarballs) - 10} more"

            if not dialogs.confirm(self, "Extract & Convert Sample Archives", f"Found {len(sample_tarballs)} sample tar.gz file(s):\n\n"
                f"{sample_list}\n\n"
                f"This will:\n"
                f"1. Extract each tar.gz to get 10X matrix files\n"
                f"2. Load and combine all samples\n"
                f"3. Save as: {output_path.name}\n\n"
                f"Continue?"):
                return

        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self._log("Extracting and converting sample archives...")

        self.extract_convert_worker = ExtractAndConvert10xWorker(str(raw_dir), str(output_path))
        run_worker(
            self.extract_convert_worker,
            on_finished=self._on_extract_convert_finished,
            on_failed=self._on_extract_convert_failed,
            on_progress=self._on_status,
            on_progress_pct=lambda v: self.progress_bar.setValue(v),
        )

    def _on_extract_convert_finished(self, payload):
        """Handle extract and convert completion."""
        output_path, message = payload
        self._log(message)
        self.progress_bar.setValue(100)
        self._refresh_local()
        self._run_auto_detect()
        self.main_window.mark_step_complete(0)
        if output_path:
            self._load_and_set_adata(output_path)

    def _on_extract_convert_failed(self, message: str):
        self._log(message)
        self.progress_bar.setValue(0)
        dialogs.warning(self, "Conversion Failed", message)

    def _find_10x_folders(self):
        """Find folders containing 10X Genomics files (matrix.mtx, genes/features, barcodes) or .h5 files."""
        raw_dir = raw_data_dir(self.project_dir)
        if not raw_dir.exists():
            return []

        found_folders = []

        # Check raw_data directly for MTX files or H5 files
        mtx_files = list(raw_dir.glob("*.mtx.gz")) + list(raw_dir.glob("*.mtx"))
        h5_files = [f for f in raw_dir.glob("*.h5") if not f.name.endswith('.h5ad')]
        if mtx_files or h5_files:
            found_folders.append(raw_dir)

        # Check subdirectories (common for GEO extractions like GSM*/filtered_feature_bc_matrix/)
        for subdir in raw_dir.iterdir():
            if subdir.is_dir():
                mtx_in_subdir = list(subdir.glob("*.mtx.gz")) + list(subdir.glob("*.mtx"))
                if mtx_in_subdir:
                    found_folders.append(subdir)
                # Also check one level deeper (e.g., GSMxxxx/filtered_feature_bc_matrix/)
                for subsubdir in subdir.iterdir():
                    if subsubdir.is_dir():
                        mtx_deep = list(subsubdir.glob("*.mtx.gz")) + list(subsubdir.glob("*.mtx"))
                        if mtx_deep:
                            found_folders.append(subsubdir)

        return found_folders

    def _find_tsv_matrices(self):
        """Find TSV/CSV gene expression matrices (common GEO format)."""
        raw_dir = raw_data_dir(self.project_dir)
        if not raw_dir.exists():
            return []

        # Look for TSV/CSV matrix files
        tsv_files = []
        for pattern in ["*matrix*.tsv.gz", "*matrix*.tsv", "*counts*.tsv.gz",
                        "*expression*.tsv.gz", "GSM*.tsv.gz", "GSM*.csv.gz",
                        "*matrix*.csv.gz", "*matrix*.csv", "*counts*.csv.gz", "*counts*.csv",
                        "*expression*.csv.gz", "*expression*.csv", "GSM*.csv",
                        "*.csv.gz", "*.csv"]:
            tsv_files.extend(raw_dir.glob(pattern))

        # Filter out metadata files (not count matrices)
        tsv_files = [f for f in tsv_files if 'metadata' not in f.name.lower()
                     and 'barcodes' not in f.name.lower()
                     and 'features' not in f.name.lower()
                     and 'genes' not in f.name.lower()]

        return list(set(tsv_files))

    def _find_matrix_txt_files(self):
        """Find dense matrix.txt.gz files (genes x cells format)."""
        raw_dir = raw_data_dir(self.project_dir)
        if not raw_dir.exists():
            return []

        # Look for matrix.txt.gz files (dense gene expression matrices)
        matrix_files = []
        for pattern in ["*matrix*.txt.gz", "*matrix*.txt", "*_matrix.txt.gz", "*_matrix.txt",
                        "GSM*_matrix.txt.gz", "GSM*_matrix.txt"]:
            matrix_files.extend(raw_dir.glob(pattern))

        return list(set(matrix_files))

    def _run_batch_convert(self, files, convert_fn, label, *,
                           sort_items=True, worker_attr='_convert_worker'):
        """Run a BatchConvertWorker with the standard progress wiring."""
        self._log(f"Found {len(files)} {label} file(s) to convert")
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        worker = BatchConvertWorker(
            files, convert_fn, label=label, sort_items=sort_items)
        setattr(self, worker_attr, worker)
        run_worker(
            worker,
            on_finished=self._on_multi_file_convert_done,
            on_failed=self._on_multi_file_convert_failed,
            on_progress=self._on_status,
            on_progress_pct=lambda v: self.progress_bar.setValue(v),
        )

    def _convert_matrix_txt_to_h5ad(self):
        """Convert dense matrix.txt.gz files to h5ad (via worker thread)."""
        if not self.project_dir:
            dialogs.warning(self, "Error", "Select a project directory first")
            return

        matrix_files = self._find_matrix_txt_files()
        if not matrix_files:
            dialogs.warning(self, "Error", "No matrix.txt.gz files found in raw_data")
            return

        self._run_batch_convert(
            matrix_files, _convert_matrix_txt, label="matrix",
            sort_items=False,  # filename order is meaningful here
            worker_attr='_matrix_worker')

    def _convert_10x_to_h5ad(self):
        """Convert extracted 10X files or TSV matrices to h5ad format (via worker threads)."""
        if not self.project_dir:
            dialogs.warning(self, "Error", "Select a project directory first")
            return

        folders = self._find_10x_folders()
        tsv_files = self._find_tsv_matrices()

        # Check for standalone H5 files (10X HDF5 format)
        raw_dir = raw_data_dir(self.project_dir)
        h5_files = [f for f in raw_dir.glob("*.h5") if not f.name.endswith('.h5ad')] if raw_dir.exists() else []

        if not folders and not tsv_files and not h5_files:
            dialogs.warning(self, "Error", "No 10X files (matrix.mtx.gz, .h5) or TSV matrices found")
            return

        # If TSV files found and no 10X data, use TSV converter
        if tsv_files and not folders and not h5_files:
            self._convert_tsv_to_h5ad(tsv_files)
            return

        # If standalone H5 files found, convert those directly
        if h5_files:
            self._convert_h5_files(h5_files)
            return

        # 10X folder conversion via worker
        self._log(f"Found {len(folders)} folder(s) with 10X data")
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)

        self._10x_worker = TenXFolderConvertWorker(folders)
        run_worker(
            self._10x_worker,
            on_finished=self._on_multi_file_convert_done,
            on_failed=self._on_multi_file_convert_failed,
            on_progress=self._on_status,
            on_progress_pct=lambda v: self.progress_bar.setValue(v),
        )

    def _convert_h5_files(self, h5_files: list):
        """Convert 10X HDF5 (.h5) files to h5ad (via worker thread)."""
        self._run_batch_convert(
            h5_files, _convert_h5, label="H5", worker_attr='_h5_worker')

    def _convert_tsv_to_h5ad(self, tsv_files: list):
        """Convert TSV/CSV expression matrices to h5ad (via worker thread)."""
        self._run_batch_convert(
            tsv_files, _convert_tsv, label="TSV", worker_attr='_tsv_worker')

    def _refresh_local(self):
        """Refresh local file state (sidebar file tree handles display)."""
        pass

    def _log(self, msg: str):
        self.log_message.emit(msg)

    # === MTX / SCP Converter Methods ===

    def _on_multi_file_convert_failed(self, message: str):
        """Shared failure handler for multi-file converters."""
        self._log(message)
        self.progress_bar.setValue(0)
        dialogs.warning(self, "Conversion Failed", message)

    def _on_multi_file_convert_done(self, payload):
        """Handle worker completion for multi-file converters (TSV, H5, matrix.txt, 10X)."""
        adatas_list, message = payload
        self._log(message)

        if not adatas_list:
            self.progress_bar.setValue(0)
            return

        processed_dir = processed_data_dir(self.project_dir)
        processed_dir.mkdir(exist_ok=True)

        try:
            if len(adatas_list) > 1:
                if dialogs.confirm(self, "Multiple Samples", f"Loaded {len(adatas_list)} samples. Merge into one h5ad file?\n\n"
                    "Yes = Merge all samples (recommended for analysis)\n"
                    "No = Save each separately"):
                    self._log(f"Merging {len(adatas_list)} samples...")
                    import anndata

                    # Make sample keys unique
                    sample_keys = []
                    seen_names = {}
                    for name, _ in adatas_list:
                        if name in seen_names:
                            seen_names[name] += 1
                            unique_name = f"{name}_{seen_names[name]}"
                        else:
                            seen_names[name] = 0
                            unique_name = name
                        sample_keys.append(unique_name)

                    merged = anndata.concat(
                        [a for _, a in adatas_list],
                        label='sample', keys=sample_keys,
                        join='outer', fill_value=0
                    )
                    merged.obs_names_make_unique()

                    output_path = processed_dir / "merged_samples.h5ad"
                    merged.write_h5ad(output_path)
                    self._log(f"Saved merged: {output_path.name} ({merged.n_obs:,} cells x {merged.n_vars:,} genes)")
                    self._load_and_set_adata(str(output_path))
                else:
                    last_path = None
                    for name, adata in adatas_list:
                        safe_name = name.replace("/", "_").replace("\\", "_")
                        output_path = processed_dir / f"{safe_name}.h5ad"
                        adata.write_h5ad(output_path)
                        self._log(f"Saved: {output_path.name}")
                        last_path = output_path
                    if last_path:
                        self._load_and_set_adata(str(last_path))
            else:
                name, adata = adatas_list[0]
                safe_name = name.replace("/", "_").replace("\\", "_")
                output_path = processed_dir / f"{safe_name}.h5ad"
                adata.write_h5ad(output_path)
                self._log(f"Saved: {output_path.name} ({adata.n_obs:,} cells x {adata.n_vars:,} genes)")
                self._load_and_set_adata(str(output_path))

            if self.progress_bar:
                self.progress_bar.setValue(100)
            self._refresh_local()
            self._run_auto_detect()
            if self.status_label:
                self.status_label.setText("Conversion complete!")
            self._tab_status_label.setText("Conversion complete!")
            self.main_window.mark_step_complete(0)

        except Exception as e:
            import traceback
            self._log(f"Merge/save error: {str(e)}")
            dialogs.warning(self, "Error", f"Failed to save: {str(e)}\n\n{traceback.format_exc()}")

    def _run_auto_detect(self):
        """Scan project directory and emit data status to sidebar + in-tab."""
        if not self.project_dir:
            # No study context: point at the Project workspace.
            self.data_status_changed.emit("No study selected")
            self._tab_status_label.setText("No study selected")
            self._loaded_frame.hide()
            self._detected_frame.hide()
            self._import_title.setText("No study selected")
            self._import_subtitle.setText(
                "Open a project and choose a study in the Project "
                "workspace.")
            self._import_section.show()
            self._update_forward_nav()
            return

        detected = self._detect_convertible_format()

        if detected is None:
            accession = Path(self.project_dir).name
            self.data_status_changed.emit(
                "No data registered — add data from the Project page."
            )
            self._tab_status_label.setText("No data registered")
            self._loaded_frame.hide()
            self._detected_frame.hide()
            # Message, not a form: intake starts on the Project page.
            self._import_title.setText(accession)
            self._import_subtitle.setText(
                "No data is registered for this study yet. Use "
                "'Add Data...' on the Project page to import a local "
                "file, download from GEO, or assemble from parts.")
            self._import_section.show()
            self._update_forward_nav()
            return

        fmt_key = detected['format_key']

        if fmt_key == 'h5ad_exists':
            file_names = ", ".join(f.name for f in detected['files'][:3])
            if len(detected['files']) > 3:
                file_names += f" (+{len(detected['files']) - 3} more)"
            self.data_status_changed.emit(f"<b>Ready:</b> {file_names}")
            self._tab_status_label.setText(f"Ready: {file_names}")
            self._detected_frame.hide()
            self.main_window.mark_step_complete(0)
            if self.main_window.current_adata is not None:
                # Loaded state: dataset overview is the hero; data
                # management lives on the Project page.
                self._loaded_frame.show()
                self._import_section.hide()
            else:
                # Load the most recent h5ad; _on_h5ad_loaded finishes the UI
                most_recent = max(detected['files'], key=lambda f: f.stat().st_mtime)
                self._load_and_set_adata(str(most_recent))
        else:
            desc = f"{detected['format_name']} — {detected['description']}"
            self.data_status_changed.emit(f"<b>Convertible:</b> {desc}")
            self._tab_status_label.setText(f"Found: {desc}")
            self._loaded_frame.hide()
            # Show the detected banner with convert button
            file_names = ", ".join(f.name for f in detected['files'][:3])
            if len(detected['files']) > 3:
                file_names += f" (+{len(detected['files']) - 3} more)"
            self._detected_label.setText(
                f"<b>Found {detected['format_name']}</b><br>"
                f"{file_names}"
            )
            self._detected_frame.show()
            self._import_section.hide()
        self._update_forward_nav()

    def on_tab_activated(self):
        self._refresh_local()
        self._run_auto_detect()
        self._update_reload_btn_visibility()
        self._refresh_dataset_combo()
        self._update_forward_nav()

    def _on_reload_raw(self):
        """Restore the originally imported file, discarding all downstream processing."""
        if self.main_window.raw_h5ad_path is None:
            dialogs.warning(self, "No Imported File", "No imported data file has been recorded.")
            return

        if not dialogs.confirm(self, "Restore Imported Data",
            "This will discard all processing (QC, clustering, annotation)\n"
            "and restore the originally imported file.\n\nContinue?"):
            return

        self.main_window.reload_raw_data()
        self._reload_btn.hide()
        self._refresh_dataset_combo()

    def _update_reload_btn_visibility(self):
        """Show restore button when raw != current (i.e. processing has occurred)."""
        raw = self.main_window.raw_h5ad_path
        current = self.main_window.current_h5ad_path
        if raw and current and raw != current:
            self._reload_btn.show()
        else:
            self._reload_btn.hide()

    def _update_forward_nav(self):
        """Continue button appears once a dataset is loaded."""
        self._continue_btn.setVisible(self.main_window.current_adata is not None)

    def _update_hero_card(self):
        """Fill the hero card and the three overview panels from the loaded dataset."""
        path = self.main_window.current_h5ad_path
        adata = self.main_window.current_adata
        if not path or adata is None:
            return
        p = Path(path)
        self._hero_name.setText(p.name)
        try:
            size_mb = p.stat().st_size / (1024 * 1024)
            size_text = (f"{size_mb / 1024:.2f} GB" if size_mb >= 1024
                         else f"{size_mb:.1f} MB")
        except OSError:
            size_text = "?"
        self._hero_stats['cells'].set_value(f"{adata.n_obs:,}")
        self._hero_stats['genes'].set_value(f"{adata.n_vars:,}")
        self._hero_stats['size'].set_value(size_text)

        # --- Data overview panel ---
        import scipy.sparse as sp
        n_obs_cols = adata.obs.shape[1]
        n_var_cols = adata.var.shape[1]
        self._overview_panel.set_value(
            'cell_meta', "Available" if n_obs_cols else "None",
            ok=bool(n_obs_cols),
            tooltip=f"{n_obs_cols} column(s): " + ", ".join(list(adata.obs.columns[:8]))
                    + ("..." if n_obs_cols > 8 else "") if n_obs_cols else None)
        self._overview_panel.set_value(
            'gene_meta', "Available" if n_var_cols else "None",
            ok=bool(n_var_cols),
            tooltip=f"{n_var_cols} column(s): " + ", ".join(list(adata.var.columns[:8]))
                    + ("..." if n_var_cols > 8 else "") if n_var_cols else None)
        # anndata 0.13+ lists X in 'layers' under the key None; count and
        # name only the named layers.
        layer_names = [k for k in adata.layers.keys() if k is not None]
        n_layers = len(layer_names)
        self._overview_panel.set_value(
            'layers', str(n_layers),
            tooltip=", ".join(layer_names) if n_layers else None)
        has_raw = adata.raw is not None
        self._overview_panel.set_value(
            'raw_counts', "Available" if has_raw else "Not stored",
            ok=True if has_raw else None)
        self._overview_panel.set_value(
            'sparse', "Yes" if sp.issparse(adata.X) else "No")

        # Classify every candidate matrix, not just X: depositors disagree
        # about where raw counts live, and using a normalised matrix for
        # pseudobulk DE fails silently.
        import numpy as np

        from kosmic.scrna.load.matrices import (
            LOG_NORMALISED, RAW_COUNTS, describe_count_matrices,
        )
        candidates = describe_count_matrices(adata)
        by_slot = {d['slot']: d for d in candidates}
        x_info = by_slot.get('X', {})
        max_val = x_info.get('max', 0.0)
        others = [d for d in candidates
                  if d['slot'] != 'X' and d['verdict'] == RAW_COUNTS]
        elsewhere = ("\n\nRaw counts are also in: "
                     + ", ".join(d['label'] for d in others)) if others else ""

        if x_info.get('verdict') == LOG_NORMALISED:
            self._overview_panel.set_value(
                'normalised', "Log-normalised",
                ok=True if others else False,
                tooltip=f"X looks log-transformed (max {max_val:.1f}).\n"
                        "Skip the Normalize step in QC -- normalising twice "
                        "corrupts the data." + elsewhere
                        + ("\nRe-import from raw_data/ to use them."
                           if others else
                           "\nNo raw counts found: pseudobulk DE needs counts."))
        elif x_info.get('verdict') == RAW_COUNTS:
            self._overview_panel.set_value(
                'normalised', "Raw counts", ok=True,
                tooltip=f"Integer values (max {max_val:,.0f}). Normalise at "
                        "the QC step before clustering." + elsewhere)
        else:
            self._overview_panel.set_value(
                'normalised', "Unclear", ok=False,
                tooltip=f"Non-integer values with a large maximum ({max_val:,.1f}) "
                        "-- possibly\nscaled or corrected counts. Check the data source.")

        # --- Data source panel ---
        path_str = str(p)
        short_path = path_str if len(path_str) <= 42 else "..." + path_str[-39:]
        self._source_panel.set_value('file_path', short_path, tooltip=path_str)
        self._source_panel.set_value('file_type', "AnnData (.h5ad)")
        source_file = Path(self.main_window.raw_h5ad_path or path)
        try:
            from datetime import datetime
            mtime = datetime.fromtimestamp(source_file.stat().st_mtime)
            self._source_panel.set_value(
                'imported', mtime.strftime("%d %b %Y, %H:%M"))
        except OSError:
            self._source_panel.set_value('imported', "?")
        import anndata as _ad
        self._source_panel.set_value('anndata_ver', _ad.__version__)

        # --- Quality checks panel (computed on the counts source) ---
        mat_source = adata.raw if has_raw else adata
        X = mat_source.X
        n_genes = X.shape[1]
        self._quality_panel.set_value(
            'dims', f"{adata.n_obs:,} × {adata.n_vars:,}", ok=True)
        try:
            total = float(X.sum())
        except (TypeError, ValueError, MemoryError):
            total = float('nan')
        self._quality_panel.set_value(
            'total_counts', f"{total:,.0f}", ok=total > 0)
        if sp.issparse(X):
            detected = int((X.getnnz(axis=0) > 0).sum())
            genes_per_cell = np.asarray(X.getnnz(axis=1)).ravel()
        else:
            X_arr = np.asarray(X)
            detected = int((X_arr > 0).any(axis=0).sum())
            genes_per_cell = (X_arr > 0).sum(axis=1)
        frac = detected / max(n_genes, 1)
        self._quality_panel.set_value(
            'detected', f"{detected:,} ({frac * 100:.1f}%)", ok=frac >= 0.5)
        median_genes = int(np.median(genes_per_cell)) if len(genes_per_cell) else 0
        self._quality_panel.set_value(
            'median_genes', f"{median_genes:,}", ok=median_genes >= 200)
        mito_mask = mat_source.var_names.str.upper().str.startswith("MT-")
        if mito_mask.any() and total > 0:
            mito_idx = np.where(np.asarray(mito_mask))[0]
            mito_total = float(X[:, mito_idx].sum())
            mito_pct = mito_total / total * 100
            self._quality_panel.set_value(
                'mito', f"{mito_pct:.1f}%", ok=mito_pct <= 20)
        else:
            self._quality_panel.set_value('mito', "No MT- genes")

    def _refresh_dataset_combo(self):
        """Refresh the hero card and list the other h5ads in processed_data/ as switch targets."""
        self._update_hero_card()

        # Clear existing cards
        while self._dataset_list.count():
            item = self._dataset_list.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        self._others_section.hide()

        if not self.project_dir:
            return

        processed_dir = processed_data_dir(self.project_dir)
        if not processed_dir.exists():
            return

        current = self.main_window.current_h5ad_path or ""
        others = [f for f in sorted(processed_dir.glob("*.h5ad"),
                                    key=lambda p: p.stat().st_mtime, reverse=True)
                  if str(f) != current]
        if not others:
            return

        for f in others:
            size_mb = f.stat().st_size / (1024 * 1024)

            # Try to read shape from h5ad metadata without loading full file
            info_text = f"{size_mb:.1f} MB"
            try:
                import h5py
                with h5py.File(str(f), 'r') as h5:
                    if 'X' in h5:
                        shape = h5['X'].attrs.get('shape', None)
                        if shape is None and hasattr(h5['X'], 'shape'):
                            shape = h5['X'].shape
                        if shape is not None and len(shape) == 2:
                            info_text = f"{int(shape[0]):,} cells x {int(shape[1]):,} genes  ({size_mb:.1f} MB)"
            except (OSError, KeyError, ValueError):
                pass

            card = QPushButton()
            card.setProperty("role", "dataset_card")
            card.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
            card.setToolTip(str(f))
            card.setText(f"  {f.name}\n  {info_text}")

            path_str = str(f)
            card.clicked.connect(lambda checked, p=path_str: self._on_dataset_card_clicked(p))
            self._dataset_list.addWidget(card)

        self._others_section.show()

    def _on_dataset_card_clicked(self, path: str):
        """User clicked a dataset card to switch to it."""
        current = self.main_window.current_h5ad_path
        if path == current:
            return

        if not dialogs.confirm(self, "Switch Dataset", f"Load {Path(path).name}?\n\n"
            "This will replace the current dataset. "
            "Downstream steps (QC, clustering, etc.) will be reset."):
            return

        # Reset downstream state then load
        self.main_window.raw_h5ad_path = path
        for tab in [self.main_window.gene_names_tab, self.main_window.inspect_tab,
                    self.main_window.qc_tab, self.main_window.cluster_tab,
                    self.main_window.filter_tab]:
            tab.reset_state()
        # Reset Gene Names..Subset; keep Set Directory (0) and Load Data (1).
        for i in range(2, len(self.main_window.WORKFLOW_STEPS)):
            self.main_window._step_completed[i] = False
        self.main_window.steps_reset.emit()

        self._load_and_set_adata(path)

