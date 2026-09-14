# GEO import page: query a GEO accession, pick files, download into the
# active study's raw_data/ (single, batch 10X, matrix.txt, metadata).
# Hosted by the Project workspace's Add Data dialog (ADR-001 Phase E).
# Downloads end by notifying the host; format detection and conversion
# happen on the study's Dataset screen.
import re

from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QCursor
from PyQt6.QtWidgets import (
    QButtonGroup, QCheckBox, QFrame, QHBoxLayout, QLabel, QLineEdit,
    QProgressBar, QPushButton, QRadioButton, QScrollArea, QStackedWidget,
    QVBoxLayout, QWidget,
)

from kosmic.gui.intake.workers import DownloadWorker, GEOQueryWorker
from kosmic.gui.shared import dialogs, run_worker
from kosmic.gui.shared.widgets import (
    CaptionLabel, HeaderLabel, IconTile, PrimaryButton, SecondaryButton,
    SecondaryLabel, SettingsGroup, StatusLabel,
)
from kosmic.paths import processed_data_dir, raw_data_dir


class _FileRow(QFrame):
    """Selectable card row that forwards clicks to its radio button."""

    def __init__(self, radio):
        super().__init__()
        self.setProperty("role", "file_row")
        self.setProperty("selected", "false")
        self.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._radio = radio

    def set_selected(self, selected: bool) -> None:
        self.setProperty("selected", "true" if selected else "false")
        style = self.style()
        if style is not None:
            style.unpolish(self)
            style.polish(self)

    def mousePressEvent(self, event):
        self._radio.setChecked(True)
        super().mousePressEvent(event)


class GeoImportPage(QWidget):
    """Download-from-GEO intake tool, decoupled from the scRNA workspace.

    Context: 'study_dir_provider' returns the target study Path (or
    None), 'log_cb' receives log lines, 'on_raw_data_changed' fires
    after files land in raw_data/ so the host can refresh.
    """

    def __init__(self, study_dir_provider, log_cb, on_raw_data_changed,
                 parent=None):
        super().__init__(parent)
        self._study_dir_provider = study_dir_provider
        self._log_cb = log_cb
        self._on_raw_changed = on_raw_data_changed

        self.geo_files = []
        self.data_files = []
        self.metadata_files = []
        self.metadata_checkboxes = []
        self.selected_file = None
        self.query_worker = None
        self.download_worker = None
        self._current_meta_worker = None
        self._current_10x_worker = None
        self._current_matrix_worker = None
        self.detected_10x_samples = {}
        self.detected_matrix_files = []
        self._metadata_download_queue = []
        self._batch_total_files = 0
        self._batch_completed_files = 0

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(6)
        self._build_ui(outer)

        footer = QHBoxLayout()
        self.status_label = SecondaryLabel("")
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
        """The target study folder (copied code reads this name)."""
        return self._study_dir_provider()

    def _log(self, msg: str):
        self._log_cb(msg)

    def _refresh_local(self):
        pass  # host owns file views

    def _run_auto_detect(self):
        # In the intake dialog, "re-detect" means: tell the host new raw
        # files exist; conversion happens on the Dataset screen.
        self._on_raw_changed()

    def _auto_convert(self):
        # Conversion happens on the study's Dataset screen; hand back.
        self._on_raw_changed()

    def _build_ui(self, outer):
        # --- Tab 2: Download from GEO ---
        geo_scroll = QScrollArea()
        geo_scroll.setWidgetResizable(True)
        geo_scroll.setFrameShape(QFrame.Shape.NoFrame)
        geo_page = QWidget()
        gl = QVBoxLayout(geo_page)
        gl.setContentsMargins(20, 16, 20, 16)
        gl.setSpacing(12)

        geo_input_layout = QHBoxLayout()
        geo_input_layout.setSpacing(10)
        geo_input_layout.addWidget(QLabel("GEO accession"))

        self.gse_input = QLineEdit()
        self.gse_input.setPlaceholderText("GSE292067")
        self.gse_input.setFixedWidth(200)
        self.gse_input.setToolTip("Enter a GEO Series accession number (e.g., GSE174574).")
        self.gse_input.returnPressed.connect(self._query_geo)
        geo_input_layout.addWidget(self.gse_input)

        self.query_btn = SecondaryButton("Query GEO")
        self.query_btn.setToolTip("Query the GEO FTP server for available supplementary files.")
        self.query_btn.clicked.connect(self._query_geo)
        geo_input_layout.addWidget(self.query_btn)

        geo_input_layout.addStretch()
        gl.addLayout(geo_input_layout)

        # Body below the search row: one stacked area that is either the
        # pre-query landing state or the fully-populated results page.
        self._geo_body = QStackedWidget()

        # Page 0: pre-query landing, centered in the available space
        geo_landing = QWidget()
        ph = QVBoxLayout(geo_landing)
        ph.setSpacing(6)
        ph.addStretch(2)
        tile_row = QHBoxLayout()
        tile_row.addStretch()
        tile_row.addWidget(IconTile(caption="GEO", size=56))
        tile_row.addStretch()
        ph.addLayout(tile_row)
        ph.addSpacing(8)
        ph_title = QLabel("Search a GEO accession")
        ph_title.setProperty("role", "header")
        ph_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        ph.addWidget(ph_title)
        ph_caption = CaptionLabel(
            "Enter a GSE accession above and press Query GEO to list the "
            "study's files.")
        ph_caption.setAlignment(Qt.AlignmentFlag.AlignCenter)
        ph.addWidget(ph_caption)
        ph.addStretch(3)

        # Page 1: results, populated before the page is shown
        self._geo_results = QWidget()
        rl = QVBoxLayout(self._geo_results)
        rl.setContentsMargins(0, 0, 0, 0)
        rl.setSpacing(12)

        # Query result hero: accession + file count
        self._geo_hero = QFrame()
        self._geo_hero.setObjectName("info_panel")
        geo_hero_row = QHBoxLayout(self._geo_hero)
        geo_hero_row.setContentsMargins(16, 12, 16, 12)
        geo_hero_row.setSpacing(16)
        geo_hero_row.addWidget(IconTile(caption="GEO", size=56))
        geo_hero_col = QVBoxLayout()
        geo_hero_col.setSpacing(2)
        self._geo_hero_title = QLabel("")
        self._geo_hero_title.setObjectName("hero_dataset_name")
        geo_hero_col.addWidget(self._geo_hero_title)
        self._geo_hero_sub = CaptionLabel("")
        geo_hero_col.addWidget(self._geo_hero_sub)
        geo_hero_row.addLayout(geo_hero_col)
        geo_hero_row.addStretch()
        rl.addWidget(self._geo_hero)

        self._choose_header = HeaderLabel("Choose data to import")
        rl.addWidget(self._choose_header)

        files_holder = QWidget()
        self._file_rows_layout = QVBoxLayout(files_holder)
        self._file_rows_layout.setContentsMargins(0, 0, 0, 0)
        self._file_rows_layout.setSpacing(8)
        rl.addWidget(files_holder)

        meta_files_group = SettingsGroup("Supplementary Metadata Files", collapsible=True, expanded=True)
        meta_info = SecondaryLabel("These files may contain cell type annotations, cluster info, or other metadata from the original authors.")
        meta_info.setWordWrap(True)
        meta_files_group.add_widget(meta_info)

        self.metadata_files_list = QVBoxLayout()
        meta_files_group.add_layout(self.metadata_files_list)

        self.no_metadata_label = SecondaryLabel("No metadata files detected")
        meta_files_group.add_widget(self.no_metadata_label)

        self.download_metadata_btn = QPushButton("Download Selected Metadata")
        self.download_metadata_btn.clicked.connect(self._download_metadata_files)
        self.download_metadata_btn.setEnabled(False)
        meta_files_group.add_widget(self.download_metadata_btn)

        rl.addWidget(meta_files_group)

        # Selected-file summary bar + primary download action
        self._selected_bar = QFrame()
        self._selected_bar.setObjectName("info_panel")
        bar_row = QHBoxLayout(self._selected_bar)
        bar_row.setContentsMargins(16, 12, 16, 12)
        bar_row.setSpacing(24)
        self._sel_tile = IconTile("file", "", size=48, icon_size=18)
        bar_row.addWidget(self._sel_tile)
        sel_col = QVBoxLayout()
        sel_col.setSpacing(2)
        sel_col.addWidget(CaptionLabel("Selected file"))
        self._sel_name = QLabel("")
        self._sel_name.setProperty("role", "file_row_name")
        sel_col.addWidget(self._sel_name)
        bar_row.addLayout(sel_col, stretch=1)
        type_col = QVBoxLayout()
        type_col.setSpacing(2)
        type_col.addWidget(CaptionLabel("Type"))
        self._sel_type = QLabel("")
        self._sel_type.setProperty("role", "kv_value")
        type_col.addWidget(self._sel_type)
        bar_row.addLayout(type_col)
        size_col = QVBoxLayout()
        size_col.setSpacing(2)
        size_col.addWidget(CaptionLabel("Size"))
        self._sel_size = QLabel("")
        self._sel_size.setProperty("role", "kv_value")
        size_col.addWidget(self._sel_size)
        bar_row.addLayout(size_col)
        self.download_btn = PrimaryButton("Download & Import")
        self.download_btn.setMinimumHeight(34)
        self.download_btn.setToolTip(
            "Download the selected file to the project's raw_data folder;\n"
            "KOSMIC then detects the format and converts it to h5ad.")
        self.download_btn.clicked.connect(self._start_download)
        self.download_btn.setEnabled(False)
        bar_row.addWidget(self.download_btn)
        rl.addWidget(self._selected_bar)

        # Batch / cancel actions (appear contextually)
        btn_layout = QHBoxLayout()

        self.download_all_10x_btn = QPushButton("Download All 10X")
        self.download_all_10x_btn.clicked.connect(self._download_all_10x_samples)
        self.download_all_10x_btn.setToolTip("Download all 10X sample files (matrix, features, barcodes)")
        self.download_all_10x_btn.hide()
        btn_layout.addWidget(self.download_all_10x_btn)

        self.download_all_matrix_btn = QPushButton("Download All Matrix TXT")
        self.download_all_matrix_btn.clicked.connect(self._download_all_matrix_files)
        self.download_all_matrix_btn.setToolTip("Download all matrix.txt.gz files from this dataset")
        self.download_all_matrix_btn.hide()
        btn_layout.addWidget(self.download_all_matrix_btn)

        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self._cancel)
        self.cancel_btn.setEnabled(False)
        btn_layout.addWidget(self.cancel_btn)

        btn_layout.addStretch()
        rl.addLayout(btn_layout)
        rl.addStretch()

        self._geo_body.addWidget(geo_landing)        # index 0
        self._geo_body.addWidget(self._geo_results)  # index 1
        gl.addWidget(self._geo_body, stretch=1)

        geo_scroll.setWidget(geo_page)

        outer.addWidget(geo_scroll, 1)

    def _query_geo(self):
        gse = self.gse_input.text().strip().upper()
        if not gse.startswith("GSE"):
            dialogs.warning(self, "Error", "Enter a valid GSE number (e.g., GSE292067)")
            return

        self.query_btn.setEnabled(False)
        if self.status_label:
            self.status_label.setText("Querying GEO...")
        if self.progress_bar:
            self.progress_bar.setRange(0, 0)  # Indeterminate
        self._tab_status_label.setText("Querying GEO...")

        self._log(f"Querying GEO for {gse}...")

        self.query_worker = GEOQueryWorker(gse)
        run_worker(
            self.query_worker,
            on_finished=self._on_query_done,
            on_failed=self._on_query_failed,
        )

    def _on_query_done(self, payload):
        files, message = payload
        self.query_btn.setEnabled(True)
        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)

        self.geo_files = files
        self._populate_files_table()
        self._geo_hero_title.setText(self.gse_input.text().strip().upper())
        self._geo_hero_sub.setText(
            f"{len(self.data_files)} file(s) found for this accession")
        # Swap landing state for the fully-populated results page in one go
        self._geo_body.setCurrentIndex(1)
        if self.status_label:
            self.status_label.setText(f"Found {len(files)} files")
        self._tab_status_label.setText(f"Found {len(files)} files")
        self._log(message)

    def _on_query_failed(self, message: str):
        self.query_btn.setEnabled(True)
        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)
        self.geo_files = []
        self.data_files = []
        self._clear_file_rows()
        self._geo_body.setCurrentIndex(0)
        # Clear metadata section
        for cb in self.metadata_checkboxes:
            cb.deleteLater()
        self.metadata_checkboxes = []
        self.metadata_files = []
        self.no_metadata_label.show()
        self.download_metadata_btn.setEnabled(False)
        if self.status_label:
            self.status_label.setText("Query failed")
        self._tab_status_label.setText(f"Query failed: {message}")
        self._log(f"Error: {message}")
        dialogs.warning(self, "Query Failed", message)

    @staticmethod
    def _fmt_size(size: float) -> str:
        if size > 1e9:
            return f"{size / 1e9:.1f} GB"
        if size > 1e6:
            return f"{size / 1e6:.1f} MB"
        return f"{size / 1e3:.1f} KB"

    @staticmethod
    def _ext_caption(fmt: str) -> str:
        """Short file-type caption for an IconTile ('RDS', 'TAR', ...)."""
        low = fmt.lower()
        if 'h5ad' in low:
            return "h5ad"
        if 'h5seurat' in low:
            return "H5S"
        if 'rds' in low:
            return "RDS"
        if 'raw' in low or 'tar' in low:
            return "TAR"
        if 'hdf5' in low:
            return "H5"
        if 'mtx' in low:
            return "MTX"
        if 'count' in low:
            return "TXT"
        return "FILE"

    @staticmethod
    def _compat_text(fmt: str) -> tuple[str, bool]:
        """One-line import outlook for a GEO file format."""
        low = fmt.lower()
        if 'h5ad' in low:
            return "✓ Ready to use: no conversion needed", True
        if 'rds' in low or 'h5seurat' in low:
            return "✓ Compatible: will convert to h5ad", True
        if 'raw' in low or 'tar' in low:
            return "✓ Raw 10X archive: will extract and convert", True
        if low in ('hdf5', 'mtx', 'count matrix'):
            return "✓ Compatible: will convert to h5ad", True
        return "Unrecognised format: may need manual conversion", False

    def _clear_file_rows(self):
        while self._file_rows_layout.count():
            item = self._file_rows_layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

    def _populate_files_table(self):
        """Rebuild the selectable file rows from self.geo_files."""
        self._clear_file_rows()
        self.button_group = QButtonGroup(self)

        # Sort by preference: h5ad > RDS > RAW.tar > others
        def sort_key(f):
            fmt = f['format'].lower()
            if 'h5ad' in fmt:
                return 0
            elif 'rds' in fmt:
                return 1
            elif 'raw' in fmt:
                return 2
            return 3

        sorted_files = sorted(self.geo_files, key=sort_key)

        # Metadata files go in their own section below
        data_files = [f for f in sorted_files if f['format'] not in ['Metadata']]

        for i, f in enumerate(data_files):
            radio = QRadioButton()
            self.button_group.addButton(radio, i)
            row = _FileRow(radio)
            h = QHBoxLayout(row)
            h.setContentsMargins(12, 8, 16, 8)
            h.setSpacing(12)
            h.addWidget(radio)
            h.addWidget(IconTile("file", self._ext_caption(f['format']),
                                 size=44, icon_size=16))
            name_col = QVBoxLayout()
            name_col.setSpacing(1)
            name_label = QLabel(f['name'])
            name_label.setProperty("role", "file_row_name")
            name_col.addWidget(name_label)
            compat_text, compat_ok = self._compat_text(f['format'])
            compat = StatusLabel(
                compat_text, state='success' if compat_ok else 'warning')
            compat.hide()
            name_col.addWidget(compat)
            row._compat_label = compat
            h.addLayout(name_col, stretch=1)
            h.addWidget(SecondaryLabel(f['format']))
            size_label = QLabel(self._fmt_size(f['size']))
            size_label.setProperty("role", "kv_value")
            h.addWidget(size_label)
            radio.toggled.connect(
                lambda checked, r=row, fi=f:
                self._on_file_radio_toggled(checked, r, fi))
            self._file_rows_layout.addWidget(row)
            if i == 0:  # Select first (best) by default
                radio.setChecked(True)

        has_files = len(data_files) > 0
        self._choose_header.setVisible(has_files)
        self._selected_bar.setVisible(has_files)
        self.download_btn.setEnabled(has_files and self.project_dir is not None)
        self.geo_files = sorted_files  # Keep all files for metadata extraction
        self.data_files = data_files  # Keep data files for selection

        # Detect 10X sample groups (files like sample1.matrix.mtx.gz, sample1.features.tsv.gz, etc.)
        self._detect_10x_sample_groups()

        # Detect matrix.txt.gz files for batch download
        self._detect_matrix_txt_samples()

        # Populate metadata files section
        self._populate_metadata_files()

    def _on_file_radio_toggled(self, checked, row, file_info):
        row.set_selected(checked)
        row._compat_label.setVisible(checked)
        if checked:
            self._sel_tile.set_caption(self._ext_caption(file_info['format']))
            self._sel_name.setText(file_info['name'])
            self._sel_type.setText(file_info['format'])
            self._sel_size.setText(self._fmt_size(file_info['size']))

    def _populate_metadata_files(self):
        """Populate the metadata files checkboxes."""
        # Clear existing checkboxes
        for cb in self.metadata_checkboxes:
            cb.deleteLater()
        self.metadata_checkboxes = []
        self.metadata_files = []

        # Find metadata files from the original query results
        for f in self.geo_files:
            if f['format'] in ['Metadata']:
                self.metadata_files.append(f)

        if self.metadata_files:
            self.no_metadata_label.hide()
            for mf in self.metadata_files:
                size = mf['size']
                if size > 1e6:
                    size_str = f"{size/1e6:.1f} MB"
                else:
                    size_str = f"{size/1e3:.1f} KB"

                cb = QCheckBox(f"{mf['name']} ({mf['format']}, {size_str})")
                cb.setProperty('file_info', mf)
                self.metadata_files_list.addWidget(cb)
                self.metadata_checkboxes.append(cb)

            self.download_metadata_btn.setEnabled(True)
        else:
            self.no_metadata_label.show()
            self.download_metadata_btn.setEnabled(False)

    def _download_metadata_files(self):
        """Download selected metadata files."""
        if not self.project_dir:
            dialogs.warning(self, "Error", "Select a project directory first")
            return

        selected = [cb.property('file_info') for cb in self.metadata_checkboxes if cb.isChecked()]
        if not selected:
            dialogs.warning(self, "No Selection", "Please check at least one metadata file to download")
            return

        raw_dir = raw_data_dir(self.project_dir)
        raw_dir.mkdir(exist_ok=True)

        self._log(f"Downloading {len(selected)} metadata file(s)...")
        self.download_metadata_btn.setEnabled(False)

        # Download files sequentially
        self._metadata_download_queue = selected.copy()
        self._download_next_metadata()

    def _download_next_metadata(self):
        """Download the next metadata file in the queue."""
        if not self._metadata_download_queue:
            self.download_metadata_btn.setEnabled(True)
            self._log("All metadata files downloaded!")
            dialogs.info(self, "Download Complete",
                "Metadata files downloaded to raw_data folder.\n\n"
                "You can apply this metadata to your h5ad file in the Convert tab.")
            return

        file_info = self._metadata_download_queue.pop(0)
        raw_dir = raw_data_dir(self.project_dir)
        output_path = raw_dir / file_info['name']

        self._log(f"Downloading: {file_info['name']}")

        self._current_meta_worker = DownloadWorker(file_info['url'], str(output_path), file_info.get('size', 0))
        run_worker(
            self._current_meta_worker,
            on_finished=self._on_metadata_download_done,
            on_failed=self._on_metadata_download_failed,
            on_progress=lambda s: self._log(f"  {s}"),
        )

    def _on_metadata_download_done(self, message: str):
        """Handle metadata file download completion (success)."""
        self._log(f"  Done: {message}")
        self._advance_metadata_queue()

    def _on_metadata_download_failed(self, message: str):
        """Handle metadata file download failure."""
        self._log(f"  Failed: {message}")
        self._advance_metadata_queue()

    def _advance_metadata_queue(self):
        if self._current_meta_worker:
            self._current_meta_worker.deleteLater()
            self._current_meta_worker = None
        QTimer.singleShot(100, self._download_next_metadata)

    def _detect_10x_sample_groups(self):
        """Detect 10X sample groups from file names (e.g., sample1.matrix.mtx.gz)."""
        self.detected_10x_samples = {}

        # Patterns for 10X files with sample prefixes
        # Matches: prefix.matrix.mtx.gz, prefix_matrix.mtx.gz, prefix.features.tsv.gz, etc.
        matrix_pattern = re.compile(r'^(.+?)[._]matrix\.mtx(?:\.gz)?$', re.IGNORECASE)
        features_pattern = re.compile(r'^(.+?)[._](?:features|genes)\.tsv(?:\.gz)?$', re.IGNORECASE)
        barcodes_pattern = re.compile(r'^(.+?)[._]barcodes\.tsv(?:\.gz)?$', re.IGNORECASE)

        # Find all matrix files and extract sample prefixes
        for f in self.data_files:
            fname = f['name']
            match = matrix_pattern.match(fname)
            if match:
                prefix = match.group(1)
                if prefix not in self.detected_10x_samples:
                    self.detected_10x_samples[prefix] = {'matrix': None, 'features': None, 'barcodes': None}
                self.detected_10x_samples[prefix]['matrix'] = f

        # Find matching features and barcodes files
        for f in self.data_files:
            fname = f['name']

            match = features_pattern.match(fname)
            if match:
                prefix = match.group(1)
                if prefix in self.detected_10x_samples:
                    self.detected_10x_samples[prefix]['features'] = f

            match = barcodes_pattern.match(fname)
            if match:
                prefix = match.group(1)
                if prefix in self.detected_10x_samples:
                    self.detected_10x_samples[prefix]['barcodes'] = f

        # Filter to only complete triplets
        complete_samples = {k: v for k, v in self.detected_10x_samples.items()
                          if v['matrix'] and v['features'] and v['barcodes']}
        self.detected_10x_samples = complete_samples

        # Show/hide the button based on detection
        n_samples = len(self.detected_10x_samples)
        if n_samples >= 2:
            n_files = n_samples * 3
            self.download_all_10x_btn.setText(f"Download All 10X ({n_samples} samples, {n_files} files)")
            self.download_all_10x_btn.show()
            self.download_all_10x_btn.setEnabled(self.project_dir is not None)
            self._log(f"Detected {n_samples} 10X sample groups: {', '.join(sorted(self.detected_10x_samples.keys())[:5])}{'...' if n_samples > 5 else ''}")
        else:
            self.download_all_10x_btn.hide()
            self.detected_10x_samples = {}

    def _detect_matrix_txt_samples(self):
        """Detect matrix.txt.gz files for batch download."""
        self.detected_matrix_files = []

        # Pattern for matrix.txt.gz files
        matrix_pattern = re.compile(r'.*matrix.*\.txt(?:\.gz)?$', re.IGNORECASE)

        for f in self.data_files:
            fname = f['name']
            if matrix_pattern.match(fname):
                self.detected_matrix_files.append(f)

        # Show/hide the button based on detection
        n_files = len(self.detected_matrix_files)
        if n_files >= 1:
            self.download_all_matrix_btn.setText(f"Download All Matrix TXT ({n_files} files)")
            self.download_all_matrix_btn.show()
            self.download_all_matrix_btn.setEnabled(self.project_dir is not None)
            self._log(f"Detected {n_files} matrix.txt files for download")
        else:
            self.download_all_matrix_btn.hide()
            self.detected_matrix_files = []

    def _download_all_matrix_files(self):
        """Download all detected matrix.txt.gz files."""
        if not self.project_dir:
            dialogs.warning(self, "Error", "Select a project directory first")
            return

        if not self.detected_matrix_files:
            dialogs.warning(self, "Error", "No matrix files detected")
            return

        n_files = len(self.detected_matrix_files)

        if not dialogs.confirm(self, "Confirm Download", f"Download {n_files} matrix.txt.gz file(s)?\n\n"
            f"Files: {', '.join([f['name'] for f in self.detected_matrix_files[:5]])}{'...' if n_files > 5 else ''}\n\n"
            "Files will be auto-converted after download."):
            return

        raw_dir = raw_data_dir(self.project_dir)
        raw_dir.mkdir(exist_ok=True)

        # Build download queue
        self._matrix_download_queue = list(self.detected_matrix_files)
        self._matrix_total_files = len(self._matrix_download_queue)
        self._matrix_completed_files = 0
        self._matrix_failed_files = []
        self._batch_total_files = self._matrix_total_files
        self._batch_completed_files = 0

        self._log(f"Starting batch download of {self._matrix_total_files} matrix files...")
        self.download_all_matrix_btn.setEnabled(False)
        self.download_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)

        self._download_next_matrix_file()

    def _download_next_matrix_file(self):
        """Download the next matrix file in the queue."""
        if not self._matrix_download_queue:
            # All done
            self.download_all_matrix_btn.setEnabled(True)
            self.download_btn.setEnabled(True)
            self.cancel_btn.setEnabled(False)

            if self._matrix_failed_files:
                self._log(f"Batch download complete: {self._matrix_completed_files} succeeded, {len(self._matrix_failed_files)} failed")
                dialogs.warning(self, "Download Complete",
                    f"Downloaded {self._matrix_completed_files} of {self._matrix_total_files} files.\n\n"
                    f"Failed: {', '.join(self._matrix_failed_files[:5])}{'...' if len(self._matrix_failed_files) > 5 else ''}")
            else:
                self._log(f"Batch download complete: all {self._matrix_completed_files} files downloaded!")
                dialogs.info(self, "Download Complete",
                    f"Successfully downloaded all {self._matrix_completed_files} matrix files!\n\n"
                    "Convert them from the study's Dataset screen (scRNA "
                    "-> Load Data).")

            self._refresh_local()
            self._auto_convert()
            return

        file_info = self._matrix_download_queue.pop(0)
        raw_dir = raw_data_dir(self.project_dir)
        output_path = raw_dir / file_info['name']

        progress_text = f"[{self._matrix_completed_files + 1}/{self._matrix_total_files}]"
        self._log(f"{progress_text} Downloading: {file_info['name']}")
        dl_text = f"Downloading {self._matrix_completed_files + 1}/{self._matrix_total_files}: {file_info['name']}"
        if self.status_label:
            self.status_label.setText(dl_text)
        self._tab_status_label.setText(dl_text)

        self._current_matrix_worker = DownloadWorker(file_info['url'], str(output_path), file_info.get('size', 0))
        self._current_matrix_filename = file_info['name']
        self._current_matrix_worker.bytes_progress.connect(self._on_progress)
        run_worker(
            self._current_matrix_worker,
            on_finished=self._on_matrix_file_download_done,
            on_failed=self._on_matrix_file_download_failed,
        )

    def _on_matrix_file_download_done(self, _message: str):
        """Handle successful completion of a single matrix file download."""
        self._matrix_completed_files += 1
        self._log(f"  ✓ {self._current_matrix_filename}")
        self._advance_matrix_queue()

    def _on_matrix_file_download_failed(self, message: str):
        """Handle failure of a single matrix file download."""
        self._matrix_failed_files.append(self._current_matrix_filename)
        self._log(f"  ✗ {self._current_matrix_filename}: {message}")
        self._advance_matrix_queue()

    def _advance_matrix_queue(self):
        self._batch_completed_files = self._matrix_completed_files
        if self._current_matrix_worker:
            self._current_matrix_worker.deleteLater()
            self._current_matrix_worker = None
        self._current_matrix_filename = None
        # 500ms delay respects NCBI rate limits
        QTimer.singleShot(500, self._download_next_matrix_file)

    def _download_all_10x_samples(self):
        """Download all detected 10X sample files."""
        if not self.project_dir:
            dialogs.warning(self, "Error", "Select a project directory first")
            return

        if not self.detected_10x_samples:
            dialogs.warning(self, "Error", "No 10X samples detected")
            return

        n_samples = len(self.detected_10x_samples)
        n_files = n_samples * 3

        if not dialogs.confirm(self, "Confirm Download", f"Download {n_files} files for {n_samples} 10X samples?\n\n"
            f"Samples: {', '.join(sorted(self.detected_10x_samples.keys())[:10])}{'...' if n_samples > 10 else ''}\n\n"
            "This may take a while for large datasets."):
            return

        raw_dir = raw_data_dir(self.project_dir)
        raw_dir.mkdir(exist_ok=True)

        # Build download queue with all files
        self._10x_download_queue = []
        for sample_name, files in sorted(self.detected_10x_samples.items()):
            for file_type in ['matrix', 'features', 'barcodes']:
                if files[file_type]:
                    self._10x_download_queue.append(files[file_type])

        self._10x_total_files = len(self._10x_download_queue)
        self._10x_completed_files = 0
        self._10x_failed_files = []
        self._batch_total_files = self._10x_total_files
        self._batch_completed_files = 0

        self._log(f"Starting batch download of {self._10x_total_files} files...")
        self.download_all_10x_btn.setEnabled(False)
        self.download_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)

        self._download_next_10x_file()

    def _download_next_10x_file(self):
        """Download the next 10X file in the queue."""
        if not self._10x_download_queue:
            # All done
            self.download_all_10x_btn.setEnabled(True)
            self.download_btn.setEnabled(True)
            self.cancel_btn.setEnabled(False)

            if self._10x_failed_files:
                self._log(f"Batch download complete: {self._10x_completed_files} succeeded, {len(self._10x_failed_files)} failed")
                dialogs.warning(self, "Download Complete",
                    f"Downloaded {self._10x_completed_files} of {self._10x_total_files} files.\n\n"
                    f"Failed: {', '.join(self._10x_failed_files[:5])}{'...' if len(self._10x_failed_files) > 5 else ''}")
            else:
                self._log(f"Batch download complete: all {self._10x_completed_files} files downloaded!")
                dialogs.info(self, "Download Complete",
                    f"Successfully downloaded all {self._10x_completed_files} files!\n\n"
                    "Convert them from the study's Dataset screen (scRNA "
                    "-> Load Data).")

            self._refresh_local()
            self._auto_convert()
            return

        file_info = self._10x_download_queue.pop(0)
        raw_dir = raw_data_dir(self.project_dir)
        output_path = raw_dir / file_info['name']

        progress_text = f"[{self._10x_completed_files + 1}/{self._10x_total_files}]"
        self._log(f"{progress_text} Downloading: {file_info['name']}")
        dl_text = f"Downloading {self._10x_completed_files + 1}/{self._10x_total_files}: {file_info['name']}"
        if self.status_label:
            self.status_label.setText(dl_text)
        self._tab_status_label.setText(dl_text)

        self._current_10x_worker = DownloadWorker(file_info['url'], str(output_path), file_info.get('size', 0))
        self._current_10x_filename = file_info['name']
        self._current_10x_worker.bytes_progress.connect(self._on_progress)
        run_worker(
            self._current_10x_worker,
            on_finished=self._on_10x_file_download_done,
            on_failed=self._on_10x_file_download_failed,
        )

    def _on_10x_file_download_done(self, _message: str):
        """Handle successful completion of a single 10X file download."""
        self._10x_completed_files += 1
        self._log(f"  ✓ {self._current_10x_filename}")
        self._advance_10x_queue()

    def _on_10x_file_download_failed(self, message: str):
        """Handle failure of a single 10X file download."""
        self._10x_failed_files.append(self._current_10x_filename)
        self._log(f"  ✗ {self._current_10x_filename}: {message}")
        self._advance_10x_queue()

    def _advance_10x_queue(self):
        self._batch_completed_files = self._10x_completed_files
        if self._current_10x_worker:
            self._current_10x_worker.deleteLater()
            self._current_10x_worker = None
        self._current_10x_filename = None
        # 500ms delay respects NCBI rate limits and lets the event loop process
        QTimer.singleShot(500, self._download_next_10x_file)

    def _get_selected_file(self):
        checked_id = self.button_group.checkedId()
        if 0 <= checked_id < len(self.data_files):
            return self.data_files[checked_id]
        return None

    def _start_download(self):
        if not self.project_dir:
            dialogs.warning(self, "Error", "Select a project directory first (Ctrl+O)")
            return

        selected = self._get_selected_file()
        if not selected:
            dialogs.warning(self, "Error", "Select a file to download")
            return

        raw_dir = raw_data_dir(self.project_dir)
        raw_dir.mkdir(exist_ok=True)
        output_path = raw_dir / selected['name']

        if output_path.exists():
            if not dialogs.confirm(self, "Exists", f"Overwrite {selected['name']}?"):
                return

        self.selected_file = selected
        self._batch_total_files = 1  # single file, not a batch
        self._batch_completed_files = 0
        self.download_worker = DownloadWorker(selected['url'], str(output_path), selected.get('size', 0))
        self.download_worker.bytes_progress.connect(self._on_progress)
        run_worker(
            self.download_worker,
            on_finished=self._on_download_done,
            on_failed=self._on_download_failed,
            on_progress=self._on_status,
        )

        # Reset and enable progress bar
        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)

        self._tab_status_label.setText(f"Downloading: {selected['name']}")

        self.download_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self._log(f"Downloading: {selected['name']}")

    def _cancel(self):
        if self.download_worker and self.download_worker.isRunning():
            self.download_worker.cancel()
        # Also cancel batch 10X download if running
        if self._current_10x_worker and self._current_10x_worker.isRunning():
            self._current_10x_worker.cancel()
            self._10x_download_queue = []  # Clear remaining queue
            self._log("Batch download cancelled")
            self.download_all_10x_btn.setEnabled(True)
            self.download_btn.setEnabled(True)
        # Also cancel batch matrix download if running
        if self._current_matrix_worker and self._current_matrix_worker.isRunning():
            self._current_matrix_worker.cancel()
            self._matrix_download_queue = []  # Clear remaining queue
            self._log("Matrix batch download cancelled")
            self.download_all_matrix_btn.setEnabled(True)
            self.download_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)

    def _on_progress(self, downloaded: int, total: int):
        if not self.progress_bar:
            return
        if total > 0:
            self.progress_bar.setRange(0, 100)
            # For batch downloads, compute cumulative progress across all files
            if self._batch_total_files > 1:
                file_pct = downloaded / total
                overall_pct = (self._batch_completed_files + file_pct) / self._batch_total_files * 100
                self.progress_bar.setValue(int(overall_pct))
            else:
                self.progress_bar.setValue(int(downloaded * 100 / total))
        else:
            # Unknown total — show indeterminate (pulsing) progress bar
            self.progress_bar.setRange(0, 0)

    def _on_status(self, status: str):
        if self.status_label:
            self.status_label.setText(status)
        self._tab_status_label.setText(status)

    def _on_download_done(self, message: str):
        self.download_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self._log(message)

        if self.progress_bar:
            self.progress_bar.setValue(100)
        self._tab_status_label.setText("Download complete")

        # An h5ad is ready to use directly: register it as the study's
        # processed dataset (raw copy stays in raw_data/).
        if self.selected_file and 'h5ad' in self.selected_file['format'].lower():
            h5ad_src = raw_data_dir(self.project_dir) / self.selected_file['name']
            if h5ad_src.exists():
                import shutil
                dest = processed_data_dir(self.project_dir)
                dest.mkdir(exist_ok=True)
                shutil.copy2(h5ad_src, dest / h5ad_src.name)
                self._log(f"Registered processed dataset: {h5ad_src.name}")
        self._on_raw_changed()

    def _on_download_failed(self, message: str):
        self.download_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self._log(message)
        if self.progress_bar:
            self.progress_bar.setValue(0)
        self._tab_status_label.setText("Download failed")
        dialogs.warning(
            self, "Download Failed",
            f"{message}\n\n"
            "Check your internet connection and the accession, then try "
            "again. Some GEO files are also mirrored as supplementary "
            "files on the study's GEO page.")

