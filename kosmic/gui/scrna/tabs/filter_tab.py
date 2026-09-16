# Subset step: extract selected cell types into a new study.
#
# The user picks a cell-type column and ticks the types to keep;
# 'FilterWorker' writes those cells to a new sibling study folder named
# '<study>_<descriptor>', clears embeddings and clusterings (they were
# computed over the full dataset), records the parent in the subset's
# provenance, and the app switches to the new study. Subsets are what
# per-cell-type DE and the meta-analysis usually run on.

import os
import re
import fnmatch
from pathlib import Path
from PyQt6.QtWidgets import (
    QHBoxLayout, QLabel,
    QPushButton, QListWidget, QListWidgetItem,
    QLineEdit, QCheckBox, QHeaderView, QAbstractItemView,
)
from PyQt6.QtCore import Qt, pyqtSignal, QTimer
import pandas as pd
from kosmic.gui.shared.theme import NoScrollComboBox
from kosmic.gui.shared.widgets import BaseWorker, Column, HintLabel, PrimaryButton, ResultsTable, SettingsGroup, SidebarPage
from kosmic.gui.shared import dialogs, run_worker


class FilterWorker(BaseWorker):
    """
    Worker for filtering cells by type.

    Emits 'finished_ok' with a tuple '(adata_filtered, summary_message)'.
    """

    def __init__(self, adata, cell_type_column, selected_types, output_path):
        super().__init__()
        self.adata = adata
        self.cell_type_column = cell_type_column
        self.selected_types = selected_types
        self.output_path = output_path

    def _run(self):
        self.progress.emit("Filtering cells...")
        self.progress_pct.emit(10)

        # Get original counts
        n_cells_original = self.adata.n_obs
        n_genes_original = self.adata.n_vars

        # Create mask for selected cell types
        mask = self.adata.obs[self.cell_type_column].isin(self.selected_types)
        n_cells_selected = mask.sum()

        self.progress.emit(
            f"Selected {n_cells_selected:,} cells from {n_cells_original:,}")
        self.progress_pct.emit(30)

        if n_cells_selected == 0:
            raise RuntimeError("No cells match the selected types!")

        # Subset the data
        adata_filtered = self.adata[mask].copy()
        self.progress_pct.emit(50)

        # The subset keeps the parent's full gene list. A gene with no
        # counts in the kept cells is a measured zero, not a missing gene:
        # HVG selection never picks it, the DE gene filter drops it per run
        # with the reason recorded, and all-zero columns cost nothing in a
        # sparse matrix. Dropping them gave every subset its own gene
        # universe, which shrank the shared-gene intersection the atlas and
        # the meta-analysis are built on.
        n_genes_filtered = adata_filtered.n_vars
        self.progress_pct.emit(70)

        # Embeddings, neighbour graph, and cluster labels are computed
        # over the full pre-subset cell set; they slice consistently but
        # are biologically inappropriate for the new cell set.
        self.progress.emit("Clearing stale embeddings and clustering results...")
        for key in ('X_pca', 'X_pca_harmony', 'X_umap', 'X_tsne'):
            if key in adata_filtered.obsm:
                del adata_filtered.obsm[key]
        for key in list(adata_filtered.obsp.keys()):
            del adata_filtered.obsp[key]
        for key in ('neighbors', 'leiden', 'umap', 'pca',
                    'rank_genes_groups'):
            if key in adata_filtered.uns:
                del adata_filtered.uns[key]
        for col in ('leiden', 'clusters', 'seurat_clusters'):
            if col in adata_filtered.obs.columns:
                del adata_filtered.obs[col]

        # Add filtering metadata
        adata_filtered.uns['filtering_info'] = {
            'original_cells': n_cells_original,
            'original_genes': n_genes_original,
            'filtered_cells': adata_filtered.n_obs,
            'filtered_genes': n_genes_filtered,
            'cell_type_column': self.cell_type_column,
            'selected_types': list(self.selected_types),
            'genes_removed': n_genes_original - n_genes_filtered
        }

        # Save filtered data with file-size progress monitoring
        import time
        import threading

        n_cells = adata_filtered.n_obs
        n_genes = adata_filtered.n_vars
        self.progress.emit(
            f"Saving {n_cells:,} cells x {n_genes:,} genes "
            f"to {Path(self.output_path).name}...")
        self.progress_pct.emit(85)
        Path(self.output_path).parent.mkdir(parents=True, exist_ok=True)

        save_error = [None]

        def _do_save():
            try:
                adata_filtered.write_h5ad(self.output_path)
            except Exception as e:
                save_error[0] = e

        save_thread = threading.Thread(target=_do_save, daemon=True)
        t0 = time.time()
        save_thread.start()

        # Poll file size while saving
        out_path = Path(self.output_path)
        while save_thread.is_alive():
            save_thread.join(timeout=2)
            if out_path.exists():
                mb = out_path.stat().st_size / (1024 * 1024)
                elapsed = time.time() - t0
                self.progress.emit(
                    f"Writing {out_path.name}... "
                    f"{mb:.0f} MB ({elapsed:.0f}s)")

        if save_error[0]:
            raise save_error[0]

        dt = time.time() - t0
        actual_mb = out_path.stat().st_size / (1024 * 1024)
        self.progress.emit(f"Saved {actual_mb:.0f} MB in {dt:.0f}s")
        self.progress_pct.emit(100)

        stats = {
            "n_cells_original": n_cells_original,
            "n_cells_filtered": adata_filtered.n_obs,
            "n_genes_original": n_genes_original,
            "n_genes_filtered": n_genes_filtered,
        }

        return (adata_filtered, stats)


class FilterTab(SidebarPage):
    """Tab for filtering cells by type."""

    help_id = "scrna/subset"
    CONTENT_MARGINS = (10, 10, 10, 10)

    # Signal to notify other tabs when filtered data is ready
    filtered_data_ready = pyqtSignal(str)  # path to filtered h5ad
    log_message = pyqtSignal(str)  # routed to OutputPanel

    def __init__(self, main_window=None):
        super().__init__(main_window)
        self.main_window = main_window
        self.adata = None
        self.h5ad_path = None
        self.cell_types = []
        self._pattern_timer = QTimer(self)
        self._pattern_timer.setSingleShot(True)
        self._pattern_timer.setInterval(400)
        self._pattern_timer.timeout.connect(self._apply_pattern_debounced)
        self.setup_ui()

    def reset_state(self):
        """Clear all cached state (called on raw data reload)."""
        self.adata = None
        self.h5ad_path = None
        self.cell_types = []
        self._last_adata_version = -1
        # Reset UI
        self.cell_type_list.clear()
        self.filter_btn.setEnabled(False)
        self.status_label.setText("Ready")

    def setup_ui(self):
        """
        Setup the UI components.

        Sidebar lives on 'self.sidebar_layout' (provided by
        'SidebarPage'). Right-side content goes on 'self.content'.
        """
        left_layout = self.sidebar_layout

        # ── Cell Types group ──
        cell_group = SettingsGroup("Cell Types", collapsible=True)
        self._cell_group = cell_group

        self.column_combo = NoScrollComboBox()
        self.column_combo.setToolTip("Select column with cell type labels.\n★ marks auto-detected columns.")
        self.column_combo.currentTextChanged.connect(self.on_column_changed)
        cell_group.add_row("Column:", self.column_combo)

        self.search_edit = QLineEdit()
        self.search_edit.setPlaceholderText("Search cell types...")
        self.search_edit.setToolTip("Supports wildcards: *Endothelial*")
        self.search_edit.textChanged.connect(self.filter_cell_types)
        cell_group.add_widget(self.search_edit)

        self.cell_type_list = QListWidget()
        self.cell_type_list.setSelectionMode(QListWidget.SelectionMode.NoSelection)
        self.cell_type_list.itemChanged.connect(self.update_selection_summary)
        cell_group.add_widget(self.cell_type_list)

        btn_layout = QHBoxLayout()
        btn_layout.setSpacing(4)
        self.select_all_btn = QPushButton("All")
        self.select_all_btn.clicked.connect(self.select_all)
        btn_layout.addWidget(self.select_all_btn)

        self.select_none_btn = QPushButton("None")
        self.select_none_btn.clicked.connect(self.select_none)
        btn_layout.addWidget(self.select_none_btn)

        self.invert_btn = QPushButton("Invert")
        self.invert_btn.clicked.connect(self.invert_selection)
        btn_layout.addWidget(self.invert_btn)
        cell_group.add_layout(btn_layout)

        left_layout.addWidget(cell_group, 1)

        # ── Quick Select group ──
        quick_group = SettingsGroup("Quick Select", collapsible=True)
        self._quick_group = quick_group

        self.pattern_edit = QLineEdit()
        self.pattern_edit.setPlaceholderText("e.g. Endothel*, *_EC")
        self.pattern_edit.textChanged.connect(self._on_pattern_changed)
        quick_group.add_widget(self.pattern_edit)

        common_row1 = QHBoxLayout()
        common_row1.setSpacing(4)
        for name, pat in [("Endo", "*ndothel*"), ("Fibro", "*ibro*"),
                          ("SMC", "*mooth*,*SMC*,*VSMC*")]:
            btn = QPushButton(name)
            btn.clicked.connect(lambda checked, p=pat: self.quick_select(p))
            common_row1.addWidget(btn)
        common_row1.addStretch()
        quick_group.add_layout(common_row1)

        left_layout.addWidget(quick_group)

        # ── New Study group ──
        # Subsetting creates a new sibling study folder under the active
        # project (e.g. GSE183852 -> GSE183852_DCMEC). The user picks a
        # short descriptor; the new study name = current accession +
        # underscore + descriptor.
        new_study_group = SettingsGroup("New Study", collapsible=True)
        self._new_study_group = new_study_group

        self.auto_name_check = QCheckBox("Auto-suggest descriptor")
        self.auto_name_check.setChecked(True)
        self.auto_name_check.stateChanged.connect(self.update_auto_filename)
        new_study_group.add_widget(self.auto_name_check)

        new_study_group.add_widget(QLabel("Descriptor:"))
        self.descriptor_edit = QLineEdit()
        self.descriptor_edit.setPlaceholderText("e.g. DCMEC")
        self.descriptor_edit.setToolTip(
            "Short label for this subset (no spaces). "
            "Combined with the active study name to make the new folder."
        )
        self.descriptor_edit.textChanged.connect(self._on_descriptor_changed)
        new_study_group.add_widget(self.descriptor_edit)

        self.new_study_preview = HintLabel("")
        self.new_study_preview.setWordWrap(True)
        new_study_group.add_widget(self.new_study_preview)

        left_layout.addWidget(new_study_group)

        self.filter_btn = PrimaryButton("Subset to New Study")
        self.filter_btn.setToolTip(
            "Create a new study folder in the project, save the subset "
            "there, and switch to it as the active study.\n"
            "Clustering results are cleared on the subset (they were "
            "computed on the full dataset) -- re-run Cluster afterwards."
        )
        self.filter_btn.setMinimumHeight(32)
        self.filter_btn.setEnabled(False)
        self.filter_btn.clicked.connect(self.run_filter)
        left_layout.addWidget(self.filter_btn)

        self.status_label = HintLabel("Ready")
        left_layout.addWidget(self.status_label)

        self.progress_bar = None  # wired to sidebar by main.py

        # ---- Right: selection summary ----
        right_layout = self.content_layout
        right_layout.setSpacing(8)

        # Data overview label
        self.data_overview_label = QLabel("No data loaded — complete the pipeline first.")
        self.data_overview_label.setWordWrap(True)
        right_layout.addWidget(self.data_overview_label)

        # Selection summary header
        self.selection_header = QLabel("No cell types selected")
        self.selection_header.setProperty("role", "selection_header")
        right_layout.addWidget(self.selection_header)

        # Selection breakdown table
        self.summary_table = ResultsTable()
        self.summary_table.set_schema([
            Column("Cell Type", "cell_type", "s"),
            Column("Cells",     "count",     ",d"),
            Column("%",         "pct",       ".1f"),
        ])
        self.summary_table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.Stretch
        )
        self.summary_table.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.ResizeMode.ResizeToContents
        )
        self.summary_table.horizontalHeader().setSectionResizeMode(
            2, QHeaderView.ResizeMode.ResizeToContents
        )
        self.summary_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.summary_table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.summary_table.verticalHeader().setVisible(False)
        right_layout.addWidget(self.summary_table, 1)

    def _update_summary(self):
        """Update the data overview label."""
        if self.adata is None:
            self.data_overview_label.setText(
                "No data loaded — complete the pipeline first."
            )
            return

        parts = []
        if self.h5ad_path:
            parts.append(os.path.basename(self.h5ad_path))
        parts.append(f"{self.adata.n_obs:,} cells, {self.adata.n_vars:,} genes")
        self.data_overview_label.setText(" · ".join(parts))

    def on_column_changed(self, text):
        """Handle cell type column change."""
        if not text or self.adata is None:
            return

        # Get actual column name (remove star prefix if present)
        col_name = self.column_combo.currentData()
        if col_name is None:
            col_name = text.replace("★ ", "")

        # Get unique cell types and their counts
        value_counts = self.adata.obs[col_name].value_counts()
        self.cell_types = list(value_counts.index)

        # Populate the list
        self.cell_type_list.clear()
        for cell_type in self.cell_types:
            count = value_counts[cell_type]
            item = QListWidgetItem(f"{cell_type} ({count:,})")
            item.setData(Qt.ItemDataRole.UserRole, cell_type)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Unchecked)
            self.cell_type_list.addItem(item)

        self.update_selection_summary()
        self.update_auto_filename()

    def filter_cell_types(self, text):
        """Filter displayed cell types based on search text."""
        if not text:
            # Show all
            for i in range(self.cell_type_list.count()):
                self.cell_type_list.item(i).setHidden(False)
            return

        # Support wildcards
        pattern = text.replace('*', '.*')
        try:
            regex = re.compile(pattern, re.IGNORECASE)
        except re.error:
            return

        for i in range(self.cell_type_list.count()):
            item = self.cell_type_list.item(i)
            cell_type = item.data(Qt.ItemDataRole.UserRole)
            matches = regex.search(cell_type) is not None
            item.setHidden(not matches)

    def select_all(self):
        """Select all visible cell types."""
        for i in range(self.cell_type_list.count()):
            item = self.cell_type_list.item(i)
            if not item.isHidden():
                item.setCheckState(Qt.CheckState.Checked)
        self.update_selection_summary()

    def select_none(self):
        """Deselect all cell types."""
        for i in range(self.cell_type_list.count()):
            item = self.cell_type_list.item(i)
            item.setCheckState(Qt.CheckState.Unchecked)
        self.update_selection_summary()

    def invert_selection(self):
        """Invert the selection."""
        for i in range(self.cell_type_list.count()):
            item = self.cell_type_list.item(i)
            if item.checkState() == Qt.CheckState.Checked:
                item.setCheckState(Qt.CheckState.Unchecked)
            else:
                item.setCheckState(Qt.CheckState.Checked)
        self.update_selection_summary()

    def apply_pattern(self):
        """Select cell types matching the pattern."""
        pattern_text = self.pattern_edit.text().strip()
        if not pattern_text:
            return

        # Split by comma for multiple patterns
        patterns = [p.strip() for p in pattern_text.split(',')]

        for i in range(self.cell_type_list.count()):
            item = self.cell_type_list.item(i)
            cell_type = item.data(Qt.ItemDataRole.UserRole)

            # Check if matches any pattern
            matches = False
            for pattern in patterns:
                if fnmatch.fnmatch(cell_type.lower(), pattern.lower()):
                    matches = True
                    break

            if matches:
                item.setCheckState(Qt.CheckState.Checked)

        self.update_selection_summary()

    def _on_pattern_changed(self, text):
        """Live-match pattern as user types (debounced)."""
        if text.strip():
            self._pattern_timer.start()
        else:
            self._pattern_timer.stop()

    def _apply_pattern_debounced(self):
        """Called after typing pause — select matching types."""
        self.select_none()
        self.apply_pattern()

    def quick_select(self, pattern):
        """Apply a quick selection pattern."""
        self.select_none()
        self.pattern_edit.setText(pattern)

    def get_selected_types(self):
        """Get list of selected cell types."""
        selected = []
        for i in range(self.cell_type_list.count()):
            item = self.cell_type_list.item(i)
            if item.checkState() == Qt.CheckState.Checked:
                selected.append(item.data(Qt.ItemDataRole.UserRole))
        return selected

    def update_selection_summary(self):
        """Update the selection summary table."""
        if self.adata is None:
            self.selection_header.setText("No data loaded")
            self.summary_table.setRowCount(0)
            return

        selected = self.get_selected_types()
        col_name = self.column_combo.currentData()
        if col_name is None:
            text = self.column_combo.currentText()
            col_name = text.replace("★ ", "") if text else None

        if not col_name:
            self.selection_header.setText("Select a cell type column")
            self.summary_table.setRowCount(0)
            return

        total_cells = self.adata.n_obs

        if not selected:
            self.selection_header.setText(
                f"No cell types selected — {total_cells:,} cells total"
            )
            self.summary_table.setRowCount(0)
            return

        # Calculate selected cell count
        value_counts = self.adata.obs[col_name].value_counts()
        selected_cells = sum(value_counts.get(t, 0) for t in selected)

        self.selection_header.setText(
            f"Keeping {selected_cells:,} of {total_cells:,} cells "
            f"({selected_cells/total_cells*100:.1f}%) — "
            f"{len(selected)} type(s)"
        )

        # Populate table
        rows = []
        for cell_type in selected:
            count = value_counts.get(cell_type, 0)
            rows.append({
                'cell_type': str(cell_type),
                'count':     int(count),
                'pct':       count / total_cells * 100 if total_cells else 0,
            })

        right_align = Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter

        def _right_align(item, row, col_i):
            if col_i in (1, 2):
                item.setTextAlignment(right_align)

        self.summary_table.set_data(pd.DataFrame(rows), decorate=_right_align)

        # Update auto filename
        if self.auto_name_check.isChecked():
            self.update_auto_filename()

    def update_auto_filename(self):
        """Auto-suggest a short descriptor based on the selected cell types."""
        if not self.auto_name_check.isChecked():
            return

        selected = self.get_selected_types()
        if not selected:
            descriptor = ""
        elif len(selected) == 1:
            descriptor = re.sub(r'[^\w]', '_', selected[0])
        elif len(selected) <= 3:
            parts = [re.sub(r'[^\w]', '_', t)[:10] for t in selected]
            descriptor = '_'.join(parts)
        else:
            descriptor = f"subset_{len(selected)}"

        self.descriptor_edit.blockSignals(True)
        self.descriptor_edit.setText(descriptor)
        self.descriptor_edit.blockSignals(False)
        self._refresh_new_study_preview()

    def _on_descriptor_changed(self, text: str):
        # Hand-edit means the user owns the descriptor; stop auto-suggesting.
        if self.auto_name_check.isChecked():
            self.auto_name_check.blockSignals(True)
            self.auto_name_check.setChecked(False)
            self.auto_name_check.blockSignals(False)
        self._refresh_new_study_preview()

    def _new_study_name(self) -> str:
        """Compose ``{active_study_or_project}_{descriptor}`` from current state."""
        descriptor = self.descriptor_edit.text().strip()
        if not descriptor:
            return ""
        app = self.window()
        base = getattr(app, 'current_accession', None)
        if not base:
            project = getattr(app, 'current_project_dir', None)
            base = project.name if project is not None else ""
        return f"{base}_{descriptor}" if base else descriptor

    def _refresh_new_study_preview(self):
        name = self._new_study_name()
        if name:
            self.new_study_preview.setText(f"Will create: {name}/")
        else:
            self.new_study_preview.setText("Pick a descriptor.")

    def run_filter(self):
        """Subset the active dataset into a new sibling study folder."""
        if self.adata is None:
            dialogs.warning(self, "Error", "No data loaded.")
            return

        selected = self.get_selected_types()
        if not selected:
            dialogs.warning(self, "Error", "No cell types selected.")
            return

        col_name = self.column_combo.currentData()
        if col_name is None:
            text = self.column_combo.currentText()
            col_name = text.replace("★ ", "") if text else None

        if not col_name:
            dialogs.warning(self, "Error", "Please select a cell type column.")
            return

        app = self.window()
        project_dir = getattr(app, 'current_project_dir', None)
        if project_dir is None:
            dialogs.warning(
                self, "No Project",
                "Open a project in the Project tab before subsetting."
            )
            return

        new_study = self._new_study_name()
        if not new_study:
            dialogs.warning(
                self, "Missing Descriptor",
                "Type a short descriptor (e.g. DCMEC) for the new study."
            )
            return

        target_dir = project_dir / new_study
        if target_dir.exists() and any(target_dir.iterdir()):
            if not dialogs.confirm(
                self, "Folder Exists",
                f"'{new_study}/' already exists and contains files. Overwrite?"
            ):
                return

        # Confirm if many cells will be removed
        mask = self.adata.obs[col_name].isin(selected)
        selected_cells = mask.sum()
        removed_pct = (1 - selected_cells / self.adata.n_obs) * 100

        if removed_pct > 90:
            if not dialogs.confirm(
                self, "Confirm Subset",
                f"This will remove {removed_pct:.1f}% of cells.\n\n"
                f"Keeping {selected_cells:,} of {self.adata.n_obs:,} cells.\n\n"
                "Continue?",
            ):
                return

        # Scaffold the new study folder via AppWindow's contract. Save
        # the subset h5ad with a name that matches the new study so the
        # auto-loader / pipeline conventions stay consistent.
        app.add_study(new_study)
        from kosmic.paths import processed_h5ad_path
        output_path = str(processed_h5ad_path(target_dir, f"{new_study}.h5ad"))
        self._pending_new_study = new_study
        # Capture the source (parent) file + selection now, before the active
        # study switches, so the subset's provenance can point back to it.
        self._subset_source_h5ad = getattr(self.main_window, 'current_h5ad_path', None)
        self._subset_col = col_name
        self._subset_types = list(selected)

        self.filter_btn.setEnabled(False)
        if self.progress_bar:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)

        self.worker = FilterWorker(
            self.adata, col_name, selected, output_path
        )
        run_worker(
            self.worker,
            on_finished=self.on_filter_finished,
            on_failed=self.on_filter_failed,
            on_progress=self._on_status,
            on_progress_pct=self.progress_bar.setValue if self.progress_bar else None,
        )

    def _on_status(self, status: str):
        self.status_label.setText(status)
        self.log_message.emit(status)

    def on_filter_finished(self, payload):
        """Switch active study to the new folder and load the subset adata into scRNA."""
        self.filter_btn.setEnabled(True)
        adata_filtered, stats = payload

        new_study = getattr(self, '_pending_new_study', None)
        self._pending_new_study = None

        from kosmic.paths import processed_h5ad_path
        output_path = None

        if new_study:
            app = self.window()
            project_dir = getattr(app, 'current_project_dir', None)
            if project_dir is not None:
                output_path = str(
                    processed_h5ad_path(project_dir / new_study, f"{new_study}.h5ad")
                )
            # Switch active study first (this clears the workspace's adata
            # via set_project_directory), then restore the subset adata so
            # the user lands in the new study with the subset already loaded.
            if hasattr(app, 'set_active_study'):
                app.set_active_study(new_study)
            ws = self.main_window
            if ws is not None and hasattr(ws, 'set_adata'):
                ws.set_adata(adata_filtered, file_path=output_path)
                # This tab's own reference would otherwise keep the parent
                # alive next to the subset until the tab is next activated.
                self.adata = adata_filtered
                # Record the subset as a provenance stage on the new study,
                # pointing back to the parent file it was derived from.
                if hasattr(ws, 'record_provenance'):
                    ws.record_provenance('subset', {
                        'cell_type_column': getattr(self, '_subset_col', None),
                        'cell_types': getattr(self, '_subset_types', []),
                        'n_cells_before': int(stats.get('n_cells_original', 0)),
                        'n_cells_after': int(stats.get('n_cells_filtered', 0)),
                        'n_genes_before': int(stats.get('n_genes_original', 0)),
                        'n_genes_after': int(stats.get('n_genes_filtered', 0)),
                    }, source=getattr(self, '_subset_source_h5ad', None))

        cell_pct = (stats["n_cells_filtered"] / stats["n_cells_original"] * 100
                    if stats["n_cells_original"] else 0)
        gene_pct = (stats["n_genes_filtered"] / stats["n_genes_original"] * 100
                    if stats["n_genes_original"] else 0)

        if new_study:
            message = (
                f"New study created: {new_study}\n"
                f"It is now the active study.\n\n"
                f"Cells: {stats['n_cells_original']:,} -> {stats['n_cells_filtered']:,} "
                f"({cell_pct:.1f}%)\n"
                f"Genes: {stats['n_genes_original']:,} -> {stats['n_genes_filtered']:,} "
                f"({gene_pct:.1f}%)\n\n"
                f"Clustering results (PCA, UMAP, t-SNE, clusters, and "
                f"per-cluster marker genes) were cleared: they were "
                f"computed on the full dataset and are not valid for the "
                f"subset. Run the Cluster step on the new study to "
                f"recompute them."
            )
        else:
            message = (
                f"Subset complete.\n\n"
                f"Cells: {stats['n_cells_original']:,} -> {stats['n_cells_filtered']:,} "
                f"({cell_pct:.1f}%)\n"
                f"Genes: {stats['n_genes_original']:,} -> {stats['n_genes_filtered']:,} "
                f"({gene_pct:.1f}%)"
            )

        dialogs.info(self, "Subset complete", message)

        if output_path:
            self.filtered_data_ready.emit(output_path)

    def on_filter_failed(self, message):
        """Handle filter failure."""
        self.filter_btn.setEnabled(True)
        dialogs.error(self, "Error", message)

    def on_tab_activated(self):
        """Called when tab becomes visible — pull latest adata from workspace."""
        ws = self.main_window
        if ws is None or ws.current_adata is None:
            return
        version = getattr(ws, '_adata_version', 0)
        if version != getattr(self, '_last_adata_version', -1):
            self._last_adata_version = version
            self.set_data(ws.current_adata, ws.current_h5ad_path)

        # Ensure button is enabled if we have data and columns
        if self.adata is not None and self.column_combo.count() > 0:
            self.filter_btn.setEnabled(True)

    def set_data(self, adata, file_path):
        """Set data from another tab."""
        self.adata = adata
        self.h5ad_path = file_path

        # Trigger column update
        self.column_combo.clear()

        # Ranked keywords — earlier entries are higher priority
        ranked_keywords = [
            'cell_type', 'celltype', 'cell type',
            'cell_annotation', 'annotation', 'cluster', 'leiden', 'louvain',
        ]

        columns = list(self.adata.obs.columns)
        priority_cols = []
        other_cols = []

        for col in columns:
            col_lower = col.lower()
            if any(kw in col_lower for kw in ranked_keywords):
                priority_cols.append(col)
            elif self.adata.obs[col].dtype == 'object' or self.adata.obs[col].dtype.name == 'category':
                other_cols.append(col)

        # Sort priority columns by keyword rank (best match first)
        def _rank(col):
            cl = col.lower()
            for i, kw in enumerate(ranked_keywords):
                if cl == kw:
                    return i  # exact match ranks highest
            for i, kw in enumerate(ranked_keywords):
                if kw in cl:
                    return i + len(ranked_keywords)
            return len(ranked_keywords) * 2

        priority_cols.sort(key=_rank)

        for col in priority_cols:
            self.column_combo.addItem(f"★ {col}", col)
        for col in other_cols:
            if col not in priority_cols:
                self.column_combo.addItem(col, col)

        self._update_summary()

        if self.auto_name_check.isChecked():
            self.update_auto_filename()

        self.filter_btn.setEnabled(True)
