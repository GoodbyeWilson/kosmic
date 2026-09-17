# Combine dialog. Modal opened from the Project page; assembles the
# selected studies' h5ads into a master at
# project_dir/_master/processed_data/master.h5ad. The master then
# appears as a synthetic study row in the project's study table.
# Label propagation back to per-study h5ads is a separate action
# triggered from the master row's kebab menu.
from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QCursor
from PyQt6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDialog, QHBoxLayout,
    QHeaderView, QLabel, QPlainTextEdit, QProgressBar, QSpinBox,
    QStackedWidget, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from kosmic.combine import concat_studies, propagate_labels
from kosmic.combine.gene_merge import (
    find_merge_asymmetry, merge_map_from_provenance,
)
from kosmic.scrna.load.gene_overlap import read_var_names
from kosmic.gui.shared.widgets import (
    BaseWorker, HeaderLabel, PrimaryButton, SecondaryButton, SecondaryLabel,
)
from kosmic.gui.shared import run_worker
from kosmic.paths import (
    MASTER_ACCESSION, list_studies, master_h5ad_path, processed_data_dir, raw_data_dir,
)


# Condition-column resolution + agnostic defaults --------------------

_PREFERRED_CONDITION_COLS = (
    '_role', 'condition', 'Condition', 'group', 'Group', 'disease',
    'Disease', 'phenotype', 'Phenotype',
)

# Only generic, tissue-agnostic terms get auto-mapped. Disease-specific
# names (DCM / NICM / ICM / AS / NF etc.) stay blank for the user to set.
_GENERIC_CONTROL = {'control', 'donor', 'normal', 'healthy', 'wt', 'ctrl'}
_GENERIC_DISEASE = {'disease', 'diseased', 'patient', 'affected'}

_ROLE_OPTIONS = ('', 'control', 'disease', 'exclude')


def _default_role_for(value: str, existing_role_col_present: bool,
                      existing_value: Optional[str] = None) -> str:
    """Return the default role for a condition-column value.

    Generic terms map to canonical roles; everything else (cardiac-
    specific or otherwise tissue-jargon) stays blank for the user.
    When the source column was already '_role', the original value
    passes through unchanged.
    """
    if existing_role_col_present and existing_value in ('control', 'disease', 'exclude'):
        return existing_value
    v = value.strip().lower()
    if v in _GENERIC_CONTROL:
        return 'control'
    if v in _GENERIC_DISEASE:
        return 'disease'
    return ''


def _pick_condition_column(obs_columns: list[str]) -> Optional[str]:
    """Pick the most likely condition column from a study's obs."""
    for c in _PREFERRED_CONDITION_COLS:
        if c in obs_columns:
            return c
    return None


def _read_obs_metadata(h5ad_path: Path) -> dict:
    """Read obs structure via h5py (no full DataFrame load).

    Categorical columns in h5ad store their unique values as a small
    'categories' string dataset alongside an int 'codes' dataset; we
    only need 'categories' to populate the role-mapping dropdowns.
    Skipping non-categorical columns avoids reading huge per-cell
    arrays just to discover values.
    """
    import h5py
    from anndata.io import read_elem

    columns: list[str] = []
    values: dict[str, list[str]] = {}
    n_obs = 0
    with h5py.File(h5ad_path, 'r') as h:
        if 'obs' not in h:
            return {'columns': [], 'values': {}, 'n_obs': 0}
        obs_group = h['obs']
        # n_obs: length of the index. read_elem handles whatever encoding
        # it was written with -- a raw '.shape[0]' breaks on pandas' now
        # common nullable-string encoding, which stores it as a Group.
        idx_key = obs_group.attrs.get('_index')
        if isinstance(idx_key, bytes):
            idx_key = idx_key.decode('utf-8')
        if idx_key and idx_key in obs_group:
            n_obs = len(read_elem(obs_group[idx_key]))

        for col, item in obs_group.items():
            if col == idx_key:
                continue
            columns.append(col)
            # Categorical: group containing 'categories' and 'codes'.
            if isinstance(item, h5py.Group) and 'categories' in item:
                try:
                    cats = item['categories'][()]
                    cats_list = [
                        c.decode('utf-8') if isinstance(c, bytes) else str(c)
                        for c in cats
                    ]
                    if 1 < len(cats_list) <= 50:
                        values[col] = sorted(cats_list)
                except (OSError, KeyError):
                    pass
            # Non-categorical small string datasets: bail; cost of
            # uniqueness across millions of cells isn't worth it for
            # this UI (the role-mapping uses categorical cols anyway).

    return {'columns': columns, 'values': values, 'n_obs': n_obs}


def _pick_source_h5ad(study_dir: Path) -> Optional[Path]:
    """Pick the h5ad a combine should read for a study.

    Preference order:
      1. The newest .h5ad in processed_data/ -- the study as KOSMIC has
         prepared it: harmonised gene names, roles, QC, annotation.
      2. The newest .h5ad in raw_data/, for a study nothing has been done
         to yet.

    This used to prefer raw_data/, reasoning that a combined object wants
    raw counts. It does, but 'raw counts' is not the same as 'the file as
    downloaded': processed_data/ holds counts too, and 'concat_studies'
    promotes '.raw' to '.X' regardless. Reading raw_data/ silently threw
    away every step of the pipeline -- and for a CellxGene deposit it
    reads the Ensembl-indexed original, which shares no gene names with
    the others, so the inner join collapses.

    The tie-break is mtime, not size: chaffin's raw download is 13 GB and
    its processed file 1.9 GB, and the small one is the right one.
    """
    for folder_fn in (processed_data_dir, raw_data_dir):
        folder = folder_fn(study_dir)
        if folder.is_dir():
            h5ads = list(folder.glob("*.h5ad"))
            if h5ads:
                return max(h5ads, key=lambda p: p.stat().st_mtime)
    return None


# Worker

class _ScanWorker(BaseWorker):
    """Background scanner that reads obs metadata for each study via h5py.

    Emits per-study progress; result is a list of
    (accession, source_path, obs_metadata_dict).
    """

    def __init__(self, project_dir: Path, parent=None):
        super().__init__(parent)
        self._project_dir = Path(project_dir)

    def _run(self):
        results = []
        accessions = [a for a in list_studies(self._project_dir)
                      if a != MASTER_ACCESSION]
        total = len(accessions)
        for i, accession in enumerate(accessions, 1):
            self.progress.emit(f"Scanning {accession} ({i}/{total})...")
            self.progress_pct.emit(int(100 * i / max(1, total)))
            study_path = self._project_dir / accession
            source = _pick_source_h5ad(study_path)
            if source is None:
                results.append((accession, None, None))
                continue
            try:
                meta = _read_obs_metadata(source)
            except (OSError, KeyError, ValueError):
                meta = {'columns': [], 'values': {}, 'n_obs': 0}
            # Gene names and the harmonisation record travel with the scan so
            # the merge-asymmetry check can be recomputed instantly whenever
            # the study selection changes, without touching disk again.
            try:
                meta['genes'] = read_var_names(source)
            except (OSError, KeyError, ValueError):
                meta['genes'] = []
            meta['merges'] = merge_map_from_provenance(source.parent)
            results.append((accession, source, meta))
        return results


class _CombineWorker(BaseWorker):
    """Builds the master h5ad on a worker thread."""

    def __init__(self, study_paths, cap_per_study, output_path,
                 role_map=None, drop_genes=None, drop_excluded=False,
                 parent=None):
        super().__init__(parent)
        self._study_paths = list(study_paths)
        self._cap = cap_per_study
        self._output_path = Path(output_path)
        self._role_map = role_map
        self._drop_genes = list(drop_genes or ())
        self._drop_excluded = bool(drop_excluded)

    def _run(self):
        import anndata as ad
        concat_studies(
            self._study_paths,
            output_path=self._output_path,
            cap_per_study=self._cap,
            role_map=self._role_map,
            drop_genes=self._drop_genes,
            drop_excluded=self._drop_excluded,
            progress_callback=lambda msg: self.progress.emit(msg),
        )
        # Quick peek for the result payload; backed mode = cheap.
        a = ad.read_h5ad(self._output_path, backed='r')
        try:
            n_cells, n_genes = int(a.n_obs), int(a.n_vars)
        finally:
            a.file.close()
        # Record what the master is, and in particular what was excluded
        # from it -- a dropped gene leaves no trace in the h5ad, so without
        # this there is no way to tell it was ever considered.
        try:
            from kosmic import provenance
            provenance.record_stage(
                self._output_path.parent, self._output_path.stem, 'combine',
                {
                    'studies': [Path(p).parent.parent.name
                                for p in self._study_paths],
                    'cap_per_study': self._cap,
                    'n_cells': n_cells,
                    'n_genes': n_genes,
                    'dropped_genes': self._drop_genes,
                    'dropped_reason': ('merged inconsistently across studies'
                                       if self._drop_genes else None),
                    'excluded_cells_dropped': self._drop_excluded,
                })
        except Exception as exc:            # provenance never blocks a build
            self.progress.emit(f"Could not record provenance: {exc}")

        return {
            'path': str(self._output_path),
            'n_cells': n_cells,
            'n_genes': n_genes,
            'dropped_genes': self._drop_genes,
        }


class _PropagateWorker(BaseWorker):
    """Projects master labels onto each full per-study h5ad."""

    def __init__(self, master_path, study_paths, label_cols, parent=None):
        super().__init__(parent)
        self._master_path = Path(master_path)
        self._study_paths = [Path(p) for p in study_paths]
        self._label_cols = list(label_cols)

    def _run(self):
        results = []
        for i, study_path in enumerate(self._study_paths, 1):
            self.progress.emit(
                f"[{i}/{len(self._study_paths)}] {study_path.parent.parent.name}")
            n = propagate_labels(
                self._master_path,
                study_path,
                label_cols=self._label_cols,
                progress_callback=lambda msg: self.progress.emit(msg),
            )
            results.append((study_path.parent.parent.name, n))
            self.progress_pct.emit(int(100 * i / len(self._study_paths)))
        return results


# Dialog

class CombineDialog(QDialog):
    """Modal that builds the project's master h5ad from selected studies."""

    help_id = "project/atlas"

    log_message = pyqtSignal(str)
    master_built = pyqtSignal(str)         # path to written master.h5ad
    open_master_requested = pyqtSignal()   # "Open in scRNA" pressed

    _STATE_CONFIGURE = 0
    _STATE_RUNNING = 1

    def __init__(self, main_window, parent=None):
        super().__init__(parent or main_window)
        self.setWindowTitle("Combine studies")
        self.setModal(True)
        self.resize(900, 800)
        self.main_window = main_window
        self._worker: Optional[BaseWorker] = None
        self._scan_worker: Optional[BaseWorker] = None
        # {accession: (gene_names, recorded merge map or None)} -- cached from
        # the scan so the asymmetry check is pure set arithmetic thereafter.
        self._gene_state: dict[str, tuple] = {}
        self._merge_asymmetry: dict = {"genes": {}, "unknown": []}
        self._build_ui()
        # Show the dialog chrome immediately; scan obs in a worker thread.
        self._start_scan()

    def _build_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        self._stack = QStackedWidget(self)
        self._stack.addWidget(self._build_configure_page())
        self._stack.addWidget(self._build_running_page())
        outer.addWidget(self._stack)

    # ---- Configure page --------------------------------------------------

    def _build_configure_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(32, 24, 32, 24)
        layout.setSpacing(12)

        layout.addWidget(HeaderLabel("Combine"))
        layout.addWidget(SecondaryLabel(
            "Concatenate the processed h5ads of selected studies into a "
            "single master dataset. The master can then be opened in the "
            "scRNA workspace for cross-study clustering and annotation."
        ))

        # Study table
        layout.addWidget(SecondaryLabel("Studies"))
        self._study_table = QTableWidget(0, 6)
        self._study_table.setHorizontalHeaderLabels(
            ["", "Study", "Source", "File", "Cells", "Will use"])
        self._study_table.setSelectionMode(
            QAbstractItemView.SelectionMode.NoSelection)
        self._study_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self._study_table.verticalHeader().setVisible(False)
        self._study_table.verticalHeader().setDefaultSectionSize(34)
        hdr = self._study_table.horizontalHeader()
        hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.Fixed)
        hdr.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        hdr.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        hdr.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        hdr.setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)
        hdr.setSectionResizeMode(5, QHeaderView.ResizeMode.ResizeToContents)
        self._study_table.setColumnWidth(0, 30)
        # ~5 rows: the variable-length gene-overlap warning below otherwise
        # squeezes this to 1-2 rows in the stretch layout.
        self._study_table.setMinimumHeight(220)
        layout.addWidget(self._study_table, 1)

        # Cap spinner
        cap_row = QHBoxLayout()
        cap_row.addWidget(QLabel("Cells per study cap:"))
        self._cap_spin = QSpinBox()
        self._cap_spin.setRange(1_000, 10_000_000)
        self._cap_spin.setSingleStep(10_000)
        self._cap_spin.setValue(50_000)
        self._cap_spin.setSuffix(" cells")
        self._cap_spin.setEnabled(False)
        self._cap_spin.valueChanged.connect(self._refresh_estimates)
        cap_row.addWidget(self._cap_spin)

        self._uncap_check = QCheckBox("Use all cells (no cap)")
        self._uncap_check.setChecked(True)
        self._uncap_check.toggled.connect(self._on_uncap_toggled)
        cap_row.addWidget(self._uncap_check)
        cap_row.addStretch()
        layout.addLayout(cap_row)

        self._estimate_label = SecondaryLabel("")
        layout.addWidget(self._estimate_label)

        # Harmonisation sums two columns into one where HGNC has merged the
        # loci, but only in studies whose annotation carried both features.
        # The result is a gene that is a sum of two features in one cohort
        # and a single feature in another -- a constant per-study offset that
        # inflates between-study heterogeneity. Dropping them is the default,
        # but the count and the list are always in view so the choice is
        # visible rather than silent.
        drop_row = QHBoxLayout()
        self._drop_merged_check = QCheckBox(
            "Drop genes merged inconsistently across studies")
        self._drop_merged_check.setChecked(True)
        self._drop_merged_check.setEnabled(False)
        self._drop_merged_check.toggled.connect(self._refresh_estimates)
        drop_row.addWidget(self._drop_merged_check)
        self._drop_details_btn = SecondaryButton("Details...")
        self._drop_details_btn.setEnabled(False)
        self._drop_details_btn.clicked.connect(self._show_merge_details)
        drop_row.addWidget(self._drop_details_btn)
        drop_row.addStretch()
        layout.addLayout(drop_row)

        self._drop_excluded_check = QCheckBox(
            "Leave out cells marked 'exclude'")
        self._drop_excluded_check.setChecked(True)
        self._drop_excluded_check.setToolTip(
            "Skip cells whose condition you set to Exclude in Inspect.\n"
            "They never enter a contrast, so carrying them into the\n"
            "master means clustering, annotating and storing them for\n"
            "nothing." + "\n\n"
            "The source studies are untouched either way -- this only\n"
            "decides what goes into the master.")
        self._drop_excluded_check.toggled.connect(self._refresh_estimates)
        layout.addWidget(self._drop_excluded_check)

        self._drop_label = SecondaryLabel("")
        self._drop_label.setWordWrap(True)
        layout.addWidget(self._drop_label)

        # Role mapping table
        layout.addWidget(SecondaryLabel(
            "Condition -> role mapping (writes a unified '_role' obs column "
            "on the master; protects disease signal during Harmony)"))
        self._role_table = QTableWidget(0, 5)
        self._role_table.setHorizontalHeaderLabels(
            ["Study", "Column", "Value", "Cells", "Role"])
        self._role_table.setSelectionMode(
            QAbstractItemView.SelectionMode.NoSelection)
        self._role_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self._role_table.verticalHeader().setVisible(False)
        self._role_table.verticalHeader().setDefaultSectionSize(34)
        rhdr = self._role_table.horizontalHeader()
        rhdr.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        rhdr.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        rhdr.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        rhdr.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        rhdr.setSectionResizeMode(4, QHeaderView.ResizeMode.Fixed)
        self._role_table.setColumnWidth(4, 110)
        self._role_table.setMaximumHeight(220)
        layout.addWidget(self._role_table)

        self._role_status = SecondaryLabel("")
        layout.addWidget(self._role_status)

        # Action row
        action_row = QHBoxLayout()
        self._refresh_btn = SecondaryButton("Refresh")
        self._refresh_btn.clicked.connect(self.refresh)
        action_row.addWidget(self._refresh_btn)

        action_row.addStretch()

        self._cancel_btn = SecondaryButton("Cancel")
        self._cancel_btn.clicked.connect(self.reject)
        action_row.addWidget(self._cancel_btn)

        self._run_btn = PrimaryButton("Build Master")
        self._run_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._run_btn.clicked.connect(self._on_run_clicked)
        action_row.addWidget(self._run_btn)

        layout.addLayout(action_row)
        return page

    # ---- Running page ----------------------------------------------------

    def _build_running_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(32, 24, 32, 24)
        layout.setSpacing(12)

        layout.addWidget(HeaderLabel("Running"))
        self._status_label = SecondaryLabel("")
        layout.addWidget(self._status_label)

        self._progress = QProgressBar()
        self._progress.setRange(0, 0)  # indeterminate
        layout.addWidget(self._progress)

        self._log = QPlainTextEdit()
        self._log.setReadOnly(True)
        layout.addWidget(self._log, 1)

        action_row = QHBoxLayout()
        action_row.addStretch()

        self._close_btn = SecondaryButton("Close")
        self._close_btn.setEnabled(False)
        self._close_btn.clicked.connect(self.accept)
        action_row.addWidget(self._close_btn)

        self._open_master_btn = PrimaryButton("Open master in scRNA")
        self._open_master_btn.setEnabled(False)
        self._open_master_btn.clicked.connect(self._emit_open_and_close)
        action_row.addWidget(self._open_master_btn)

        layout.addLayout(action_row)
        return page

    def _emit_open_and_close(self) -> None:
        self.open_master_requested.emit()
        self.accept()

    # ---- Public API ------------------------------------------------------

    def _start_scan(self) -> None:
        """Kick off the obs-metadata scan in a worker thread."""
        self._stack.setCurrentIndex(self._STATE_CONFIGURE)
        self._study_table.setRowCount(0)
        self._study_sources: dict[str, Path] = {}
        self._study_obs: dict[str, dict] = {}
        self._role_state: dict[str, dict] = {}
        self._run_btn.setEnabled(False)
        project_dir = self.main_window.current_project_dir
        if project_dir is None:
            self._refresh_role_table()
            self._estimate_label.setText("No project open.")
            return

        self._estimate_label.setText("Scanning studies...")
        self._role_status.setText("")

        self._scan_worker = _ScanWorker(project_dir, parent=self)
        run_worker(
            self._scan_worker,
            on_finished=self._on_scan_done,
            on_failed=self._on_scan_failed,
            on_progress=lambda msg: self._estimate_label.setText(msg),
        )

    def refresh(self) -> None:
        """Public hook so callers (or the dialog itself) can restart the scan."""
        self._start_scan()

    def _on_scan_failed(self, message: str) -> None:
        self._estimate_label.setText(f"Scan failed: {message}")

    def _on_scan_done(self, results) -> None:
        """Populate the study table + role state once the scan finishes."""
        for accession, source, meta in results:
            usable = source is not None
            if usable:
                self._study_sources[accession] = source
                self._study_obs[accession] = meta or {
                    'columns': [], 'values': {}, 'n_obs': 0}
                cols = self._study_obs[accession]['columns']
                values = self._study_obs[accession]['values']
                n_cells = self._study_obs[accession]['n_obs']
                col = _pick_condition_column(cols)
                mapping: dict[str, str] = {}
                if col is not None:
                    is_role_col = (col == '_role')
                    for v in values.get(col, []):
                        mapping[v] = _default_role_for(
                            v, existing_role_col_present=is_role_col,
                            existing_value=v if is_role_col else None)
                self._role_state[accession] = {'col': col, 'mapping': mapping}
                self._gene_state[accession] = (
                    self._study_obs[accession].get('genes') or [],
                    self._study_obs[accession].get('merges'))
                source_folder = source.parent.name
                file_name = source.name
            else:
                source_folder = "(none)"
                file_name = ""
                n_cells = 0

            row = self._study_table.rowCount()
            self._study_table.insertRow(row)
            cb = QCheckBox()
            cb.setChecked(usable)
            cb.setEnabled(usable)
            cb.stateChanged.connect(self._refresh_estimates)
            cb.stateChanged.connect(self._refresh_role_table)
            self._study_table.setCellWidget(row, 0, cb)
            name_item = QTableWidgetItem(accession)
            if not usable:
                name_item.setToolTip(
                    "No .h5ad found in raw_data/ or processed_data/")
            self._study_table.setItem(row, 1, name_item)
            self._study_table.setItem(row, 2, QTableWidgetItem(source_folder))
            self._study_table.setItem(row, 3, QTableWidgetItem(file_name))
            self._study_table.setItem(
                row, 4,
                QTableWidgetItem(f"{n_cells:,}" if n_cells else "-"))
            self._study_table.setItem(row, 5, QTableWidgetItem(""))

        self._refresh_estimates()
        self._refresh_role_table()
        self._run_btn.setEnabled(bool(results))

    def _selected_accessions(self) -> list[str]:
        out = []
        for row in range(self._study_table.rowCount()):
            cb = self._study_table.cellWidget(row, 0)
            if cb and cb.isChecked():
                out.append(self._study_table.item(row, 1).text())
        return out

    def _refresh_role_table(self) -> None:
        """Rebuild the role-mapping table from the currently-selected studies."""
        self._role_table.setRowCount(0)
        unmapped = 0
        total = 0
        for accession in self._selected_accessions():
            meta = self._study_obs.get(accession, {})
            state = self._role_state.get(accession)
            if not state or not meta:
                continue
            col = state['col']
            mapping = state['mapping']
            values = meta.get('values', {})

            # Column-selector dropdown shows all short-categorical columns,
            # so user can override the auto-pick.
            col_options = [c for c in meta['columns'] if c in values]
            if col is not None and col not in col_options:
                col_options = [col] + col_options

            if not col or not mapping:
                # Show one placeholder row so user can pick a column.
                row = self._role_table.rowCount()
                self._role_table.insertRow(row)
                self._role_table.setItem(row, 0, QTableWidgetItem(accession))
                self._role_table.setCellWidget(
                    row, 1, self._make_col_combo(accession, col_options, col))
                self._role_table.setItem(
                    row, 2, QTableWidgetItem("(pick a column)"))
                self._role_table.setItem(row, 3, QTableWidgetItem(""))
                self._role_table.setItem(row, 4, QTableWidgetItem(""))
                unmapped += 1
                total += 1
                continue

            for v in values.get(col, []):
                row = self._role_table.rowCount()
                self._role_table.insertRow(row)
                self._role_table.setItem(row, 0, QTableWidgetItem(accession))
                self._role_table.setCellWidget(
                    row, 1, self._make_col_combo(accession, col_options, col))
                self._role_table.setItem(row, 2, QTableWidgetItem(v))
                self._role_table.setItem(row, 3, QTableWidgetItem(""))
                self._role_table.setCellWidget(
                    row, 4,
                    self._make_role_combo(accession, v, mapping.get(v, '')))
                if not mapping.get(v):
                    unmapped += 1
                total += 1

        if total == 0:
            self._role_status.setText(
                "No studies selected with discoverable condition values.")
        elif unmapped:
            self._role_status.setText(
                f"{unmapped} of {total} value(s) unmapped. Set a role for "
                f"every row before building.")
        else:
            self._role_status.setText(
                f"All {total} condition value(s) mapped.")

    def _make_col_combo(self, accession: str, options: list[str],
                        current: Optional[str]) -> QComboBox:
        combo = QComboBox()
        combo.addItems(options)
        if current is not None and current in options:
            combo.setCurrentText(current)
        combo.currentTextChanged.connect(
            lambda new_col, a=accession: self._on_col_changed(a, new_col))
        return combo

    def _make_role_combo(self, accession: str, value: str,
                         current: str) -> QComboBox:
        combo = QComboBox()
        combo.addItems(_ROLE_OPTIONS)
        combo.setCurrentText(current or '')
        combo.currentTextChanged.connect(
            lambda new_role, a=accession, v=value:
            self._on_role_changed(a, v, new_role))
        return combo

    def _on_col_changed(self, accession: str, new_col: str) -> None:
        state = self._role_state.get(accession)
        meta = self._study_obs.get(accession, {})
        if not state or not meta:
            return
        state['col'] = new_col
        is_role_col = (new_col == '_role')
        new_mapping: dict[str, str] = {}
        for v in meta.get('values', {}).get(new_col, []):
            new_mapping[v] = _default_role_for(
                v, existing_role_col_present=is_role_col,
                existing_value=v if is_role_col else None)
        state['mapping'] = new_mapping
        self._refresh_role_table()

    def _on_role_changed(self, accession: str, value: str,
                         new_role: str) -> None:
        state = self._role_state.get(accession)
        if not state:
            return
        state['mapping'][value] = new_role
        # Re-evaluate the unmapped count without rebuilding all widgets.
        unmapped = 0
        total = 0
        for accession in self._selected_accessions():
            s = self._role_state.get(accession)
            if not s or not s['col']:
                continue
            for v, r in s['mapping'].items():
                total += 1
                if not r:
                    unmapped += 1
        if total == 0:
            return
        if unmapped:
            self._role_status.setText(
                f"{unmapped} of {total} value(s) unmapped. Set a role for "
                f"every row before building.")
        else:
            self._role_status.setText(
                f"All {total} condition value(s) mapped.")

    # ---- Internals -------------------------------------------------------

    def _on_uncap_toggled(self, on: bool) -> None:
        self._cap_spin.setEnabled(not on)
        self._refresh_estimates()

    def _refresh_estimates(self) -> None:
        cap = None if self._uncap_check.isChecked() else self._cap_spin.value()
        total_master_cells = 0
        n_selected = 0
        for row in range(self._study_table.rowCount()):
            cb = self._study_table.cellWidget(row, 0)
            n_cells_text = self._study_table.item(row, 4).text().replace(',', '')
            try:
                n_cells = int(n_cells_text)
            except ValueError:
                n_cells = 0
            if cb and cb.isChecked():
                will_use = n_cells if cap is None else min(cap, n_cells)
                total_master_cells += will_use
                n_selected += 1
                self._study_table.item(row, 5).setText(f"{will_use:,}")
            else:
                self._study_table.item(row, 5).setText("")

        cap_note = "all cells" if cap is None else f"cap {cap:,}"
        self._estimate_label.setText(
            f"Master will contain ~{total_master_cells:,} cells from "
            f"{n_selected} stud{'y' if n_selected == 1 else 'ies'} ({cap_note}).")

        self._refresh_merge_asymmetry()

    def _refresh_merge_asymmetry(self) -> None:
        """Recheck which genes the selected studies merged inconsistently.

        Depends on the selection: deselecting the one study that is out of
        step removes the asymmetry entirely, so this is recomputed rather
        than cached against the project.
        """
        selected = {a: self._gene_state[a] for a in self._selected_accessions()
                    if a in self._gene_state}
        self._merge_asymmetry = find_merge_asymmetry(selected)
        genes = self._merge_asymmetry["genes"]
        unknown = self._merge_asymmetry["unknown"]

        has_genes = bool(genes)
        self._drop_merged_check.setEnabled(has_genes)
        self._drop_details_btn.setEnabled(has_genes)

        parts = []
        if has_genes:
            parts.append(
                f"{len(genes):,} gene(s) are a sum of several features in "
                f"some of these studies and a single feature in others, "
                f"because HGNC has merged those loci. Within a study that is "
                f"harmless; across studies it shifts the gene's counts by a "
                f"constant amount per cohort.")
        if unknown:
            parts.append(
                f"Gene names have not been harmonised in: "
                f"{', '.join(unknown)}. Those studies cannot be checked, so "
                f"the count above may be incomplete.")
        self._drop_label.setText(" ".join(parts))
        self._drop_label.setVisible(bool(parts))

    def _show_merge_details(self) -> None:
        """List the affected genes and how many features each study summed."""
        genes = self._merge_asymmetry["genes"]
        if not genes:
            return
        studies = sorted(self._selected_accessions())
        header = "gene".ljust(18) + "  ".join(s[:12].rjust(12) for s in studies)
        rows = [header, "-" * len(header)]
        for gene in sorted(genes):
            counts = genes[gene]
            rows.append(gene.ljust(18) + "  ".join(
                str(counts.get(s, "-")).rjust(12) for s in studies))
        from kosmic.gui.shared import dialogs
        dialogs.info(
            self, "Genes merged inconsistently",
            f"{len(genes):,} gene(s) affected. The number under each study "
            f"is how many source features were summed into that gene."
            + "\n\n"
            + "\n".join(rows))

    def _selected_paths(self) -> list[Path]:
        """Return the resolved source h5ad path per selected study."""
        paths = []
        for row in range(self._study_table.rowCount()):
            cb = self._study_table.cellWidget(row, 0)
            if cb and cb.isChecked():
                accession = self._study_table.item(row, 1).text()
                source = self._study_sources.get(accession)
                if source is not None:
                    paths.append(source)
        return paths

    def _on_run_clicked(self) -> None:
        paths = self._selected_paths()
        if not paths:
            return

        # Assemble role_map and validate every value is mapped.
        role_map: dict[str, tuple[str, dict[str, str]]] = {}
        missing: list[str] = []
        for accession in self._selected_accessions():
            state = self._role_state.get(accession)
            if not state or not state.get('col'):
                missing.append(f"{accession}: no condition column picked")
                continue
            mapping = state['mapping']
            blanks = [v for v, r in mapping.items() if not r]
            if blanks:
                missing.append(
                    f"{accession}: unmapped value(s) -- {', '.join(blanks)}")
                continue
            role_map[accession] = (state['col'], dict(mapping))

        if missing:
            from kosmic.gui.shared import dialogs
            dialogs.warning(self, "Role mapping incomplete",
                            "\n".join(missing))
            return

        cap = None if self._uncap_check.isChecked() else self._cap_spin.value()
        project_dir = self.main_window.current_project_dir
        output = master_h5ad_path(project_dir)

        self._status_label.setText("Building master h5ad...")
        self._log.clear()
        self._progress.setRange(0, 0)
        self._stack.setCurrentIndex(self._STATE_RUNNING)
        self._open_master_btn.setEnabled(False)
        self._close_btn.setEnabled(False)

        drop_genes = (sorted(self._merge_asymmetry["genes"])
                      if self._drop_merged_check.isChecked() else [])
        if drop_genes:
            self._append_log(
                f"Dropping {len(drop_genes):,} gene(s) merged inconsistently "
                f"across studies.")

        self._worker = _CombineWorker(
            paths, cap, output, role_map=role_map, drop_genes=drop_genes,
            drop_excluded=self._drop_excluded_check.isChecked(), parent=self)
        run_worker(
            self._worker,
            on_finished=self._on_build_finished,
            on_failed=self._on_failed,
            on_progress=self._append_log,
        )

    def _on_build_finished(self, result: dict) -> None:
        self._progress.setRange(0, 100)
        self._progress.setValue(100)
        dropped = result.get('dropped_genes') or []
        note = (f", {len(dropped):,} inconsistently-merged gene(s) excluded"
                if dropped else "")
        self._status_label.setText(
            f"Master built: {result['n_cells']:,} cells x "
            f"{result['n_genes']:,} genes{note}")
        self._open_master_btn.setEnabled(True)
        self._close_btn.setEnabled(True)
        self.master_built.emit(result['path'])

    def _on_failed(self, message: str) -> None:
        self._progress.setRange(0, 100)
        self._progress.setValue(0)
        self._status_label.setText("Failed.")
        self._append_log(f"ERROR: {message}")
        self._close_btn.setEnabled(True)

    def _append_log(self, msg: str) -> None:
        self._log.appendPlainText(msg)
        self.log_message.emit(msg)
