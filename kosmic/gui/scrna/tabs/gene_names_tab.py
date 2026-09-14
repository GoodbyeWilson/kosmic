# Gene Names Tab
# Harmonises gene names to current HGNC-approved symbols.
# Shows a summary and detailed mapping table so the user can review
# before accepting.

from PyQt6.QtWidgets import (
    QHBoxLayout, QLabel, QFrame, QGridLayout, QPushButton,
    QTableWidget, QTableWidgetItem, QHeaderView, QLineEdit,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QCursor, QFont

import pandas as pd

from kosmic.gui.shared.widgets import (
    BaseWorker, CaptionLabel, Column, PrimaryButton, ResultsTable,
    SecondaryButton, SecondaryLabel, SectionHeader, SidebarPage,
)
from kosmic.gui.shared import dialogs, run_worker


# ---------------------------------------------------------------------------
# Worker thread
# ---------------------------------------------------------------------------

class _ApplyWorker(BaseWorker):
    """
    Apply harmonisation to adata and save, in a background thread.

    Emits 'finished_ok' with the harmonised adata, or 'None' if the
    harmonisation or save step failed (errors are logged via 'progress').
    """

    def __init__(self, adata, lookup, reasons, save_path=None, parent=None):
        super().__init__(parent)
        self._adata = adata
        self._lookup = lookup
        self._reasons = reasons
        self._save_path = save_path

    def _run(self):
        try:
            from kosmic.scrna.inspect.gene_names import harmonise_adata
            adata_new, _ = harmonise_adata(
                self._adata, lookup=self._lookup, reasons=self._reasons,
                progress_callback=lambda msg: self.progress.emit(msg),
            )
            # Make sure we have a concrete AnnData (not a view) before writing.
            try:
                if getattr(adata_new, 'is_view', False):
                    adata_new = adata_new.copy()
            except Exception:
                pass
            if self._save_path:
                self.progress.emit(f"Saving h5ad to {self._save_path} ...")
                try:
                    adata_new.write_h5ad(self._save_path)
                    import os
                    size = os.path.getsize(self._save_path)
                    self.progress.emit(f"Saved {size/1e6:.1f} MB to {self._save_path}")
                except Exception as save_err:
                    import traceback
                    self.progress.emit(f"SAVE FAILED ({self._save_path}): {save_err}")
                    self.progress.emit(traceback.format_exc())
                    raise
            else:
                self.progress.emit("WARNING: no save path set — harmonisation is in-memory only")
            self.progress.emit("Done")
            return adata_new
        except Exception as e:
            import traceback
            self.progress.emit(f"Error: {e}")
            self.progress.emit(traceback.format_exc())
            return None


class _HarmoniseWorker(BaseWorker):
    """
    Analyse gene names against HGNC — lightweight, no adata copy.

    Emits 'finished_ok' with a tuple '(report, (lookup, reasons))' on
    success or '(None, error_str)' on failure so the caller's single
    handler can branch on 'report is None'.
    """

    def __init__(self, var_names, parent=None):
        super().__init__(parent)
        self._var_names = var_names

    def _run(self):
        try:
            from kosmic.scrna.inspect.gene_names import (
                HarmonisationReport, count_locus_groups, load_hgnc_lookup,
                load_locus_groups,
            )
            self.progress.emit("Loading HGNC symbol table...")
            lookup, reasons = load_hgnc_lookup()
            locus_groups = load_locus_groups()

            self.progress.emit("Matching gene names...")
            # Build report without touching adata
            var_names = self._var_names
            report = HarmonisationReport(total_genes=len(var_names))
            mappings = []
            for name in var_names:
                if name in lookup:
                    mappings.append((name, lookup[name], reasons.get(name, 'alias')))
            report.renamed = len(mappings)
            report.mappings = mappings
            report.unchanged = report.total_genes - report.renamed

            # Check for duplicates
            from collections import Counter
            new_names = [lookup.get(n, n) for n in var_names]
            dupes = {n: c for n, c in Counter(new_names).items() if c > 1}
            merge_info = []
            if dupes:
                name_to_positions = {}
                for i, n in enumerate(new_names):
                    name_to_positions.setdefault(n, []).append(i)
                for n, positions in name_to_positions.items():
                    if len(positions) > 1:
                        merge_info.append((n, [var_names[p] for p in positions]))
            report.merged_genes = merge_info
            report.duplicates_merged = len(merge_info)

            self.progress.emit("Classifying gene names...")
            report.locus_before = count_locus_groups(
                var_names, locus_groups=locus_groups)
            report.locus_after = count_locus_groups(
                var_names, locus_groups=locus_groups, lookup=lookup)

            # The sidebar counts these; without the names there is no way to
            # see *which* genes match nothing in HGNC, and those are exactly
            # the ones that vanish from a cross-study join.
            report.unrecognised = sorted(
                {lookup.get(n, n) for n in var_names} - set(locus_groups))

            return (report, (lookup, reasons, locus_groups))
        except Exception as e:
            return (None, str(e))


# ---------------------------------------------------------------------------
# Gene Names Tab
# ---------------------------------------------------------------------------

class GeneNamesTab(SidebarPage):
    """Tab for reviewing and applying HGNC gene name harmonisation."""

    help_id = "scrna/gene_names"
    CONTENT_MARGINS = (16, 16, 16, 16)
    CONTENT_SPACING = 12

    log_message = pyqtSignal(str)

    #: Reason key -> table label. 'merged' marks a column that keeps its
    #: name but has another column's counts added into it.
    _REASON_LABELS = {
        'previous_symbol': 'Previous symbol',
        'alias': 'Alias',
        'mitochondrial': 'Mitochondrial',
        'merged': 'Merged',
        'unrecognised': 'Unrecognised',
    }

    #: Chip key -> base label; the row count is appended at runtime.
    _CHIP_LABELS = {
        None: 'All',
        'previous_symbol': 'Previous symbol',
        'alias': 'Alias',
        'mitochondrial': 'Mitochondrial',
        '_merged': 'Merged',
        '_unrecognised': 'Unrecognised',
    }


    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self.progress_bar = None  # wired to sidebar progress bar by AppWindow
        self._worker = None
        self._apply_worker = None
        self._pending_adata = None
        self._report = None
        self._lookup_data = None
        self._locus_groups = None
        self._harmonised = False
        self._applying = False  # guard against re-entrant set_adata
        self._last_adata_id = None  # track which adata we've already analysed
        self._setup_ui()

    def reset_state(self):
        """Clear all cached analysis state (called on raw data reload)."""
        # Kill running workers
        for w in (self._worker, self._apply_worker):
            if w is not None and w.isRunning():
                w.quit()
                w.wait(2000)
        self._worker = None
        self._apply_worker = None
        self._pending_adata = None
        self._report = None
        self._lookup_data = None
        self._locus_groups = None
        self._harmonised = False
        self._applying = False
        self._last_adata_id = None
        # Reset UI
        self._table.setRowCount(0)
        self._apply_btn.setEnabled(False)
        self._set_status("")
        self._update_stats()

    def _setup_ui(self):
        # The QC tab uses an accordion because it has four stages to page
        # between. Harmonisation is one step, so there is nothing to collapse
        # for: the stats sit in a plain panel that is always open, with the
        # actions directly beneath them.
        side = self.sidebar_layout

        side.addWidget(SectionHeader("GENE NAMES"))

        stats_box = QFrame()
        stats_box.setObjectName("stage_card")
        stats_grid = QGridLayout(stats_box)
        stats_grid.setContentsMargins(12, 10, 12, 10)
        stats_grid.setHorizontalSpacing(8)
        stats_grid.setVerticalSpacing(4)
        stats_grid.setColumnStretch(0, 1)

        # 'Unrecognised' is the row that predicts trouble: a retired
        # clone-based name (RP11-34P13.3) or a bare Ensembl ID matches
        # nothing in another study, so those genes drop out of any
        # cross-study join.
        self._stat_rows = {}
        rows = (('total', 'Total genes'), ('renamed', 'To rename'),
                ('merged', 'To merge'), ('protein_coding', 'Protein-coding'),
                ('unrecognised', 'Unrecognised'))
        for row, (key, caption) in enumerate(rows):
            name = SecondaryLabel(caption)
            value = QLabel("—")
            value.setFont(QFont("Segoe UI", 10, QFont.Weight.Bold))
            value.setAlignment(Qt.AlignmentFlag.AlignRight
                               | Qt.AlignmentFlag.AlignVCenter)
            stats_grid.addWidget(name, row, 0)
            stats_grid.addWidget(value, row, 1)
            self._stat_rows[key] = value
        side.addWidget(stats_box)

        self._apply_btn = PrimaryButton("Apply Harmonisation")
        self._apply_btn.clicked.connect(self._on_apply)
        self._apply_btn.setEnabled(False)
        side.addWidget(self._apply_btn)

        self._skip_btn = SecondaryButton("Skip (keep original names)")
        self._skip_btn.clicked.connect(self._on_skip)
        self._skip_btn.setEnabled(False)
        side.addWidget(self._skip_btn)

        # Genes that harmonisation would build by summing two or more source
        # columns. Within this study the sum is correct -- HGNC says the loci
        # are one gene -- but it is also what later makes the gene mean
        # something different here than in a study whose annotation already
        # had them merged, so it is called out rather than left buried in the
        # mapping table.
        self._merge_note = SecondaryLabel("")
        self._merge_note.setWordWrap(True)
        self._merge_note.hide()
        side.addWidget(self._merge_note)

        self._status_label = SecondaryLabel("")
        self._status_label.setWordWrap(True)
        self._status_label.setProperty("role", "padded_status")
        self._status_label.hide()
        side.addWidget(self._status_label)

        side.addStretch()

        self._update_btn = SecondaryButton("Update HGNC Table")
        self._update_btn.setToolTip(
            "Download the latest HGNC symbol table from genenames.org")
        self._update_btn.clicked.connect(self._on_update_hgnc)
        side.addWidget(self._update_btn)


        # ---- Content: description, filters, mapping table ----
        # The counts live in the sidebar; repeating them here as a strip of
        # big numbers said the same thing twice, so the content pane is just
        # the detail the sidebar summarises.
        layout = self.content_layout

        header = QLabel("Gene Name Harmonisation")
        header.setFont(QFont("Segoe UI", 14, QFont.Weight.Bold))
        layout.addWidget(header)

        desc = SecondaryLabel(
            "Maps gene names to current HGNC-approved symbols, so that a "
            "gene renamed between genome annotation builds still matches "
            "the same gene in every other study."
        )
        desc.setWordWrap(True)
        layout.addWidget(desc)

        layout.addSpacing(10)

        search_row = QHBoxLayout()
        search_row.setSpacing(8)
        self._search_edit = QLineEdit()
        self._search_edit.setPlaceholderText(
            "Search gene name (e.g. AARS, CDK6)...")
        self._search_edit.setClearButtonEnabled(True)
        self._search_edit.setMaximumWidth(360)
        self._search_edit.textChanged.connect(self._filter_table)
        search_row.addWidget(self._search_edit)
        self._match_count_label = CaptionLabel("")
        search_row.addWidget(self._match_count_label)
        search_row.addStretch()
        layout.addLayout(search_row)

        layout.addSpacing(8)

        # Chips rather than a dropdown: the counts are the point. Seeing
        # 486 alias against 412 previous-symbol says at a glance how much of
        # this run rests on HGNC's looser category, which a collapsed
        # dropdown hides.
        chip_row = QHBoxLayout()
        chip_row.setSpacing(6)
        self._chips = {}
        chips = (
            (None, "All",
             "Every name this run would change."),
            ('previous_symbol', "Previous symbol",
             "HGNC's most authoritative mapping: the name was formally "
             "retired in favour of the new one."),
            ('alias', "Alias",
             "A looser HGNC mapping. Worth reading before you accept the "
             "run."),
            ('mitochondrial', "Mitochondrial",
             "Legacy mtDNA names (ND1, CYTB) mapped explicitly to MT-*."),
            ('_merged', "Merged",
             "Columns whose counts are added together because HGNC has "
             "merged the loci."),
            ('_unrecognised', "Unrecognised",
             "Names HGNC does not know at all -- retired clone-based names "
             "and bare Ensembl IDs. Nothing is renamed; they are listed "
             "because they match nothing in another study either, so they "
             "drop out of any cross-study join."),
        )
        for key, label, tip in chips:
            chip = QPushButton(label)
            chip.setProperty("role", "filter_chip")
            chip.setCheckable(True)
            chip.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
            chip.setToolTip(tip)
            chip.clicked.connect(
                lambda _c=False, k=key: self._set_chip_filter(k))
            chip_row.addWidget(chip)
            self._chips[key] = chip
        self._chips[None].setChecked(True)
        self._chip_filter = None
        chip_row.addStretch()
        layout.addLayout(chip_row)

        layout.addSpacing(10)

        self._table = ResultsTable()
        self._table.set_schema([
            Column("Current name", "old", "s",
                   tooltip="The name as it appears in this dataset."),
            Column("HGNC symbol", "new", "s",
                   tooltip="The approved symbol it will be renamed to."),
            Column("Change", "reason", "s",
                   tooltip="Why HGNC maps this name: a previous symbol, an "
                           "alias, or the mitochondrial special case."),
            Column("Merged with", "merged_with", "s",
                   tooltip="Other columns whose counts are added into the "
                           "same gene. Blank when nothing is summed."),
            Column("Details", "details", "s",
                   tooltip="What kind of gene HGNC says this is."),
        ])
        # The two name columns carry the content, so they take the width.
        # 'Merged with' is populated on a small minority of rows -- letting
        # it stretch left the real columns crushed against a half-empty pane.
        hdr = self._table.horizontalHeader()
        for col in (0, 1, 3, 4):
            hdr.setSectionResizeMode(col, QHeaderView.ResizeMode.Stretch)
        hdr.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        hdr.setHighlightSections(False)
        # Row numbers are relabelled by _renumber_rows so they read 1..n
        # down the visible rows rather than carrying positional gaps.
        vhdr = self._table.verticalHeader()
        vhdr.setVisible(True)
        vhdr.setFixedWidth(52)
        hdr.sortIndicatorChanged.connect(
            lambda _c, _o: self._renumber_rows())
        # Flat rows, like the Project table. Alternating bands fight the
        # hover and selection colours and make a long list harder to read
        # across, not easier.
        self._table.setAlternatingRowColors(False)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSortingEnabled(True)
        self._table.setSelectionBehavior(
            QTableWidget.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(
            QTableWidget.SelectionMode.SingleSelection)
        vhdr.setDefaultSectionSize(32)
        vhdr.setSectionResizeMode(QHeaderView.ResizeMode.Fixed)
        # Per-pixel scrolling: the default jumps a whole row at a time,
        # which reads as a stutter on a list this long. Columns stretch to
        # fit, so the horizontal bar is never needed.
        self._table.setVerticalScrollMode(
            QTableWidget.ScrollMode.ScrollPerPixel)
        self._table.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        layout.addWidget(self._table, 1)

    def _update_stats(self):
        """Fill the sidebar stats panel from the report.

        Protein-coding and unrecognised are shown as they will stand *after*
        harmonising, with the change in brackets until it has been applied.
        """
        report = self._report
        if report is None:
            for value in self._stat_rows.values():
                value.setText("—")
            self._merge_note.hide()
            return

        done = self._harmonised
        self._stat_rows['total'].setText(f"{report.total_genes:,}")
        suffix = " (done)" if done else ""
        self._stat_rows['renamed'].setText(f"{report.renamed:,}{suffix}")
        self._stat_rows['merged'].setText(
            f"{report.duplicates_merged:,}{suffix}")

        after, before = report.locus_after, report.locus_before
        for key, group in (('protein_coding', 'protein-coding gene'),
                           ('unrecognised', 'unrecognised')):
            n = after.get(group)
            if n is None:
                self._stat_rows[key].setText("—")
                continue
            delta = n - before.get(group, n)
            self._stat_rows[key].setText(
                f"{n:,}" if (done or not delta) else f"{n:,} ({delta:+,})")

        n_merged = report.duplicates_merged
        if n_merged:
            verb = "were" if done else "will be"
            self._merge_note.setText(
                f"{n_merged:,} gene(s) {verb} built by adding two or more "
                f"columns together, because HGNC has merged those loci. "
                f"Correct here; check it at Combine time, where studies "
                f"differ in how many columns they contribute.")
            self._merge_note.show()
        else:
            self._merge_note.hide()

    def _set_status(self, text: str) -> None:
        """Set the sidebar status line, hiding it when there is nothing
        to say. The headline state ("N genes. Check for outdated symbols")
        belongs to the workflow STATUS pane; this is action feedback only.
        """
        self._status_label.setText(text)
        self._status_label.setVisible(bool(text))

    def _recorded_harmonisation(self):
        """The study's last recorded gene_names stage, or None.

        Re-deriving "has this been done?" from the data does not work: a
        harmonised study and an Ensembl-indexed one both yield zero
        matches, and those mean opposite things. The provenance sidecar
        says which, and unlike the in-memory flag it survives a restart.
        """
        h5ad = getattr(self.main_window, 'current_h5ad_path', None)
        if not h5ad:
            return None
        try:
            from pathlib import Path

            from kosmic.provenance import load
            record = load(Path(h5ad).parent)
        except Exception:
            return None
        if not record:
            return None
        for stage in reversed(record.get('stages', [])):
            if stage.get('stage') == 'gene_names':
                return stage
        return None

    def _apply_recorded_state(self, stage) -> None:
        """Present a study that was already harmonised, as history.

        The analysis still runs -- it is what supplies the locus groups for
        the Details column and the list of unrecognised names, both of
        which describe the data as it stands. But its rename list is empty,
        because harmonising leaves nothing further to rename. Folding the
        recorded stage back in means the chips and the table show what was
        actually done rather than an empty table, which is indistinguishable
        from "we found nothing".
        """
        params = stage.get('params', {})
        when = (stage.get('timestamp') or '')[:16].replace('T', ' ')

        renames = params.get('renames')
        report = self._report
        report.mappings = [tuple(entry) for entry in (renames or [])]
        report.merged_genes = [
            (target, list(sources))
            for target, sources in (params.get('merged') or {}).items()]
        report.duplicates_merged = len(report.merged_genes)
        report.renamed = (len(report.mappings) if renames is not None
                          else params.get('renamed', 0))

        self._harmonised = True
        self._apply_btn.setEnabled(False)
        self._skip_btn.setEnabled(False)

        stamp = f" on {when}" if when else ""
        if renames is None:
            # Harmonised by a build that recorded only the counts. The
            # merges are still listed; the individual renames are gone.
            self._set_status(
                f"Harmonised{stamp}: {params.get('renamed', 0):,} renamed, "
                f"{report.duplicates_merged:,} merged. Only the merges were "
                "recorded in enough detail to list here — the full "
                "rename list is kept for runs from this version onwards.")
        else:
            self._set_status(
                f"Harmonised{stamp}. The table shows the changes that "
                "were applied.")

    def on_tab_activated(self):
        """Called when the tab becomes visible."""
        adata = self.main_window.current_adata
        if adata is None:
            self._set_status("No data loaded. Load data first.")
            return
        if self._harmonised:
            return

        if self._worker is not None and self._worker.isRunning():
            return
        # Only run if this is genuinely new data we haven't seen
        adata_id = id(adata)
        if adata_id == self._last_adata_id:
            return
        self._last_adata_id = adata_id
        self._run_harmonisation(adata)

    def set_adata(self, adata):
        """Receive adata broadcast — reset state for new data."""
        # Don't reset if we are the source of this update
        if self._applying:
            return
        self._harmonised = False
        self._report = None
        self._pending_adata = None
        self._last_adata_id = None
        self._table.setRowCount(0)
        self._set_status("")

    # ------------------------------------------------------------------
    # Harmonisation
    # ------------------------------------------------------------------

    def _set_sidebar_busy(self, busy):
        """Set sidebar progress bar to indeterminate or reset."""
        if self.progress_bar is not None:
            if busy:
                self.progress_bar.setRange(0, 0)
            else:
                self.progress_bar.setRange(0, 100)
                self.progress_bar.setValue(0)

    def _run_harmonisation(self, adata):
        """Start the harmonisation worker (lightweight — only reads var_names)."""
        self._set_sidebar_busy(True)
        self._set_status("Analysing gene names...")
        self._apply_btn.setEnabled(False)
        self._skip_btn.setEnabled(False)

        self._worker = _HarmoniseWorker(list(adata.var_names))
        run_worker(
            self._worker,
            on_finished=self._on_finished,
            on_progress=self._on_progress,
        )

    def _on_progress(self, msg):
        self._set_status(msg)
        self.log_message.emit(msg)

    def _on_finished(self, payload):
        report_or_error, lookup_data = payload
        self._set_sidebar_busy(False)

        if report_or_error is None:
            # Error — lookup_data contains the error string
            self._set_status(f"Error: {lookup_data}")
            self.log_message.emit(f"Gene name harmonisation failed: {lookup_data}")
            return

        self._report = report_or_error
        # (lookup, reasons) drive Apply; the locus groups fill the table's
        # Details column without reloading the 45k-row table per repaint.
        lookup, reasons, locus_groups = lookup_data
        self._lookup_data = (lookup, reasons)
        self._locus_groups = locus_groups

        # Update summary stats

        # A study harmonised in an earlier session has nothing left to
        # rename, so the fresh analysis comes back empty. Fold in what was
        # recorded then, so the chips and table show the applied changes.
        recorded = self._recorded_harmonisation()
        if recorded is not None:
            self._apply_recorded_state(recorded)
            self.main_window.mark_step_complete(1)

        # Populate table
        self._populate_table()
        self._update_stats()

        # Enable buttons
        if self._harmonised:
            # Already applied in an earlier session: the rename counts now
            # describe history, not work outstanding, so leave the buttons
            # disabled and the recorded status in place.
            self.log_message.emit(
                f"Gene names: already harmonised "
                f"({self._report.renamed} renamed, "
                f"{self._report.duplicates_merged} merged)")
            return

        self._apply_btn.setEnabled(self._report.renamed > 0)
        self._skip_btn.setEnabled(True)

        if self._report.renamed == 0:
            # Zero renames on an Ensembl-ID dataset means nothing matched,
            # not that everything is current -- say so instead of a false
            # all-clear.
            adata = getattr(self.main_window, 'current_adata', None)
            looks_ensembl = False
            if adata is not None and adata.n_vars:
                sample_names = adata.var_names[:200]
                n_ens = int(sample_names.str.upper()
                            .str.startswith("ENSG").sum())
                looks_ensembl = n_ens > len(sample_names) * 0.5
            if looks_ensembl:
                self._set_status(
                    "Gene names look like Ensembl IDs (ENSG...), which the "
                    "HGNC table cannot match. Re-import the data so symbols "
                    "are used, or continue knowing downstream marker lists "
                    "and pathways will not match.")
                self.log_message.emit(
                    "Gene names: var_names look like Ensembl IDs -- "
                    "harmonisation skipped, downstream symbol matching "
                    "will fail")
            else:
                self._set_status(
                    "All gene names are already up to date with HGNC.")
                self.log_message.emit(
                    "Gene names: all already current HGNC symbols")
            self._harmonised = True
            self._update_stats()
            self.main_window.mark_step_complete(1)
        else:
            # The workflow STATUS pane already says "N genes. Check for
            # outdated symbols"; the chips carry the breakdown.
            self.log_message.emit(
                f"Gene names: {self._report.renamed} can be renamed, "
                f"{self._report.duplicates_merged} duplicates to merge"
            )

    def _set_chip_filter(self, key):
        """Switch the active chip and refilter."""
        self._chip_filter = key
        for chip_key, chip in self._chips.items():
            chip.setChecked(chip_key == key)
        self._populate_table()

    def _refresh_chip_counts(self):
        """Put each category's row count on its chip.

        A category with nothing in it is disabled rather than hidden, so
        the absence is itself readable -- 'Mitochondrial (0)' says the
        legacy mtDNA names were not present, not that we forgot to look.
        """
        report = self._report
        if report is None:
            for key, chip in self._chips.items():
                chip.setText(self._CHIP_LABELS[key])
                chip.setEnabled(False)
            return

        counts = {key: 0 for key in self._chips}
        renamed_from = set()
        for old, _new, reason in report.mappings:
            renamed_from.add(old)
            if reason in counts:
                counts[reason] += 1
        # 'Merged' counts every source column folded into another gene;
        # those that were not also renamed are extra rows in the All view,
        # so 'All' has to include them or the chip disagrees with the table.
        merge_sources = [source for _t, sources in report.merged_genes
                         for source in sources]
        counts['_merged'] = len(merge_sources)
        counts[None] = len(report.mappings) + sum(
            1 for source in merge_sources if source not in renamed_from)
        counts['_unrecognised'] = len(report.unrecognised)

        for key, chip in self._chips.items():
            chip.setText(f"{self._CHIP_LABELS[key]} ({counts[key]:,})")
            # 'All' stays live even at zero so there is always a way back.
            chip.setEnabled(key is None or counts[key] > 0)

    def _populate_table(self):
        """Fill the mapping table for the active chip.

        One row per column whose identity changes: every rename, plus any
        column that keeps its name but absorbs another (the primary of a
        merge). Merges were previously computed here and then dropped on
        the floor, which hid exactly the rows worth looking at.
        """
        report = self._report
        self._refresh_chip_counts()
        if report is None:
            self._table.setRowCount(0)
            self._match_count_label.setText("")
            return

        groups = self._locus_groups or {}

        if self._chip_filter == '_unrecognised':
            # These are not renames -- there is no target -- so the table
            # shows the name alone and says why it is here.
            entries = [{'old': name, 'new': '', 'reason': 'unrecognised',
                        'merged_with': '',
                        'details': 'Not an HGNC symbol'}
                       for name in report.unrecognised]
        else:
            merge_sources = {target: list(sources)
                             for target, sources in report.merged_genes}
            entries = []
            renamed_from = set()
            for old, new, reason in report.mappings:
                renamed_from.add(old)
                entries.append({
                    'old': old, 'new': new, 'reason': reason,
                    'merged_with': ", ".join(
                        s for s in merge_sources.get(new, []) if s != old),
                    'details': groups.get(new, ''),
                })
            # A merge always has a column that kept its name; it has no
            # mapping entry, but its counts change, so it belongs here.
            for target, sources in merge_sources.items():
                for source in sources:
                    if source in renamed_from:
                        continue
                    entries.append({
                        'old': source, 'new': target, 'reason': 'merged',
                        'merged_with': ", ".join(
                            s for s in sources if s != source),
                        'details': groups.get(target, ''),
                    })
            if self._chip_filter == '_merged':
                entries = [e for e in entries if e['merged_with']
                           or e['reason'] == 'merged']
            elif self._chip_filter is not None:
                entries = [e for e in entries
                           if e['reason'] == self._chip_filter]

        self._all_entries = entries
        df = pd.DataFrame({
            'old': [e['old'] for e in entries],
            'new': [e['new'] for e in entries],
            'reason': [self._REASON_LABELS.get(e['reason'], e['reason'])
                       for e in entries],
            'merged_with': [e['merged_with'] for e in entries],
            'details': [e['details'] for e in entries],
        })
        if not df.empty:
            df = df.sort_values('old', key=lambda c: c.str.upper())
        self._table.set_data(df)
        self._table.sortItems(0, Qt.SortOrder.AscendingOrder)
        self._filter_table(self._search_edit.text())

    def _filter_table(self, text):
        """Apply the search box on top of the active chip."""
        text = text.strip().upper()
        total = self._table.rowCount()
        visible = 0
        for row in range(total):
            if text:
                old_item = self._table.item(row, 0)
                new_item = self._table.item(row, 1)
                keeps = (text in (old_item.text() if old_item else '').upper()
                         or text in (new_item.text() if new_item else '').upper())
            else:
                keeps = True
            self._table.setRowHidden(row, not keeps)
            visible += keeps

        self._renumber_rows()
        noun = "row" if visible == 1 else "rows"
        self._match_count_label.setText(
            f"{visible:,} of {total:,} {noun}" if text
            else f"{total:,} {noun}")

    def _renumber_rows(self):
        """Number the visible rows 1..n down the screen.

        Qt's own row numbers are positional, so filtering leaves gaps
        (1, 2, 5, 9...) and sorting leaves them out of order. Relabelling
        after every filter and sort keeps the count meaning "the nth row
        you can see".
        """
        header = self._table.verticalHeader()
        if not header.isVisible():
            return
        n = 0
        for row in range(self._table.rowCount()):
            if self._table.isRowHidden(row):
                self._table.setVerticalHeaderItem(row, QTableWidgetItem(""))
                continue
            n += 1
            item = QTableWidgetItem(str(n))
            item.setTextAlignment(Qt.AlignmentFlag.AlignRight
                                  | Qt.AlignmentFlag.AlignVCenter)
            self._table.setVerticalHeaderItem(row, item)

    # ------------------------------------------------------------------
    # Actions
    # ------------------------------------------------------------------

    def _on_apply(self):
        """Apply harmonisation — modify adata in a background worker."""
        adata = self.main_window.current_adata
        if adata is None or self._report is None:
            return

        self._apply_btn.setEnabled(False)
        self._skip_btn.setEnabled(False)
        self._set_sidebar_busy(True)
        self._set_status("Applying harmonisation...")
        self.log_message.emit("Applying gene name harmonisation...")

        lookup, reasons = self._lookup_data
        save_path = self.main_window.current_h5ad_path
        self._apply_worker = _ApplyWorker(adata, lookup, reasons, save_path=save_path)
        run_worker(
            self._apply_worker,
            on_finished=self._on_apply_finished,
            on_progress=self._on_progress,
        )

    def _on_apply_finished(self, adata_new):
        """Handle completed apply worker."""
        self._set_sidebar_busy(False)

        if adata_new is None:
            self._set_status("Error applying harmonisation.")
            self._apply_btn.setEnabled(True)
            self._skip_btn.setEnabled(True)
            return

        self._applying = True
        self.main_window.set_adata(adata_new)
        self._applying = False
        self._harmonised = True

        n = self._report.renamed
        m = self._report.duplicates_merged

        # Also update the raw data file so reloads keep harmonised names
        raw_path = getattr(self.main_window, 'raw_h5ad_path', None)
        if raw_path:
            try:
                adata_new.write_h5ad(raw_path)
                self.log_message.emit(f"Raw data updated with harmonised gene names: {raw_path}")
            except Exception as e:
                self.log_message.emit(f"Warning: could not update raw file: {e}")

        self._set_status(
            f"Applied: {n} gene(s) renamed, {m} duplicate(s) merged. "
            f"Now {adata_new.n_vars:,} genes. Saved."
        )
        self.log_message.emit(f"Gene names harmonised: {n} renamed, {m} merged → {adata_new.n_vars:,} genes")

        self._record_provenance(adata_new)
        self._update_stats()

        # Mark step complete
        self.main_window.mark_step_complete(1)

    def _record_provenance(self, adata_new):
        """Write what harmonisation did into the study's sidecar.

        The merge map is the load-bearing part. Once two columns have been
        summed there is nothing in the h5ad to say it happened -- the study
        just has a gene called HNRNPU like everyone else -- so this record
        is the only way the Combine step can later tell that the same gene
        was built from two features here and one elsewhere.
        """
        report = self._report
        self.main_window.record_provenance(
            'gene_names',
            {
                'total_genes_before': report.total_genes,
                'total_genes_after': int(adata_new.n_vars),
                'renamed': report.renamed,
                # The individual renames, not just the count: once applied
                # they leave no trace in the data, so this list is the only
                # way the tab can still show what it did on a later launch.
                'renames': [[old, new, reason]
                            for old, new, reason in report.mappings],
                'merged': {target: list(sources)
                           for target, sources in report.merged_genes},
                'protein_coding_after':
                    report.locus_after.get('protein-coding gene'),
                'unrecognised_after': report.locus_after.get('unrecognised'),
            },
        )

    def _on_skip(self):
        """Skip harmonisation — keep original names."""
        self._harmonised = True
        self._apply_btn.setEnabled(False)
        self._skip_btn.setEnabled(False)
        self._set_status("Skipped — keeping original gene names.")
        self._update_stats()
        self.log_message.emit("Gene name harmonisation skipped")
        self.main_window.mark_step_complete(1)

    def _on_update_hgnc(self):
        """Download the latest HGNC table."""
        try:
            from kosmic.scrna.inspect.gene_names import download_hgnc_table
            path = download_hgnc_table()
            dialogs.info(self, "HGNC Table Updated",
                                    f"Downloaded latest HGNC symbol table to:\n{path}")
            self.log_message.emit(f"HGNC table updated: {path}")
            # Reset so it re-runs with new table
            if self.main_window.current_adata is not None:
                self._harmonised = False
                self._run_harmonisation(self.main_window.current_adata)
        except Exception as e:
            dialogs.warning(self, "Download Failed",
                                f"Could not download HGNC table:\n{e}")

    def refresh_theme(self):
        """Re-apply theme -- palette-driven; no custom action needed."""
        pass
