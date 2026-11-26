# The Methods Comparison workspace's one page (ADR-017).
#
# Sidebar: the counts source and cell type, the studies, the methods and
# the call rule, and Run comparison. Content: the merged design and its
# warnings, where each method's result came from, a summary per method and
# the pairwise agreement.
#
# Choosing a selection checks its merged design in a worker (the atlas's
# obs are read) and shows the saved comparison for it, if there is one.
# Running asks first when the atlas's saved merged run does not match the
# design: Resolve (fit again), Run anyway (use it and record the
# differences) or Cancel. Changing the rule after a run re-cuts the
# comparison from its saved table, without fitting anything, and saves it.

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox, QHeaderView, QListWidget,
    QFrame, QListWidgetItem, QScrollArea, QSpinBox, QTableWidgetItem, QVBoxLayout, QWidget,
)

from kosmic.comparison.agreement import RULE_LABELS, RULES, Rule, summary
from kosmic.comparison.run import ALL_METHODS
from kosmic.comparison.standard import METHOD_LABELS
from kosmic.gui.shared import dialogs, run_worker
from kosmic.gui.shared.widgets import (
    CopyableTableWidget, HintLabel, PrimaryButton, SectionHeader, SettingsGroup,
    SidebarPage, StatusLabel,
)

_RULE_NAMES = {'significance': 'Adjusted p-value', 'top_n': 'Top N by rank',
               'significance_and_direction': 'Adjusted p-value, same direction'}

_SOURCE_TEXT = {'saved': 'Saved run', 'pooled': 'Pooled again', 'fitted': 'Fitted'}

RESOLVE, RUN_ANYWAY, CANCEL = 'Resolve', 'Run anyway', 'Cancel'


def _table(headers) -> CopyableTableWidget:
    t = CopyableTableWidget(0, len(headers))
    t.setHorizontalHeaderLabels(headers)
    t.verticalHeader().setVisible(False)
    t.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    t.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
    t.setSizeAdjustPolicy(QAbstractItemView.SizeAdjustPolicy.AdjustToContents)
    return t


def _fill(table, rows, first_number: int = 1) -> None:
    """Fill ``table``; columns from ``first_number`` on are right-aligned numbers."""
    table.setRowCount(len(rows))
    for r, row in enumerate(rows):
        for c, value in enumerate(row):
            item = QTableWidgetItem(value)
            if c >= first_number:
                item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            table.setItem(r, c, item)
    table.resizeColumnsToContents()
    # As tall as its rows, so the two tables sit one under the other.
    height = table.horizontalHeader().height() + 2 * table.frameWidth()
    height += sum(table.rowHeight(r) for r in range(table.rowCount()))
    table.setFixedHeight(height)


def _num(value, fmt='{:,}') -> str:
    try:
        if value != value:   # NaN
            return '–'
    except TypeError:
        return '–'
    return fmt.format(value)


def mismatch_text(checked) -> str:
    """The dialog text for a saved merged run whose design differs."""
    lines = [f"The studies' design: {checked.design.formula}",
             f"The saved merged run: {(checked.saved.record or {}).get('design') or 'not recorded'}",
             "", "Differences:"]
    lines += [f"- {d}" for d in checked.saved.differences]
    lines += ["", "Resolve fits merged pseudobulk again with the studies' design and "
                  "replaces the saved merged run for this counts source and cell type.",
              "Run anyway uses the saved run as it is and records the differences "
              "with the comparison."]
    if not checked.design.can_fit:
        lines += ["", "Resolve is not available: the atlas has no column for "
                      + ', '.join(checked.design.unresolved) + '.']
    return '\n'.join(lines)


class ComparePage(SidebarPage):
    """Merged pseudobulk against the meta-analysis methods."""

    help_id = "comparison/index"

    def __init__(self, workspace):
        super().__init__()
        self.workspace = workspace
        self.project = None
        self.selections: list = []
        self.checked = None
        self._check_worker = None
        self._run_worker = None
        # Every worker started, until it finishes: a QThread destroyed while
        # running aborts the process, so a check replaced by a newer one is
        # kept here rather than dropped.
        self._workers: list = []
        self._build_sidebar()
        self._build_content()
        self._refresh_enabled()

    # -- layout ----------------------------------------------------------------

    def _build_sidebar(self):
        sel = SettingsGroup("Selection")
        self.selection_combo = QComboBox()
        self.selection_combo.currentIndexChanged.connect(self._on_selection_changed)
        sel.add_widget(self.selection_combo)
        self.study_list = QListWidget()
        self.study_list.setMinimumHeight(110)   # about four studies; the sidebar scrolls
        self.study_list.setMaximumHeight(130)
        self.study_list.itemChanged.connect(self._on_studies_changed)
        sel.add_widget(self.study_list)
        self.sidebar_layout.addWidget(sel)

        meth = SettingsGroup("Methods")
        self.method_checks = {}
        for mid in ALL_METHODS:
            box = QCheckBox(METHOD_LABELS[mid])
            box.setChecked(True)
            box.toggled.connect(self._on_rule_changed)
            meth.add_widget(box)
            self.method_checks[mid] = box
        self.sidebar_layout.addWidget(meth)

        rule = SettingsGroup("Calls")
        self.rule_combo = QComboBox()
        for key in RULES:
            self.rule_combo.addItem(_RULE_NAMES[key], key)
            self.rule_combo.setItemData(self.rule_combo.count() - 1, RULE_LABELS[key],
                                        Qt.ItemDataRole.ToolTipRole)
        self.rule_combo.currentIndexChanged.connect(self._on_rule_changed)
        rule.add_widget(self.rule_combo)
        self.alpha_spin = QDoubleSpinBox()
        self.alpha_spin.setDecimals(3)
        self.alpha_spin.setRange(0.001, 0.5)
        self.alpha_spin.setSingleStep(0.01)
        self.alpha_spin.setValue(0.05)
        self.alpha_spin.valueChanged.connect(self._on_rule_changed)
        rule.add_row("Adjusted p below", self.alpha_spin)
        self.top_spin = QSpinBox()
        self.top_spin.setRange(1, 100000)
        self.top_spin.setValue(100)
        self.top_spin.valueChanged.connect(self._on_rule_changed)
        rule.add_row("Top N genes", self.top_spin)
        self.sidebar_layout.addWidget(rule)

        self.run_btn = PrimaryButton("Run comparison")
        self.run_btn.clicked.connect(self._on_run)
        self.sidebar_layout.addWidget(self.run_btn)
        self.status = StatusLabel("")
        self.status.setWordWrap(True)
        self.sidebar_layout.addWidget(self.status)
        self.sidebar_layout.addStretch()
        self._sync_rule_inputs()

    def _build_content(self):
        # The tables are as tall as their rows, so the content is taller
        # than a short window (or one with the output panel raised); it
        # scrolls rather than being squeezed until its widgets overlap.
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setProperty("role", "thin_scroll")
        inner = QWidget()
        lay = QVBoxLayout(inner)
        lay.setContentsMargins(0, 0, 8, 0)
        lay.setAlignment(Qt.AlignmentFlag.AlignTop)
        scroll.setWidget(inner)
        self.content_layout.addWidget(scroll)
        self.intro = HintLabel("")
        lay.addWidget(self.intro)
        lay.addWidget(SectionHeader("Merged pseudobulk design"))
        self.design_label = HintLabel("")
        self.design_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        lay.addWidget(self.design_label)
        self.warnings_label = StatusLabel("", state='warning')
        self.warnings_label.setWordWrap(True)
        lay.addWidget(self.warnings_label)
        lay.addWidget(SectionHeader("Genes called"))
        self.rule_label = HintLabel("")
        lay.addWidget(self.rule_label)
        self.summary_table = _table(["Method", "Result from", "Genes compared", "Called",
                                     "Up", "Down"])
        lay.addWidget(self.summary_table)
        lay.addWidget(SectionHeader("Agreement between methods"))
        self.pairwise_table = _table(["Methods", "Both call", "Jaccard",
                                      "Same direction", "Rank correlation"])
        lay.addWidget(self.pairwise_table)
        self.notes_label = HintLabel("")
        lay.addWidget(self.notes_label)
        self._show_bundle(None)

    # -- state -----------------------------------------------------------------

    def selection(self):
        i = self.selection_combo.currentIndex()
        if i < 0 or i >= len(self.selections):
            return None
        chosen = [self.study_list.item(r).data(Qt.ItemDataRole.UserRole)
                  for r in range(self.study_list.count())
                  if self.study_list.item(r).checkState() == Qt.CheckState.Checked]
        return self.selections[i].restricted(chosen)

    def methods(self) -> list:
        return [m for m, box in self.method_checks.items() if box.isChecked()]

    def rule(self) -> Rule:
        return Rule(self.rule_combo.currentData(), alpha=float(self.alpha_spin.value()),
                    top_n=int(self.top_spin.value()))

    def _busy(self) -> bool:
        return self._run_worker is not None and self._run_worker.isRunning()

    def _start(self, worker, **callbacks):
        self._workers.append(worker)
        worker.finished.connect(lambda w=worker: self._workers.remove(w)
                                if w in self._workers else None)
        return run_worker(worker, **callbacks)

    def _refresh_enabled(self):
        sel = self.selection()
        ok = (sel is not None and len(sel.runs) >= 2 and len(self.methods()) >= 2
              and self.checked is not None and not self._busy())
        self.run_btn.setEnabled(ok)

    def _sync_rule_inputs(self):
        top = self.rule_combo.currentData() == 'top_n'
        self.alpha_spin.setEnabled(not top)
        self.top_spin.setEnabled(top)

    # -- project and selection ---------------------------------------------------

    def set_project(self, project):
        """Show ``project``'s selections, keeping the one chosen when the
        project and its selections are unchanged."""
        from kosmic.comparison.run import available_selections

        selections = available_selections(project) if project is not None else []
        if project == self.project and selections == self.selections:
            return
        current = self.selection()
        self.project = project
        self.selections = selections
        self.selection_combo.blockSignals(True)
        self.selection_combo.clear()
        for s in self.selections:
            self.selection_combo.addItem(s.label)
            self.selection_combo.setItemData(self.selection_combo.count() - 1, s.label,
                                             Qt.ItemDataRole.ToolTipRole)
        paths = [s.path for s in self.selections]
        if current is not None and current.path in paths:
            self.selection_combo.setCurrentIndex(paths.index(current.path))
        self.selection_combo.blockSignals(False)
        if project is None:
            self.intro.setText("Open a project to compare merged pseudobulk with "
                               "meta-analysis.")
        elif not self.selections:
            self.intro.setText("No cell type has a saved DESeq2 run in two or more studies. "
                               "Run DE per cell type in each study first.")
        else:
            self.intro.setText("Compare the genes merged pseudobulk on the shared atlas "
                               "calls with those the meta-analysis methods call, for one "
                               "counts source and cell type.")
        self._on_selection_changed()

    def _on_selection_changed(self, *_):
        i = self.selection_combo.currentIndex()
        self.study_list.blockSignals(True)
        self.study_list.clear()
        if 0 <= i < len(self.selections):
            for study in self.selections[i].studies:
                item = QListWidgetItem(study)
                item.setData(Qt.ItemDataRole.UserRole, study)
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(Qt.CheckState.Checked)
                self.study_list.addItem(item)
        self.study_list.blockSignals(False)
        self.checked = None
        self._show_design(None)
        saved = self._saved_bundle()
        self.workspace.set_bundle(saved)
        self._show_bundle(saved)
        self._start_check()

    def _on_studies_changed(self, *_):
        """Another set of studies has its own design check."""
        self.checked = None
        self._show_design(None)
        self._start_check()

    def _saved_bundle(self):
        from kosmic.comparison.bundle import load
        from kosmic.paths import comparison_output_dir
        sel = self.selection()
        if self.project is None or sel is None:
            return None
        return load(comparison_output_dir(self.project, sel.path))

    def _start_check(self):
        sel = self.selection()
        if self.project is None or sel is None:
            self._refresh_enabled()
            return
        self.design_label.setText("Checking the merged design...")
        key = (sel.path, sel.studies)
        self._check_worker = self._start(
            _check_worker(self.project, sel),
            on_finished=lambda checked, key=key: self._on_checked(checked, key),
            on_failed=self._on_failed)
        self._refresh_enabled()

    def _on_checked(self, checked, key=None):
        sel = self.selection()
        if key is not None and (sel is None or key != (sel.path, sel.studies)):
            return      # a check of a selection no longer shown
        self.checked = checked
        self._show_design(checked)
        self._refresh_enabled()

    # -- running ---------------------------------------------------------------

    def _on_run(self):
        sel, checked = self.selection(), self.checked
        if sel is None or checked is None:
            return
        action = 'auto'
        if 'pb_merged' in self.methods() and checked.needs_choice:
            options = ([RESOLVE] if checked.design.can_fit else []) + [RUN_ANYWAY, CANCEL]
            answer = dialogs.choose(self, "Merged design differs", mismatch_text(checked),
                                    options, escape=CANCEL)
            if answer in (None, CANCEL):
                return
            action = 'refit' if answer == RESOLVE else 'saved'
        elif 'pb_merged' in self.methods() and not checked.design.can_fit \
                and checked.saved is None:
            dialogs.warning(self, "Merged pseudobulk cannot be fitted",
                            "\n".join(checked.design.warnings))
            return
        self.status.set_state('info')
        self.status.setText("Running...")
        self.workspace.status_message.emit("Running the methods comparison...")
        self._run_worker = self._start(
            _comparison_worker(self.project, sel, self.methods(), self.rule(), action,
                               checked),
            on_finished=self._on_run_finished, on_failed=self._on_failed,
            on_progress=self.workspace.log_message.emit)
        self._refresh_enabled()

    def _on_run_finished(self, bundle):
        self.workspace.set_bundle(bundle)
        self._show_bundle(bundle)
        self.status.set_state('success')
        self.status.setText("Comparison saved.")
        self.workspace.status_message.emit("Methods comparison saved.")
        if self.checked is not None and bundle.params.get('merged_action') == 'refit':
            self._start_check()     # the saved merged run changed
        self._refresh_enabled()

    def _on_failed(self, message):
        first = message.splitlines()[0] if message else "Failed."
        self.status.set_state('error')
        self.status.setText(first)
        self.workspace.log_message.emit(message)
        self._refresh_enabled()

    def _on_rule_changed(self, *_):
        self._sync_rule_inputs()
        self._refresh_enabled()
        bundle = self.workspace.bundle
        methods = [m for m in self.methods() if bundle is not None and m in bundle.methods]
        if bundle is None or len(methods) < 2 or self.project is None:
            return
        recut = bundle.recut(self.rule(), methods)
        recut.save(self.project, bundle.params.get('selection'))
        self.workspace.set_bundle(recut)
        self._show_bundle(recut)

    # -- display -----------------------------------------------------------------

    def _show_design(self, checked):
        if checked is None:
            self.design_label.setText("")
            self.warnings_label.setText("")
            self.warnings_label.setVisible(False)
            return
        d, saved = checked.design, checked.saved
        lines = [f"Design: {d.formula}"]
        mapped = [f"{k} → {v}" for k, v in d.mapping.items() if v and v != k]
        if mapped:
            lines.append("Per-study covariates on the atlas: " + ', '.join(mapped))
        if saved is None:
            lines.append("No saved merged run for this selection: it will be fitted.")
        elif saved.matches:
            lines.append(f"Saved merged run: {saved.path.name} (design matches).")
        else:
            lines.append(f"Saved merged run: {saved.path.name} (design differs; you will be "
                         "asked whether to fit it again).")
        if d.filters is None:
            lines.append("The per-study runs did not record their gene filters, so they "
                         "cannot be compared with the merged run's, and a new merged fit "
                         "uses KOSMIC's default gene filters, applied within each study.")
        self.design_label.setText('\n'.join(lines))
        self.warnings_label.setText('\n'.join(d.warnings))
        self.warnings_label.setVisible(bool(d.warnings))

    def _show_bundle(self, bundle):
        if bundle is None:
            self.rule_label.setText("No comparison has been run for this selection.")
            _fill(self.summary_table, [])
            _fill(self.pairwise_table, [])
            self.notes_label.setText("")
            return
        sources = bundle.params.get('sources') or {}
        self.rule_label.setText(f"{bundle.n_universe:,} genes tested by every method; "
                                f"called: {bundle.rule.text}.")
        _fill(self.summary_table, [
            [METHOD_LABELS[r.method], _SOURCE_TEXT.get(
                (sources.get(r.method) or {}).get('source'), '–'),
             _num(r.n_tested), _num(r.n_called), _num(r.n_up), _num(r.n_down)]
            for r in summary(bundle.results).itertuples()], first_number=2)
        _fill(self.pairwise_table, [
            [f"{METHOD_LABELS[r.method_a]} – {METHOD_LABELS[r.method_b]}", _num(r.n_both),
             _num(r.jaccard, '{:.2f}'), _num(r.direction_concordance, '{:.0%}'),
             _num(r.spearman, '{:.2f}')]
            for r in bundle.pairwise.itertuples()])
        notes = [f"{METHOD_LABELS.get(m, m)}: {s.get('note', '')}" for m, s in sources.items()]
        if bundle.params.get('design_differences'):
            notes.append("Merged pseudobulk was used although its design differs: "
                         + ' '.join(bundle.params['design_differences']))
        self.notes_label.setText('\n'.join(notes))


def _check_worker(project, selection):
    from kosmic.gui.comparison.workers import CheckWorker
    return CheckWorker(project, selection)


def _comparison_worker(*args):
    from kosmic.gui.comparison.workers import ComparisonWorker
    return ComparisonWorker(*args)
