# The DE workspace's run bar (ADR-016): which gene DE run the steps after
# gene DE use.
#
# A study analysed per cell type has one gene DE run per cell type, and
# may have a whole-study run as well. The bar sits above the workspace's
# pages and is shown on Gene DE, Pathway DE, Enrichment and Methods. Its box
# lists the runs: the whole-study run first, then each cell type with its
# significant genes and donors; a cell type the batch run skipped, or whose
# run failed, is listed but cannot be chosen, with the reason. The line
# under the box gives the chosen run's method, counts, donors per arm and
# significant genes under its rule.
#
# The bar holds no state of its own: the workspace keeps the runs and the
# selection ('DEWorkspace.cell_type_runs', 'DEWorkspace.run_index') and
# fills the bar from them; choosing an entry asks the workspace to select
# that run ('run_chosen').
from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QFrame, QHBoxLayout, QLabel, QVBoxLayout

from kosmic.gui.shared.theme import NoScrollComboBox, get_color
from kosmic.gui.shared.widgets import SecondaryLabel

#: The box's entry for a run made on the whole study.
WHOLE_STUDY = "Whole study"


def is_whole_study(run) -> bool:
    """Whether 'run' (a 'kosmic.de.batch.CellTypeRun') tested the whole study."""
    return not run.cell_type


def run_name(run) -> str:
    return WHOLE_STUDY if is_whole_study(run) else run.cell_type


def entry_text(run) -> str:
    """The box's text for 'run': its name, then its significant genes and
    donors, or why it was not run."""
    if run.status != 'ok':
        reason = f": {run.message}" if run.message else ""
        return f"{run_name(run)}  ({run.status}{reason})"
    donors = f", {run.n_samples} donors" if run.n_samples else ""
    return f"{run_name(run)}  ({run.n_significant:,} significant{donors})"


_METHOD_TEXT = {'deseq2': 'DESeq2', 'welch_cpm': 'Welch t-test (CPM)',
                'welch_cpm_eb': 'Welch t-test (CPM, moderated)'}


def detail_text(run, *, method: str, counts: Optional[str],
                rule: Optional[tuple]) -> str:
    """The line under the box: method, counts, donors per arm and the
    number of significant genes under the run's rule (fdr, lfc, assumed)."""
    from kosmic.meta_analysis.io import counts_text
    from kosmic.paths import ORIGINAL_COUNTS_FOLDER

    parts = [_METHOD_TEXT.get(method, method)]
    if counts:
        text = ("original counts" if counts == ORIGINAL_COUNTS_FOLDER
                else counts_text(counts))
        parts.append(text[0].upper() + text[1:])
    if run.n_donors_control or run.n_donors_disease:
        parts.append(f"{run.n_donors_control} control, {run.n_donors_disease} disease donors")
    elif run.n_samples:
        parts.append(f"{run.n_samples} donors")
    sig = f"{run.n_significant:,} significant"
    if rule:
        fdr, lfc, assumed = rule
        sig += f" (FDR < {fdr:g}" + (f", |log2FC| > {lfc:g}" if lfc else "")
        sig += ", assumed)" if assumed else ")"
    parts.append(sig)
    return " · ".join(parts)


class RunBar(QFrame):
    """'Results for [run]' and the chosen run's details, above the DE pages."""

    run_chosen = pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("de_run_bar")
        outer = QVBoxLayout(self)
        outer.setContentsMargins(12, 8, 12, 8)
        outer.setSpacing(4)

        row = QHBoxLayout()
        row.setSpacing(8)
        title = QLabel("Results for")
        title.setStyleSheet("font-weight: 600;")
        row.addWidget(title)
        self.combo = NoScrollComboBox()
        self.combo.setMinimumWidth(320)
        self.combo.setToolTip(
            "The gene DE run that Gene DE, Pathway DE, Enrichment and Methods "
            "show. Changing it changes every page.")
        self.combo.currentIndexChanged.connect(self._on_index_changed)
        row.addWidget(self.combo, 1)
        row.addStretch(1)
        outer.addLayout(row)

        self.detail = SecondaryLabel("")
        self.detail.setWordWrap(True)
        outer.addWidget(self.detail)
        self.refresh_theme()

    def set_runs(self, runs, current: int) -> None:
        """List 'runs' and show entry 'current'; no 'run_chosen' is emitted."""
        self.combo.blockSignals(True)
        self.combo.clear()
        model = self.combo.model()
        for i, run in enumerate(runs):
            self.combo.addItem(entry_text(run), i)
            if run.status != 'ok':
                model.item(i).setEnabled(False)
        if 0 <= current < len(runs):
            self.combo.setCurrentIndex(current)
        self.combo.blockSignals(False)

    def set_detail(self, text: str) -> None:
        self.detail.setText(text)

    def _on_index_changed(self, index: int) -> None:
        if index >= 0:
            self.run_chosen.emit(index)

    def refresh_theme(self) -> None:
        self.setStyleSheet(
            f"QFrame#de_run_bar {{ background: {get_color('bg_secondary')}; "
            f"border-bottom: 1px solid {get_color('border')}; }}")
