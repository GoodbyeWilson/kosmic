# Meta-Analysis Settings Dialog
# Configures which pooling methods are available in the meta-analysis
# workflow, grouped by statistical family. Accessed from View menu.

from __future__ import annotations

from typing import Dict, List, Tuple

from PyQt6.QtCore import QSettings
from PyQt6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QGroupBox,
    QVBoxLayout,
)
from kosmic.gui.shared import dialogs
from kosmic.gui.shared.widgets import SecondaryLabel

# (family_key, family_label, [(method_key, method_label), ...])
METHOD_FAMILIES: List[Tuple[str, str, List[Tuple[str, str]]]] = [
    ('effect_size', 'Effect-size', [
        ('dl',   'DL'),
        ('reml', 'REML'),
    ]),
    ('rank', 'Rank-based', [
        ('sumrank', 'SumRank'),
        ('gwop',    'gwOP'),
    ]),
    ('pvalue', 'P-value', [
        ('fisher',   'Fisher'),
        ('stouffer', 'Stouffer'),
    ]),
]

# All method keys across all families
ALL_METHOD_KEYS = [k for _, _, methods in METHOD_FAMILIES
                   for k, _ in methods]

# Default enabled set: one best-performing method per family,
# determined by robust-frontier score on held-out LOO direction
#   Effect-size:       REML (score 23, most meta-discoveries)
#   Rank-based:        SumRank (score 25, best direction replication)
#   P-value combining: Stouffer (score 23, direction-aware)
# Others are available in View > Meta-Analysis Settings.
DEFAULT_ENABLED = {
    'reml',
    'sumrank',
    'stouffer',
}

_SETTINGS_KEY = 'meta_analysis/enabled_methods'
_DE_METHOD_KEY = 'meta_analysis/de_method'

DE_METHODS = ['DESeq2', 'Welch CPM', 'Welch Raw']
DEFAULT_DE_METHOD = 'DESeq2'


def load_enabled_methods(settings: QSettings) -> set:
    """Load enabled method keys from QSettings, or return defaults."""
    saved = settings.value(_SETTINGS_KEY)
    if saved is not None and isinstance(saved, list):
        return set(saved) & set(ALL_METHOD_KEYS)
    return set(DEFAULT_ENABLED)


def save_enabled_methods(settings: QSettings, enabled: set):
    """Persist enabled method keys to QSettings."""
    settings.setValue(_SETTINGS_KEY, sorted(enabled))


def load_de_method(settings: QSettings) -> str:
    saved = settings.value(_DE_METHOD_KEY)
    if saved and saved in DE_METHODS:
        return saved
    return DEFAULT_DE_METHOD


def save_de_method(settings: QSettings, method: str):
    settings.setValue(_DE_METHOD_KEY, method)


class MASettingsDialog(QDialog):
    """Dialog for enabling/disabling meta-analysis methods by family."""

    def __init__(self, settings: QSettings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Meta-Analysis Settings")
        self.setMinimumWidth(420)
        self._settings = settings
        self._checks: Dict[str, QCheckBox] = {}

        enabled = load_enabled_methods(settings)

        layout = QVBoxLayout(self)

        hint = SecondaryLabel(
            "Enable or disable methods available in the meta-analysis "
            "workflow. Methods in the same family use the same "
            "statistical framework; cross-family agreement is "
            "scientifically meaningful.")
        hint.setWordWrap(True)
        layout.addWidget(hint)

        for fam_key, fam_label, methods in METHOD_FAMILIES:
            group = QGroupBox(fam_label)
            g_lay = QVBoxLayout(group)
            for m_key, m_label in methods:
                cb = QCheckBox(m_label)
                cb.setChecked(m_key in enabled)
                self._checks[m_key] = cb
                g_lay.addWidget(cb)
            layout.addWidget(group)

        # DE method
        de_group = QGroupBox("DE Method")
        de_lay = QVBoxLayout(de_group)
        self._de_combo = QComboBox()
        self._de_combo.addItems(DE_METHODS)
        self._de_combo.setCurrentText(load_de_method(settings))
        self._de_combo.setToolTip(
            "Which DE method's results to load for meta-analysis.\n"
            "DESeq2: median-of-ratios (recommended).\n"
            "Welch CPM: CPM-normalized Welch t-test.\n"
            "Welch Raw: mean counts, log2(mean+1).")
        de_lay.addWidget(self._de_combo)
        layout.addWidget(de_group)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _accept(self):
        enabled = {k for k, cb in self._checks.items() if cb.isChecked()}
        if not enabled:
            dialogs.warning(
                self, "No methods",
                "At least one method must be enabled.")
            return
        save_enabled_methods(self._settings, enabled)
        save_de_method(self._settings, self._de_combo.currentText())
        self.accept()
