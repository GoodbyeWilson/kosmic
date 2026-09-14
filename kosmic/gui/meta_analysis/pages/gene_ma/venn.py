# Venn tab: area-proportional Venn + intersection summary.
#
# 'populate' computes per-method sig sets and emits
# 'log_message' for each summary line.
from __future__ import annotations

from typing import Optional

import pandas as pd
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QVBoxLayout, QWidget

from kosmic import DEFAULT_FDR
from kosmic.gui.meta_analysis.venn_widget import VennWidget
from kosmic.gui.shared.widgets import SecondaryLabel


class VennTab(QWidget):
    """Area-proportional Venn diagram + intersection-size summary."""

    log_message = pyqtSignal(str)

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)

        self.venn_widget = VennWidget()
        lay.addWidget(self.venn_widget, stretch=1)

        self.text_label = SecondaryLabel("")
        self.text_label.setWordWrap(True)
        lay.addWidget(self.text_label)

    # -- public API ----------------------------------------------------

    def clear(self) -> None:
        try:
            self.venn_widget.set_data({})
            self.text_label.setText("")
        except Exception:
            # The Venn widget can be in an inconsistent state during
            # workspace tear-down (no harm leaving it as-is).
            pass

    def populate(
        self,
        meta_df: Optional[pd.DataFrame],
        method_keys: Optional[list],
    ) -> None:
        """
        Render the Venn + summary, emitting 'log_message'
        for each line the parent should echo into its own log.
        """
        if meta_df is None or not method_keys:
            return

        df = meta_df
        sig_sets: dict = {}
        for key in method_keys:
            fdr_col = f'fdr_{key}'
            if fdr_col not in df.columns:
                sig_sets[key] = set()
                continue
            fdr = df[fdr_col].values.astype(float)
            sig_sets[key] = set(df.loc[fdr < DEFAULT_FDR, 'names'].values)

        self.venn_widget.set_data(sig_sets)

        n_methods = len(method_keys)
        if n_methods == 1 and sig_sets:
            key = method_keys[0]
            n = len(sig_sets[key])
            lines = [f"{key}: {n} significant genes"]
            self.log_message.emit(f"{key}: {n} significant genes")
        else:
            consensus = (set.intersection(*sig_sets.values())
                         if sig_sets else set())
            lines = [f"All {n_methods} families: {len(consensus)} genes"]
            for key in method_keys:
                lines.append(f"  {key}: {len(sig_sets[key])} significant")
            self.log_message.emit(
                f"{len(consensus)} genes in all {n_methods} families")
        self.text_label.setText('\n'.join(lines))


__all__ = ["VennTab"]
