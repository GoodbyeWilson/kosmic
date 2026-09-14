# DE Analysis — Gene Set Selection Page
# Wraps PathwayExplorer as a standalone workflow step.

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QHBoxLayout

from kosmic.gui.shared.pathway_explorer import PathwayExplorer
from kosmic.gui.shared.widgets import (
    PrimaryButton, SecondaryLabel, SimplePage,
)


class GeneSetPage(SimplePage):
    """Select Gene Sets page: choose which pathways to analyse."""

    help_id = "de/gene_set"
    CONTENT_MARGINS = (0, 0, 0, 0)
    CONTENT_SPACING = 0

    def __init__(self, workspace):
        super().__init__()
        self.ws = workspace
        self.pathway_explorer = PathwayExplorer(workspace=self.ws)
        self.pathway_explorer.pathways_changed.connect(self._on_pathways_changed)
        self.body_layout.addWidget(self.pathway_explorer)

        # Confirmation strip: summary + explicit confirm to tick the step.
        confirm_row = QHBoxLayout()
        confirm_row.setContentsMargins(16, 8, 16, 12)
        self._summary_label = SecondaryLabel("No gene set loaded")
        confirm_row.addWidget(self._summary_label, 1, Qt.AlignmentFlag.AlignVCenter)
        self._confirm_btn = PrimaryButton("Confirm Gene Sets")
        self._confirm_btn.setMinimumHeight(32)
        self._confirm_btn.clicked.connect(self._on_confirm_clicked)
        confirm_row.addWidget(self._confirm_btn, 0, Qt.AlignmentFlag.AlignVCenter)
        self.body_layout.addLayout(confirm_row)

        self._refresh_summary()

    def on_activated(self):
        """Called when page becomes visible. Update coverage if adata available."""
        if self.ws.current_adata is not None:
            self.pathway_explorer.set_var_names(self.ws.current_adata.var_names)
        # The explorer auto-loads a default gene set without emitting
        # 'pathways_changed', so pull its current selection into the
        # workspace so the Confirm button reflects what's on screen.
        if not self.ws.pathway_gene_sets:
            current = self.pathway_explorer.get_flat_pathways()
            if current:
                self._on_pathways_changed(current)
        self._refresh_summary()

    def _on_pathways_changed(self, flat_pathways):
        """Called when the Pathway Explorer's active pathways change."""
        from kosmic.reference.pathways import derive_hierarchy, build_parent_gene_sets

        self.ws.pathway_gene_sets = flat_pathways

        hierarchy = derive_hierarchy(flat_pathways)
        if hierarchy:
            self.ws.current_hierarchy = hierarchy
            self.ws.parent_gene_sets = build_parent_gene_sets(flat_pathways, hierarchy)
        else:
            self.ws.current_hierarchy = None
            self.ws.parent_gene_sets = None

        n_pw = len(flat_pathways)
        n_genes = len({g for genes in flat_pathways.values() for g in genes})
        self.ws.geneset_label = f"{n_pw} pathways ({n_genes} genes)"
        self.ws.log_message.emit(f"Pathways updated: {n_pw} pathways, {n_genes} genes")
        self._refresh_summary()

    def _refresh_summary(self):
        gs = self.ws.pathway_gene_sets or {}
        if not gs:
            self._summary_label.setText("No gene set loaded")
            self._confirm_btn.setEnabled(False)
            return
        n_pw = len(gs)
        n_genes = len({g for genes in gs.values() for g in genes})
        self._summary_label.setText(
            f"{n_pw} pathways selected ({n_genes} unique genes)")
        self._confirm_btn.setEnabled(True)

    def _on_confirm_clicked(self):
        if not self.ws.pathway_gene_sets:
            return
        self.ws.mark_step_complete(2)  # Step 2 = Select Gene Sets (scoring mode)
        # Advance to Gene DE (step 3) -- it must run before the pathway analysis.
        if hasattr(self.ws, 'switch_tab'):
            self.ws.switch_tab(3)
