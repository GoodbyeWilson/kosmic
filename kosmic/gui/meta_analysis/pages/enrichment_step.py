# Discovery step 4: pathway enrichment of the pooled gene list.
#
# This IS the DE workspace's EnrichmentPage, subclassed rather than
# reimplemented. The meta side previously had its own cut-down copy --
# ORA only, seven hardcoded libraries, fixed thresholds -- while DE had
# ORA plus preranked GSEA, eleven libraries and configurable cutoffs.
# Two UIs over one shared worker is how they drifted apart, so there is
# now one UI and one implementation.
#
# EnrichmentPage reads four things off its "workspace": the DE table, a
# progress bar, log/status signals, and a writable gsea_library_name.
# '_MetaEnrichmentContext' supplies those from the pooled table.
from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import QObject, pyqtSignal
from PyQt6.QtWidgets import QWidget

from kosmic.gui.de_analysis.pages.enrichment_page import EnrichmentPage


class _MetaEnrichmentContext(QObject):
    """Presents the pooled meta-analysis table as a DE result.

    The column names differ by history, not by meaning: pooling calls
    the adjusted p-value 'fdr' and the per-gene p-value 'pvals_pooled',
    where DE calls them 'pvals_adj' and 'pvals'.
    """

    log_message = pyqtSignal(str)
    status_message = pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.de_results = None
        self.gene_de_genome_wide = None
        self.progress_bar = None
        # EnrichmentPage writes a GSEA run into the project so Figure
        # Export can reload it. Without this the meta arm's GSEA results
        # were silently not persisted -- the guard there is a plain
        # getattr(ws, 'project_dir', None).
        self.project_dir = None
        self.gsea_library_name = None
        self.enrichment_results = None

    def set_meta_results(self, meta_df):
        """Take the pooled table, or None to clear."""
        if meta_df is None or meta_df.empty:
            self.de_results = None
            self.gene_de_genome_wide = None
            return

        df = meta_df.copy()
        if 'fdr' in df.columns:
            df['pvals_adj'] = df['fdr']
        if 'pvals' not in df.columns and 'pvals_pooled' in df.columns:
            df['pvals'] = df['pvals_pooled']
        self.de_results = df
        # Pooling is genome-wide by construction, so the same table is
        # the unbiased ranking GSEA wants -- there is no narrower
        # hypothesis-mode subset to fall back from.
        self.gene_de_genome_wide = df


class MetaEnrichmentPage(EnrichmentPage):
    """The DE enrichment page, fed by the pooled meta-analysis table."""

    help_id = "meta/enrichment"

    # Rank preranked GSEA by the pooled z-statistic rather than the
    # pooled log2FC: z divides the effect by its pooled SE, so a large
    # effect resting on one noisy study does not outrank a
    # well-estimated smaller one.
    rank_column = 'z_stat'

    def __init__(self, parent: Optional[QWidget] = None):
        self._ctx = _MetaEnrichmentContext()
        self._project_folder = None
        super().__init__(self._ctx)

    def set_project_folder(self, folder):
        """Where the run's provenance sidecar and outputs live."""
        self._project_folder = folder
        self._ctx.project_dir = folder

    @property
    def context(self):
        """The stand-in workspace, for signal wiring and progress bar."""
        return self._ctx

    @property
    def log_message(self):
        return self._ctx.log_message

    @property
    def progress_bar(self):
        return self._ctx.progress_bar

    @progress_bar.setter
    def progress_bar(self, bar):
        self._ctx.progress_bar = bar

    def set_meta_results(self, meta_df):
        """Hand over the pooled table after a discovery run."""
        self._ctx.set_meta_results(meta_df)
        self._update_gene_count()

    # -- provenance ----------------------------------------------------
    #
    # Recorded here rather than in the DE page so that page keeps its
    # own behaviour: the meta run writes into the project-level sidecar
    # that the Methods step renders, which the DE workspace has no part
    # in. Both finish handlers funnel through _record_enrichment.

    def _on_enrichr_finished(self, payload):
        super()._on_enrichr_finished(payload)
        self._record_enrichment('ora')

    def _on_gsea_finished(self, payload):
        super()._on_gsea_finished(payload)
        self._record_enrichment('gsea')

    def _record_enrichment(self, method):
        from kosmic.meta_analysis.io import record_meta_stage
        df = self._results_df
        n_terms = 0 if df is None else len(df)
        n_sig = 0
        if df is not None and 'FDR' in getattr(df, 'columns', []):
            n_sig = int((df['FDR'] < self._p_thresh.value()).sum())
        study, background, _ = self._get_sig_genes()
        record_meta_stage(self._project_folder, 'meta_enrichment', {
            'method':      method,
            'library':     self._library_combo.currentText(),
            'direction':   self._direction_combo.currentData(),
            'rank_by':     self.rank_column if method == 'gsea' else None,
            'fdr':         self._p_thresh.value(),
            'gene_fdr':    self._fdr_spin.value(),
            'gene_lfc':    self._lfc_spin.value(),
            'min_genes':   self._min_genes.value(),
            'n_genes_in':  len(study or []),
            'n_background': len(background or []),
            'n_terms':     n_terms,
            'n_terms_sig': n_sig,
        })
