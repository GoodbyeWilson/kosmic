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

from pathlib import Path
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
        self._selection = None
        super().__init__(self._ctx)

    def set_project_folder(self, folder, selection=None):
        """Where the run's provenance sidecar and outputs live.

        'selection' is the output folder of the pooling run the gene list
        came from (ADR-007).
        """
        self._project_folder = folder
        self._selection = selection

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

    def _on_go_finished(self, payload):
        super()._on_go_finished(payload)
        self._record_enrichment('go')

    def _on_enrichr_finished(self, payload):
        super()._on_enrichr_finished(payload)
        self._record_enrichment('ora')

    def _on_gsea_finished(self, payload):
        super()._on_gsea_finished(payload)
        self._record_enrichment('gsea')

    def _gsea_results_path(self):
        """The selection's meta-analysis folder, one file per library.

        The DE page's location is a study folder; the project root is not
        one, and a single file name let each library's run replace the
        last while the record named the newest.
        """
        if not self._project_folder:
            return None
        from kosmic.paths import meta_output_dir
        from kosmic.de.fgsea import GSEA_LIBRARY_RESULTS_FILE
        stem = Path(GSEA_LIBRARY_RESULTS_FILE).stem
        return (meta_output_dir(self._project_folder, self._selection)
                / f"{stem}_{self._gsea_library_id}.csv")

    def _record_enrichment(self, method):
        """One entry per method and library in the selection's record:
        rerunning a library replaces its entry, another library adds one."""
        from kosmic.meta_analysis.io import record_meta_stage
        params = self._enrichment_params(method)
        if method == 'gsea':
            params['library'] = self.ws.gsea_library_name
        results_file = None
        if method == 'gsea' and params['n_terms']:
            path = self._gsea_results_path()
            if path is not None:
                results_file = path.relative_to(self._project_folder).as_posix()
        params['results_file'] = results_file
        record_meta_stage(
            self._project_folder, 'meta_enrichment', params,
            selection=self._selection,
            replace_matching={'method': method, 'library': params['library']})
