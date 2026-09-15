# Background workers for the Differential Expression workspace.
#
# Every long computation the DE pages start runs here, off the GUI
# thread, as a 'BaseWorker' launched through 'run_worker'. Each worker
# is a thin wrapper that calls into the analysis package ('kosmic.de')
# and reports progress and a result back to the page; none of them
# contains analysis logic of its own.
#
# - 'MetabolicDEWorker'   -- pseudobulk gene-level DE (Welch or DESeq2)
# - 'PathwayScoringWorker' -- pathway-level DE on pseudobulk scores
# - 'PathwayDataWorker'   -- loads the gene sets a pathway run needs
# - 'CellTypeBatchWorker' -- the per-cell-type batch run ('kosmic.de.batch')
# - 'GOEnrichmentWorker', 'EnrichrORAWorker' -- enrichment of a gene list
# - 'SpilloverWorker'     -- the ambient-spillover diagnostic for one pathway
# - 'MetaExportWorker'    -- writes the DE result + pseudobulk where the
#                            Meta-Analysis workspace will find them

from pathlib import Path
import numpy as np
import pandas as pd

from kosmic.gui.shared.theme import get_font_sizes
from kosmic.gui.shared.widgets import BaseWorker
from kosmic.paths import de_result_path, de_stats_dir, pseudobulk_path
from kosmic import DEFAULT_FDR


class MetabolicDEWorker(BaseWorker):
    """
    Worker thread for pseudobulk DE analysis (Welch's t-test or DESeq2).

    Emits 'finished_ok' with a tuple
    '(de_results, significant_genes, pathway_coverage, sample_df, message)'.
    """

    def __init__(self, adata, sample_col: str, condition_col: str,
                 control_label: str, disease_label: str,
                 pathway_gene_sets: dict, min_cells: int = 10,
                 min_expressing_samples: int = 2, de_method: str = 'ttest',
                 full_genome: bool = False, moderate: bool = False,
                 unit: str = 'sample', counts_layer=None, detection_min_pct=None,
                 detection_min_donor_frac=None, detection_study_col=None,
                 fdr_genes=None, deseq2_independent_filter=None,
                 deseq2_cooks_filter=None, filter_min_count=None,
                 filter_min_samples=None, covariates=None, min_counts=0):
        super().__init__()
        self.adata = adata
        self.sample_col = sample_col
        self.condition_col = condition_col
        self.control_label = control_label
        self.disease_label = disease_label
        self.pathway_gene_sets = pathway_gene_sets
        self.min_cells = min_cells
        self.min_counts = min_counts
        self.min_expressing_samples = min_expressing_samples
        self.de_method = de_method
        self.full_genome = full_genome
        self.moderate = moderate
        self.unit = unit
        self.counts_layer = counts_layer
        self.detection_min_pct = detection_min_pct
        self.detection_min_donor_frac = detection_min_donor_frac
        self.detection_study_col = detection_study_col
        self.fdr_genes = fdr_genes
        self.deseq2_independent_filter = deseq2_independent_filter
        self.deseq2_cooks_filter = deseq2_cooks_filter
        self.filter_min_count = filter_min_count
        self.filter_min_samples = filter_min_samples
        self.covariates = list(covariates or ())

    def _run(self):
        from kosmic.de.de_analysis import run_de_pipeline

        mode = "full-genome" if self.full_genome else "pathway-only"
        engine = ("cell-level Wilcoxon" if self.unit == 'cell'
                  else f"method={self.de_method}")
        adjust = (f", adjusting for {', '.join(self.covariates)}"
                  if self.covariates else "")
        self.progress.emit(
            f"Running DE pipeline ({mode}, {engine}{adjust})...")
        if self.fdr_genes:
            self.progress.emit(
                f"Hypothesis mode: FDR corrected over {len(self.fdr_genes):,} "
                f"committed genes; the fit stays genome-wide")
        if self.detection_study_col:
            self.progress.emit(
                f"Detection stratified by '{self.detection_study_col}' -- a "
                f"gene must clear the threshold in every one of them")
        self.progress_pct.emit(10)

        de_results, significant_genes, pathway_coverage, sample_df = run_de_pipeline(
            self.adata,
            sample_col=self.sample_col,
            condition_col=self.condition_col,
            pathway_gene_sets=self.pathway_gene_sets,
            min_cells=self.min_cells,
            min_counts=self.min_counts,
            min_expressing_samples=self.min_expressing_samples,
            de_method=self.de_method,
            full_genome=self.full_genome,
            moderate=self.moderate,
            unit=self.unit,
            counts_layer=self.counts_layer,
            fdr_genes=self.fdr_genes,
            covariates=self.covariates,
            **({} if self.deseq2_independent_filter is None
               else {'deseq2_independent_filter': self.deseq2_independent_filter}),
            **({} if self.deseq2_cooks_filter is None
               else {'deseq2_cooks_filter': self.deseq2_cooks_filter}),
            **({} if self.detection_min_pct is None
               else {'detection_min_pct': self.detection_min_pct}),
            **({} if self.detection_min_donor_frac is None
               else {'detection_min_donor_frac': self.detection_min_donor_frac}),
            detection_study_col=self.detection_study_col,
            **({} if self.filter_min_count is None
               else {'filter_min_count': self.filter_min_count}),
            **({} if self.filter_min_samples is None
               else {'filter_min_samples': self.filter_min_samples}),
            progress_callback=lambda msg: self.progress.emit(msg),
        )

        self.progress_pct.emit(90)

        # Log coverage summary
        for pathway_name, cov in pathway_coverage.items():
            self.progress.emit(
                f"  {pathway_name}: {cov['available_genes']}/{cov['total_genes']} genes "
                f"({cov['coverage_pct']:.1f}%)"
            )

        if not sample_df.empty:
            for condition in [self.control_label, self.disease_label]:
                n = (sample_df['condition'] == condition).sum()
                self.progress.emit(f"  {condition}: {n} samples")

        if de_results.empty:
            raise RuntimeError("No genes passed filtering")

        self.progress_pct.emit(100)
        n_sig = len(significant_genes)
        # Two criteria, and only the first carries an error rate: the FDR
        # tests "this gene changed", not "it changed by more than 0.25".
        # Saying "changed" rather than "significant" keeps the sentence
        # to what the statistics support.
        self.progress.emit(
            f"Analysis complete: {n_sig} genes changed "
            "(FDR<0.05), of which the estimated |log2FC| exceeds 0.25")

        return (de_results, significant_genes, pathway_coverage, sample_df,
                "Analysis complete")










class PathwayScoringWorker(BaseWorker):
    """
    Worker thread for pathway-level DE using pseudobulk scores.

    Computes per-sample pathway scores from raw counts, then runs the
    same run_pseudobulk_de() used for gene DE. Output format is identical
    to gene DE (names, logfoldchanges, se, pvals) so it feeds directly
    into meta-analysis.

    Emits 'finished_ok' with a tuple '(pathway_de, pathway_de_copy, message)'.
    """

    def __init__(self, adata, sample_col: str, condition_col: str,
                 control_cond: str, disease_cond: str,
                 pathway_gene_sets: dict, output_dir: Path,
                 level_label: str = '', method: str = 'score_genes',
                 gene_de_results=None, rho: float = 0.1, counts_layer=None,
                 gene_detection_pct=None, min_coverage=None):
        super().__init__()
        self.adata = adata
        self.sample_col = sample_col
        self.condition_col = condition_col
        self.control_cond = control_cond
        self.disease_cond = disease_cond
        self.pathway_gene_sets = pathway_gene_sets
        self.output_dir = output_dir
        self.level_label = level_label
        self.method = method
        self.gene_de_results = gene_de_results
        self.rho = rho
        self.counts_layer = counts_layer
        self.gene_detection_pct = gene_detection_pct
        self.min_coverage = min_coverage

    def _coverage_kwargs(self):
        """Only override compute_pathway_scores defaults when the GUI set a value."""
        kw = {}
        if self.gene_detection_pct is not None:
            kw['gene_detection_pct'] = self.gene_detection_pct
        if self.min_coverage is not None:
            kw['min_coverage'] = self.min_coverage
        return kw

    def _run(self):
        if self.method == 'deseq2':
            return self._run_deseq2(normalization='deseq2')
        if self.method == 'deseq2_vst':
            return self._run_deseq2(normalization='vst')
        if self.method == 'deseq2_vst_zscore':
            return self._run_deseq2(normalization='vst', standardize_genes=True)
        if self.method == 'pseudobulk_cpm':
            return self._run_pseudobulk_cpm()
        return self._run_cell_level_scoring()

    def _run_deseq2(self, normalization='deseq2', standardize_genes=False):
        """Pseudobulk sum -> DESeq2 normalisation -> per-donor pathway mean
        -> plain Welch t-test. Effect = average per-gene log2FC (poolable).
        'normalization' is 'deseq2' (size-factor log2) or 'vst'."""
        from kosmic.de.de_analysis import (
            compute_pathway_scores, run_pseudobulk_de, annotate_de_results,
        )

        if standardize_genes:
            label = 'VST gene-standardised'
        elif normalization == 'vst':
            label = 'VST'
        else:
            label = 'DESeq2 log2FC'
        self.progress.emit(f"Computing pathway scores ({label}, per donor)...")

        score_matrix, sample_df, pathway_names = compute_pathway_scores(
            self.adata, self.pathway_gene_sets,
            self.sample_col, self.condition_col,
            normalization=normalization, counts_layer=self.counts_layer,
            standardize_genes=standardize_genes,
            **self._coverage_kwargs(),
        )

        if score_matrix.size == 0 or not pathway_names:
            raise RuntimeError("No pathways with sufficient gene coverage")

        self._log_dropped(set(pathway_names))
        self.progress.emit(f"Testing {len(pathway_names)} pathways...")

        pathway_de = run_pseudobulk_de(
            score_matrix, sample_df, pathway_names,
            pre_transformed=True, moderate=False,
        )

        if pathway_de.empty:
            raise RuntimeError("No pathways passed filtering")

        pathway_coverage = self._build_pathway_coverage()
        pathway_de, _ = annotate_de_results(pathway_de, pathway_coverage)

        return self._save_and_emit(pathway_de, pathway_names, scores={
            'matrix': score_matrix, 'sample_df': sample_df, 'names': pathway_names},
            coverage_report=self._coverage_report())

    def _run_pseudobulk_cpm(self):
        """Pseudobulk sum → CPM (full-genome lib sizes) → log2 → t-test."""
        from kosmic.de.de_analysis import (
            compute_pathway_scores, run_pseudobulk_de, annotate_de_results,
        )

        self.progress.emit("Computing pathway scores (pseudobulk CPM-log2)...")

        score_matrix, sample_df, pathway_names = compute_pathway_scores(
            self.adata, self.pathway_gene_sets,
            self.sample_col, self.condition_col,
            counts_layer=self.counts_layer,
            **self._coverage_kwargs(),
        )

        if score_matrix.size == 0 or not pathway_names:
            raise RuntimeError("No pathways with sufficient gene coverage")

        self._log_dropped(set(pathway_names))

        self.progress.emit(f"Testing {len(pathway_names)} pathways...")

        pathway_de = run_pseudobulk_de(
            score_matrix, sample_df, pathway_names,
            pre_transformed=True, moderate=False,
        )

        if pathway_de.empty:
            raise RuntimeError("No pathways passed filtering")

        pathway_coverage = self._build_pathway_coverage()
        pathway_de, _ = annotate_de_results(pathway_de, pathway_coverage)

        return self._save_and_emit(pathway_de, pathway_names, scores={
            'matrix': score_matrix, 'sample_df': sample_df, 'names': pathway_names},
            coverage_report=self._coverage_report())

    def _run_cell_level_scoring(self):
        """Cell-level scoring (score_genes, mean, zscore, singscore) → sample aggregation → t-test."""
        from kosmic.de.pathway_scoring import score_pathways, test_pathway_scores
        from kosmic.de.de_analysis import annotate_de_results

        self.progress.emit(f"Computing per-cell pathway scores ({self.method})...")

        score_adata, available_pathways = score_pathways(
            self.adata, self.pathway_gene_sets, method=self.method,
        )

        if not available_pathways:
            raise RuntimeError("No pathways with sufficient gene coverage")

        self._log_dropped(set(available_pathways))

        self.progress.emit(f"Testing {len(available_pathways)} pathways...")

        stats_df = test_pathway_scores(
            score_adata, available_pathways,
            self.sample_col, self.condition_col,
            self.control_cond, self.disease_cond,
        )

        if stats_df.empty:
            raise RuntimeError("No pathways passed filtering")

        # Convert to standard DE format for the Pathway Activity plot
        pathway_de = pd.DataFrame({
            'names': stats_df['Pathway'],
            'logfoldchanges': stats_df['Effect_Size'],
            'pvals': stats_df['P_Value'],
            'pvals_adj': stats_df['P_Adjusted'] if 'P_Adjusted' in stats_df.columns else stats_df['P_Value'],
        })

        pathway_coverage = self._build_pathway_coverage()
        pathway_de, _ = annotate_de_results(pathway_de, pathway_coverage)

        return self._save_and_emit(pathway_de, list(available_pathways.keys()))

    def _log_dropped(self, kept_names):
        """Surface pathways that were silently filtered out by the scorer."""
        dropped = sorted(set(self.pathway_gene_sets) - set(kept_names))
        if not dropped:
            return
        preview = ", ".join(dropped[:5])
        if len(dropped) > 5:
            preview += f", ... (+{len(dropped) - 5} more)"
        self.progress.emit(
            f"Skipped {len(dropped)} pathways with insufficient coverage: {preview}"
        )

    def _build_pathway_coverage(self):
        from kosmic.scrna.inspect.detection import (
            detect_species, format_gene_for_species,
        )

        pathway_coverage = {}
        var_names_seq = (self.adata.raw.var_names if self.adata.raw is not None
                         else self.adata.var_names)
        var_names = set(var_names_seq)
        # Convert human pathway symbols to the dataset's species namespace,
        # matching what kosmic.de.de_analysis.prepare_gene_coverage does so
        # the scorer and the gene-DE engine agree on coverage.
        species = detect_species(list(var_names_seq))
        if species == 'mouse':
            def _match(g):
                return format_gene_for_species(g, 'mouse')
        else:
            def _match(g):
                return g

        for pw_name, genes in self.pathway_gene_sets.items():
            available = [g for g in genes if _match(g) in var_names]
            pathway_coverage[pw_name] = {
                'genes': available,
                'total': len(genes),
                'coverage': len(available) / len(genes) if genes else 0,
            }
        return pathway_coverage

    def _coverage_report(self):
        """Per-pathway detection/coverage table matching the scored set."""
        from kosmic.de.de_analysis import pathway_coverage_report
        return pathway_coverage_report(
            self.adata, self.pathway_gene_sets,
            counts_layer=self.counts_layer, **self._coverage_kwargs())

    def _save_and_emit(self, pathway_de, pathway_names, scores=None,
                       coverage_report=None):
        n_sig = (pathway_de['pvals_adj'] < DEFAULT_FDR).sum() if 'pvals_adj' in pathway_de.columns else 0
        self.progress.emit(
            f"Pathway DE complete: {len(pathway_de)} pathways, {n_sig} significant"
        )

        scoring_dir = self.output_dir / 'pathway_scoring'
        scoring_dir.mkdir(parents=True, exist_ok=True)
        pathway_de.to_csv(
            scoring_dir / f'{self.level_label}pathway_de_results.csv', index=False
        )
        if coverage_report is not None and not coverage_report.empty:
            coverage_report.to_csv(
                scoring_dir / f'{self.level_label}pathway_coverage.csv', index=False)

        # Persist the per-donor score matrix so figure export reloads the exact
        # values this run produced (any scoring method) rather than recomputing.
        if scores is not None:
            from kosmic.de.pathway_scoring import (
                save_per_donor_scores, PER_DONOR_SCORES_FILE,
            )
            save_per_donor_scores(
                scores, scoring_dir / f'{self.level_label}{PER_DONOR_SCORES_FILE}')

        # 'scores' carries the per-donor score matrix so the plot can reuse
        # it instead of recomputing (None for methods that don't produce one);
        # 'coverage_report' is the per-pathway detection/coverage table.
        return (pathway_de, scores, coverage_report, "Pathway DE complete")









class MetaExportWorker(BaseWorker):
    """Export the chosen study's DE result + pseudobulk for meta-analysis.

    The primary DE is already computed, so this writes it to the meta-analysis
    CSV layout and saves the per-donor pseudobulk (read back by the consensus
    case-control permutation). It does not recompute other DE methods. Emits
    'finished_ok' with the completion message.
    """

    def __init__(self, adata, sample_col, condition_col, min_cells,
                 primary_method, primary_moderate, primary_de_results,
                 output_dir, gse_id, min_counts=0):
        super().__init__()
        self.adata = adata
        self.sample_col = sample_col
        self.condition_col = condition_col
        self.min_cells = min_cells
        self.min_counts = min_counts
        self.primary_method = primary_method
        self.primary_moderate = primary_moderate
        self.primary_de_results = primary_de_results
        self.output_dir = Path(output_dir)
        self.gse_id = gse_id

    def _run(self):
        from kosmic.de.batch import method_label, pseudobulk_frame
        from kosmic.de.de_analysis import create_pseudobulk

        stats_dir = de_stats_dir(self.output_dir)
        stats_dir.mkdir(parents=True, exist_ok=True)
        self.progress_pct.emit(0)

        # Export the chosen method's results (already computed) in the layout
        # the meta-analysis loads.
        label = method_label(self.primary_method, self.primary_moderate)
        if self.primary_de_results is not None and not self.primary_de_results.empty:
            out_path = de_result_path(self.output_dir, self.gse_id, label)
            self.primary_de_results.to_csv(out_path, index=False)
            self.progress.emit(f"[{label}] Exported DE results for meta-analysis")

        # Per-donor pseudobulk for the consensus case-control permutation.
        self.progress_pct.emit(40)
        try:
            if (self.primary_de_results is not None
                    and 'names' in self.primary_de_results.columns):
                gene_names = [str(g) for g in self.primary_de_results['names'].tolist()]
            else:
                gene_names = [str(g) for g in self.adata.var_names]
            self.progress.emit(f"Creating pseudobulk ({len(gene_names)} genes)...")
            pb_matrix, pb_sample_df, pb_genes = create_pseudobulk(
                self.adata, gene_names, self.sample_col, self.condition_col,
                min_cells=self.min_cells, min_counts=self.min_counts,
                aggregate='sum')

            # Drop the excluded arm. The DE dropped it, so leaving it here
            # hands the consensus case-control permutation and
            # leave-one-out donors the contrast never used -- 11 of 25 on
            # GSE292067. The per-cell-type runner already excludes them,
            # so this also stops the two exports disagreeing.
            if pb_matrix.size > 0 and 'role' in pb_sample_df.columns:
                keep = pb_sample_df['role'].astype(str).str.lower().isin(
                    ('disease', 'control')).values
                if keep.any() and not keep.all():
                    n_dropped = int((~keep).sum())
                    pb_matrix = pb_matrix[keep]
                    pb_sample_df = pb_sample_df[keep].reset_index(drop=True)
                    self.progress.emit(
                        f"Excluded {n_dropped} donor(s) marked 'exclude' from "
                        f"the pseudobulk export")

            if pb_matrix.size > 0:
                # Shared with the per-cell-type batch runner so both write
                # the same schema -- including the 'role' column the
                # consensus case-control permutation reads back.
                pb_df = pseudobulk_frame(pb_matrix, pb_sample_df, pb_genes)
                pb_path = pseudobulk_path(self.output_dir, self.gse_id)
                pb_df.to_csv(pb_path)
                self.progress.emit(
                    f"Saved pseudobulk: {pb_path.name} "
                    f"({pb_matrix.shape[0]} samples x {pb_matrix.shape[1]} genes)")
        except Exception as e:
            self.progress.emit(f"Warning: could not save pseudobulk: {e}")
            import traceback
            traceback.print_exc()

        self.progress_pct.emit(100)
        return "Meta-analysis export complete"


class CellTypeBatchWorker(BaseWorker):
    """Run the DE pipeline once per cell type, writing per-type results.

    The single-run worker above is for one contrast on whatever is
    loaded. This one loops cell types, because a per-cell-type
    meta-vs-mega design needs the same contrast run dozens of times and
    doing that by hand -- subset, run, export, reload -- is where the
    mistakes come from.

    All the statistics live in 'kosmic.de.batch'; this is the thread
    wrapper. Emits 'finished_ok' with '(runs, summary_df, message)'.
    """

    def __init__(self, adata, cell_type_col, sample_col, condition_col,
                 output_dir, accession, cell_types=None,
                 de_method='deseq2', moderate=False, covariates=None,
                 min_cells=10, pathway_gene_sets=None, fdr_genes=None,
                 min_counts=0, **pipeline_kwargs):
        super().__init__()
        self.adata = adata
        self.cell_type_col = cell_type_col
        self.sample_col = sample_col
        self.condition_col = condition_col
        self.output_dir = Path(output_dir)
        self.accession = accession
        self.cell_types = list(cell_types) if cell_types else None
        self.de_method = de_method
        self.moderate = moderate
        self.covariates = list(covariates or ())
        self.min_cells = min_cells
        self.min_counts = min_counts
        self.pathway_gene_sets = pathway_gene_sets
        self.fdr_genes = fdr_genes
        self.pipeline_kwargs = pipeline_kwargs
        self._cancelled = False

    def cancel(self):
        """Stop after the cell type currently running finishes."""
        self._cancelled = True

    def _run(self):
        from kosmic.de.batch import batch_summary, run_de_by_cell_type

        n_types = len(self.cell_types) if self.cell_types else 0
        adjust = (f", adjusting for {', '.join(self.covariates)}"
                  if self.covariates else "")
        self.progress.emit(
            f"Batch DE over {n_types or 'all'} cell types "
            f"(method={self.de_method}{adjust})...")

        done = {'n': 0}

        def _progress(msg):
            self.progress.emit(msg)
            # Coarse progress: each cell type's headline message advances
            # the bar. Finer granularity would need the pipeline to
            # report its own fraction.
            if n_types and msg.startswith('['):
                head = msg.split(']', 1)[0]
                if ':' not in head:
                    done['n'] = min(done['n'] + 1, n_types)
                    self.progress_pct.emit(int(100 * done['n'] / n_types))

        runs = run_de_by_cell_type(
            self.adata, self.cell_type_col, self.sample_col,
            self.condition_col, self.output_dir, self.accession,
            cell_types=self.cell_types, de_method=self.de_method,
            moderate=self.moderate, covariates=self.covariates,
            min_cells=self.min_cells, min_counts=self.min_counts,
            pathway_gene_sets=self.pathway_gene_sets,
            fdr_genes=self.fdr_genes,
            progress_callback=_progress,
            cancelled=lambda: self._cancelled,
            **self.pipeline_kwargs)

        summary = batch_summary(runs)
        n_ok = int((summary['status'] == 'ok').sum()) if len(summary) else 0
        n_failed = int((summary['status'] == 'failed').sum()) if len(summary) else 0

        summary_path = de_stats_dir(self.output_dir) / (
            f"{self.accession}_by_cell_type_summary.csv")
        try:
            summary_path.parent.mkdir(parents=True, exist_ok=True)
            summary.to_csv(summary_path, index=False)
            self.progress.emit(f"Saved summary: {summary_path.name}")
        except Exception as exc:  # noqa: BLE001
            self.progress.emit(f"Warning: could not save summary: {exc}")

        self.progress_pct.emit(100)
        message = (f"Batch complete: {n_ok} cell type(s) analysed"
                   + (f", {n_failed} failed" if n_failed else "")
                   + (" (cancelled)" if self._cancelled else ""))
        return runs, summary, message


class GOEnrichmentWorker(BaseWorker):
    """
    Worker thread for GO term enrichment (elim algorithm).

    Wraps src.gene_sets.go_enrichment.run_enrichment().
    Emits 'finished_ok' with '(results_df, message)'.
    """

    def __init__(self, study_genes, background_genes,
                 ontology='BP', method='elim', alpha=DEFAULT_FDR,
                 min_genes=5):
        super().__init__()
        self.study_genes = study_genes
        self.background_genes = background_genes
        self.ontology = ontology
        self.method = method
        self.alpha = alpha
        self.min_genes = min_genes

    def _run(self):
        from kosmic.de.go_enrichment import run_enrichment

        df = run_enrichment(
            study_genes=self.study_genes,
            background_genes=self.background_genes,
            ontology=self.ontology,
            method=self.method,
            alpha=self.alpha,
            min_genes=self.min_genes,
            progress_callback=self._on_progress,
        )

        n_sig = (df['P_value'] < self.alpha).sum() if not df.empty else 0
        message = (
            f"GO enrichment complete: {n_sig} significant terms "
            f"from {len(df)} tested ({self.ontology}, {self.method})")
        return (df, message)

    def _on_progress(self, msg):
        self.progress.emit(msg)


class EnrichrORAWorker(BaseWorker):
    """
    Worker thread for over-representation analysis using Enrichr libraries.

    Downloads the gene set library from Enrichr, then runs Fisher's exact test
    for each term against the study gene list. BH correction applied.
    Emits 'finished_ok' with '(results_df, message)'.
    """

    def __init__(self, study_genes, background_genes, library_name,
                 alpha=DEFAULT_FDR, min_genes=5):
        super().__init__()
        self.study_genes = set(study_genes)
        self.background_genes = set(background_genes)
        self.library_name = library_name
        self.alpha = alpha
        self.min_genes = min_genes

    def _run(self):
        from kosmic.de.ora import run_ora
        from kosmic.reference.pathways.enrichr import fetch_library

        self.progress.emit(f"Downloading {self.library_name} from Enrichr...")
        library = fetch_library(self.library_name)
        if not library:
            raise RuntimeError(f"Empty library: {self.library_name}")

        n_study = len(set(self.study_genes) & set(self.background_genes))
        self.progress.emit(
            f"Testing {len(library)} terms ({n_study} study genes)..."
        )

        def _on_progress(i, total):
            self.progress_pct.emit(int(i / max(total, 1) * 100))

        df = run_ora(
            self.study_genes, self.background_genes, library,
            min_genes=self.min_genes, progress_callback=_on_progress,
        )

        if df.empty:
            return (df, "No enriched terms found.")

        n_sig = int((df['FDR'] < self.alpha).sum())
        message = (
            f"Enrichment complete: {n_sig} significant terms from "
            f"{len(df)} tested ({self.library_name})")
        return (df, message)


class PathwayDataWorker(BaseWorker):
    """
    Compute pseudobulk expression and pathway score for one pathway.

    Used by PathwayDEPage for pathway detail plots.

    Emits 'finished_ok' with a dict of results, or 'None' if nothing
    could be computed.
    """

    def __init__(self, adata, genes, sample_col, condition_col,
                 disease_label, control_label, pathway_name,
                 method='pseudobulk_cpm', counts_layer=None,
                 precomputed_scores=None, allowed_genes=None, parent=None):
        super().__init__(parent)
        self._adata = adata
        self._genes = genes
        self._sample_col = sample_col
        self._condition_col = condition_col
        self._disease_label = disease_label
        self._control_label = control_label
        self._pathway_name = pathway_name
        self._method = method
        self._counts_layer = counts_layer
        # (scores, sample_df) already computed by the pathway run; when set,
        # the per-donor score is reused rather than recomputed.
        self._precomputed_scores = precomputed_scores
        # DE-tested gene set; the per-gene heatmap restricts to it so it matches
        # the gene DE table (and the exported heatmap).
        self._allowed_genes = allowed_genes

    def _run(self):
        try:
            from kosmic.de.de_analysis import create_pseudobulk
            from scipy import sparse

            adata = self._adata
            genes = self._genes

            pb_matrix, sample_df, genes_used = create_pseudobulk(
                adata, genes, self._sample_col, self._condition_col,
                aggregate='mean', counts_layer=self._counts_layer,
            )
            if len(genes_used) == 0 or pb_matrix.size == 0:
                return None

            # Per-donor score: reuse the one the pathway run already computed;
            # only recompute when it wasn't retained (e.g. results loaded from
            # disk), which is when the expensive full-genome pseudobulk is
            # rebuilt.
            if self._precomputed_scores is not None:
                pathway_scores, score_sdf = self._precomputed_scores
            else:
                from kosmic.de.pathway_scoring import per_donor_pathway_score
                pathway_scores, score_sdf = per_donor_pathway_score(
                    adata, self._pathway_name, genes,
                    self._sample_col, self._condition_col, self._method,
                    counts_layer=self._counts_layer)
            if pathway_scores.size == 0:
                # Method couldn't score it (coverage) -> mean-expression fallback.
                pathway_scores = pb_matrix.mean(axis=1)
                score_sdf = sample_df
            disease_mask = (score_sdf['condition'] == self._disease_label).values
            control_mask = (score_sdf['condition'] == self._control_label).values

            # Per-gene summary uses pb_matrix (mean expression) with its own
            # condition masks aligned to pb_matrix rows.
            gene_dmask = (sample_df['condition'] == self._disease_label).values
            gene_cmask = (sample_df['condition'] == self._control_label).values
            gene_summary = []
            for i, gene in enumerate(genes_used):
                d_vals = pb_matrix[gene_dmask, i]
                c_vals = pb_matrix[gene_cmask, i]
                d_mean = d_vals.mean() if len(d_vals) else 0
                c_mean = c_vals.mean() if len(c_vals) else 0
                d_sem = d_vals.std(ddof=1) / np.sqrt(len(d_vals)) if len(d_vals) > 1 else 0
                c_sem = c_vals.std(ddof=1) / np.sqrt(len(c_vals)) if len(c_vals) > 1 else 0

                pct = 0.0
                if gene in adata.var_names:
                    idx = list(adata.var_names).index(gene)
                    col = adata.X[:, idx]
                    if sparse.issparse(col):
                        col = col.toarray().flatten()
                    pct = (col > 0).mean() * 100

                gene_summary.append({
                    'gene': gene,
                    'disease_mean': d_mean, 'control_mean': c_mean,
                    'disease_sem': d_sem, 'control_sem': c_sem,
                    'disease_samples': d_vals.tolist(),
                    'control_samples': c_vals.tolist(),
                    'log2fc': np.log2((d_mean + 0.01) / (c_mean + 0.01)),
                    'pct_expressing': pct,
                })

            # Heatmap
            heatmap_bytes = self._render_heatmap(adata, genes_used)

            return {
                'pathway_name': self._pathway_name,
                'method': self._method,
                'genes_used': genes_used,
                'pb_matrix': pb_matrix,
                'sample_df': sample_df,
                'pathway_scores': pathway_scores,
                'gene_summary': gene_summary,
                'heatmap_bytes': heatmap_bytes,
                'disease_label': self._disease_label,
                'control_label': self._control_label,
                'disease_mask': disease_mask,
                'control_mask': control_mask,
            }
        except Exception:
            import traceback
            traceback.print_exc()
            return None

    def _render_heatmap(self, adata, genes_used):
        """Render heatmap using create_gene_heatmap."""
        from kosmic.visualisation.de.gene_heatmap import prepare_heatmap_data, create_gene_heatmap
        from kosmic.gui.shared.theme import get_color
        import io

        pw_sets = {self._pathway_name: genes_used}

        try:
            sample_df, available = prepare_heatmap_data(
                adata, pw_sets,
                self._sample_col, self._condition_col,
                allowed_genes=self._allowed_genes,
            )
        except Exception:
            return None

        pw_genes = available.get(self._pathway_name)
        if not pw_genes or sample_df is None or sample_df.empty:
            return None

        theme_colors = {
            'bg_primary': get_color('bg_primary'),
            'fg_primary': get_color('fg_primary'),
            'fg_secondary': get_color('fg_secondary'),
        }

        fig = create_gene_heatmap(
            sample_df, pw_genes, self._pathway_name,
            self._control_label, self._disease_label,
            theme_colors=theme_colors,
            font_sizes=get_font_sizes(),
        )
        if fig is None:
            return None

        import matplotlib.pyplot as plt
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=130, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()


class SpilloverWorker(BaseWorker):
    """Run the ambient-spillover diagnostic for one pathway off-thread.

    Thin wrapper around ``kosmic.de.spillover.pathway_spillover``. Emits
    'finished_ok' with its result dict, or raises (handled via on_failed).
    """

    def __init__(self, adata, pathway_name, genes, sample_col, condition_col,
                 disease_label, control_label, contam_key,
                 method='pseudobulk_cpm', counts_layer=None, parent=None):
        super().__init__(parent)
        self._adata = adata
        self._pathway_name = pathway_name
        self._genes = genes
        self._sample_col = sample_col
        self._condition_col = condition_col
        self._disease_label = disease_label
        self._control_label = control_label
        self._contam_key = contam_key
        self._method = method
        self._counts_layer = counts_layer

    def _run(self):
        from kosmic.de.spillover import pathway_spillover
        self.progress.emit(
            f"Testing '{self._pathway_name}' against {self._contam_key}...")
        return pathway_spillover(
            self._adata, self._pathway_name, self._genes,
            self._sample_col, self._condition_col,
            self._disease_label, self._control_label, self._contam_key,
            method=self._method, counts_layer=self._counts_layer,
        )


class FgseaWorker(BaseWorker):
    """Run preranked GSEA (fgsea) on the DESeq2 gene ranks off-thread.

    Emits 'finished_ok' with the fgsea result DataFrame (names, nes, es,
    pvals, pvals_adj, leading_edge).
    """

    def __init__(self, gene_de_results, pathway_gene_sets,
                 rank_by='logfoldchanges', parent=None):
        super().__init__(parent)
        self._gene_de = gene_de_results
        self._gene_sets = pathway_gene_sets
        self._rank_by = rank_by

    def _run(self):
        from kosmic.de.fgsea import run_fgsea
        self.progress.emit("Running preranked GSEA (fgsea)...")
        return run_fgsea(self._gene_de, self._gene_sets, rank_by=self._rank_by,
                         return_details=True)


class LibraryGseaWorker(BaseWorker):
    """Preranked GSEA of a fetched Enrichr library against the ranked gene DE.

    Fetches the library's gene sets (network) then runs fgsea, so the
    unbiased 'which pathways are enriched genome-wide' question is answered
    with the same rank-based method as the a-priori panel. Emits
    'finished_ok' with '(summary_df, message)'.
    """

    def __init__(self, gene_de_results, library_name,
                 rank_by='logfoldchanges', parent=None):
        super().__init__(parent)
        self._gene_de = gene_de_results
        self._library_name = library_name
        self._rank_by = rank_by

    def _run(self):
        import gseapy as gp
        from kosmic.de.fgsea import run_fgsea

        self.progress.emit(f"Fetching {self._library_name}...")
        gene_sets = gp.get_library(name=self._library_name, organism='Human')
        if not gene_sets:
            raise RuntimeError(f"Could not fetch library '{self._library_name}'.")

        self.progress.emit(
            f"Preranked GSEA over {len(gene_sets):,} gene sets...")
        res, details = run_fgsea(
            self._gene_de, gene_sets, rank_by=self._rank_by,
            return_details=True)
        # Keep running-enrichment curves only for the most-enriched sets (what a
        # mountain-plot figure would ever show), so a big library does not pin
        # hundreds of full-length curves in memory.
        if details and not res.empty:
            keep = set(res.reindex(res['nes'].abs().sort_values(
                ascending=False).index)['names'].astype(str).head(40))
            details = {k: v for k, v in details.items() if k in keep}
        n_sig = int((res['pvals_adj'] < DEFAULT_FDR).sum()) if not res.empty else 0
        return res, details, (f"GSEA complete -- {len(res):,} gene sets, "
                              f"{n_sig} at FDR < {DEFAULT_FDR}")
