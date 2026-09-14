# Background workers for the Discovery Meta-Analysis page and the
# Validation step.
#
# 'GeneMAWorker' runs the selected pooling methods over the studies'
# DE results, with one shared case/control permutation calibration
# across methods. 'ReproducibilityWorker' computes cross-dataset
# reproducibility; 'ConsensusLOOWorker' and 'LOOBenchmarkWorker' run
# the leave-one-study-out replication for the consensus and for each
# constituent method. All call into 'kosmic.meta_analysis' and carry no
# analysis logic themselves.
from __future__ import annotations

from PyQt6.QtCore import pyqtSignal

from kosmic import DEFAULT_FDR
from kosmic.gui.shared.widgets import BaseWorker


class ReproducibilityWorker(BaseWorker):
    """Run cross-dataset reproducibility analysis off the GUI thread."""

    def __init__(self, datasets, labels):
        super().__init__()
        self.datasets = datasets
        self.labels = labels

    def _run(self):
        from kosmic.meta_analysis.reproducibility import run_reproducibility
        return run_reproducibility(
            self.datasets, self.labels, threshold=DEFAULT_FDR,
            progress_cb=lambda m: self.progress.emit(m))


class GeneMAWorker(BaseWorker):
    """
    Run N pooling methods with shared CC permutation.

    Emits 'finished_ok' with a dict '{meta_df, method_dfs, method_keys}'.
    """

    def __init__(self, datasets, labels, params):
        super().__init__()
        self.datasets = datasets
        self.labels = labels
        self.params = params

    def _run(self):
        method_keys = self.params['method_keys']
        calibration = self.params.get('calibration', 'analytical')
        min_st = self.params.get('min_studies', 3)
        n_cc_perms = self.params.get('cc_n_perms', 1000)
        n_methods = len(method_keys)

        # Step 1: Analytical pooling for each method
        from kosmic.meta_analysis.pooling_dispatch import (
            get_analytical_pool_fn)
        analytical_dfs = {}
        for i, key in enumerate(method_keys):
            self.progress.emit(
                f"Analytical pooling {i+1}/{n_methods}: {key}...")
            pool_fn = get_analytical_pool_fn(key, min_st)
            meta_df = pool_fn(self.datasets, min_studies=min_st)
            analytical_dfs[key] = meta_df

        # Step 2: CC permutation calibration (if selected)
        if calibration == 'cc_perm':
            pb_paths = self.params.get('pseudobulk_paths', [])
            if not pb_paths:
                raise RuntimeError(
                    "No pseudobulk data. Run DE with "
                    "'Save pseudobulk' first.")

            pb_data = self._load_pseudobulk(pb_paths)
            if not pb_data:
                raise RuntimeError("Failed to load pseudobulk data.")

            from kosmic.meta_analysis.cc_permutation import (
                consensus_cc_permutation)

            cc_de_method = self.params.get('de_method', 'deseq2_fast')
            self.progress.emit(
                f"CC permutation (shared DE, {cc_de_method})...")
            method_dfs = consensus_cc_permutation(
                pb_data, analytical_dfs,
                n_perms=n_cc_perms, seed=42,
                de_method=cc_de_method,
                min_studies=min_st,
                progress_cb=lambda cur, tot: self.progress.emit(
                    f"CC permutation {cur}/{tot}"),
                status_cb=lambda msg: self.progress.emit(msg),
            )
            self.progress.emit("CC permutation complete")
        else:
            method_dfs = analytical_dfs

        # Step 3: Merge (computes per-method FDR + consensus FDR =
        # max of per-method FDR, so a gene passes consensus iff it
        # passes every individual method's FDR threshold).
        self.progress.emit("Merging consensus results...")
        from kosmic.meta_analysis.pooling_dispatch import (
            merge_consensus_methods)
        consensus_df = merge_consensus_methods(method_keys, method_dfs)

        self.progress.emit(
            f"Pooling complete: {len(consensus_df)} genes")

        return {
            'meta_df': consensus_df,
            'method_dfs': method_dfs,
            'method_keys': method_keys,
        }

    def _load_pseudobulk(self, pb_paths):
        from kosmic.meta_analysis.io import load_pseudobulk_set
        return load_pseudobulk_set(
            pb_paths, labels=self.labels, progress_cb=self.progress.emit)


class ConsensusLOOWorker(BaseWorker):
    """
    Held-out replication LOO across all studies.

    Re-runs the consensus pipeline once per fold (K folds), holding
    out one study each time, then tests whether each fold's consensus
    genes are independently replicated by the held-out study.
    """

    fold_progress = pyqtSignal(int, int, str)  # (fold_1based, n_folds, study_name)

    def __init__(self, datasets, labels, params,
                 replication_mode='strict', pb_paths=None):
        super().__init__()
        self.datasets = datasets
        self.labels = labels
        self.params = params
        self.replication_mode = replication_mode
        self.pb_paths = pb_paths or []

    def _run(self):
        from kosmic.meta_analysis.consensus_loo import run_consensus_loo

        K = len(self.datasets)
        if K < 3:
            raise RuntimeError(
                f"LOO validation requires at least 3 studies; got K={K}.")

        # Load pseudobulk if calibration is CC perm
        pseudobulk_data = None
        if self.params.get('calibration') == 'cc_perm':
            if not self.pb_paths:
                raise RuntimeError(
                    "CC perm calibration requires pseudobulk CSVs. "
                    "Run a normal consensus first so the pseudobulk "
                    "paths are cached.")
            self.progress.emit(
                f"Loading pseudobulk for {len(self.pb_paths)} studies...")
            pseudobulk_data = self._load_pseudobulk(self.pb_paths)
            if not pseudobulk_data:
                raise RuntimeError("Failed to load pseudobulk data.")

        self.progress.emit(
            f"LOO: K={K} folds, mode={self.replication_mode}, "
            f"calibration={self.params.get('calibration', 'analytical')}")

        def _fold_cb(j, n, name):
            self.fold_progress.emit(j, n, name)
            self.progress.emit(f"Fold {j}/{n}: held out {name}")

        def _stage_cb(msg):
            self.progress.emit(msg)

        result = run_consensus_loo(
            self.datasets, self.labels, self.params,
            replication_mode=self.replication_mode,
            pseudobulk_data=pseudobulk_data,
            fold_progress_cb=_fold_cb,
            stage_progress_cb=_stage_cb,
        )

        self.progress.emit(
            f"LOO complete: avg replication "
            f"{result['avg_replication_pct']:.1f}% "
            f"(range {result['min_replication_pct']:.1f}-"
            f"{result['max_replication_pct']:.1f}%)")
        return result

    def _load_pseudobulk(self, pb_paths):
        """
        Load pseudobulk for every fold. Returns 'None' on the first
        per-path failure so the caller can fail-fast (LOO cannot proceed
        with partial data because every fold needs every study).
        """
        from kosmic.meta_analysis.io import load_pseudobulk_set
        try:
            return load_pseudobulk_set(pb_paths, labels=self.labels)
        except Exception as e:  # noqa: BLE001
            self.progress.emit(f"  pseudobulk load failed: {e}")
            return None


class LOOBenchmarkWorker(BaseWorker):
    """
    Run held-out replication LOO for each individual constituent
    method plus the consensus, to answer the question 'is the
    consensus actually better than its parts?'.

    Uses analytical calibration for all configs so the total runtime
    stays tractable; CC perm across N methods x K folds would be
    hours. Analytical on K=9 x 4 configs is ~1-2 minutes.
    """

    config_progress = pyqtSignal(int, int, str)

    def __init__(self, datasets, labels, method_keys, base_params,
                 replication_mode='strict',
                 level='current',
                 per_method_top50=None):
        super().__init__()
        self.datasets = datasets
        self.labels = labels
        self.method_keys = list(method_keys)
        self.base_params = dict(base_params)
        self.replication_mode = replication_mode
        self.level = level
        self.per_method_top50 = per_method_top50 or {}

    def _run(self):
        from kosmic.meta_analysis.consensus_loo import (
            enumerate_method_configurations,
            run_consensus_loo,
        )

        K = len(self.datasets)
        if K < 3:
            raise RuntimeError(
                f"LOO benchmark requires at least 3 studies; got K={K}.")

        # Build the list of configs to benchmark depending on level
        if self.level == 'current':
            if not self.method_keys:
                raise RuntimeError("No methods selected for benchmarking.")
            # Existing behaviour: each chosen method alone + the
            # full consensus of the chosen set.
            configs = []
            for m in self.method_keys:
                configs.append({
                    'label':       f"{m} alone",
                    'method_keys': [m],
                    'n_methods':   1,
                    'states':      None,
                })
            if len(self.method_keys) > 1:
                configs.append({
                    'label': f"Consensus ({'+'.join(self.method_keys)})",
                    'method_keys': list(self.method_keys),
                    'n_methods':   len(self.method_keys),
                    'states':      None,
                })
        else:
            # Combinatorial sweep level
            configs = enumerate_method_configurations(
                self.level,
                user_top50_per_method=self.per_method_top50)
            if not configs:
                raise RuntimeError(
                    f"Enumeration produced no configs for level={self.level}.")

        n_total = len(configs)
        self.progress.emit(
            f"Benchmark ({self.level}): K={K} studies, "
            f"{n_total} configs (analytical calibration, "
            f"mode={self.replication_mode})")

        # Force analytical calibration for benchmark runs
        params = dict(self.base_params)
        params['calibration'] = 'analytical'

        results = []
        for i, cfg in enumerate(configs, 1):
            self.config_progress.emit(i, n_total, cfg['label'])
            self.progress.emit(
                f"Config {i}/{n_total}: {cfg['label']}")

            def _fold_cb(j, nk, name, _i=i, _n=n_total, _l=cfg['label']):
                self.progress.emit(
                    f"  config {_i}/{_n} ({_l})  "
                    f"fold {j}/{nk}: held out {name}")

            params_c = dict(params)
            params_c['method_keys'] = list(cfg['method_keys'])
            loo_res = run_consensus_loo(
                self.datasets, self.labels, params_c,
                replication_mode=self.replication_mode,
                pseudobulk_data=None,
                fold_progress_cb=_fold_cb,
            )
            results.append({
                'config_label': cfg['label'],
                'method_keys':  list(cfg['method_keys']),
                'n_methods':    cfg['n_methods'],
                'states':       cfg.get('states'),
                'loo_result':   loo_res,
            })

        self.progress.emit(
            f"Benchmark complete: {len(results)} configs.")
        return results


__all__ = [
    "ConsensusLOOWorker",
    "GeneMAWorker",
    "LOOBenchmarkWorker",
    "ReproducibilityWorker",
]
