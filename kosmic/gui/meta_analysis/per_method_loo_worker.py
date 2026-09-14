# Shared worker that runs held-out K-fold LOO independently for each of
# a set of pooling methods, returning per-method direction-replication
# rates plus optional held-out case/control AUC (when pseudobulk is
# available). Reuses 'run_consensus_loo' with method_keys=[single_method]
# so each per-method LOO is the same machinery the consensus page uses.
#
# Direction replication runs analytical-only (no CC perm) since the
# threshold lives in the fold's consensus build, not the replication
# test. AUC is computed only when pseudobulk is provided.

from __future__ import annotations

from typing import List, Optional

from kosmic.gui.shared.widgets import BaseWorker
from kosmic import DEFAULT_FDR


class PerMethodLOOWorker(BaseWorker):
    """
    Run held-out LOO for each method independently.

    Emits finished_ok with a dict:
        {
            'per_method': {key: loo_result_dict, ...},
            'method_keys': list of keys with valid results,
        }
    """

    def __init__(self, datasets, labels, method_keys, params,
                 pseudobulk_paths: Optional[List[str]] = None,
                 parent=None):
        super().__init__(parent)
        self.datasets = datasets
        self.labels = labels
        self.method_keys = method_keys
        self.params = params
        self.pseudobulk_paths = pseudobulk_paths or []

    def _run(self):
        from kosmic.meta_analysis.consensus_loo import run_consensus_loo

        # Load pseudobulk once if we have paths -- this enables AUC.
        pb_data = None
        if self.pseudobulk_paths:
            pb_data = self._load_pseudobulk(self.pseudobulk_paths)
            if not pb_data:
                self.progress.emit(
                    "pseudobulk load failed; AUC will be unavailable")

        per_method = {}
        for i, key in enumerate(self.method_keys):
            self.progress.emit(
                f"LOO {i+1}/{len(self.method_keys)}: {key}...")
            params_single = dict(self.params)
            params_single['method_keys'] = [key]
            try:
                result = run_consensus_loo(
                    self.datasets, self.labels, params_single,
                    pseudobulk_data=pb_data,
                    consensus_fdr_threshold=self.params.get(
                        'consensus_fdr_threshold', DEFAULT_FDR),
                    replication_threshold=self.params.get(
                        'replication_threshold', DEFAULT_FDR),
                    top_n=self.params.get('top_n', 1000),
                    replication_mode='direction_only',
                    # Stash per-fold per-method gene sets + minimal
                    # held-out data so Consensus Evaluation can
                    # compute direction-replication for arbitrary
                    # combinations without re-running LOO.
                    stash_fold_details=True,
                )
            except Exception as e:
                import traceback
                self.progress.emit(
                    f"LOO failed for {key}: {e}\n"
                    f"{traceback.format_exc()}")
                continue
            per_method[key] = result

        if not per_method:
            raise RuntimeError("LOO produced no results for any method.")

        return {
            'per_method': per_method,
            'method_keys': list(per_method.keys()),
            'pseudobulk_available': pb_data is not None,
        }

    def _load_pseudobulk(self, pb_paths):
        """
        Skip per-path failures so a single unloadable study does not
        kill the whole per-method LOO run.
        """
        from kosmic.meta_analysis.io import load_pseudobulk_set
        return load_pseudobulk_set(
            pb_paths, labels=self.labels,
            progress_cb=self.progress.emit, skip_errors=True)
