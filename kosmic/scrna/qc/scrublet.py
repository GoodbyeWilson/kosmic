# Doublet detection via Scrublet (Wolock, Lopez & Klein 2019).
#
# Uses scanpy's port of the algorithm ('scanpy.pp.scrublet') rather than
# the standalone 'scrublet' package. Same method, same outputs -- on
# planted cross-type doublets the two agree on 99.5% of calls with AUC
# 1.0 each -- but the standalone package hard-depends on 'annoy', which
# ships no Windows wheel and so needs a C++ compiler at install time.
# scanpy's version only touches annoy when asked for approximate
# neighbours, which this never does.
#
# Doublets form within one droplet library, so each sample is scored on
# its own: pooled, Scrublet simulates doublets from nuclei of different
# donors, which cannot occur. The loop is written here rather than using
# scanpy's 'batch_key' so that a sample too small to score is skipped
# and reported instead of failing the whole run.
#
# Scrublet's automatic threshold is the dip between the two modes of the
# simulated-doublet scores. A sample without real doublets (an
# already-cleaned deposit, or one dominated by a single cell type) has
# no such dip, and the cut then lands in noise: on a clean single-type
# sample it called 26% of cells. The threshold is accepted only when
# most simulated doublets score above it; otherwise the sample gets no
# calls. A cut at the expected rate is available, but only on request.
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import pandas as pd

#: Fewest scored cells a sample needs; below this Scrublet's PCA and
#: neighbour graph are not meaningful.
MIN_CELLS_PER_SAMPLE = 100

#: Fraction of simulated doublets that must score above the automatic
#: threshold for it to be trusted (Scrublet's detectable doublet fraction).
MIN_DETECTABLE_FRACTION = 0.5


def _score_one(counts, expected_doublet_rate, random_state):
    """Scrublet on one sample's counts: (scores, calls, threshold, detectable)."""
    import anndata as ad
    import scanpy as sc

    sub = ad.AnnData(X=counts)
    n_comps = int(min(30, sub.n_obs - 1, sub.n_vars - 1))
    sc.pp.scrublet(
        sub,
        expected_doublet_rate=expected_doublet_rate,
        sim_doublet_ratio=2.0,
        n_neighbors=None,
        n_prin_comps=n_comps,
        use_approx_neighbors=False,   # keep annoy out of it
        random_state=random_state,
        verbose=False,
    )
    uns = sub.uns.get('scrublet', {})
    threshold = uns.get('threshold')
    sim = np.asarray(uns.get('doublet_scores_sim', []), dtype=float)
    detectable = (float((sim > threshold).mean())
                  if threshold is not None and np.isfinite(threshold) and sim.size
                  else 0.0)
    return (sub.obs['doublet_score'].to_numpy(dtype=float),
            sub.obs['predicted_doublet'].to_numpy().astype(bool),
            threshold, detectable)


def run_scrublet(
    adata,
    sample_col: Optional[str] = None,
    expected_doublet_rate: float = 0.06,
    min_counts: int = 2,
    force_expected_rate: bool = False,
    progress_callback=None,
) -> Tuple:
    """Run Scrublet doublet detection on each sample of an AnnData object.

    Parameters
    ----------
    adata : anndata.AnnData
        Data to check for doublets. Scored on its raw counts.
    sample_col : str, optional
        obs column identifying the droplet library (sample). None scores
        the whole object as one sample.
    expected_doublet_rate : float
        Expected fraction of doublets, passed to Scrublet.
    min_counts : int
        Minimum UMI counts for a cell to be scored. Cells below this
        get a score of 0 and are never called doublets.
    force_expected_rate : bool
        Where a sample's automatic threshold is rejected, call the top
        'expected_doublet_rate' of its scores instead of calling none.
    progress_callback : callable(msg), optional

    Returns
    -------
    tuple of (anndata.AnnData, dict)
        adata with 'doublet_score' and 'predicted_doublet' in obs. The dict
        has n_doublets, n_total, pct_doublets, doublet_scores, and
        'per_sample': one row per sample with its threshold, detectable
        fraction, doublets called and status ('called', 'no clear
        threshold', 'forced cut' or 'too few cells').
    """
    import scipy.sparse as sp

    from kosmic.scrna.counts import count_source
    counts, _, _ = count_source(adata)
    counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)

    totals = np.asarray(counts.sum(axis=1)).ravel()
    scorable = totals >= min_counts
    if sample_col is None:
        samples = np.zeros(adata.n_obs, dtype=object)
        samples[:] = 'all cells'
    else:
        samples = adata.obs[sample_col].astype(str).to_numpy()

    scores = np.zeros(adata.n_obs, dtype=float)
    predicted = np.zeros(adata.n_obs, dtype=bool)
    rows = []
    unique = list(pd.unique(samples))
    for i, s in enumerate(unique):
        idx = np.flatnonzero((samples == s) & scorable)
        if progress_callback:
            progress_callback(f"Scrublet: sample {i + 1}/{len(unique)} '{s}' "
                              f"({idx.size:,} cells)...")
        row = {'sample': s, 'n_cells': int((samples == s).sum()),
               'n_scored': int(idx.size), 'threshold': np.nan,
               'detectable_fraction': np.nan, 'n_doublets': 0}
        if idx.size < MIN_CELLS_PER_SAMPLE:
            rows.append({**row, 'status': 'too few cells'})
            continue

        s_scores, s_calls, threshold, detectable = _score_one(
            counts[idx].astype(np.float32), expected_doublet_rate, 42)
        scores[idx] = s_scores
        row['detectable_fraction'] = round(detectable, 3)
        if (threshold is not None and np.isfinite(threshold)
                and detectable >= MIN_DETECTABLE_FRACTION):
            predicted[idx] = s_calls
            row.update(threshold=float(threshold), status='called')
        elif force_expected_rate:
            # Scores take few distinct values, so ties at the cut are kept
            # rather than dropped; the count can exceed the expected rate.
            cut = float(np.percentile(s_scores, 100 * (1 - expected_doublet_rate)))
            predicted[idx] = s_scores >= cut
            row.update(threshold=cut, status='forced cut')
        else:
            row['status'] = 'no clear threshold'
        row['n_doublets'] = int(predicted[idx].sum())
        rows.append(row)

    adata.obs['doublet_score'] = scores
    adata.obs['predicted_doublet'] = predicted

    n_doublets = int(predicted.sum())
    n_total = adata.n_obs
    results = {
        'n_doublets': n_doublets,
        'n_total': n_total,
        'pct_doublets': (n_doublets / n_total) * 100 if n_total else 0.0,
        'doublet_scores': scores,
        'per_sample': pd.DataFrame(rows),
    }
    return adata, results
