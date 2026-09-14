# Provenance-based staleness check for figure-export pages: has the underlying
# data changed since the DE / pathway-scoring results a figure renders were
# computed? Compares the data token each stage consumed against a live
# fingerprint of the current adata, so a figure warns when it reflects an old run
# (e.g. QC / decontX / role assignment changed after DE ran).
from __future__ import annotations

from pathlib import Path
from typing import Optional


def de_results_stale_note(de_ws) -> Optional[str]:
    """Warning text if the loaded DE / pathway results predate the current data,
    else None. Best-effort: returns None on any missing input or error."""
    if de_ws is None:
        return None
    h5ad = getattr(de_ws, 'h5ad_path', None)
    adata = getattr(de_ws, 'current_adata', None)
    if not h5ad or adata is None:
        return None
    try:
        from kosmic import provenance
        live = provenance.compute_fingerprint(adata, h5ad_path=h5ad)['token']
        study_dir = Path(h5ad).parent
        for stage in ('gene_de', 'pathway_de'):
            st = provenance.last_stage(study_dir, stage)
            consumed = st.get('consumed_token') if st else None
            if consumed and consumed != live:
                return ("Results may be out of date: the data changed since "
                        "these were computed. Re-run DE / pathway scoring in the "
                        "DE window to refresh this figure.")
    except Exception:
        return None
    return None
