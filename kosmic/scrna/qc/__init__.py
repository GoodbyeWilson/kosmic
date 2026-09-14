# QC tab: filter pipeline, normalisation, doublet detection.

from kosmic.scrna.qc.filter import (
    detect_mitochondrial_genes,
    fix_nan_in_sparse,
    run_qc_pipeline,
)
from kosmic.scrna.qc.mad_thresholds import compute_mad_thresholds
from kosmic.scrna.qc.normalize import normalize_adata
from kosmic.scrna.qc.scrublet import run_scrublet
from kosmic.scrna.qc.soupx import run_soupx

__all__ = [
    "detect_mitochondrial_genes",
    "fix_nan_in_sparse",
    "run_qc_pipeline",
    "compute_mad_thresholds",
    "normalize_adata",
    "run_scrublet",
    "run_soupx",
]
