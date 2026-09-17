"""KOSMIC -- pure analysis logic. No GUI framework imports allowed in
this package.

Configuration values live in ``config.toml`` at the repo root and are
loaded once at import time. To change a default, edit the TOML and
restart Python. Consumers do ``from kosmic import DEFAULT_FDR`` etc.
"""
from __future__ import annotations

import os
import tomllib
from pathlib import Path

_CONFIG_PATH = Path(__file__).resolve().parent.parent / "config.toml"
with _CONFIG_PATH.open("rb") as _f:
    _cfg = tomllib.load(_f)

# Significance thresholds
DEFAULT_FDR                  = _cfg["significance"]["fdr"]
DEFAULT_LFC_THRESHOLD        = _cfg["significance"]["lfc_threshold"]

# Meta-analysis
MIN_STUDIES                  = _cfg["meta_analysis"]["min_studies"]
CC_N_PERMS                   = _cfg["meta_analysis"]["cc_n_perms"]

# Pathway scoring
PATHWAY_GENE_DETECTION_PCT   = _cfg["pathway_scoring"]["gene_detection_pct"]
PATHWAY_MIN_GENES            = _cfg["pathway_scoring"]["min_genes"]
PATHWAY_MIN_COVERAGE         = _cfg["pathway_scoring"]["min_coverage"]
PATHWAY_BORDERLINE_COVERAGE  = _cfg["pathway_scoring"]["borderline_coverage"]
PATHWAY_VIF_RHO              = _cfg["pathway_scoring"]["vif_rho"]
PATHWAY_N_CTRL_GENES         = _cfg["pathway_scoring"]["n_ctrl_genes"]

# GWAS overlap
GWAS_LP_THRESHOLD            = _cfg["gwas"]["lp_threshold"]
GWAS_WINDOW_KB               = _cfg["gwas"]["window_kb"]

# Help: documentation site (None = offline browser only)
HELP_SITE_URL                = (os.environ.get("KOSMIC_HELP_SITE_URL")
                                or _cfg.get("help", {}).get("site_url") or None)

# DE pipeline
DE_MIN_CELLS                 = _cfg["de"]["min_cells"]
DE_MIN_COUNTS                = _cfg["de"].get("min_counts", 0)
DE_MIN_EXPRESSING_SAMPLES    = _cfg["de"]["min_expressing_samples"]
DE_FILTER_MIN_COUNT          = _cfg["de"]["filter_min_count"]
DE_FILTER_MIN_TOTAL_COUNT    = _cfg["de"]["filter_min_total_count"]
DE_FILTER_MIN_SAMPLES        = _cfg["de"]["filter_min_samples"]
DE_MAX_SE                    = _cfg["de"]["max_se"]
DE_SE_FLOOR                  = _cfg["de"]["se_floor"]
DE_DETECTION_MIN_PCT         = _cfg["de"]["detection_min_pct"]
DE_DETECTION_ON              = _cfg["de"]["detection_on"]
DE_DETECTION_MIN_DONOR_FRAC  = _cfg["de"].get("detection_min_donor_frac", 0.0)
DE_DESEQ2_INDEPENDENT_FILTER = _cfg["de"]["deseq2_independent_filter"]
DE_DESEQ2_COOKS_FILTER       = _cfg["de"]["deseq2_cooks_filter"]

# DE independent filtering
DE_ENRICH_IF_DEFAULT_PCT     = _cfg["de"]["independent_filtering"]["enrichment_default_pct"]
DE_ENRICH_IF_ON_BY_DEFAULT   = _cfg["de"]["independent_filtering"]["enrichment_on_by_default"]

# Meta-analysis independent filtering
META_IF_DEFAULT_PCT          = _cfg["meta_analysis"]["independent_filtering"]["default_pct"]
META_IF_ON_BY_DEFAULT        = _cfg["meta_analysis"]["independent_filtering"]["on_by_default"]
META_IF_FILTER_STAT          = _cfg["meta_analysis"]["independent_filtering"]["filter_stat"]

# scRNA processing
MARKER_MAX_CELLS_PER_CLUSTER = _cfg.get("scrna", {}).get("marker_max_cells_per_cluster", 5000)

# UI / typography
DEFAULT_BASE_FONT_SIZE       = _cfg["ui"]["base_font_size"]
UI_FONT_SCALE                = _cfg["ui"].get("font_scale", 1.0)
DEFAULT_EXPORT_DPI           = _cfg["ui"]["export_dpi"]
UI_EMBEDDING_MAX_POINTS      = _cfg["ui"].get("embedding_max_points", 100_000)
