# Differential expression engines.
#
# - 'pathway_scoring.py' -- pathway scoring methods (CAMERA-style + VIF)
# - 'de_analysis.py'     -- pseudobulk DE engine (DESeq2 / Welch)
#
# 'deseq2_fast.py' lives in 'kosmic/meta_analysis/' alongside its only
# live caller (cc_permutation.py).

from kosmic.de.de_analysis import (
    prepare_gene_coverage,
    create_pseudobulk,
    run_pseudobulk_de,
    filter_by_expression,
    annotate_de_results,
    run_de_pipeline,
)

__all__ = [
    "prepare_gene_coverage",
    "create_pseudobulk",
    "run_pseudobulk_de",
    "filter_by_expression",
    "annotate_de_results",
    "run_de_pipeline",
]
