# Annotation: marker scoring, reference-based assignment, CellTypist.

from kosmic.scrna.annotate.score import (
    DEFAULT_MARKERS,
    collapse_celltypist_label,
    score_marker_genes,
    assign_cell_types_per_cluster,
    run_ora_per_cluster,
    annotate_clusters,
    run_celltypist,
)
from kosmic.scrna.annotate.reference import (
    REFERENCES_DIR,
    BUILTIN_REFERENCES,
    list_available_references,
    load_reference,
    annotate_by_reference,
    run_reference_annotation,
    build_reference_from_h5ad,
)
from kosmic.scrna.annotate.marker_genes import (
    compute_cluster_marker_genes,
    get_top_marker_genes,
)

__all__ = [
    "DEFAULT_MARKERS",
    "collapse_celltypist_label",
    "score_marker_genes",
    "assign_cell_types_per_cluster",
    "run_ora_per_cluster",
    "annotate_clusters",
    "run_celltypist",
    "REFERENCES_DIR",
    "BUILTIN_REFERENCES",
    "list_available_references",
    "load_reference",
    "annotate_by_reference",
    "run_reference_annotation",
    "build_reference_from_h5ad",
    "compute_cluster_marker_genes",
    "get_top_marker_genes",
]
