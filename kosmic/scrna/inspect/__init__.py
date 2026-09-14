# Inspect tab utilities: column detection, role assignment, gene-name harmonisation.

from kosmic.scrna.inspect.detection import detect_species, format_gene_for_species
from kosmic.scrna.inspect.batch import (
    detect_batch_column,
    detect_condition_column,
    summarise_samples,
    STANDARDISED_OBS_NAMES,
    USER_OBS_PATTERNS,
    COMPUTED_OBS_PATTERNS,
)
from kosmic.scrna.inspect.roles import (
    parse_role_map,
    roles_to_labels,
    resolve_roles,
    CONTROL_VALUE_PATTERNS,
    DISEASE_VALUE_PATTERNS,
)
from kosmic.scrna.inspect.gene_names import (
    HarmonisationReport,
    load_hgnc_lookup,
    harmonise_adata,
    download_hgnc_table,
)
from kosmic.scrna.inspect.gene_expression import (
    analyze_gene_expression,
    compute_dotplot_stats,
)
from kosmic.scrna.inspect.transcriptome_size import test_transcriptome_size_adata

__all__ = [
    "detect_species",
    "format_gene_for_species",
    "detect_batch_column",
    "detect_condition_column",
    "summarise_samples",
    "STANDARDISED_OBS_NAMES",
    "USER_OBS_PATTERNS",
    "COMPUTED_OBS_PATTERNS",
    "parse_role_map",
    "roles_to_labels",
    "resolve_roles",
    "CONTROL_VALUE_PATTERNS",
    "DISEASE_VALUE_PATTERNS",
    "HarmonisationReport",
    "load_hgnc_lookup",
    "harmonise_adata",
    "download_hgnc_table",
    "analyze_gene_expression",
    "compute_dotplot_stats",
    "test_transcriptome_size_adata",
]
