"""Visualisation package. Sub-packages: scrna, de, meta, shared."""

from kosmic.meta_analysis.io import load_de_results
from kosmic.visualisation.scrna.umap import (
    create_embedding_plot,
    create_highlight_plot,
    create_overview_figure,
    create_scatter_plots,
)
from kosmic.visualisation.scrna.clustree import create_clustree_figure
from kosmic.visualisation.scrna.variable_genes import (
    compute_variable_genes,
    create_variable_gene_plot,
)
from kosmic.visualisation.de.volcano import create_volcano_plot
from kosmic.visualisation.de.gene_heatmap import (
    create_gene_heatmap,
    prepare_heatmap_data,
)
from kosmic.visualisation.de.bubble_plot import create_bubble_plot
from kosmic.visualisation.de.pathway_plots import create_pathway_bar_plot
from kosmic.visualisation.de.top_de_plots import (
    create_pathway_summary_chart,
    create_gene_bar_chart,
)
from kosmic.visualisation.de.dataset_summary import (
    compute_dataset_summary,
    create_summary_figure,
)
from kosmic.visualisation.de.dotplots import (
    pathway_scores_to_sample_df,
    create_patient_dotplot,
    create_gene_dotplot,
    create_expression_dotplot,
)
from kosmic.visualisation.de.enrichment import create_enrichment_bar_plot
from kosmic.visualisation.meta.heatmap_forest import create_multi_dataset_figure
from kosmic.visualisation.meta.forest_plot import create_family_forest_plot


# Shared typography scale used by every plot module.
_TYPE_RATIOS = {
    'title': 1.4,
    'suptitle': 1.5,
    'axis_label': 1.3,
    'tick': 0.9,
    'legend': 0.9,
    'annotation': 0.8,
    'colorbar_label': 1.0,
    'colorbar_tick': 0.85,
}

_DEFAULT_BASE = 11  # pt


def default_font_sizes(base: int = _DEFAULT_BASE) -> dict:
    """Return the canonical font-size dict from a *base* size in pt.

    Keys: title, suptitle, axis_label, tick, legend, annotation,
    colorbar_label, colorbar_tick, base.
    """
    sizes = {name: round(base * ratio, 1) for name, ratio in _TYPE_RATIOS.items()}
    sizes['base'] = base
    return sizes


__all__ = [
    "default_font_sizes",
    "load_de_results",
    "create_embedding_plot",
    "create_highlight_plot",
    "create_overview_figure",
    "create_scatter_plots",
    "create_clustree_figure",
    "compute_variable_genes",
    "create_variable_gene_plot",
    "create_volcano_plot",
    "create_gene_heatmap",
    "prepare_heatmap_data",
    "create_bubble_plot",
    "create_pathway_bar_plot",
    "create_pathway_summary_chart",
    "create_gene_bar_chart",
    "compute_dataset_summary",
    "create_summary_figure",
    "pathway_scores_to_sample_df",
    "create_patient_dotplot",
    "create_gene_dotplot",
    "create_expression_dotplot",
    "create_enrichment_bar_plot",
    "create_multi_dataset_figure",
    "create_family_forest_plot",
]
