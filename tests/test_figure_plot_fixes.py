"""Figures: enrichment significance labelling, GSEA mountain FDR text and
the highly variable genes plot."""
from __future__ import annotations

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from kosmic.visualisation.de.enrichment import create_enrichment_bar_plot
from kosmic.visualisation.de.gsea_mountain import _format_fdr
from kosmic.visualisation.scrna.hvg_plot import create_hvg_plot


def _terms(fdr):
    return pd.DataFrame({'Term': [f'term {i}' for i in range(len(fdr))],
                         'P_value': fdr, 'FDR': fdr})


def test_enrichment_axis_names_the_statistic():
    fig = create_enrichment_bar_plot(_terms([0.001, 0.01]),
                                     significance_label='elim P (uncorrected)')
    assert fig.axes[0].get_xlabel() == '-log10 elim P (uncorrected)'
    plt.close(fig)


def test_enrichment_says_when_nothing_passes():
    fig = create_enrichment_bar_plot(_terms([0.2, 0.4]))
    assert fig.axes[0].get_title().startswith('No terms at FDR < 0.05')
    plt.close(fig)


def test_zero_fdr_is_reported_as_a_bound():
    assert _format_fdr(0.0) == 'FDR < 0.001'
    assert _format_fdr(0.0123) == 'FDR = 0.012'


@pytest.mark.parametrize('flavor, col', [('seurat_v3', 'variances_norm'),
                                         ('seurat', 'dispersions_norm')])
def test_hvg_plot_both_flavours(flavor, col):
    var = pd.DataFrame({'means': [0.01, 0.1, 1.0, 0.0], col: [1.0, 3.0, 0.5, 1.0],
                        'highly_variable': [False, True, False, False]})
    fig = create_hvg_plot(var, flavor)
    assert fig.axes[0].get_xscale() == 'log'
    assert '1 of 4' in fig.axes[0].get_title()
    plt.close(fig)


def test_hvg_plot_needs_its_columns():
    assert create_hvg_plot(pd.DataFrame({'means': np.ones(3)})) is None


def test_volcano_labels_do_not_overlap():
    from kosmic.visualisation.de.volcano import create_volcano_plot
    # Twenty significant genes packed into a small region: naive placement
    # stacks their labels on top of each other.
    rng = np.random.default_rng(0)
    n = 20
    de = pd.DataFrame({
        'names': [f'GENE{i}' for i in range(n)],
        'logfoldchanges': 1.0 + rng.normal(0, 0.05, n),
        'pvals_adj': 10.0 ** -(10 + rng.normal(0, 0.2, n)),
    })
    fig = create_volcano_plot(de, 'DCM', 'NF', max_labels=n)
    renderer = fig.canvas.get_renderer()
    boxes = [t.get_window_extent(renderer) for t in fig.axes[0].texts]
    assert boxes
    for i, a in enumerate(boxes):
        for b in boxes[i + 1:]:
            assert not a.overlaps(b)
    plt.close(fig)


def test_embedding_labels_on_plot():
    from kosmic.visualisation.scrna.umap import create_embedding_plot
    rng = np.random.default_rng(1)
    coords = np.vstack([rng.normal(c, 0.3, (200, 2)) for c in ((0, 0), (0.2, 0), (5, 5))])
    labels = np.repeat(['Alpha cell', 'Beta cell', 'Gamma cell'], 200)
    fig = create_embedding_plot(coords, labels, labels_on_plot=True)
    ax = fig.axes[0]
    assert ax.get_legend() is None
    renderer = fig.canvas.get_renderer()
    boxes = [t.get_window_extent(renderer) for t in ax.texts]
    assert len(boxes) == 3
    for i, a in enumerate(boxes):
        for b in boxes[i + 1:]:
            assert not a.overlaps(b)
    plt.close(fig)


def test_patient_dotplot_shows_pathway_de_fdr():
    from kosmic.visualisation.de.dotplots import create_patient_dotplot
    rng = np.random.default_rng(2)
    sdf = pd.DataFrame({'condition': ['NF'] * 6 + ['DCM'] * 6,
                        'TCA_raw': rng.normal(0, 1, 12), 'OXPHOS_raw': rng.normal(0, 1, 12)})
    fig = create_patient_dotplot(sdf, ['TCA', 'OXPHOS'], 'DCM', 'NF',
                                 fdr_by_pathway={'TCA': 0.003})
    texts = [t.get_text() for ax in fig.axes for t in ax.texts]
    assert texts == ['FDR = 0.003']   # OXPHOS has no Pathway DE result: no label
    plt.close(fig)


def test_enrichment_dot_plot():
    from kosmic.visualisation.de.enrichment import create_enrichment_dot_plot, top_terms
    df = pd.DataFrame({'Term': ['GO:0006099 tricarboxylic acid cycle', 'GO:1 b', 'GO:2 c'],
                       'FDR': [0.001, 0.01, 0.5], 'Fold_Enrichment': [3.5, 1.4, 1.1],
                       'Gene_Count': [14, 120, 30]})
    terms, note = top_terms(df)
    assert list(terms['Term']) == ['tricarboxylic acid cycle', 'b'] and note is None
    fig = create_enrichment_dot_plot(df, significance_label='elim P (uncorrected)')
    ax = fig.axes[0]
    assert ax.get_xlabel() == 'Fold enrichment'
    assert ax.get_yticklabels()[0].get_text() == 'tricarboxylic acid cycle'   # top row
    assert 'elim P' in fig.axes[1].get_ylabel()
    plt.close(fig)


def test_dataset_summary_gates_on_fdr_unless_lfc_set():
    from kosmic.visualisation.de.dataset_summary import compute_dataset_summary
    de = pd.DataFrame({'names': ['A', 'B', 'C'], 'logfoldchanges': [-0.1, -1.0, 0.5],
                       'pvals_adj': [0.01, 0.01, 0.5]})
    sets = {'P': ['A', 'B', 'C']}
    row = compute_dataset_summary(de, sets, {}).iloc[0]
    assert row['sig_genes_down'] == 2                       # FDR alone
    row = compute_dataset_summary(de, sets, {}, lfc_threshold=0.25).iloc[0]
    assert row['sig_genes_down'] == 1                       # and |log2FC| > 0.25
