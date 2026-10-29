"""A small study folder laid out as KOSMIC writes one, for the tests of
Figure Export reading a study's saved data (ADR-014).

'make_study_folder' writes processed_data/<name>.h5ad, the study manifest
and, as the stages that read them arrive, the saved DE outputs. The h5ad
comes in the two forms real studies have:

- 'log': X log-normalised, the counts in layers['counts'] (and a
  'decontX_counts' layer), no uns['log1p'] -- as the DCM studies, so the
  reader must decide from the values.
- 'counts': the counts in X and no layers -- as a combined atlas.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

#: Genes named after real markers, so the bundled marker panels find them.
MARKER_GENES = ['PTPRC', 'CD3E', 'CD14', 'LYZ', 'PECAM1', 'VWF', 'COL1A1',
                'DCN', 'TNNT2', 'MYH7', 'ACTA2', 'MYH11']
N_CELLS, N_FILLER = 240, 48


def make_adata(form: str = 'log', seed: int = 0):
    """The study's AnnData in memory, before it is written."""
    import anndata as ad
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    genes = MARKER_GENES + [f'G{i}' for i in range(N_FILLER)]
    counts = rng.negative_binomial(2, 0.3, size=(N_CELLS, len(genes))).astype(np.float32)
    counts[rng.random(counts.shape) < 0.6] = 0
    counts[0, 0] = 45                     # a count of 20 or more early on
    counts[5] = 0                         # a cell with no counts
    samples = np.repeat([f's{i}' for i in range(1, 7)], N_CELLS // 6)
    condition = np.where(np.isin(samples, ['s1', 's2', 's3']), 'Control', 'Disease')
    obs = pd.DataFrame({
        'sample': pd.Categorical(samples),
        'condition': pd.Categorical(condition),
        '_role': pd.Categorical(np.where(condition == 'Control', 'control', 'disease')),
        'cell_type': pd.Categorical(rng.choice(['Endothelial', 'Fibroblast', 'Myeloid'],
                                               N_CELLS)),
        'n_genes_by_counts': (counts > 0).sum(axis=1),
        'total_counts': counts.sum(axis=1),
        'pct_counts_mt': rng.random(N_CELLS) * 5,
    }, index=[f'cell{i}' for i in range(N_CELLS)])
    var = pd.DataFrame({'highly_variable': rng.random(len(genes)) < 0.3,
                        'means': counts.mean(axis=0),
                        'variances_norm': rng.random(len(genes)) * 3}, index=genes)
    csr = sp.csr_matrix(counts)
    if form == 'log':
        totals = counts.sum(axis=1, keepdims=True)
        lognorm = np.log1p(np.divide(counts * 1e4, totals, out=np.zeros_like(counts),
                                     where=totals > 0))
        a = ad.AnnData(X=sp.csr_matrix(lognorm.astype(np.float32)), obs=obs, var=var)
        a.layers['counts'] = csr
        a.layers['decontX_counts'] = sp.csr_matrix(np.floor(counts * 0.9))
    elif form == 'counts':
        a = ad.AnnData(X=csr, obs=obs, var=var)
    else:
        raise ValueError(form)
    a.obsm['X_umap'] = rng.normal(size=(N_CELLS, 2)).astype(np.float32)
    a.obsm['X_pca'] = rng.normal(size=(N_CELLS, 10)).astype(np.float32)
    a.obsp['connectivities'] = sp.random(N_CELLS, N_CELLS, density=0.02, format='csr',
                                         random_state=seed)
    a.uns['pca'] = {'variance_ratio': np.linspace(0.3, 0.01, 10)}
    a.uns['hvg'] = {'flavor': 'seurat_v3'}
    a.uns['role_map'] = {'Control': 'control', 'Disease': 'disease'}
    return a


def de_table(genes, seed: int = 0) -> pd.DataFrame:
    """A gene DE table with the columns the DE steps write."""
    rng = np.random.default_rng(seed)
    lfc = rng.normal(size=len(genes))
    p = rng.random(len(genes)) ** 3
    return pd.DataFrame({
        'names': list(genes), 'logfoldchanges': lfc, 'abs_logfoldchange': np.abs(lfc),
        'pvals': p, 'pvals_adj': np.minimum(p * 4, 1.0), 'se': rng.random(len(genes)) + 0.1,
        'disease_samples': 3, 'control_samples': 3})


def write_de_outputs(study: Path) -> Path:
    """Saved DE outputs named as on the real projects: per-cell-type runs
    (one with the study name in lower case, as GSE174574 has), a whole-study
    run of another method, the significant and pseudobulk files beside
    them, a batch summary whose 'de_path' column points to another
    computer, and the pathway scoring outputs at both levels."""
    from kosmic.de.fgsea import FGSEA_RESULTS_FILE, GSEA_LIBRARY_RESULTS_FILE
    from kosmic.paths import de_pathway_scoring_dir, de_stats_dir

    name = Path(study).name
    stats, scoring = de_stats_dir(study), de_pathway_scoring_dir(study)
    genes = MARKER_GENES + [f'G{i}' for i in range(N_FILLER)]
    for i, acc in enumerate([f'{name}_Endothelial', f'{name}_Fibroblast',
                             f'{name.lower()}_Myeloid_NN']):
        de_table(genes, seed=i).to_csv(stats / f'{acc}_DE_deseq2.csv', index=False)
        de_table(genes[:5], seed=i).to_csv(stats / f'{acc}_significant_genes.csv',
                                           index=False)
        pd.DataFrame({'G0': [1, 2]}).to_csv(stats / f'{acc}_pseudobulk.csv')
    de_table(genes, seed=9).to_csv(stats / f'{name}_DE_welch_cpm_eb.csv', index=False)
    pd.DataFrame({
        'cell_type': ['Endothelial', 'Fibroblast'],
        'accession': [f'{name}_Endothelial', f'{name}_Fibroblast'],
        'status': ['ok', 'ok'],
        'de_path': [rf'C:\Users\someone\{name}\results\x_DE_deseq2.csv'] * 2,
    }).to_csv(stats / f'{name}_by_cell_type_summary.csv', index=False)
    pw = pd.DataFrame({'names': ['Glycolysis', 'OXPHOS'], 'logfoldchanges': [0.5, -0.4],
                       'pvals': [0.01, 0.2], 'pvals_adj': [0.02, 0.3],
                       'Effect_Size': [0.8, -0.5]})
    pw.to_csv(scoring / 'pathway_de_results.csv', index=False)
    pw.iloc[:1].to_csv(scoring / 'parent_pathway_de_results.csv', index=False)
    pd.DataFrame({'sample': ['s1', 's4'], 'condition': ['Control', 'Disease'],
                  'Glycolysis': [0.1, 0.9], 'OXPHOS': [0.4, 0.2]}).to_csv(
        scoring / 'pathway_scores_per_donor.csv', index=False)
    gsea = pd.DataFrame({'names': ['HALLMARK_A', 'HALLMARK_B'], 'nes': [1.8, -1.2],
                         'es': [0.5, -0.3], 'pvals': [0.001, 0.04],
                         'pvals_adj': [0.01, 0.2], 'leading_edge': ['VWF;PECAM1', 'LYZ']})
    gsea.to_csv(scoring / FGSEA_RESULTS_FILE, index=False)
    gsea.to_csv(scoring / GSEA_LIBRARY_RESULTS_FILE, index=False)
    return study


def make_study_folder(root: Path, name: str = 'StudyA', form: str = 'log') -> Path:
    """Write a study folder under 'root' and return it."""
    from kosmic.manifest import write_manifest
    from kosmic.paths import processed_data_dir, scaffold_study

    study = Path(root) / name
    scaffold_study(study)
    adata = make_adata(form)
    h5ad = processed_data_dir(study) / f'{name}.h5ad'
    adata.write_h5ad(h5ad)
    write_manifest(
        study,
        dataset={'file': h5ad.name, 'n_cells': adata.n_obs, 'n_genes': adata.n_vars},
        semantics={'sample_column': 'sample', 'condition_column': 'condition',
                   'cell_type_column': 'cell_type',
                   'role_map': {'Control': 'control', 'Disease': 'disease'}})
    return study


def figure_workspace(*, scrna_adata=None, de=None, meta=None, study=None,
                     name: str = 'HeldStudy'):
    """A stand-in for FigureExportWorkspace, for page tests: the pages read
    their study through a real StudySource.

    'scrna_adata' is held by a stand-in scRNA workspace and 'de' (an object
    with the DE workspace's attribute names) is the DE workspace, both on
    'study'; 'study' defaults to an empty temporary folder called 'name'.
    'meta' is the Meta-Analysis workspace.
    """
    import tempfile
    from types import SimpleNamespace

    from kosmic.gui.figure_export.shared.study_source import StudySource
    if study is None:
        study = Path(tempfile.mkdtemp(prefix='kosmic_fig_')) / name
        study.mkdir()
    study = Path(study)
    scrna = (SimpleNamespace(current_adata=scrna_adata, current_project_dir=str(study))
             if scrna_adata is not None else None)
    if de is not None:
        de.project_dir = study
    src = StudySource(lambda: scrna, lambda: de)
    src.set_study(study)
    return SimpleNamespace(scrna_ws=scrna, de_ws=de, meta_ws=meta, source=src,
                           data_version=lambda: 0)
