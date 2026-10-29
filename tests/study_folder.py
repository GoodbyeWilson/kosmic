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
                        'means': counts.mean(axis=0)}, index=genes)
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
