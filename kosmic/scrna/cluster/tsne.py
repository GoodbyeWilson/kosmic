# t-SNE embedding with openTSNE / scanpy fallback.
#
# Prefers 'openTSNE' (FFT-accelerated, much faster on >10k cells) and
# falls back to 'sc.tl.tsne' (scikit-learn under the hood) when openTSNE
# isn't installed. Reads 'X_pca_harmony' if present, else 'X_pca'.
from __future__ import annotations

import numpy as np


def compute_tsne(adata, perplexity: float = 100, random_state: int = 0):
    """Compute t-SNE; write 'adata.obsm['X_tsne']'.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have 'X_pca' (or 'X_pca_harmony') in 'obsm'.
    perplexity : float
    random_state : int

    Returns
    -------
    anndata.AnnData
        Modified in place.
    """
    use_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'
    pca_data = adata.obsm[use_rep]
    learning_rate = max(200, adata.n_obs / 12)

    try:
        from openTSNE import TSNE as OpenTSNE
        tsne = OpenTSNE(
            perplexity=perplexity,
            learning_rate=learning_rate,
            n_iter=750,
            n_jobs=-1,
            initialization='pca',
            random_state=random_state,
        )
        embedding = tsne.fit(pca_data)
        adata.obsm['X_tsne'] = np.array(embedding)
    except ImportError:
        import scanpy as sc
        sc.tl.tsne(
            adata,
            use_rep=use_rep,
            perplexity=perplexity,
            learning_rate=learning_rate,
            random_state=random_state,
        )

    return adata
