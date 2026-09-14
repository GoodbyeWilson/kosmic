# Gene-group expression summary.
#
# Given an AnnData and a list of genes (a 'gene group' such as ion
# channels), summarise how that group is represented across the dataset:
# which genes are detected, their mean expression and fraction of cells
# expressing, and -- when an obs grouping is supplied -- the per-group
# mean / percent matrices that drive a dot plot.

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scipy.sparse as sp


@dataclass
class GeneGroupSummary:
    """Result of summarising a gene group against a dataset."""

    present: list          # genes found in the dataset (input order, deduped)
    missing: list          # requested genes not in the dataset
    group_col: str = None  # obs column used for the per-group breakdown
    groups: list = field(default_factory=list)        # category labels
    per_gene: pd.DataFrame = field(default_factory=pd.DataFrame)    # index=gene
    mean_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)  # genes x groups
    pct_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)   # genes x groups


def _as_1d(arr):
    """Flatten a numpy / np.matrix result of a sparse reduction to 1d."""
    return np.asarray(arr).ravel()


def _col_mean(X):
    """Per-column mean of a cells x genes matrix (sparse or dense)."""
    if X.shape[0] == 0:
        return np.full(X.shape[1], np.nan)
    return _as_1d(X.mean(axis=0))


def _col_pct(X, threshold):
    """Per-column fraction of rows with value above 'threshold' (0-1)."""
    n = X.shape[0]
    if n == 0:
        return np.full(X.shape[1], np.nan)
    if sp.issparse(X):
        if threshold == 0:
            nz = _as_1d((X > 0).sum(axis=0))
        else:
            nz = _as_1d((X > threshold).sum(axis=0))
        return nz / n
    return _as_1d((X > threshold).sum(axis=0)) / n


def summarise_gene_group(adata, genes, group_col=None, layer=None,
                         expr_threshold=0.0):
    """Summarise expression of a gene group across a dataset.

    Parameters
    ----------
    adata : AnnData
        Cells x genes. Uses 'adata.X' unless 'layer' is given.
    genes : sequence of str
        The gene group. Order is preserved; duplicates and genes absent
        from the dataset are dropped from 'present' (the latter recorded
        in 'missing').
    group_col : str, optional
        An 'adata.obs' column to break expression down by (e.g. 'cell_type'
        or 'leiden'). If None, only the dataset-wide per-gene stats are
        computed and the matrices are empty.
    layer : str, optional
        Layer to read instead of 'adata.X'.
    expr_threshold : float
        A cell counts as 'expressing' a gene when its value exceeds this
        (default 0).

    Returns
    -------
    GeneGroupSummary
    """
    var_set = set(map(str, adata.var_names))
    seen = set()
    present, missing = [], []
    for g in genes:
        g = str(g)
        if g in seen:
            continue
        seen.add(g)
        (present if g in var_set else missing).append(g)

    if not present:
        return GeneGroupSummary(present=[], missing=missing, group_col=group_col)

    sub = adata[:, present]
    X = sub.layers[layer] if layer is not None else sub.X

    mean_all = _col_mean(X)
    pct_all = _col_pct(X, expr_threshold)

    per_gene = pd.DataFrame(
        {'mean_expr': mean_all, 'pct_expressing': pct_all},
        index=pd.Index(present, name='gene'),
    )

    groups, mean_matrix, pct_matrix = [], pd.DataFrame(), pd.DataFrame()
    if group_col is not None:
        if group_col not in adata.obs.columns:
            raise KeyError(f"group_col '{group_col}' not in adata.obs")
        col = adata.obs[group_col]
        if isinstance(col.dtype, pd.CategoricalDtype):
            groups = [g for g in col.cat.categories if (col == g).any()]
        else:
            groups = sorted(map(str, col.dropna().unique()))
        codes = col.astype(str).values

        means, pcts = {}, {}
        for g in groups:
            mask = codes == str(g)
            Xg = X[mask]
            means[g] = _col_mean(Xg)
            pcts[g] = _col_pct(Xg, expr_threshold)
        mean_matrix = pd.DataFrame(means, index=pd.Index(present, name='gene'))
        pct_matrix = pd.DataFrame(pcts, index=pd.Index(present, name='gene'))

        # Peak group = where each gene's mean expression is highest.
        if len(groups):
            per_gene['peak_group'] = mean_matrix.idxmax(axis=1)

    return GeneGroupSummary(
        present=present, missing=missing, group_col=group_col,
        groups=list(groups), per_gene=per_gene,
        mean_matrix=mean_matrix, pct_matrix=pct_matrix,
    )


def marker_panel(markers, cell_types=None, max_per_type=None):
    """Assemble an ordered marker panel for a cell-type sense check.

    A dot plot only reads as a check when the genes stay grouped by the
    cell type they are supposed to mark: annotated types run along one
    axis, their markers along the other, and a correct annotation shows a
    diagonal. Sorting the genes alphabetically, or letting a gene shared
    by two types appear twice, destroys that.

    Parameters
    ----------
    markers : dict
        ``{cell_type: [gene, ...]}``, e.g. ``DEFAULT_MARKERS`` or the
        result of a PanglaoDB / CellMarker2 lookup.
    cell_types : sequence of str, optional
        Which types to include, in the order given. Defaults to every
        key, sorted. Names not in ``markers`` are ignored.
    max_per_type : int, optional
        Keep at most this many markers per type. Panels from a database
        can run to dozens per type, which is unreadable as a plot.

    Returns
    -------
    (genes, owner)
        ``genes`` is the ordered, de-duplicated gene list. ``owner`` maps
        each gene to the first cell type that claimed it -- a gene listed
        under two types is plotted once, under the earlier one, so the
        panel has one row per gene.
    """
    if cell_types is None:
        cell_types = sorted(markers)

    genes, owner = [], {}
    for cell_type in cell_types:
        entries = markers.get(cell_type)
        if not entries:
            continue
        kept = 0
        for gene in entries:
            gene = str(gene)
            if gene in owner:
                continue          # already claimed by an earlier type
            owner[gene] = cell_type
            genes.append(gene)
            kept += 1
            if max_per_type and kept >= max_per_type:
                break
    return genes, owner


def compare_annotations(adata, left: str = 'cell_type',
                        right: str = 'cell_type_atlas'):
    """Compare two label columns cell by cell.

    Built for independent per-study annotation against labels propagated
    from an atlas. Agreement between them is a result worth reporting;
    disagreement is worth understanding before building on either.

    Cells missing a label on either side are excluded from the rate --
    a cell the atlas never saw (because its role is 'exclude') is not a
    disagreement, and counting it as one would understate agreement in
    proportion to how much you excluded.

    Returns
    -------
    dict
        ``{'n_compared', 'n_agree', 'agreement', 'n_unlabelled',
        'crosstab', 'disagreements'}``. ``crosstab`` is left x right
        counts; ``disagreements`` lists ``(left, right, n)`` worst first.
        ``agreement`` is ``None`` when nothing is comparable.
    """
    import pandas as pd

    empty = {'n_compared': 0, 'n_agree': 0, 'agreement': None,
             'n_unlabelled': 0, 'crosstab': pd.DataFrame(),
             'disagreements': []}
    if left not in adata.obs.columns or right not in adata.obs.columns:
        return empty

    a = adata.obs[left].astype('object')
    b = adata.obs[right].astype('object')
    both = a.notna() & b.notna()
    n_unlabelled = int((~both).sum())
    if not both.any():
        return {**empty, 'n_unlabelled': n_unlabelled}

    a, b = a[both].astype(str), b[both].astype(str)
    agree = (a == b)
    crosstab = pd.crosstab(a, b)

    disagreements = []
    for (la, rb), n in crosstab.stack().items():
        if la != rb and n:
            disagreements.append((str(la), str(rb), int(n)))
    disagreements.sort(key=lambda t: -t[2])

    return {
        'n_compared': int(both.sum()),
        'n_agree': int(agree.sum()),
        'agreement': float(agree.mean()),
        'n_unlabelled': n_unlabelled,
        'crosstab': crosstab,
        'disagreements': disagreements,
    }


def match_marker_types(data_types, marker_types):
    """Map each observed cell type to the marker set that describes it.

    Annotation vocabularies rarely match a marker database word for
    word: CellTypist's heart model says "Ventricular Cardiomyocyte"
    where the curated set says "Cardiomyocyte", and "LYVE1+ Macrophage"
    where it says "Macrophage". Picking markers by the database's names
    would mean picking on an axis nobody thinks in.

    Matching is on whole words, not substrings. Raw substrings are far
    too loose -- "T_cell" normalises to "tcell", which is contained in
    "mastcell", so Mast Cell silently acquired T-cell markers. A marker
    set matches when its words are all present in the observed name (or
    the reverse); the most specific match wins, so "Cardiomyocyte" beats
    a hypothetical "Myocyte".

    Returns
    -------
    dict
        ``{observed_type: marker_type}`` for those that matched. Types
        with no match are absent, so a caller can say which of them
        cannot be checked.
    """
    import re

    def words(s):
        return {w for w in re.split(r'[^a-z0-9]+', str(s).lower()) if w}

    marker_words = {m: words(m) for m in marker_types}
    out = {}
    for observed in data_types:
        seen = words(observed)
        if not seen:
            continue
        candidates = [
            m for m, mw in marker_words.items()
            if mw and (mw <= seen or seen <= mw)
        ]
        if candidates:
            out[observed] = max(candidates,
                                key=lambda m: len(marker_words[m] & seen))
    return out
