# Subclustering one population of a study (ADR-008).
#
# The Cluster step embeds and clusters the whole study at one resolution,
# and the genes and components chosen over every cell are the ones that
# separate the major cell types, not the states within one of them. Focus
# mode takes the cells of one or more groups of an obs column (Leiden
# clusters, cell types), gives them their own highly variable genes, PCA,
# Harmony, neighbour graph, UMAP and Leiden clustering, and lets them be
# annotated on their own.
#
# Applying the result writes two things onto the full study, matched by
# cell name: 'leiden_sub', each cell's parent cluster with the focused
# cells split as '<groups>.<subcluster>', and the 'cell_type' labels the
# focused cells were given. The study's 'leiden' and embedding are never
# changed. 'apply_focus' returns the parameters of the 'subcluster'
# provenance stage.
#
# The population is a working subset (see workset.py): its own copy of
# the focused rows of X, log-normalised in place when X holds counts, and
# fresh obs and var frames. Its 'leiden' column is overwritten by the
# subclusters, so each cell's parent group is kept in 'Focus.parent'.
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np
import pandas as pd

# No anndata at import: the Cluster step imports this module while the
# window is built, and a worker started by then may be importing anndata.
from kosmic.scrna.annotate.labels import cluster_column, set_cell_type
from kosmic.scrna.inspect.roles import included_mask

#: The obs column subcluster ids are written to.
SUBCLUSTER_COLUMN = 'leiden_sub'

#: A population smaller than this is not given an embedding of its own.
MIN_CELLS = 50

#: uns key of the hand relabels made before the population was last
#: re-clustered, each with the resolution whose subcluster it names.
EARLIER_EDITS_KEY = 'manual_relabels_earlier'

#: Settings used where neither the caller nor the study gives one.
DEFAULT_SETTINGS = {
    'n_hvg': 2000, 'n_pcs': 30, 'harmony': False, 'harmony_batch_key': None,
    'harmony_merge_across': None, 'scale': False, 'n_neighbors': 15,
    'resolution': 0.5,
}

Progress = Optional[Callable[[str], None]]


def natural_key(label) -> list:
    """Sort key that puts '2' before '10' and '7.2' before '7.10'."""
    parts = re.split(r'(\d+)', str(label))
    return [int(p) if i % 2 else p.lower() for i, p in enumerate(parts)]


@dataclass
class Focus:
    """One population opened in focus mode.

    'adata' is the working subset and 'parent' each of its cells' group
    in 'column', indexed by cell name. 'settings' records what the
    population was embedded and clustered with.
    """

    column: str
    groups: list
    adata: object
    parent: pd.Series
    n_excluded: int = 0
    settings: dict = field(default_factory=dict)

    @property
    def prefix(self) -> str:
        """The groups as they lead each subcluster id: '7', '7+12'."""
        return '+'.join(self.groups)


def order_groups(values: pd.Series, groups: Sequence) -> list:
    """'groups' as strings, without repeats, in the column's order.

    The order of a categorical column is its category order; any other
    column is put in natural order.
    """
    wanted = list(dict.fromkeys(str(g) for g in groups))
    if isinstance(values.dtype, pd.CategoricalDtype):
        rank = {str(c): i for i, c in enumerate(values.cat.categories)}
        return sorted(wanted, key=lambda g: (rank.get(g, len(rank)), natural_key(g)))
    return sorted(wanted, key=natural_key)


def select_population(adata, column: str, groups: Sequence, *,
                      honour_roles: bool = True) -> Focus:
    """The cells whose 'obs[column]' is one of 'groups', as a Focus.

    Cells marked '_role == exclude' are left out when 'honour_roles' is
    set, as in the PCA and cluster steps, and counted in 'n_excluded'.
    Raises ValueError when fewer than MIN_CELLS cells remain.
    """
    from kosmic.scrna.cluster.workset import working_subset

    if column not in adata.obs.columns:
        raise ValueError(f"'{column}' is not a column of this study.")
    values = adata.obs[column]
    groups = order_groups(values, groups)
    in_groups = (values.notna() & values.astype(str).isin(groups)).to_numpy()
    keep = included_mask(adata) if honour_roles else None
    mask = in_groups if keep is None else in_groups & keep
    n = int(mask.sum())
    if n < MIN_CELLS:
        raise ValueError(
            f"{n} cells have {column} = {', '.join(groups)}; at least "
            f"{MIN_CELLS} are needed to subcluster them.")
    sub = working_subset(adata, mask)
    if SUBCLUSTER_COLUMN in sub.obs.columns:
        del sub.obs[SUBCLUSTER_COLUMN]
    parent = pd.Series(values.astype(str).to_numpy()[mask],
                       index=sub.obs_names, name=column)
    return Focus(column, groups, sub, parent,
                 n_excluded=int((in_groups & ~mask).sum()))


def embed_and_cluster(focus: Focus, settings: dict,
                      progress: Progress = None) -> dict:
    """Give the population its own embedding and clustering, in place.

    The steps are those of the PCA and cluster steps, on the subset:
    highly variable genes (within the Harmony batch key when Harmony is
    on), log-normalisation when X holds counts, PCA, Harmony, the
    neighbour graph, UMAP, Leiden, and marker genes ranked between the
    subclusters. With 'scale' the PCA runs on scaled genes, but X itself
    is never scaled, because marker genes and annotation read it.

    'settings' takes the keys of DEFAULT_SETTINGS. Returns the settings
    actually used, which are also stored on 'focus.settings'.
    """
    import scanpy as sc

    from kosmic.scrna.cluster.harmony import run_harmony
    from kosmic.scrna.cluster.hvg import _looks_like_counts, find_hvg
    from kosmic.scrna.cluster.neighbors import compute_neighbors

    say = progress or (lambda _msg: None)
    sub = focus.adata
    s = {**DEFAULT_SETTINGS, **settings}

    batch = s['harmony_batch_key'] if s['harmony'] else None
    if batch and (batch not in sub.obs.columns or sub.obs[batch].nunique() < 2):
        say(f"These cells have fewer than two values of '{batch}', so "
            "Harmony is not run.")
        batch = None

    counts = _looks_like_counts(sub)
    say(f"Selecting up to {int(s['n_hvg']):,} highly variable genes"
        + (f" within '{batch}'..." if batch else "..."))
    find_hvg(sub, n_hvg=int(s['n_hvg']), subset=False, batch_key=batch)
    hvg_batch = sub.uns.get('hvg', {}).get('batch_key', batch) or None
    if batch and not hvg_batch:
        say(f"Some values of '{batch}' have too few of these cells to select "
            "variable genes within each, so they were selected across all "
            "cells.")
    if counts:
        say("Log-normalising these cells' counts (10,000 per cell)...")
        if not np.issubdtype(sub.X.dtype, np.floating):
            sub.X = sub.X.astype(np.float32)
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)

    hvg = sub.var['highly_variable'].to_numpy()
    n_pcs = int(min(s['n_pcs'], sub.n_obs - 1, int(hvg.sum()) - 1))
    say(f"Running PCA ({n_pcs} components)...")
    if s['scale']:
        work = sub[:, hvg].copy()
        sc.pp.scale(work, max_value=10, zero_center=True)
        sc.tl.pca(work, n_comps=n_pcs, zero_center=True)
        sub.obsm['X_pca'] = work.obsm['X_pca']
        sub.uns['pca'] = work.uns['pca']
    else:
        sc.tl.pca(sub, n_comps=n_pcs, zero_center=True)
    sub.obsm.pop('X_pca_harmony', None)
    merge = (s['harmony_merge_across'] or None) if batch else None
    if batch:
        say(f"Running Harmony on '{batch}'...")
        run_harmony(sub, batch_key=batch, condition_key=merge)

    n_neighbors = int(min(s['n_neighbors'], sub.n_obs - 1))
    say(f"Computing {n_neighbors} neighbours ({n_pcs} PCs)...")
    compute_neighbors(sub, n_neighbors=n_neighbors, n_pcs=n_pcs, force=True)
    say("Computing UMAP...")
    sc.tl.umap(sub)

    focus.settings = {
        'n_hvg': int(hvg.sum()),
        'hvg_flavor': str(sub.uns.get('hvg', {}).get('flavor', '')),
        'hvg_batch_key': hvg_batch,
        'scaled': bool(s['scale']),
        'n_pcs': n_pcs,
        'harmony': bool(batch),
        'harmony_batch_key': batch,
        'harmony_merge_across': merge,
        'n_neighbors': n_neighbors,
    }
    recluster(focus, float(s['resolution']), progress)
    return focus.settings


def recluster(focus: Focus, resolution: float, progress: Progress = None) -> int:
    """Run Leiden on the population's graph at 'resolution'.

    The embedding and graph are kept, as when the cluster step is re-run
    at another resolution, and each cell keeps its label. Annotation
    details that name the previous subclusters are dropped; hand relabels
    move to 'uns[EARLIER_EDITS_KEY]' with the resolution whose subclusters
    they name. Returns the number of subclusters.
    """
    from kosmic.scrna.annotate.marker_genes import compute_cluster_marker_genes
    from kosmic.scrna.cluster.leiden import cluster_leiden

    say = progress or (lambda _msg: None)
    sub = focus.adata
    say(f"Leiden clustering (resolution={resolution})...")
    cluster_leiden(sub, resolution=resolution)
    edits = sub.uns.pop('manual_relabels', None)
    if edits:
        sub.uns[EARLIER_EDITS_KEY] = list(sub.uns.get(EARLIER_EDITS_KEY, [])) + [
            {**{k: (v.item() if hasattr(v, 'item') else v) for k, v in zip(edits, row)},
             'resolution': focus.settings.get('resolution')}
            for row in zip(*edits.values())]
    for key in ('cluster_annotation_details', 'cluster_type_source',
                'rank_genes_groups'):
        sub.uns.pop(key, None)
    n = int(sub.obs['leiden'].nunique())
    say(f"Found {n} subclusters; ranking their marker genes...")
    try:
        compute_cluster_marker_genes(sub, cluster_col='leiden')
    except ValueError as exc:   # a single subcluster has nothing to rank against
        say(f"Marker genes not ranked: {exc}")
    focus.settings['resolution'] = round(float(resolution), 3)
    return n


def subcluster_ids(focus: Focus) -> pd.Series:
    """Each cell's subcluster id, '<groups>.<subcluster>'."""
    ids = focus.prefix + '.' + focus.adata.obs['leiden'].astype(str)
    return ids.rename(SUBCLUSTER_COLUMN)


def apply_focus(adata, focus: Focus, *, annotation: Optional[dict] = None,
                manual_edits: Optional[list] = None) -> dict:
    """Write the population's subclusters and labels onto the full study.

    Cells are matched by name; a focused cell no longer in the study is
    skipped and counted. 'leiden_sub' starts as a copy of the study's
    cluster column the first time; the focused cells then get their
    subcluster ids. Each focused cell whose 'cell_type' changed in focus
    mode gets the new label; the rest keep theirs. The focused cells'
    'cell_type_fine', when focus mode gave them one, is written too.

    'annotation' is the settings of the automated annotation run in
    focus mode, as the 'annotate' stage records them, and 'manual_edits'
    the hand relabels, as 'annotate_manual' records them; relabels made
    before a re-clustering are added from 'uns[EARLIER_EDITS_KEY]'. Each
    relabel names its subcluster in full and the resolution it was made
    at. Returns the parameters of the 'subcluster' provenance stage.
    """
    sub = focus.adata
    ids = subcluster_ids(focus).to_numpy(dtype=object)
    rows = adata.obs_names.get_indexer(sub.obs_names)
    found = rows >= 0
    if not found.any():
        raise ValueError("None of the subclustered cells are in the study.")
    hit = rows[found]

    if SUBCLUSTER_COLUMN in adata.obs.columns:
        current = adata.obs[SUBCLUSTER_COLUMN]
    else:
        base = cluster_column(adata.obs)
        current = (adata.obs[base] if base
                   else pd.Series(np.nan, index=adata.obs_names))
    current = current.astype(object).to_numpy()
    values = current.astype(str).astype(object)
    values[pd.isna(current)] = np.nan
    values[hit] = ids[found]
    present = sorted(pd.unique(values[pd.notna(values)]), key=natural_key)
    adata.obs[SUBCLUSTER_COLUMN] = pd.Categorical(values, categories=present)

    n_relabelled = 0
    if 'cell_type' in sub.obs.columns:
        new = sub.obs['cell_type'].astype(object).to_numpy()[found]
        old = (adata.obs['cell_type'].astype(object).to_numpy()[hit]
               if 'cell_type' in adata.obs.columns
               else np.full(len(hit), np.nan, dtype=object))
        changed = pd.notna(new) & (new.astype(str) != old.astype(str))
        for label in pd.unique(new[changed]):
            mask = np.zeros(adata.n_obs, dtype=bool)
            mask[hit[changed & (new == label)]] = True
            n_relabelled += set_cell_type(adata.obs, mask, str(label))

    # The fine labels a CellTypist pass in focus mode gave, so the study's
    # 'cell_type_fine' does not keep the parent cluster's winning label.
    if 'cell_type_fine' in sub.obs.columns:
        # An explicit copy: with pandas copy-on-write, to_numpy() can
        # return a read-only view of the column.
        fine = (adata.obs['cell_type_fine'].to_numpy(dtype=object, copy=True)
                if 'cell_type_fine' in adata.obs.columns
                else np.full(adata.n_obs, np.nan, dtype=object))
        fine[hit] = sub.obs['cell_type_fine'].astype(object).to_numpy()[found]
        adata.obs['cell_type_fine'] = pd.Categorical(fine)

    parent = focus.parent.to_numpy()
    labels = (sub.obs['cell_type'].astype(object).to_numpy()
              if 'cell_type' in sub.obs.columns else None)
    subclusters = {}
    for sid in sorted(set(ids), key=natural_key):
        in_sub = ids == sid
        entry = {}
        if labels is not None:
            counts = pd.Series(labels[in_sub]).dropna().astype(str).value_counts()
            entry['label'] = counts.index[0] if len(counts) else ''
        entry['n_cells'] = int(in_sub.sum())
        if len(focus.groups) > 1:
            entry['from'] = {g: int(n) for g, n in
                             pd.Series(parent[in_sub]).value_counts().items()}
        subclusters[sid] = entry

    edits = [*sub.uns.get(EARLIER_EDITS_KEY, []),
             *({**e, 'resolution': focus.settings.get('resolution')}
               for e in manual_edits or [])]
    params = {
        'column': focus.column,
        'groups': list(focus.groups),
        'n_cells_by_group': {g: int((parent == g).sum()) for g in focus.groups},
        'n_cells': int(sub.n_obs),
        'n_excluded': int(focus.n_excluded),
        **focus.settings,
        'n_subclusters': len(subclusters),
        'annotation': dict(annotation) if annotation else {
            'method': 'manual' if edits else 'none'},
    }
    if edits:
        params['manual_edits'] = [
            {**e, 'cluster': f"{focus.prefix}.{e.get('cluster')}"} for e in edits]
    params['n_cells_relabelled'] = int(n_relabelled)
    if not found.all():
        params['n_cells_not_in_study'] = int((~found).sum())
    params['subclusters'] = subclusters
    return params


def clear_subclusters(adata) -> bool:
    """Drop 'leiden_sub', whose ids mean nothing once 'leiden' is recomputed."""
    if SUBCLUSTER_COLUMN not in adata.obs.columns:
        return False
    del adata.obs[SUBCLUSTER_COLUMN]
    return True


def focus_defaults(adata, record: Optional[dict] = None) -> dict:
    """Settings to start a focus from: what the study was embedded with.

    Taken from the study's last 'embedding' and 'cluster' provenance
    stages ('record' is the loaded record), then from what the study
    itself holds (its HVGs, Harmony embedding, neighbour graph and
    Leiden parameters), then from DEFAULT_SETTINGS.
    """
    from kosmic.provenance import last_stage

    def params(stage):
        entry = last_stage(record, stage) if record else None
        return (entry or {}).get('params') or {}

    emb, clu = params('embedding'), params('cluster')
    graph = (adata.uns.get('neighbors') or {}).get('params') or {}
    leiden = (adata.uns.get('leiden') or {}).get('params') or {}
    n_hvg = emb.get('n_hvg') or (int(adata.var['highly_variable'].sum())
                                 if 'highly_variable' in adata.var.columns else 0)
    return {
        'n_hvg': int(n_hvg or DEFAULT_SETTINGS['n_hvg']),
        'n_pcs': int(graph.get('n_pcs') or DEFAULT_SETTINGS['n_pcs']),
        'harmony': bool(emb.get('harmony', 'X_pca_harmony' in adata.obsm)),
        'harmony_batch_key': emb.get('harmony_batch_key') or None,
        'harmony_merge_across': emb.get('harmony_merge_across') or None,
        'scale': bool(emb.get('scaled', False)),
        'n_neighbors': int(graph.get('n_neighbors') or DEFAULT_SETTINGS['n_neighbors']),
        'resolution': float(clu.get('resolution') or leiden.get('resolution')
                            or DEFAULT_SETTINGS['resolution']),
    }
