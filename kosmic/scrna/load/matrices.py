# Candidate expression-matrix discovery and selection.
#
# An h5ad can carry the same cells' expression several times over: 'X', any
# number of 'layers', and '.raw'. Which one holds *raw integer counts* varies
# by depositor and is not knowable from the file name:
#
#   - CellxGene deposits normalise 'X' and put counts in '.raw'.
#   - Broad SCP's cardiac deposit has CellBender counts in 'X' and the
#     uncorrected CellRanger counts in "layers['cellranger_raw']".
#   - Seurat exports usually have counts in 'X' alone.
#
# Feeding log-normalised values to a count model (DESeq2 pseudobulk) produces
# no error, just wrong answers -- so the choice has to be surfaced rather than
# assumed. This module describes every candidate and names the best one; the
# GUI presents that and lets the user override.
from __future__ import annotations

from typing import Optional

import numpy as np

# Below this maximum, non-integer data is almost certainly log-transformed
# rather than a scaled count (log1p of a 10x UMI count rarely exceeds ~12).
_LOG_MAX = 20.0

RAW_COUNTS = "raw counts"
LOG_NORMALISED = "log-normalised"
UNCLEAR = "unclear"


def _sample_values(matrix, n_sample: int) -> tuple[np.ndarray, float]:
    """Return up to 'n_sample' stored values and the matrix maximum."""
    # 'issparse', not 'hasattr(matrix, "data")': a dense ndarray also has
    # a '.data' attribute, but it is a memoryview of the underlying buffer
    # with no '.size' or '.max()', so the hasattr test sent every dense
    # matrix down the sparse branch and raised.
    import scipy.sparse as sp

    if sp.issparse(matrix):
        data = matrix.data
        if data.size == 0:
            return np.empty(0), 0.0
        return np.asarray(data[:n_sample]), float(data.max())
    arr = np.asarray(matrix)
    if arr.size == 0:
        return np.empty(0), 0.0
    return arr.ravel()[:n_sample], float(arr.max())


def _verdict(values: np.ndarray, max_val: float, log1p_recorded: bool) -> str:
    """Classify a matrix from sampled values."""
    if values.size == 0:
        return UNCLEAR
    is_integer = bool(np.allclose(values, np.round(values)))
    if is_integer:
        return RAW_COUNTS
    if log1p_recorded or max_val < _LOG_MAX:
        return LOG_NORMALISED
    return UNCLEAR


def describe_count_matrices(adata, n_sample: int = 2000) -> list[dict]:
    """Describe every candidate expression matrix in ``adata``.

    Parameters
    ----------
    adata : anndata.AnnData
        Loaded object. Backed objects work but '.raw' may be unavailable.
    n_sample : int
        Stored values sampled per matrix for the integer test.

    Returns
    -------
    list of dict
        One entry per candidate, ``X`` first, with keys:
        ``slot`` (``'X'`` / ``'raw'`` / ``'layers:<name>'`` -- the token
        :func:`promote_matrix` accepts), ``label`` (human-readable),
        ``shape``, ``is_integer``, ``max``, ``n_stored`` and ``verdict``
        (one of :data:`RAW_COUNTS`, :data:`LOG_NORMALISED`, :data:`UNCLEAR`).
    """
    log1p_recorded = "log1p" in getattr(adata, "uns", {})
    out: list[dict] = []

    def add(slot: str, label: str, matrix, shape):
        if matrix is None:
            return
        values, max_val = _sample_values(matrix, n_sample)
        out.append({
            "slot": slot,
            "label": label,
            "shape": tuple(shape),
            "is_integer": (bool(np.allclose(values, np.round(values)))
                           if values.size else False),
            "max": max_val,
            "n_stored": int(matrix.nnz) if hasattr(matrix, "nnz") else int(np.size(matrix)),
            "verdict": _verdict(values, max_val, log1p_recorded),
        })

    add("X", "X", adata.X, adata.shape)
    for key in getattr(adata, "layers", {}) or {}:
        add(f"layers:{key}", f"layers['{key}']", adata.layers[key], adata.shape)
    raw = getattr(adata, "raw", None)
    if raw is not None:
        add("raw", "raw.X", raw.X, raw.shape)
    return out


def best_counts_slot(adata, n_sample: int = 2000) -> Optional[str]:
    """Return the slot token most likely to hold raw counts, or ``None``.

    Prefers ``X`` when it already holds counts so the common case is a no-op.
    Otherwise returns the first other candidate classified as raw counts.
    """
    described = describe_count_matrices(adata, n_sample=n_sample)
    by_slot = {d["slot"]: d for d in described}
    if by_slot.get("X", {}).get("verdict") == RAW_COUNTS:
        return "X"
    for d in described:
        if d["verdict"] == RAW_COUNTS:
            return d["slot"]
    return None


def describe_count_matrices_path(path, n_sample: int = 2000) -> list[dict]:
    """Describe an h5ad's candidate matrices without loading any of them.

    Same output as :func:`describe_count_matrices`, read straight from the
    file. Use this before deciding which matrix to load: a big object can
    hold several copies of the expression data (reichart's CellxGene h5ad
    carries 1.35 billion non-zeros twice over), and loading them all to pick
    one exhausts memory.

    ``max`` is taken from the sampled values rather than the whole matrix,
    which is enough to separate counts from log-transformed data.
    """
    import h5py

    def _sample(group) -> tuple[np.ndarray, float]:
        node = group["data"] if isinstance(group, h5py.Group) and "data" in group else group
        if node.size == 0:
            return np.empty(0), 0.0
        values = np.asarray(node[:n_sample]).ravel()
        nz = values[values != 0]
        if nz.size == 0:
            return np.empty(0), 0.0
        return nz, float(nz.max())

    def _shape(group, fallback) -> tuple:
        shape = group.attrs.get("shape") if hasattr(group, "attrs") else None
        return tuple(int(x) for x in shape) if shape is not None else fallback

    out: list[dict] = []
    with h5py.File(path, "r") as f:
        log1p_recorded = "log1p" in (f["uns"].keys() if "uns" in f else [])

        def add(slot: str, label: str, node, fallback_shape=()):
            if node is None:
                return
            values, max_val = _sample(node)
            out.append({
                "slot": slot,
                "label": label,
                "shape": _shape(node, fallback_shape),
                "is_integer": (bool(np.allclose(values, np.round(values)))
                               if values.size else False),
                "max": max_val,
                "n_stored": int(node["data"].size) if (
                    isinstance(node, h5py.Group) and "data" in node) else 0,
                "verdict": _verdict(values, max_val, log1p_recorded),
            })

        add("X", "X", f.get("X"))
        for key in (f["layers"].keys() if "layers" in f else []):
            add(f"layers:{key}", f"layers['{key}']", f[f"layers/{key}"])
        if "raw" in f and "X" in f["raw"]:
            add("raw", "raw.X", f["raw/X"])
    return out


def read_counts(path, slot: str):
    """Load only the matrix at ``slot`` from an h5ad, as an AnnData.

    Peak memory is one matrix rather than all of them. ``obs`` always comes
    from the top level; ``var`` follows the matrix (``raw`` carries its own,
    which may cover a different gene set).
    """
    import anndata as ad
    import h5py
    from anndata.io import read_elem

    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
        if slot == "X":
            X, var = read_elem(f["X"]), read_elem(f["var"])
        elif slot == "raw":
            if "raw" not in f or "X" not in f["raw"]:
                raise ValueError(f"{path} has no raw.X")
            X = read_elem(f["raw/X"])
            var = (read_elem(f["raw/var"]) if "var" in f["raw"]
                   else read_elem(f["var"]))
        elif slot.startswith("layers:"):
            key = slot.split(":", 1)[1]
            if "layers" not in f or key not in f["layers"]:
                raise ValueError(f"{path} has no layer {key!r}")
            X, var = read_elem(f[f"layers/{key}"]), read_elem(f["var"])
        else:
            raise ValueError(
                f"Unrecognised slot {slot!r}; expected 'X', 'raw' or "
                "'layers:<name>'")

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    # Studies are combined by joining on var_names, so the identifier has to
    # be the same kind everywhere. CellxGene indexes on Ensembl IDs; most
    # other sources index on symbols. Left alone, an inner join across them
    # matches nothing at all.
    from kosmic.scrna.load.converters import normalise_gene_index
    return normalise_gene_index(adata)


def promote_matrix(adata, slot: str, copy: bool = True):
    """Return ``adata`` with the matrix at ``slot`` moved into ``X``.

    Parameters
    ----------
    adata : anndata.AnnData
    slot : str
        ``'X'`` (no-op), ``'raw'``, or ``'layers:<name>'``.
    copy : bool
        Operate on a copy (default) rather than in place.

    Raises
    ------
    ValueError
        If ``slot`` is not a recognised token or is absent from ``adata``.

    Notes
    -----
    Promoting ``'raw'`` rebuilds the object from ``adata.raw``, whose gene set
    may differ from ``adata.var_names``; ``obs`` is carried across but
    ``obsm``/``varm`` are dropped because they are indexed on the old var axis.
    """
    if slot == "X":
        return adata.copy() if copy else adata

    if slot == "raw":
        if getattr(adata, "raw", None) is None:
            raise ValueError("adata has no .raw to promote")
        promoted = adata.raw.to_adata()
        promoted.obs = adata.obs.copy()
        promoted.uns = dict(adata.uns)
        promoted.uns.pop("log1p", None)      # X is no longer log-transformed
        return promoted

    if slot.startswith("layers:"):
        key = slot.split(":", 1)[1]
        if key not in (adata.layers or {}):
            raise ValueError(f"adata has no layer {key!r}")
        promoted = adata.copy() if copy else adata
        promoted.X = promoted.layers[key].copy()
        del promoted.layers[key]
        promoted.uns.pop("log1p", None)
        return promoted

    raise ValueError(
        f"Unrecognised slot {slot!r}; expected 'X', 'raw' or 'layers:<name>'")
